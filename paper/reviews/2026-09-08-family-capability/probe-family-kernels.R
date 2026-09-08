#!/usr/bin/env Rscript

# Probe the vectorised distribution and score interfaces required by the
# longitudinal likelihood. Run audit-gamlss-families.R first.

if (!requireNamespace("devtools", quietly = TRUE)) stop("Install devtools.")
devtools::load_all(quiet = TRUE)

root <- file.path("paper", "reviews", "2026-09-08-family-capability")
audit_path <- file.path(root, "gamlss-family-audit.csv")
if (!file.exists(audit_path)) stop("Run audit-gamlss-families.R first.")
audit <- utils::read.csv(audit_path, stringsAsFactors = FALSE)
audit <- audit[!duplicated(audit$family) & (
  grepl("^candidate_", audit$proposed_route) |
    audit$proposed_route == "requires_binomial_denominator"
), ]
namespace <- asNamespace("gamlss.dist")

evaluate_defaults <- function(fun, parameters) {
  defaults <- list()
  environment <- new.env(parent = namespace)
  for (parameter in parameters) {
    value <- tryCatch(eval(formals(fun)[[parameter]], environment), error = identity)
    if (inherits(value, "error") || !is.numeric(value) || length(value) != 1L ||
        !is.finite(value)) return(NULL)
    defaults[[parameter]] <- value
    assign(parameter, value, envir = environment)
  }
  defaults
}

row_args <- function(args, index) {
  lapply(args, function(value) if (length(value) == 1L) value else value[[index]])
}

probe <- function(constructor, family) {
  warnings <- character()
  result <- withCallingHandlers(tryCatch({
    object <- do.call(get(constructor, envir = namespace), list())
    object <- .normalise_margin_dist_links(object)
    parameters <- names(object$parameters)
    dfun <- get(paste0("d", family), envir = namespace)
    pfun <- get(paste0("p", family), envir = namespace)
    qfun <- get(paste0("q", family), envir = namespace)
    defaults <- evaluate_defaults(qfun, parameters)
    if (is.null(defaults)) stop("non-numeric parameter default")
    fixed <- attr(object, "fixed_unlinked_values")
    if (length(fixed)) defaults <- c(defaults, fixed)
    if (family %in% .gl_binomial_families()) fixed <- c(fixed, list(bd = 10))

    n <- 5L
    offsets <- seq(-0.1, 0.1, length.out = n)
    parameters_varying <- lapply(parameters, function(parameter) {
      linkfun <- object[[paste0(parameter, ".linkfun")]]
      linkinv <- object[[paste0(parameter, ".linkinv")]]
      value <- defaults[[parameter]]
      as.numeric(linkinv(linkfun(value) + offsets))
    })
    names(parameters_varying) <- parameters
    probabilities <- seq(0.2, 0.8, length.out = n)
    quantile_args <- c(list(p = probabilities), parameters_varying, fixed)
    response <- do.call(qfun, quantile_args)
    if (length(response) != n || any(!is.finite(response))) {
      stop("quantile function did not return finite vector output")
    }
    if (!.gl_capability_family_response_matches(response, family, object)) {
      stop("quantile output failed family y.valid()")
    }

    density_args <- c(list(x = response), parameters_varying, fixed)
    probability_args <- c(list(q = response), parameters_varying, fixed)
    density_vector <- do.call(dfun, density_args)
    probability_vector <- do.call(pfun, probability_args)
    density_scalar <- vapply(seq_len(n), function(ii) {
      do.call(dfun, row_args(density_args, ii))
    }, numeric(1))
    probability_scalar <- vapply(seq_len(n), function(ii) {
      do.call(pfun, row_args(probability_args, ii))
    }, numeric(1))

    derivative_suffix <- c(mu = "m", sigma = "d", nu = "v", tau = "t")
    derivative_finite <- vapply(parameters, function(parameter) {
      derivative <- object[[paste0("dld", derivative_suffix[[parameter]])]]
      if (!is.function(derivative)) return(FALSE)
      derivative_args <- c(list(y = response), parameters_varying, fixed)
      needed <- intersect(names(derivative_args), names(formals(derivative)))
      value <- do.call(derivative, derivative_args[needed])
      length(value) == n && all(is.finite(value))
    }, logical(1))

    list(
      success = TRUE,
      parameters = paste(parameters, collapse = ","),
      vector_density_max_abs = max(abs(density_vector - density_scalar)),
      vector_probability_max_abs = max(abs(probability_vector - probability_scalar)),
      all_derivatives_finite = all(derivative_finite),
      failed_derivatives = paste(parameters[!derivative_finite], collapse = ","),
      failure = ""
    )
  }, error = function(e) list(
    success = FALSE,
    parameters = "",
    vector_density_max_abs = NA_real_,
    vector_probability_max_abs = NA_real_,
    all_derivatives_finite = FALSE,
    failed_derivatives = "",
    failure = conditionMessage(e)
  )), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  data.frame(
    constructor = constructor,
    family = family,
    success = result$success,
    parameters = result$parameters,
    vector_density_max_abs = result$vector_density_max_abs,
    vector_probability_max_abs = result$vector_probability_max_abs,
    all_derivatives_finite = result$all_derivatives_finite,
    failed_derivatives = result$failed_derivatives,
    warnings = paste(unique(warnings), collapse = " | "),
    failure = result$failure,
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, Map(probe, audit$constructor, audit$family))
output <- file.path(root, "family-kernel-probe.csv")
utils::write.csv(results, output, row.names = FALSE, na = "")
cat("Wrote", nrow(results), "family probes to", output, "\n")
cat("Successful:", sum(results$success), "\n")
cat("Vector discrepancies:", sum(
  results$vector_density_max_abs > 1e-10 |
    results$vector_probability_max_abs > 1e-10,
  na.rm = TRUE
), "\n")
cat("Derivative failures:", sum(!results$all_derivatives_finite), "\n")
print(results[!results$success | !results$all_derivatives_finite |
  results$vector_density_max_abs > 1e-10 |
  results$vector_probability_max_abs > 1e-10, ], row.names = FALSE)
