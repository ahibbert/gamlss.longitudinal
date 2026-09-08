#!/usr/bin/env Rscript

# Run one bounded optimizer pass for every immediately eligible pure family.
# This is a compatibility probe, not a convergence or inferential validation.

if (!requireNamespace("devtools", quietly = TRUE)) stop("Install devtools.")
devtools::load_all(quiet = TRUE)

root <- file.path("paper", "reviews", "2026-09-08-family-capability")
audit <- utils::read.csv(
  file.path(root, "gamlss-family-audit.csv"), stringsAsFactors = FALSE
)
audit <- audit[!duplicated(audit$family) & (
  grepl("^candidate_", audit$proposed_route) |
    audit$proposed_route == "requires_binomial_denominator"
), ]
namespace <- asNamespace("gamlss.dist")

default_parameters <- function(fun, parameters) {
  environment <- new.env(parent = namespace)
  out <- list()
  for (parameter in parameters) {
    value <- tryCatch(eval(formals(fun)[[parameter]], environment), error = identity)
    if (inherits(value, "error") || !is.numeric(value) || length(value) != 1L ||
        !is.finite(value)) return(NULL)
    out[[parameter]] <- value
    assign(parameter, value, envir = environment)
  }
  out
}

control <- gamlss_longitudinal_control(
  outer_tol = 1e-8,
  max_outer_iter = 1L,
  max_elapsed_sec = 10,
  stop_on_convergence = FALSE,
  rs = list(
    max_inner_iter = 1L,
    max_steps = 1L,
    warm_start_joint = FALSE,
    use_backtracking = FALSE,
    update_lambda = FALSE,
    discrete_score_method = "analytical"
  )
)

probe <- function(constructor, family) {
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  result <- withCallingHandlers(tryCatch({
    object <- do.call(get(constructor, envir = namespace), list())
    if (family %in% .gl_binomial_families()) {
      object <- .gl_attach_margin_args_spec(object, list(bd = 10))
    }
    object <- .normalise_margin_dist_links(object)
    spec <- .gl_capability_margin_spec(family)
    if (is.null(spec) || !identical(spec$status, "supported")) {
      stop("family is not supported by the current capability registry")
    }
    copula <- spec$compatible_copulas[[1L]]
    parameters <- names(object$parameters)
    qfun <- get(paste0("q", family), envir = namespace)
    defaults <- default_parameters(qfun, parameters)
    if (is.null(defaults)) stop("non-numeric linked parameter default")
    if (identical(family, "PIG2")) {
      defaults$mu <- 5
      defaults$sigma <- 0.1
    }
    if (identical(family, "BB")) defaults$sigma <- 0.1
    fixed <- attr(object, "fixed_unlinked_values")
    if (length(fixed)) defaults <- c(defaults, fixed)

    data <- expand.grid(time = seq_len(4L), subject = seq_len(20L))
    probabilities <- ((data$subject * 7L + data$time * 11L) %% 89L + 0.5) / 90
    quantile_args <- c(list(p = probabilities), defaults)
    if (family %in% .gl_binomial_families()) quantile_args$bd <- 10
    data$response <- do.call(qfun, quantile_args)
    if (any(!is.finite(data$response))) stop("non-finite generated response")

    fit <- gamlss_longitudinal(
      dataset = data,
      margin_dist = object,
      copula_dist = copula,
      time_var = "time",
      subject_var = "subject",
      mu.formula = response ~ factor(time),
      sigma.formula = ~ 1,
      nu.formula = ~ 1,
      tau.formula = ~ 1,
      theta.formula = ~ 1,
      optimizer_control = control,
      compute_vcov = FALSE,
      verbose = 0
    )
    if (!inherits(fit, "gamlss.longitudinal")) stop("fit object was not returned")
    fitted_loglik <- tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_)
    if (!is.finite(fitted_loglik)) stop("fit did not return a finite joint log-likelihood")
    list(success = TRUE, copula = copula, logLik = fitted_loglik, failure = "")
  }, error = function(e) list(
    success = FALSE, copula = NA_character_, logLik = NA_real_,
    failure = conditionMessage(e)
  )),
  warning = function(w) {
    if (!inherits(w, "gamlss.longitudinal_nonconvergence_warning")) {
      warnings <<- c(warnings, conditionMessage(w))
    }
    invokeRestart("muffleWarning")
  })
  data.frame(
    constructor = constructor,
    family = family,
    copula = result$copula,
    success = result$success,
    finite_logLik = is.finite(result$logLik),
    logLik = result$logLik,
    elapsed_sec = proc.time()[["elapsed"]] - started,
    warnings = paste(unique(warnings), collapse = " | "),
    failure = result$failure,
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, Map(probe, audit$constructor, audit$family))
output <- file.path(root, "family-fit-probe.csv")
utils::write.csv(results, output, row.names = FALSE, na = "")
cat("Wrote", nrow(results), "bounded family fits to", output, "\n")
cat("Successful:", sum(results$success), "\n")
cat("Total elapsed seconds:", round(sum(results$elapsed_sec), 2), "\n")
print(results[!results$success | nzchar(results$warnings), ], row.names = FALSE)
