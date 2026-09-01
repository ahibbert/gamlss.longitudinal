#' Supported longitudinal model capabilities
#'
#' Return the versioned, conservative capability registry used by
#' [gamlss_longitudinal()] and the distribution-selection helpers. The registry
#' describes capabilities that are exercised by end-to-end package or paper
#' workflows. The presence of a family in `gamlss.dist`, or of `d`, `p`, and `q`
#' functions, is not by itself a support promise.
#'
#' @param component Registry component to return: supported margin/copula
#'   combinations (`"routes"`), marginal families (`"margins"`), or copulas
#'   (`"copulas"`).
#'
#' @return A documentation-ready data frame. Every row contains the registry
#'   version so exported tables retain their provenance.
#' @export
longitudinal_capabilities <- function(component = c("routes", "margins", "copulas")) {
  component <- match.arg(component)
  switch(component,
    routes = .gl_capability_route_table(),
    margins = .gl_capability_margin_table(),
    copulas = .gl_capability_copula_table()
  )
}

.gl_capability_registry_version <- function() "2026.1"

.gl_capability_all_copulas <- function() c("N", "C", "F", "G", "J", "t")

.gl_capability_margin_specs <- function() {
  all_copulas <- .gl_capability_all_copulas()
  count_note <- paste(
    "One homogeneous count family per fit; exact discrete-discrete rectangle",
    "likelihood; automated omnibus calibration requires randomized PIT work."
  )
  list(
    NO = list(
      status = "supported", family_type = "continuous", response_domain = "real",
      response_requirement = "finite real values", parameters = c("mu", "sigma"),
      likelihood_route = "continuous_density", compatible_copulas = all_copulas,
      hessian = "analytical_and_numerical", randomized_pit = FALSE,
      diagnostics = "continuous_pit_and_pair_diagnostics",
      limitations = "One homogeneous marginal family per fit.", paper_route = FALSE
    ),
    GA = list(
      status = "supported", family_type = "continuous", response_domain = "positive_real",
      response_requirement = "finite values strictly greater than zero", parameters = c("mu", "sigma"),
      likelihood_route = "continuous_density", compatible_copulas = all_copulas,
      hessian = "analytical_and_numerical", randomized_pit = FALSE,
      diagnostics = "continuous_pit_and_pair_diagnostics",
      limitations = "One homogeneous positive continuous family per fit.", paper_route = FALSE
    ),
    GG = list(
      status = "supported", family_type = "continuous", response_domain = "positive_real",
      response_requirement = "finite values strictly greater than zero", parameters = c("mu", "sigma", "nu"),
      likelihood_route = "continuous_density", compatible_copulas = c("N", "C"),
      hessian = "analytical_and_numerical", randomized_pit = FALSE,
      diagnostics = "continuous_pit_and_pair_diagnostics",
      limitations = "Validated with Gaussian and Clayton copulas only.", paper_route = FALSE
    ),
    BCPE = list(
      status = "supported", family_type = "continuous", response_domain = "positive_real",
      response_requirement = "finite values strictly greater than zero", parameters = c("mu", "sigma", "nu", "tau"),
      likelihood_route = "continuous_density", compatible_copulas = c("N", "t"),
      hessian = "analytical_and_numerical", randomized_pit = FALSE,
      diagnostics = "continuous_pit_and_pair_diagnostics",
      limitations = "Validated with Gaussian and Student-t copulas; BCPE/t is the public paper route.", paper_route = TRUE
    ),
    LOGNO = list(
      status = "supported", family_type = "continuous", response_domain = "positive_real",
      response_requirement = "finite values strictly greater than zero", parameters = c("mu", "sigma"),
      likelihood_route = "continuous_density", compatible_copulas = "N",
      hessian = "analytical_and_numerical", randomized_pit = FALSE,
      diagnostics = "continuous_pit_and_pair_diagnostics",
      limitations = "Validated with the Gaussian copula only.", paper_route = FALSE
    ),
    PO = list(
      status = "supported", family_type = "discrete", response_domain = "count",
      response_requirement = "finite non-negative integers", parameters = "mu",
      likelihood_route = "exact_discrete_rectangle", compatible_copulas = c("N", "C"),
      hessian = "analytical_and_numerical", randomized_pit = TRUE,
      diagnostics = "randomized_pit_opt_in_and_descriptive_pair_diagnostics",
      limitations = count_note, paper_route = FALSE
    ),
    NBI = list(
      status = "supported", family_type = "discrete", response_domain = "count",
      response_requirement = "finite non-negative integers", parameters = c("mu", "sigma"),
      likelihood_route = "exact_discrete_rectangle", compatible_copulas = c("N", "C"),
      hessian = "analytical_and_numerical", randomized_pit = TRUE,
      diagnostics = "randomized_pit_opt_in_and_descriptive_pair_diagnostics",
      limitations = paste(count_note, "NBI/Clayton is the public paper route."), paper_route = TRUE
    ),
    DEL = list(
      status = "supported", family_type = "discrete", response_domain = "count",
      response_requirement = "finite non-negative integers", parameters = c("mu", "sigma", "nu"),
      likelihood_route = "exact_discrete_rectangle", compatible_copulas = c("N", "C"),
      hessian = "analytical_and_numerical", randomized_pit = TRUE,
      diagnostics = "randomized_pit_opt_in_and_descriptive_pair_diagnostics",
      limitations = count_note, paper_route = FALSE
    ),
    ZIP = list(
      status = "supported", family_type = "discrete", response_domain = "count",
      response_requirement = "finite non-negative integers", parameters = c("mu", "sigma"),
      likelihood_route = "exact_discrete_rectangle", compatible_copulas = "N",
      hessian = "analytical_and_numerical", randomized_pit = TRUE,
      diagnostics = "randomized_pit_opt_in_and_descriptive_pair_diagnostics",
      limitations = paste(count_note, "Validated with the Gaussian copula only."), paper_route = FALSE
    ),
    ZAP = list(
      status = "supported", family_type = "discrete", response_domain = "count",
      response_requirement = "finite non-negative integers", parameters = c("mu", "sigma"),
      likelihood_route = "exact_discrete_rectangle", compatible_copulas = "N",
      hessian = "analytical_and_numerical", randomized_pit = TRUE,
      diagnostics = "randomized_pit_opt_in_and_descriptive_pair_diagnostics",
      limitations = paste(count_note, "Validated with the Gaussian copula only."), paper_route = FALSE
    ),
    ZINBI = list(
      status = "supported", family_type = "discrete", response_domain = "count",
      response_requirement = "finite non-negative integers", parameters = c("mu", "sigma", "nu"),
      likelihood_route = "exact_discrete_rectangle", compatible_copulas = "N",
      hessian = "analytical_and_numerical", randomized_pit = TRUE,
      diagnostics = "randomized_pit_opt_in_and_descriptive_pair_diagnostics",
      limitations = paste(count_note, "Validated with the Gaussian copula only."), paper_route = FALSE
    )
  )
}

.gl_capability_denied_margin_specs <- function() {
  binomial_reason <- paste(
    "Bounded/binomial responses require an explicit denominator contract;",
    "the longitudinal likelihood must not assume Bernoulli denominators."
  )
  multinomial_reason <- "Ordinal and multinomial outcomes require a family-specific response and likelihood route."
  specs <- lapply(c("BI", "BB", "DBI", "ZABB", "ZABI", "ZIBB", "ZIBI"), function(family) {
    list(
      status = "unsupported", family_type = "bounded_binomial", response_domain = "bounded_count",
      response_requirement = "response plus an explicit denominator", parameters = character(),
      likelihood_route = "unsupported", compatible_copulas = character(), hessian = "not_available",
      randomized_pit = FALSE, diagnostics = "not_available", limitations = binomial_reason,
      paper_route = FALSE
    )
  })
  names(specs) <- c("BI", "BB", "DBI", "ZABB", "ZABI", "ZIBB", "ZIBI")
  ordinal <- lapply(c("MN3", "MN4", "MN5"), function(family) {
    list(
      status = "unsupported", family_type = "ordinal_multinomial", response_domain = "category",
      response_requirement = "family-specific categorical response", parameters = character(),
      likelihood_route = "unsupported", compatible_copulas = character(), hessian = "not_available",
      randomized_pit = FALSE, diagnostics = "not_available", limitations = multinomial_reason,
      paper_route = FALSE
    )
  })
  names(ordinal) <- c("MN3", "MN4", "MN5")
  c(specs, ordinal)
}

.gl_capability_copula_specs <- function() {
  list(
    N = list(
      status = "supported", parameters = "theta", links = "theta=fisher_z",
      inverse_links = "theta=tanh(eta)", parameter_domain = "-1 < theta < 1",
      limitations = "Numerical evaluation is bounded away from +/-1."
    ),
    C = list(
      status = "supported", parameters = "theta", links = "theta=log",
      inverse_links = "theta=exp(eta)", parameter_domain = "theta > 0",
      limitations = "Positive dependence only; independence is approached at the boundary."
    ),
    F = list(
      status = "supported", parameters = "theta", links = "theta=identity",
      inverse_links = "theta=eta", parameter_domain = "theta is real; theta=0 is independence",
      limitations = "Near-independence evaluation uses dedicated numerical handling."
    ),
    G = list(
      status = "supported", parameters = "theta", links = "theta=log(theta-1)",
      inverse_links = "theta=1+exp(eta)", parameter_domain = "1 <= theta <= 17 (implemented bound)",
      limitations = "Positive dependence only; the implemented inverse link caps theta at 17."
    ),
    J = list(
      status = "supported", parameters = "theta", links = "theta=log(theta-1)",
      inverse_links = "theta=1+exp(eta)", parameter_domain = "theta >= 1",
      limitations = "Positive dependence only; independence is at the boundary."
    ),
    t = list(
      status = "supported", parameters = c("theta", "zeta"),
      links = "theta=fisher_z; zeta=log(df-2)",
      inverse_links = "theta=tanh(eta); df=2+exp(zeta)",
      parameter_domain = "-1 < theta < 1; df > 2",
      limitations = "Both correlation and degrees of freedom are fitted."
    )
  )
}

.gl_capability_margin_code <- function(margin_dist) {
  if (is.character(margin_dist)) {
    if (length(margin_dist) == 1L && !is.na(margin_dist) && nzchar(margin_dist)) {
      return(as.character(margin_dist))
    }
    return(NA_character_)
  }
  if (inherits(margin_dist, c("margin_selection", "margin_screen"))) {
    selected <- best_fit_family(margin_dist)
    if (is.null(selected)) return(NA_character_)
    margin_dist <- selected
  }
  if (is.list(margin_dist) && !is.null(margin_dist$family)) {
    code <- as.character(margin_dist$family[[1L]])
    if (length(code) == 1L && !is.na(code) && nzchar(code)) return(code)
  }
  NA_character_
}

.gl_capability_margin_spec <- function(family) {
  family <- as.character(family)[1L]
  supported <- .gl_capability_margin_specs()
  if (!is.na(family) && family %in% names(supported)) return(supported[[family]])
  denied <- .gl_capability_denied_margin_specs()
  if (!is.na(family) && family %in% names(denied)) return(denied[[family]])
  NULL
}

.gl_capability_copula_spec <- function(copula) {
  copula <- as.character(copula)[1L]
  specs <- .gl_capability_copula_specs()
  if (!is.na(copula) && copula %in% names(specs)) specs[[copula]] else NULL
}

.gl_capability_family_object <- function(family) {
  family_fun <- tryCatch(
    get(family, envir = asNamespace("gamlss.dist"), mode = "function", inherits = FALSE),
    error = function(e) NULL
  )
  if (is.null(family_fun)) return(NULL)
  tryCatch(do.call(family_fun, list()), error = function(e) NULL)
}

.gl_capability_default_link_metadata <- function(family, parameters) {
  object <- .gl_capability_family_object(family)
  if (is.null(object) || length(parameters) == 0L) {
    return(list(links = "", inverse_links = ""))
  }
  links <- vapply(parameters, function(parameter) {
    value <- object[[paste0(parameter, ".link")]]
    if (is.null(value)) "unavailable" else as.character(value)[1L]
  }, character(1))
  inverse <- vapply(parameters, function(parameter) {
    if (is.function(object[[paste0(parameter, ".linkinv")]])) {
      paste0(parameter, ".linkinv")
    } else {
      "unavailable"
    }
  }, character(1))
  list(
    links = paste(paste(parameters, links, sep = "="), collapse = "; "),
    inverse_links = paste(inverse, collapse = "; ")
  )
}

.gl_capability_margin_table <- function() {
  specs <- c(.gl_capability_margin_specs(), .gl_capability_denied_margin_specs())
  rows <- lapply(names(specs), function(family) {
    spec <- specs[[family]]
    parameters <- spec$parameters
    if (length(parameters) == 0L) {
      object <- .gl_capability_family_object(family)
      parameters <- if (is.null(object)) character() else names(object$parameters)
    }
    link_meta <- .gl_capability_default_link_metadata(family, parameters)
    data.frame(
      registry_version = .gl_capability_registry_version(),
      family = family,
      status = spec$status,
      family_type = spec$family_type,
      response_domain = spec$response_domain,
      response_requirement = spec$response_requirement,
      parameters = paste(parameters, collapse = ","),
      default_links = link_meta$links,
      inverse_links = link_meta$inverse_links,
      likelihood_route = spec$likelihood_route,
      compatible_copulas = paste(spec$compatible_copulas, collapse = ","),
      hessian = spec$hessian,
      randomized_pit = isTRUE(spec$randomized_pit),
      diagnostics = spec$diagnostics,
      paper_route = isTRUE(spec$paper_route),
      limitations = spec$limitations,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.gl_capability_copula_table <- function() {
  specs <- .gl_capability_copula_specs()
  rows <- lapply(names(specs), function(copula) {
    spec <- specs[[copula]]
    data.frame(
      registry_version = .gl_capability_registry_version(),
      copula = copula,
      status = spec$status,
      parameters = paste(spec$parameters, collapse = ","),
      links = spec$links,
      inverse_links = spec$inverse_links,
      parameter_domain = spec$parameter_domain,
      limitations = spec$limitations,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.gl_capability_route_table <- function() {
  margin_specs <- .gl_capability_margin_specs()
  copula_specs <- .gl_capability_copula_specs()
  rows <- list()
  for (family in names(margin_specs)) {
    margin <- margin_specs[[family]]
    for (copula in margin$compatible_copulas) {
      copula_spec <- copula_specs[[copula]]
      rows[[length(rows) + 1L]] <- data.frame(
        registry_version = .gl_capability_registry_version(),
        margin_family = family,
        family_type = margin$family_type,
        response_domain = margin$response_domain,
        likelihood_route = margin$likelihood_route,
        copula = copula,
        copula_parameters = paste(copula_spec$parameters, collapse = ","),
        hessian = margin$hessian,
        randomized_pit = isTRUE(margin$randomized_pit),
        diagnostics = margin$diagnostics,
        paper_route = (identical(family, "BCPE") && identical(copula, "t")) ||
          (identical(family, "NBI") && identical(copula, "C")),
        limitations = margin$limitations,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.gl_capability_stop <- function(class, message, ...) {
  condition <- structure(
    c(list(message = message, call = NULL), list(...)),
    class = c(class, "gamlss_longitudinal_capability_error", "error", "condition")
  )
  stop(condition)
}

.gl_capability_multiple_margin_input <- function(margin_dist) {
  is.list(margin_dist) && is.null(margin_dist$family) && length(margin_dist) > 0L &&
    all(vapply(margin_dist, function(x) is.list(x) && !is.null(x$family), logical(1)))
}

.gl_validate_capability_response <- function(response, family, spec) {
  observed <- response[!is.na(response)]
  if (length(observed) == 0L || any(!is.finite(observed))) return(invisible(TRUE))
  valid <- switch(spec$response_domain,
    real = TRUE,
    positive_real = all(observed > 0),
    count = all(observed >= 0) && all(abs(observed - round(observed)) <= sqrt(.Machine$double.eps)),
    FALSE
  )
  if (!isTRUE(valid)) {
    .gl_capability_stop(
      "gamlss_longitudinal_response_domain_error",
      paste0(
        "Margin family '", family, "' requires ", spec$response_requirement,
        ". The observed response does not satisfy this registered domain."
      ),
      margin_family = family,
      response_domain = spec$response_domain
    )
  }
  invisible(TRUE)
}

.gl_validate_capability_route <- function(margin_dist, copula_dist, response = NULL, context = "fit") {
  if (.gl_capability_multiple_margin_input(margin_dist) ||
    (is.character(margin_dist) && length(margin_dist) > 1L)) {
    .gl_capability_stop(
      "gamlss_longitudinal_mixed_margin_error",
      "A fit must use one homogeneous marginal family across all visits; mixed continuous/discrete or family-by-visit inputs are unsupported."
    )
  }
  family <- .gl_capability_margin_code(margin_dist)
  if (is.na(family)) {
    .gl_capability_stop(
      "gamlss_longitudinal_unsupported_margin_error",
      "'margin_dist' must be one registered gamlss.dist family object. Use longitudinal_capabilities('margins') to inspect the allowlist."
    )
  }
  spec <- .gl_capability_margin_spec(family)
  supplied_type <- if (is.list(margin_dist)) tolower(paste(as.character(margin_dist$type), collapse = " ")) else ""
  if (grepl("mixed", supplied_type, fixed = TRUE)) {
    .gl_capability_stop(
      "gamlss_longitudinal_mixed_margin_error",
      paste0("Margin family '", family, "' has mixed continuous/discrete support, which has no validated longitudinal pair likelihood."),
      margin_family = family
    )
  }
  if (is.null(spec) || !identical(spec$status, "supported")) {
    reason <- if (is.null(spec)) {
      "it is not in the conservative end-to-end allowlist"
    } else {
      spec$limitations
    }
    .gl_capability_stop(
      "gamlss_longitudinal_unsupported_margin_error",
      paste0(
        "Margin family '", family, "' is not supported for longitudinal fitting because ", reason,
        ". Use longitudinal_capabilities('margins') for supported families."
      ),
      margin_family = family
    )
  }
  if (is.list(margin_dist)) {
    supplied_parameters <- names(margin_dist$parameters)
    if (!setequal(supplied_parameters, spec$parameters)) {
      .gl_capability_stop(
        "gamlss_longitudinal_unsupported_margin_parameter_error",
        paste0(
          "Margin family '", family, "' must expose the registered parameters ",
          paste(spec$parameters, collapse = ", "), "."
        ),
        margin_family = family
      )
    }
    linked <- vapply(spec$parameters, .margin_parameter_has_link, logical(1), margin_dist = margin_dist)
    if (!all(linked)) {
      .gl_capability_stop(
        "gamlss_longitudinal_unsupported_link_error",
        paste0(
          "Margin family '", family, "' lacks link, inverse-link, or derivative functions for: ",
          paste(spec$parameters[!linked], collapse = ", "), "."
        ),
        margin_family = family,
        parameters = spec$parameters[!linked]
      )
    }
  }
  copula <- as.character(copula_dist)
  if (length(copula) != 1L || is.na(copula) || is.null(.gl_capability_copula_spec(copula))) {
    .gl_capability_stop(
      "gamlss_longitudinal_unsupported_copula_error",
      paste0(
        "Copula family '", paste(copula, collapse = ","),
        "' is unsupported. Use one of: ", paste(.gl_capability_all_copulas(), collapse = ", "), "."
      )
    )
  }
  if (!copula %in% spec$compatible_copulas) {
    .gl_capability_stop(
      "gamlss_longitudinal_unsupported_route_error",
      paste0(
        "The ", family, "/", copula, " margin-copula route is not in the tested allowlist. ",
        "Registered copulas for ", family, " are: ", paste(spec$compatible_copulas, collapse = ", "), "."
      ),
      margin_family = family,
      copula = copula,
      context = context
    )
  }
  for (prefix in c("d", "p", "q")) {
    function_name <- paste0(prefix, family)
    if (!exists(function_name, envir = asNamespace("gamlss.dist"), mode = "function", inherits = FALSE)) {
      .gl_capability_stop(
        "gamlss_longitudinal_missing_family_function_error",
        paste0("Registered family '", family, "' is unavailable because gamlss.dist::", function_name, "() was not found."),
        margin_family = family,
        function_name = function_name
      )
    }
  }
  if (!is.null(response)) .gl_validate_capability_response(response, family, spec)
  invisible(list(
    registry_version = .gl_capability_registry_version(),
    margin_family = family,
    copula = copula,
    family_type = spec$family_type,
    likelihood_route = spec$likelihood_route,
    diagnostics = spec$diagnostics
  ))
}

.gl_preflight_fit_capabilities <- function(dataset, margin_dist, copula_dist, mu.formula) {
  response <- NULL
  response_vars <- tryCatch(all.vars(stats::as.formula(mu.formula)), error = function(e) character())
  if (length(response_vars) > 0L && is.data.frame(dataset) && response_vars[[1L]] %in% names(dataset)) {
    response <- dataset[[response_vars[[1L]]]]
  }
  .gl_validate_capability_route(
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    response = response,
    context = "fit"
  )
}

.gl_capability_likelihood_route <- function(margin_dist) {
  family <- .gl_capability_margin_code(margin_dist)
  spec <- .gl_capability_margin_spec(family)
  if (is.null(spec) || !identical(spec$status, "supported")) return(NA_character_)
  spec$likelihood_route
}

.gl_capability_route_supported <- function(family, copula) {
  spec <- .gl_capability_margin_spec(family)
  !is.null(spec) && identical(spec$status, "supported") && copula %in% spec$compatible_copulas
}
