#' Audit a fitted longitudinal GAMLSS-copula model specification
#'
#' `model_spec()` returns a compact, printable audit trail for a fitted model:
#' formulas, distributions, links, optimisation settings, missing-response
#' counts, likelihood components, and variance-covariance metadata.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#'
#' @return A list of class `gamlss_longitudinal_model_spec`.
#' @export
model_spec <- function(object) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }

  family <- object$margin_dist
  link_names <- names(family)[grepl("\\.link$", names(family))]
  margin_links <- data.frame(
    parameter = sub("\\.link$", "", link_names),
    link = as.character(unlist(family[link_names], use.names = FALSE)),
    stringsAsFactors = FALSE
  )
  if (nrow(margin_links) == 0L) {
    margin_links <- data.frame(parameter = character(0), link = character(0))
  }

  response <- object$response
  model_selection <- tryCatch(summary(object, include_vcov = FALSE)$fit$model_selection, error = function(e) NULL)
  dep <- tryCatch(copula_time_summary(object), error = function(e) NULL)

  out <- list(
    variables = list(
      response = object$response_var %||% "response",
      subject = object$subject_var %||% "subject",
      time = object$time_var %||% "time"
    ),
    distributions = list(
      margin = as.character(object$margin_dist$family[1]),
      margin_label = as.character(object$margin_dist$family[2] %||% object$margin_dist$family[1]),
      copula = as.character(object$copula_dist)
    ),
    formulas = object$formulas,
    margin_links = margin_links,
    optimisation = list(
      method = object$optim_method %||% object$convergence$method %||% NA_character_,
      converged = isTRUE(object$convergence$converged),
      stop_reason = object$convergence$stop_reason %||% NA_character_,
      outer_iterations = object$convergence$outer_iterations %||% NA_integer_,
      max_outer_iter = object$convergence$max_outer_iter %||% NA_integer_
    ),
    missingness = list(
      n_rows = length(response),
      n_missing_response = sum(is.na(response)),
      n_nonfinite_response = sum(!is.na(response) & !is.finite(response))
    ),
    likelihood = model_selection,
    vcov = c(object$vcov_meta %||% list(precomputed = FALSE), list(
      hessian_diagnostics = object$vcov$hessian_diagnostics %||% NULL
    )),
    dependence = if (!is.null(dep)) dep$time_summary else NULL
  )
  class(out) <- "gamlss_longitudinal_model_spec"
  out
}

#' @export
print.gamlss_longitudinal_model_spec <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nGAMLSS Longitudinal Model Specification\n")
  cat("--------------------------------------\n")
  cat("Response:", x$variables$response, " | Subject:", x$variables$subject, " | Time:", x$variables$time, "\n")
  cat("Margin:", x$distributions$margin, " | Copula:", x$distributions$copula, "\n")
  cat("Optimisation:", x$optimisation$method, " | Converged:", if (isTRUE(x$optimisation$converged)) "yes" else "no", "\n")
  cat("Missing responses:", x$missingness$n_missing_response, "of", x$missingness$n_rows, "\n")
  if (!is.null(x$vcov)) {
    cat("VCOV:", x$vcov$method_used %||% x$vcov$method %||% "not precomputed", "\n")
  }

  if (!is.null(x$formulas) && length(x$formulas) > 0L) {
    cat("\nFormulas\n")
    print(vapply(x$formulas, function(f) paste(deparse(f), collapse = " "), character(1)))
  }
  if (!is.null(x$margin_links) && nrow(x$margin_links) > 0L) {
    cat("\nMargin Links\n")
    print(x$margin_links, row.names = FALSE)
  }
  if (!is.null(x$likelihood)) {
    cat("\nLikelihood And Information Criteria\n")
    print(round(x$likelihood, digits))
  }
  if (!is.null(x$dependence) && nrow(x$dependence) > 0L) {
    cat("\nDependence By Time (includes Kendall's tau)\n")
    print(x$dependence, digits = digits, row.names = FALSE)
  }
  invisible(x)
}
