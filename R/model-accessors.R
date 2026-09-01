#' Log-likelihood for a fitted longitudinal GAMLSS-copula model
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return A scalar object of class `logLik` for the fitted joint likelihood.
#'   Its `df` attribute is the joint effective degrees of freedom and its
#'   `nobs` attribute is the total number of observed responses.
#' @export
logLik.gamlss.longitudinal <- function(object, ...) {
  value <- .gl_joint_loglik(object)
  if (!is.finite(value)) {
    stop("The fitted object does not contain a finite joint log-likelihood.", call. = FALSE)
  }
  structure(
    value,
    class = "logLik",
    df = .gl_model_edf(object),
    nobs = stats::nobs(object),
    likelihood_contract = object$likelihood_contract %||% list(objective = "ordinary")
  )
}

#' Information criteria for fitted longitudinal GAMLSS-copula models
#'
#' AIC and BIC use the fitted joint likelihood, joint effective degrees of
#' freedom, and the total number of observed responses. For
#' `missingness = "segment"`, they describe the explicit segmented likelihood
#' in which observations separated by gaps are independent.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param ... Additional fitted models for `AIC()`; passed to [stats::AIC()].
#' @param k Penalty multiplier.
#' @return A numeric criterion, or the standard comparison table when multiple
#'   models are supplied to `AIC()`.
#' @export
AIC.gamlss.longitudinal <- function(object, ..., k = 2) {
  .gl_require_converged_fit(object, "AIC model comparison")
  others <- list(...)
  if (length(others)) {
    invisible(lapply(others, .gl_require_converged_fit, operation = "AIC model comparison"))
  }
  NextMethod("AIC")
}

#' @rdname AIC.gamlss.longitudinal
#' @export
BIC.gamlss.longitudinal <- function(object, ...) {
  .gl_require_converged_fit(object, "BIC model comparison")
  others <- list(...)
  if (length(others)) {
    invisible(lapply(others, .gl_require_converged_fit, operation = "BIC model comparison"))
  }
  NextMethod("BIC")
}

#' Coefficients for a fitted longitudinal GAMLSS-copula model
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return Named coefficient vector.
#' @export
coef.gamlss.longitudinal <- function(object, ...) {
  return(object$par)
}

#' Access components of a fitted longitudinal GAMLSS-copula model
#'
#' These S3 methods expose the standard regression components expected by
#' model-auditing workflows: formulas, terms, observation counts, model frames,
#' fitted values, and residuals. They return the expanded internal data by
#' default because `gamlss_longitudinal()` represents structurally missing
#' subject-time combinations explicitly.
#'
#' @param x,object,formula A fitted `gamlss.longitudinal` object.
#' @param parameter Distributional parameter to extract. Defaults to `"mu"`.
#' @param internal Logical; return internally translated formulas when `TRUE`.
#' @param type Observation-count, model-frame, or residual type.
#' @param finite Logical; restrict fitted values or residuals to finite observed
#'   responses and finite fitted parameters.
#' @param ... Additional arguments reserved for future methods.
#'
#' @return A formula, terms object, integer count, data frame, or numeric vector.
#' @name gamlss_longitudinal_accessors
NULL

#' @rdname gamlss_longitudinal_accessors
#' @export
formula.gamlss.longitudinal <- function(x, parameter = c("mu", "sigma", "nu", "tau", "theta", "zeta"), internal = FALSE, ...) {
  parameter <- match.arg(parameter)
  formulas <- if (isTRUE(internal)) x$formulas_int else x$formulas
  fml <- formulas[[parameter]]
  if (inherits(fml, "formula")) {
    return(fml)
  }
  stats::as.formula(fml)
}

#' @rdname gamlss_longitudinal_accessors
#' @export
terms.gamlss.longitudinal <- function(x, parameter = c("mu", "sigma", "nu", "tau", "theta", "zeta"), internal = FALSE, ...) {
  stats::terms(formula.gamlss.longitudinal(x, parameter = parameter, internal = internal, ...))
}

#' @rdname gamlss_longitudinal_accessors
#' @export
nobs.gamlss.longitudinal <- function(object, type = c("observed", "expanded", "submitted"), ...) {
  type <- match.arg(type)
  if (identical(type, "submitted") && !is.null(object$dataset_original)) {
    return(nrow(object$dataset_original))
  }
  if (identical(type, "expanded") && !is.null(object$dataset)) {
    return(nrow(object$dataset))
  }
  sum(is.finite(object$response))
}

#' @rdname gamlss_longitudinal_accessors
#' @export
model.frame.gamlss.longitudinal <- function(formula, type = c("expanded", "observed", "submitted"), ...) {
  object <- formula
  type <- match.arg(type)
  out <- if (identical(type, "submitted") && !is.null(object$dataset_original)) {
    as.data.frame(object$dataset_original, stringsAsFactors = FALSE)
  } else {
    as.data.frame(object$dataset, stringsAsFactors = FALSE)
  }
  if (identical(type, "observed") && "response" %in% names(out)) {
    out <- out[is.finite(out$response), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}
