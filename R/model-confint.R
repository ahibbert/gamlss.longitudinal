#' Conditional Wald confidence intervals for fixed coefficients

#'

#' @param object A fitted `gamlss.longitudinal` object.

#' @param parm Optional coefficient names or numeric indices.

#' @param level Confidence level.

#' @param method Variance-covariance method passed to [vcov.gamlss.longitudinal()].

#' @param ... Additional arguments passed to [vcov.gamlss.longitudinal()].

#'

#' @return A matrix with lower and upper confidence limits and an
#'   `inference_contract` attribute. Fixed-smooth covariance and
#'   smoothing-parameter uncertainty are excluded.

#' @importFrom stats confint

#' @export

confint.gamlss.longitudinal <- function(
    object,
    parm = NULL,
    level = 0.95,
    method = "analytical",
    ...) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }

  estimates <- stats::coef(object)

  vc <- .resolve_vcov(object, extra_args = list(method = method, ...))
  .gl_require_available_inference(vc)

  se <- vc$se$overall

  if (is.null(names(se))) {
    names(se) <- names(estimates)
  }

  idx <- if (is.null(parm)) {
    seq_along(estimates)
  } else if (is.character(parm)) {
    match(parm, names(estimates))
  } else {
    parm
  }

  if (any(is.na(idx)) || any(idx < 1L) || any(idx > length(estimates))) {
    stop("'parm' contains unknown coefficient names or indices.", call. = FALSE)
  }

  z <- stats::qnorm((1 + level) / 2)

  est <- estimates[idx]

  se_use <- se[names(est)]

  out <- cbind(
    lower = as.numeric(est) - z * as.numeric(se_use),
    upper = as.numeric(est) + z * as.numeric(se_use)
  )

  colnames(out) <- c(
    paste0(round((1 - level) / 2 * 100, 1), " %"),
    paste0(round((1 + level) / 2 * 100, 1), " %")
  )

  rownames(out) <- names(est)

  contract <- .gl_inference_contract(
    "confint_fixed",
    coefficient_names = names(est),
    method = vc$method %||% method,
    validity_status = vc$hessian_diagnostics$status %||% "not_recorded"
  )
  contract$covariance_contract <- vc$inference_contract %||%
    .gl_fixed_inference_contract(vc, coefficient_names = names(estimates))
  .gl_attach_inference_contract(out, contract)
}
