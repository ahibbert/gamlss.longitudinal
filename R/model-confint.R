#' Confidence intervals for fixed coefficients

#'

#' @param object A fitted `gamlss.longitudinal` object.

#' @param parm Optional coefficient names or numeric indices.

#' @param level Confidence level.

#' @param method Variance-covariance method passed to [vcov.gamlss.longitudinal()].

#' @param ... Additional arguments passed to [vcov.gamlss.longitudinal()].

#'

#' @return A matrix with lower and upper confidence limits.

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

  out
}
