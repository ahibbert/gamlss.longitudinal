#' Summarise a fitted longitudinal GAMLSS-copula model
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param include_vcov Logical; include coefficient uncertainty summaries when
#'   possible.
#' @param numderiv Logical; use the numerical Hessian path when computing
#'   variance-covariance output.
#' @param ci_level Confidence level for coefficient intervals.
#' @param vcov_method Optional variance-covariance method passed as `method` to
#'   `vcov.gamlss.longitudinal()`.
#' @param ... Additional arguments passed to `vcov.gamlss.longitudinal()`.
#'
#' @return A `summary.gamlss.longitudinal` object with model metadata,
#'   coefficient summaries, smooth-term summaries, and optional vcov output.
#'   For a returned nonconverged fit, `include_vcov = FALSE` still provides
#'   point estimates and convergence diagnostics, while AIC and BIC are marked
#'   provisional and must not be used for model selection. Requesting vcov or
#'   other inference remains an error under the central convergence gate.
#'   For a fit created with `missingness = "segment"`, the summary records the
#'   number of gaps and segments and labels AIC and BIC as available under the
#'   explicit assumption that different observed segments are independent.
#'   `fit$inference_contract`, `fit$smooth_inference_contract`, and the object
#'   `inference_contract` attribute record methods, coefficient blocks,
#'   conditioning, omitted uncertainty, diagnostics, fallback, and validity.
#' @export

summary.gamlss.longitudinal <- function(
    object,
    include_vcov = TRUE,
    numderiv = FALSE,
    ci_level = 0.95,
    vcov_method = NULL,
    ...) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be of class 'gamlss.longitudinal'.")
  }

  model_info <- .gl_summary_model_info(object)
  fit_criteria <- .gl_summary_fit_criteria(object, n_obs = model_info$n_obs)

  vcov_out <- NULL
  if (isTRUE(include_vcov)) {
    vcov_args <- list(...)
    if (!is.null(vcov_method) && is.null(vcov_args$method)) {
      vcov_args$method <- vcov_method
    }
    vcov_out <- .resolve_vcov(
      object = object,
      numderiv = numderiv,
      extra_args = vcov_args
    )
  }

  coef_tbl <- .gl_summary_coefficient_table(object, vcov_out = vcov_out)
  smooth_terms <- .gl_summary_smooth_terms(object)

  .gl_build_summary_object(
    object = object,
    model_info = model_info,
    fit_criteria = fit_criteria,
    coef_tbl = coef_tbl,
    smooth_terms = smooth_terms,
    vcov_out = vcov_out,
    include_vcov = include_vcov,
    numderiv = numderiv,
    ci_level = ci_level
  )
}
