#' Assemble a summary.gamlss.longitudinal object
#'
#' @noRd
.gl_build_summary_object <- function(
    object,
    model_info,
    fit_criteria,
    coef_tbl,
    smooth_terms,
    vcov_out,
    include_vcov,
    numderiv,
    ci_level) {
  out <- list(
    model = model_info,
    convergence = object$convergence,
    fit = list(
      logLik = fit_criteria$loglik_joint,
      AIC = as.numeric(fit_criteria$aic_vec["joint"]),
      BIC = as.numeric(fit_criteria$bic_vec["joint"]),
      ci_level = ci_level,
      vcov_included = isTRUE(include_vcov),
      vcov_numderiv = isTRUE(numderiv),
      vcov_method = vcov_out$method %||% object$vcov_meta$method_used %||% object$vcov_meta$method %||% NA_character_,
      vcov_method_requested = vcov_out$method_requested %||% object$vcov_meta$method %||% NA_character_,
      hessian_diagnostics = vcov_out$hessian_diagnostics %||% NULL,
      criteria_status = if (identical(object$convergence$converged, FALSE)) {
        "provisional_nonconverged"
      } else {
        "available"
      },
      stop_reason = object$convergence$stop_reason %||% NA_character_,
      model_selection = fit_criteria$model_selection
    ),
    smooth_terms = smooth_terms,
    coefficients = coef_tbl,
    vcov = vcov_out
  )
  class(out) <- "summary.gamlss.longitudinal"
  out
}
