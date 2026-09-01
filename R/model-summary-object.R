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
      inference_contract = vcov_out$inference_contract %||%
        if (isTRUE(include_vcov) && !is.null(vcov_out)) {
          .gl_fixed_inference_contract(vcov_out, coefficient_names = names(object$par))
        } else NULL,
      smooth_inference_contract = vcov_out$smooth_inference_contract %||%
        .gl_inference_contract(
          "smooth_penalized_conditional",
          coefficient_names = .gl_smooth_coefficient_names(object),
          validity_status = if (length(.gl_smooth_coefficient_names(object))) "approximate" else "not_applicable"
        ),
      model_selection = fit_criteria$model_selection
    ),
    smooth_terms = smooth_terms,
    coefficients = coef_tbl,
    vcov = vcov_out
  )
  class(out) <- "summary.gamlss.longitudinal"
  attr(out, "inference_contract") <- out$fit$inference_contract
  out
}
