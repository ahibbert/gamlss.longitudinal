#' Attach smooth covariance output and build vcov return object
#'
#' @noRd
.gl_vcov_build_result <- function(object, eta_inv, response, vcov_final, se_final,
                                  method_used, method_requested,
                                  hessian_diagnostics) {
  smooth_vcov_out <- .gl_compute_smooth_vcov(object, eta_inv, response)
  smooth_vcov_list <- smooth_vcov_out$smooth_vcov
  smooth_se_list <- smooth_vcov_out$smooth_se

  list(
    vcov = list(
      overall = vcov_final,
      smooth_vcov = smooth_vcov_list,
      smooth_se = smooth_se_list
    ),
    se = list(
      overall = se_final,
      smooth_se = smooth_se_list
    ),
    method = method_used,
    method_requested = method_requested,
    hessian_diagnostics = hessian_diagnostics %||% NULL
  )
}
