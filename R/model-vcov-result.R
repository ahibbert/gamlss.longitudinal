#' Attach smooth covariance output and build vcov return object
#'
#' @noRd
.gl_vcov_build_result <- function(object, eta_inv, response, vcov_final, se_final,
                                  method_used, method_requested,
                                  hessian_diagnostics) {
  smooth_vcov_out <- .gl_compute_smooth_vcov(object, eta_inv, response)
  smooth_vcov_list <- smooth_vcov_out$smooth_vcov
  smooth_se_list <- smooth_vcov_out$smooth_se

  out <- list(
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
  .gl_enrich_vcov_contract(out, object)
}

#' Format detailed covariance output for the standard `vcov()` interface
#'
#' @noRd
.gl_vcov_format_result <- function(details, return_details = FALSE) {
  if (isTRUE(return_details)) {
    return(details)
  }
  out <- details$vcov$overall
  attr(out, "gamlss_longitudinal_details") <- details
  attr(out, "method") <- details$method %||% NA_character_
  attr(out, "method_requested") <- details$method_requested %||% NA_character_
  attr(out, "hessian_diagnostics") <- details$hessian_diagnostics %||% NULL
  attr(out, "inference_contract") <- details$inference_contract %||%
    attr(details$vcov$overall, "inference_contract")
  class(out) <- unique(c("gamlss_longitudinal_vcov", class(out)))
  out
}

#' Backward-compatible access to detailed covariance components
#'
#' @param x A matrix returned by `vcov.gamlss.longitudinal()`.
#' @param name Detailed component name.
#' @return The requested detailed component.
#' @export
`$.gamlss_longitudinal_vcov` <- function(x, name) {
  details <- attr(x, "gamlss_longitudinal_details", exact = TRUE)
  details[[name]]
}
