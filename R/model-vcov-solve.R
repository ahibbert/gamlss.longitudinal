#' Solve missing fixed-effect vcov output from a Hessian
#'
#' @noRd
.gl_vcov_solve_if_needed <- function(vcov_path, method, d2_mat, response,
                                     parameter_names = colnames(vcov_path$hessian_nd),
                                     inference = inference_control("standard"),
                                     gradient = NULL) {
  vcov_final <- vcov_path$vcov_final
  se_final <- vcov_path$se_final
  hessian_diagnostics <- vcov_path$hessian_diagnostics

  if (method %in% c("numderiv", "analytical", "analytical_only")) {
    if (is.null(vcov_final) || is.null(se_final)) {
      vcov_solved <- .gl_solve_hessian_vcov(
        vcov_path$hessian_nd,
        parameter_names = parameter_names,
        inference = inference,
        source = vcov_path$method_used,
        fallback = vcov_path$hessian_fallback %||% NULL,
        gradient = gradient,
        reference_hessian = vcov_path$reference_hessian %||% NULL
      )
      vcov_final <- vcov_solved$vcov
      se_final <- vcov_solved$se
      hessian_diagnostics <- vcov_solved$hessian_diagnostics
    }
  } else {
    vcov_solved <- .gl_solve_hessian_vcov(
      d2_mat * length(response),
      parameter_names = parameter_names,
      inference = inference,
      source = "legacy",
      gradient = gradient
    )
    vcov_final <- vcov_solved$vcov
    se_final <- sqrt(diag(vcov_final))
    hessian_diagnostics <- vcov_solved$hessian_diagnostics
  }

  list(
    vcov_final = vcov_final,
    se_final = se_final,
    hessian_diagnostics = hessian_diagnostics
  )
}
