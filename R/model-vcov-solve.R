#' Solve missing fixed-effect vcov output from a Hessian
#'
#' @noRd
.gl_vcov_solve_if_needed <- function(vcov_path, method, d2_mat, response) {
  vcov_final <- vcov_path$vcov_final
  se_final <- vcov_path$se_final
  hessian_diagnostics <- vcov_path$hessian_diagnostics

  if (method %in% c("numderiv", "analytical", "analytical_only")) {
    if (is.null(vcov_final) || is.null(se_final)) {
      vcov_solved <- .gl_solve_hessian_vcov(vcov_path$hessian_nd)
      vcov_final <- vcov_solved$vcov
      se_final <- vcov_solved$se
      hessian_diagnostics <- vcov_solved$hessian_diagnostics
    }
  } else {
    vcov_final <- -(solve((d2_mat))) / (length(response))
    se_final <- sqrt(abs(diag(vcov_final)))
  }

  list(
    vcov_final = vcov_final,
    se_final = se_final,
    hessian_diagnostics = hessian_diagnostics
  )
}
