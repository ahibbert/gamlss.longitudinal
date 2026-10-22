#' Compatibility entry point for the formerly Gaussian-only covariance patch
#'
#' The comprehensive correction is now shared with all continuous copulas and
#' CG. The original implementation and its evidence remain in history.
#' @noRd
.gl_covariance_gaussian_pair_derivatives <- function(u1, u2, fam_num, par1_e,
                                                    par2_e, pair_ok, has_zeta) {
  .copula_hessian_pair_derivatives(u1, u2, fam_num, par1_e, par2_e, pair_ok, has_zeta)
}

#' Assemble covariance using the same corrected Hessian consumed by CG
#' @noRd
.gl_covariance_analytical_hessian <- function(object, progress = interactive(),
                                             h = 1e-4) {
  calc_analytical_hessian(object, progress = progress, h = h)
}
