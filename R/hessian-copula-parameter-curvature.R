#' @noRd
.copula_hessian_parameter_curvature <- function(derivs, par1_e, par2_e,
                                                pair_ok, has_zeta) {
  dr_pair <- 1 - par1_e^2
  d2_linkinv_par <- -2 * par1_e * (1 - par1_e^2)
  cop_d2l_theta_pair <- derivs$d2ldth2_pair * dr_pair^2 +
    derivs$dldth * d2_linkinv_par
  cop_d2l_theta_pair <- .copula_hessian_zero_invalid(cop_d2l_theta_pair, pair_ok)

  cop_d2l_zeta_pair <- NULL
  cop_d2l_thetazeta_pair <- NULL
  if (has_zeta) {
    zeta_dr_pair <- pmax(par2_e - 2, 0)
    zeta_d2_linkinv_pair <- zeta_dr_pair
    cop_d2l_zeta_pair <- derivs$zeta$d2ldz2_pair * zeta_dr_pair^2 +
      derivs$zeta$dldz_pair * zeta_d2_linkinv_pair
    cop_d2l_zeta_pair <- .copula_hessian_zero_invalid(cop_d2l_zeta_pair, pair_ok)
    cop_d2l_thetazeta_pair <- derivs$zeta$d2ldthdz_pair
  }

  list(
    cop_d2l_theta_pair = cop_d2l_theta_pair,
    cop_d2l_zeta_pair = cop_d2l_zeta_pair,
    cop_d2l_thetazeta_pair = cop_d2l_thetazeta_pair
  )
}
