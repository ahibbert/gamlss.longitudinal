#' @noRd
.copula_hessian_parameter_curvature <- function(derivs, par1_e, par2_e,
                                                pair_ok, has_zeta, copula_dist = "N") {
  # Dependence parameters use different links. Gaussian/t use Fisher-z;
  # Clayton uses log, Frank identity, and Gumbel/Joe log(theta-1).
  dr_pair <- switch(copula_dist,
    N = 1 - par1_e^2, t = 1 - par1_e^2,
    C = par1_e, F = rep(1, length(par1_e)),
    G = ifelse(par1_e >= 17, 0, par1_e - 1), J = par1_e - 1)
  d2_linkinv_par <- switch(copula_dist,
    N = -2 * par1_e * (1 - par1_e^2), t = -2 * par1_e * (1 - par1_e^2),
    C = par1_e, F = rep(0, length(par1_e)),
    G = ifelse(par1_e >= 17, 0, par1_e - 1), J = par1_e - 1)
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
