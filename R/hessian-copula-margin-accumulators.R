#' Initialize margin-side copula Hessian accumulators
#'
#' @noRd
.copula_hessian_margin_accumulators <- function(margin_pars, n_obs, n_pairs, has_zeta) {
  cop_d2l_margin <- vector("list", length(margin_pars))
  names(cop_d2l_margin) <- margin_pars
  for (pn in margin_pars) {
    cop_d2l_margin[[pn]] <- vector("list", length(margin_pars))
    names(cop_d2l_margin[[pn]]) <- margin_pars
    for (pn2 in margin_pars) {
      cop_d2l_margin[[pn]][[pn2]] <- numeric(n_obs)
    }
  }

  cop_d1l_margin <- setNames(lapply(margin_pars, function(p) numeric(n_obs)), margin_pars)
  cop_d2l_margin_theta_u1 <- setNames(lapply(margin_pars, function(p) numeric(n_pairs)), margin_pars)
  cop_d2l_margin_theta_u2 <- setNames(lapply(margin_pars, function(p) numeric(n_pairs)), margin_pars)
  cop_d2l_margin_zeta_u1 <- if (has_zeta) {
    setNames(lapply(margin_pars, function(p) numeric(n_pairs)), margin_pars)
  } else {
    NULL
  }
  cop_d2l_margin_zeta_u2 <- if (has_zeta) {
    setNames(lapply(margin_pars, function(p) numeric(n_pairs)), margin_pars)
  } else {
    NULL
  }

  list(
    cop_d1l_margin = cop_d1l_margin,
    cop_d2l_margin = cop_d2l_margin,
    cop_d2l_margin_theta_u1 = cop_d2l_margin_theta_u1,
    cop_d2l_margin_theta_u2 = cop_d2l_margin_theta_u2,
    cop_d2l_margin_zeta_u1 = cop_d2l_margin_zeta_u1,
    cop_d2l_margin_zeta_u2 = cop_d2l_margin_zeta_u2
  )
}
