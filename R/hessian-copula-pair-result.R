#' Pair-level copula derivative helpers for Hessian assembly
#'
#' These helpers keep `.calc_copula_hessian_contributions()` readable by
#' separating pair indexing, copula derivative evaluation, and link-scale
#' curvature calculations.
#'
#' @noRd
.copula_hessian_empty_result <- function(n_obs, margin_pars, has_zeta) {
  empty <- rep(0, n_obs)

  cop_d2l_margin <- vector("list", length(margin_pars))
  names(cop_d2l_margin) <- margin_pars
  for (pn1 in margin_pars) {
    cop_d2l_margin[[pn1]] <- vector("list", length(margin_pars))
    names(cop_d2l_margin[[pn1]]) <- margin_pars
    for (pn2 in margin_pars) cop_d2l_margin[[pn1]][[pn2]] <- empty
  }

  list(
    cop_d1l_margin = setNames(lapply(margin_pars, function(p) empty), margin_pars),
    cop_d2l_margin = cop_d2l_margin,
    cop_d2l_theta = empty,
    cop_d2l_zeta = if (has_zeta) empty else NULL,
    cop_d2l_thetazeta = if (has_zeta) empty else NULL,
    cop_d2l_margin_theta = setNames(lapply(margin_pars, function(p) empty), margin_pars),
    cop_d2l_margin_zeta = if (has_zeta) {
      setNames(lapply(margin_pars, function(p) empty), margin_pars)
    } else {
      NULL
    }
  )
}

#' @noRd
.copula_hessian_zero_invalid <- function(x, pair_ok) {
  x[!pair_ok] <- 0
  x[!is.finite(x)] <- 0
  x
}
