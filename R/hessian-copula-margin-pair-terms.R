#' Replace non-finite Hessian contribution values with zero
#'
#' @noRd
.copula_hessian_finite_or_zero <- function(x) {
  x[!is.finite(x)] <- 0
  x
}

#' Pair-level first derivative contribution of log copula to one margin parameter
#'
#' @noRd
.copula_hessian_margin_d1_pair <- function(derivs, k, dFi1, dFi2, cv) {
  list(
    i1 = .copula_hessian_finite_or_zero(derivs$dcdu1[k] * dFi1 / cv),
    i2 = .copula_hessian_finite_or_zero(derivs$dcdu2[k] * dFi2 / cv)
  )
}

#' Pair-level diagonal second derivative contribution for one margin parameter
#'
#' @noRd
.copula_hessian_margin_diag_pair <- function(derivs, k, dFi1, dFi2, d2Fi1, d2Fi2, cv) {
  list(
    i1 = .copula_hessian_finite_or_zero(
      (derivs$d2cdu1_2[k] * dFi1^2 + derivs$dcdu1[k] * d2Fi1) / cv -
        (derivs$dcdu1[k] * dFi1 / cv)^2
    ),
    i2 = .copula_hessian_finite_or_zero(
      (derivs$d2cdu2_2[k] * dFi2^2 + derivs$dcdu2[k] * d2Fi2) / cv -
        (derivs$dcdu2[k] * dFi2 / cv)^2
    )
  )
}

#' Pair-level cross derivative contribution between two margin parameters
#'
#' @noRd
.copula_hessian_margin_cross_pair <- function(
    derivs, k, dFi1, dFi2, dFj1, dFj2, d2Fij_cross1, d2Fij_cross2, cv) {
  list(
    i1 = .copula_hessian_finite_or_zero(
      (derivs$d2cdu1_2[k] * dFi1 * dFj1 + derivs$dcdu1[k] * d2Fij_cross1) / cv -
        (derivs$dcdu1[k] * dFi1 / cv) * (derivs$dcdu1[k] * dFj1 / cv)
    ),
    i2 = .copula_hessian_finite_or_zero(
      (derivs$d2cdu2_2[k] * dFi2 * dFj2 + derivs$dcdu2[k] * d2Fij_cross2) / cv -
        (derivs$dcdu2[k] * dFi2 / cv) * (derivs$dcdu2[k] * dFj2 / cv)
    )
  )
}

#' Pair-level derivative contribution between a margin parameter and theta
#'
#' @noRd
.copula_hessian_margin_theta_pair <- function(derivs, k, dFi1, dFi2, cv) {
  list(
    i1 = .copula_hessian_finite_or_zero(
      derivs$d2cdthu1[k] * dFi1 / cv -
        (derivs$dcdu1[k] * dFi1 / cv) * (derivs$dcdth[k] / cv)
    ),
    i2 = .copula_hessian_finite_or_zero(
      derivs$d2cdthu2[k] * dFi2 / cv -
        (derivs$dcdu2[k] * dFi2 / cv) * (derivs$dcdth[k] / cv)
    )
  )
}

#' Pair-level derivative contribution between a margin parameter and zeta
#'
#' @noRd
.copula_hessian_margin_zeta_pair <- function(derivs, k, dFi1, dFi2, cv) {
  list(
    i1 = .copula_hessian_finite_or_zero(
      derivs$zeta$d2cdzu1[k] * dFi1 / cv -
        (derivs$dcdu1[k] * dFi1 / cv) * (derivs$zeta$dcdz[k] / cv)
    ),
    i2 = .copula_hessian_finite_or_zero(
      derivs$zeta$d2cdzu2[k] * dFi2 / cv -
        (derivs$dcdu2[k] * dFi2 / cv) * (derivs$zeta$dcdz[k] / cv)
    )
  )
}
