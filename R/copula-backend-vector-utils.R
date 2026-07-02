#' Clamp Gaussian copula correlation away from singular boundaries
#'
#' @noRd
.copula_gaussian_rho <- function(par) {
  pmin(pmax(as.numeric(par), -0.999999), 0.999999)
}

#' Clamp pseudo-observations to the open unit interval
#'
#' @noRd
.copula_clamp01 <- function(u) {
  pmin(pmax(as.numeric(u), 1e-12), 1 - 1e-12)
}

#' Recycle copula vector inputs to common length
#'
#' @noRd
.copula_recycle <- function(...) {
  args <- list(...)
  n <- max(vapply(args, length, integer(1)))
  lapply(args, rep, length.out = n)
}

#' Zero derivative helper for independence-limit copula paths
#'
#' @noRd
.copula_indep_deriv <- function(u1, deriv, log = FALSE) {
  rep(0, length(u1))
}
