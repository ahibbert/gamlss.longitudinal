.copula_gaussian_pdf <- function(u1, u2, par) {
  u1 <- .copula_clamp01(u1)

  u2 <- .copula_clamp01(u2)

  rho <- .copula_gaussian_rho(par)

  z1 <- stats::qnorm(u1)

  z2 <- stats::qnorm(u2)

  denom <- 1 - rho^2

  exp((2 * rho * z1 * z2 - rho^2 * (z1^2 + z2^2)) / (2 * denom)) / sqrt(denom)
}

.copula_gaussian_hfunc1 <- function(u1, u2, par) {
  u1 <- .copula_clamp01(u1)

  u2 <- .copula_clamp01(u2)

  rho <- .copula_gaussian_rho(par)

  stats::pnorm((stats::qnorm(u2) - rho * stats::qnorm(u1)) / sqrt(1 - rho^2))
}
