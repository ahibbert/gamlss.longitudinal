.copula_gumbel_theta <- function(par) {
  pmax(as.numeric(par), 1)
}

.copula_gumbel_parts <- function(u1, u2, par) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_gumbel_theta(par))

  u1 <- vals[[1]]

  u2 <- vals[[2]]

  theta <- vals[[3]]

  x <- -log(u1)

  y <- -log(u2)

  a <- x^theta + y^theta

  s <- a^(1 / theta)

  list(u1 = u1, u2 = u2, theta = theta, x = x, y = y, a = a, s = s, cdf = exp(-s))
}

.copula_gumbel_cdf <- function(u1, u2, par) {
  .copula_gumbel_parts(u1, u2, par)$cdf
}

.copula_gumbel_pdf <- function(u1, u2, par) {
  p <- .copula_gumbel_parts(u1, u2, par)

  out <- p$cdf * (p$x * p$y)^(p$theta - 1) *

    p$s^(1 - 2 * p$theta) * (p$s + p$theta - 1) / (p$u1 * p$u2)

  out[!is.finite(out)] <- 0

  out
}

.copula_gumbel_hfunc1 <- function(u1, u2, par) {
  p <- .copula_gumbel_parts(u1, u2, par)

  out <- p$cdf * p$s^(1 - p$theta) * p$x^(p$theta - 1) / p$u1

  .copula_clamp01(out)
}
