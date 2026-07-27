.copula_clayton_theta <- function(par) {
  pmax(as.numeric(par), 0)
}

.copula_clayton_s <- function(u1, u2, theta) {
  u1^(-theta) + u2^(-theta) - 1
}

.copula_clayton_cdf <- function(u1, u2, par) {
  u1 <- .copula_clamp01(u1)

  u2 <- .copula_clamp01(u2)

  theta <- .copula_clayton_theta(par)

  out <- rep(NA_real_, length.out = max(length(u1), length(u2), length(theta)))

  u1 <- rep(u1, length.out = length(out))

  u2 <- rep(u2, length.out = length(out))

  theta <- rep(theta, length.out = length(out))

  out <- u1 * u2

  dep <- theta > 1e-10

  if (any(dep)) {
    s <- .copula_clayton_s(u1[dep], u2[dep], theta[dep])

    out[dep] <- s^(-1 / theta[dep])
  }

  out
}

.copula_clayton_pdf <- function(u1, u2, par) {
  u1 <- .copula_clamp01(u1)

  u2 <- .copula_clamp01(u2)

  theta <- .copula_clayton_theta(par)

  out <- rep(1, length.out = max(length(u1), length(u2), length(theta)))

  u1 <- rep(u1, length.out = length(out))

  u2 <- rep(u2, length.out = length(out))

  theta <- rep(theta, length.out = length(out))

  dep <- theta > 1e-10

  if (any(dep)) {
    s <- .copula_clayton_s(u1[dep], u2[dep], theta[dep])

    out[dep] <- (theta[dep] + 1) *

      (u1[dep] * u2[dep])^(-theta[dep] - 1) *

      s^(-2 - 1 / theta[dep])
  }

  out
}

.copula_clayton_hfunc1 <- function(u1, u2, par) {
  u1 <- .copula_clamp01(u1)

  u2 <- .copula_clamp01(u2)

  theta <- .copula_clayton_theta(par)

  out <- rep(NA_real_, length.out = max(length(u1), length(u2), length(theta)))

  u1 <- rep(u1, length.out = length(out))

  u2 <- rep(u2, length.out = length(out))

  theta <- rep(theta, length.out = length(out))

  indep <- theta <= 1e-10

  out[indep] <- u2[indep]

  dep <- !indep

  if (any(dep)) {
    s <- .copula_clayton_s(u1[dep], u2[dep], theta[dep])

    out[dep] <- u1[dep]^(-theta[dep] - 1) * s^(-1 / theta[dep] - 1)
  }

  .copula_clamp01(out)
}
