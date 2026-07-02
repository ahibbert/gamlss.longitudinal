.copula_frank_par <- function(par) {
  as.numeric(par)
}


.copula_frank_cdf <- function(u1, u2, par) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_frank_par(par))

  u1 <- vals[[1]]

  u2 <- vals[[2]]

  theta <- vals[[3]]

  out <- u1 * u2

  dep <- abs(theta) > 1e-8

  if (any(dep)) {
    et <- exp(-theta[dep])

    eu <- exp(-theta[dep] * u1[dep])

    ev <- exp(-theta[dep] * u2[dep])

    out[dep] <- -log1p((eu - 1) * (ev - 1) / (et - 1)) / theta[dep]
  }

  .copula_clamp01(out)
}


.copula_frank_pdf <- function(u1, u2, par) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_frank_par(par))

  u1 <- vals[[1]]

  u2 <- vals[[2]]

  theta <- vals[[3]]

  out <- rep(1, length(theta))

  dep <- abs(theta) > 1e-8

  if (any(dep)) {
    et <- exp(-theta[dep])

    eu <- exp(-theta[dep] * u1[dep])

    ev <- exp(-theta[dep] * u2[dep])

    den <- et - 1 + (eu - 1) * (ev - 1)

    out[dep] <- -theta[dep] * (et - 1) * eu * ev / den^2
  }

  out[!is.finite(out)] <- 0

  out
}


.copula_frank_hfunc1 <- function(u1, u2, par) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_frank_par(par))

  u1 <- vals[[1]]

  u2 <- vals[[2]]

  theta <- vals[[3]]

  out <- u2

  dep <- abs(theta) > 1e-8

  if (any(dep)) {
    et <- exp(-theta[dep])

    eu <- exp(-theta[dep] * u1[dep])

    ev <- exp(-theta[dep] * u2[dep])

    out[dep] <- eu * (ev - 1) / (et - 1 + (eu - 1) * (ev - 1))
  }

  .copula_clamp01(out)
}
