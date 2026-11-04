.copula_gaussian_deriv <- function(u1, u2, par, deriv, log = FALSE) {
  u1 <- .copula_clamp01(u1)

  u2 <- .copula_clamp01(u2)

  rho <- .copula_gaussian_rho(par)

  n <- max(length(u1), length(u2), length(rho))

  u1 <- rep(u1, length.out = n)

  u2 <- rep(u2, length.out = n)

  rho <- rep(rho, length.out = n)

  z1 <- stats::qnorm(u1)

  z2 <- stats::qnorm(u2)

  denom <- 1 - rho^2

  density <- .copula_gaussian_pdf(u1, u2, rho)

  out <- switch(deriv,
    u1 = {
      dlog_du1 <- ((rho * z2 - rho^2 * z1) / denom) / stats::dnorm(z1)

      density * dlog_du1
    },
    u2 = {
      dlog_du2 <- ((rho * z1 - rho^2 * z2) / denom) / stats::dnorm(z2)

      density * dlog_du2
    },
    par = {
      sum_z2 <- z1^2 + z2^2

      numerator <- 2 * rho * z1 * z2 - rho^2 * sum_z2

      numerator_dr <- 2 * z1 * z2 - 2 * rho * sum_z2

      dlog_drho <- rho / denom +

        (numerator_dr * denom + 2 * rho * numerator) / (2 * denom^2)

      if (isTRUE(log)) dlog_drho else density * dlog_drho
    },
    stop("Unsupported Gaussian derivative: ", deriv, call. = FALSE)
  )

  out[!is.finite(out)] <- 0

  out
}
