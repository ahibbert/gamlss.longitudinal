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


  out <- switch(

    deriv,

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


.copula_gaussian_deriv2 <- function(u1, u2, par, deriv) {

  u1 <- .copula_clamp01(u1)

  u2 <- .copula_clamp01(u2)

  rho <- .copula_gaussian_rho(par)

  n <- max(length(u1), length(u2), length(rho))

  u1 <- rep(u1, length.out = n)

  u2 <- rep(u2, length.out = n)

  rho <- rep(rho, length.out = n)


  z1 <- stats::qnorm(u1)

  z2 <- stats::qnorm(u2)

  phi1 <- stats::dnorm(z1)

  phi2 <- stats::dnorm(z2)

  denom <- 1 - rho^2

  density <- .copula_gaussian_pdf(u1, u2, rho)


  out <- switch(

    deriv,

    u1 = {

      dlog_du1 <- ((rho * z2 - rho^2 * z1) / denom) / phi1

      d2log_du1 <- (((rho * z2 - rho^2 * z1) / denom) * z1 -

        rho^2 / denom) / phi1^2

      density * (dlog_du1^2 + d2log_du1)

    },

    u2 = {

      dlog_du2 <- ((rho * z1 - rho^2 * z2) / denom) / phi2

      d2log_du2 <- (((rho * z1 - rho^2 * z2) / denom) * z2 -

        rho^2 / denom) / phi2^2

      density * (dlog_du2^2 + d2log_du2)

    },

    par = {

      sum_z2 <- z1^2 + z2^2

      numerator <- 2 * rho * z1 * z2 - rho^2 * sum_z2

      numerator_dr <- 2 * z1 * z2 - 2 * rho * sum_z2

      numerator_d2r <- -2 * sum_z2

      q <- numerator_dr * denom + 2 * rho * numerator

      q_dr <- numerator_d2r * denom + 2 * numerator

      dlog_drho <- rho / denom + q / (2 * denom^2)

      d2log_drho <- (1 + rho^2) / denom^2 +

        q_dr / (2 * denom^2) +

        2 * rho * q / denom^3

      density * (dlog_drho^2 + d2log_drho)

    },

    stop("Unsupported Gaussian second derivative: ", deriv, call. = FALSE)

  )


  out[!is.finite(out)] <- 0

  out

}


