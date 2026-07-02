.copula_clayton_deriv <- function(u1, u2, par, deriv, log = FALSE) {
  u1 <- .copula_clamp01(u1)

  u2 <- .copula_clamp01(u2)

  theta <- .copula_clayton_theta(par)

  n <- max(length(u1), length(u2), length(theta))

  u1 <- rep(u1, length.out = n)

  u2 <- rep(u2, length.out = n)

  theta <- rep(theta, length.out = n)


  out <- numeric(n)

  indep <- theta <= 1e-8

  if (any(indep)) {
    out[indep] <- switch(deriv,
      u1 = .copula_indep_deriv(u1[indep], deriv, log),
      u2 = .copula_indep_deriv(u1[indep], deriv, log),
      par = .copula_one_sided_par_deriv(

        .copula_clayton_pdf,
        u1[indep],
        u2[indep],
        par0 = 0,
        h = 1e-5,
        log = log
      ),
      stop("Unsupported Clayton derivative: ", deriv, call. = FALSE)
    )
  }

  if (all(indep)) {
    out[!is.finite(out)] <- 0

    return(out)
  }


  dep <- !indep

  u1_dep <- u1[dep]

  u2_dep <- u2[dep]

  theta_dep <- theta[dep]


  s <- .copula_clayton_s(u1_dep, u2_dep, theta_dep)

  density <- .copula_clayton_pdf(u1_dep, u2_dep, theta_dep)


  out[dep] <- switch(deriv,
    u1 = {
      dlog_du1 <- -(theta_dep + 1) / u1_dep + (2 * theta_dep + 1) * u1_dep^(-theta_dep - 1) / s

      density * dlog_du1
    },
    u2 = {
      dlog_du2 <- -(theta_dep + 1) / u2_dep + (2 * theta_dep + 1) * u2_dep^(-theta_dep - 1) / s

      density * dlog_du2
    },
    par = {
      log_u1 <- log(u1_dep)

      log_u2 <- log(u2_dep)

      a <- -2 - 1 / theta_dep

      ds_dtheta <- -log_u1 * u1_dep^(-theta_dep) - log_u2 * u2_dep^(-theta_dep)

      dlog_dtheta <- 1 / (theta_dep + 1) -

        (log_u1 + log_u2) +

        log(s) / theta_dep^2 +

        a * ds_dtheta / s

      if (isTRUE(log)) dlog_dtheta else density * dlog_dtheta
    },
    stop("Unsupported Clayton derivative: ", deriv, call. = FALSE)
  )


  out[!is.finite(out)] <- 0

  out
}
