.copula_frank_deriv <- function(u1, u2, par, deriv, log = FALSE) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_frank_par(par))

  u1 <- vals[[1]]

  u2 <- vals[[2]]

  theta <- vals[[3]]


  indep <- abs(theta) <= 1e-8

  if (any(indep)) {
    out <- numeric(length(theta))

    out[indep] <- switch(deriv,
      u1 = rep(0, sum(indep)),
      u2 = rep(0, sum(indep)),
      par = .copula_central_par_deriv(

        .copula_frank_pdf,
        u1[indep],
        u2[indep],
        par0 = 0,
        h = 1e-5,
        log = log
      ),
      stop("Unsupported Frank derivative: ", deriv, call. = FALSE)
    )

    if (all(indep)) {
      out[!is.finite(out)] <- 0

      return(out)
    }
  } else {
    out <- numeric(length(theta))
  }


  dep <- !indep

  theta_dep <- theta[dep]

  u1_dep <- u1[dep]

  u2_dep <- u2[dep]

  et <- exp(-theta_dep)

  eu <- exp(-theta_dep * u1_dep)

  ev <- exp(-theta_dep * u2_dep)

  a <- et - 1

  den <- a + (eu - 1) * (ev - 1)

  density <- .copula_frank_pdf(u1_dep, u2_dep, theta_dep)


  score <- switch(deriv,
    u1 = -theta_dep + 2 * theta_dep * eu * (ev - 1) / den,
    u2 = -theta_dep + 2 * theta_dep * ev * (eu - 1) / den,
    par = {
      da <- -et

      dden <- da - u1_dep * eu * (ev - 1) - u2_dep * ev * (eu - 1)

      1 / theta_dep + da / a - (u1_dep + u2_dep) - 2 * dden / den
    },
    stop("Unsupported Frank derivative: ", deriv, call. = FALSE)
  )


  out[dep] <- if (isTRUE(log)) score else density * score

  out[!is.finite(out)] <- 0

  out
}
