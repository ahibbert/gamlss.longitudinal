.copula_frank_deriv2 <- function(u1, u2, par, deriv) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_frank_par(par))

  u1 <- vals[[1]]

  u2 <- vals[[2]]

  theta <- vals[[3]]


  indep <- abs(theta) <= 1e-6

  if (any(indep)) {
    out <- numeric(length(theta))

    out[indep] <- switch(deriv,
      u1 = rep(0, sum(indep)),
      u2 = rep(0, sum(indep)),
      par = .copula_central_par_deriv2(

        .copula_frank_pdf,
        u1[indep],
        u2[indep],
        par0 = 0,
        h = 1e-4
      ),
      stop("Unsupported Frank second derivative: ", deriv, call. = FALSE)
    )

    if (all(indep)) {
      out[!is.finite(out)] <- 0

      return(out)
    }
  } else {
    out <- numeric(length(theta))
  }


  dep <- !indep

  u1_dep <- u1[dep]

  u2_dep <- u2[dep]

  theta_dep <- theta[dep]


  out[dep] <- switch(deriv,
    u1 = {
      h <- pmin(1e-5, 0.25 * u1_dep, 0.25 * (1 - u1_dep))

      (

        .copula_frank_deriv(u1_dep + h, u2_dep, theta_dep, deriv = "u1") -

          .copula_frank_deriv(u1_dep - h, u2_dep, theta_dep, deriv = "u1")

      ) / (2 * h)
    },
    u2 = {
      h <- pmin(1e-5, 0.25 * u2_dep, 0.25 * (1 - u2_dep))

      (

        .copula_frank_deriv(u1_dep, u2_dep + h, theta_dep, deriv = "u2") -

          .copula_frank_deriv(u1_dep, u2_dep - h, theta_dep, deriv = "u2")

      ) / (2 * h)
    },
    par = {
      h <- pmin(1e-4, 0.25 * abs(theta_dep))

      (

        .copula_frank_deriv(u1_dep, u2_dep, theta_dep + h, deriv = "par") -

          .copula_frank_deriv(u1_dep, u2_dep, theta_dep - h, deriv = "par")

      ) / (2 * h)
    },
    stop("Unsupported Frank second derivative: ", deriv, call. = FALSE)
  )


  out[!is.finite(out)] <- 0

  out
}
