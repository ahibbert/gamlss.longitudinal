.copula_joe_deriv2 <- function(u1, u2, par, deriv) {
  p <- .copula_joe_parts(u1, u2, par)


  indep <- p$theta <= 1 + 1e-6

  if (any(indep)) {
    out <- numeric(length(p$theta))

    out[indep] <- switch(deriv,
      u1 = rep(0, sum(indep)),
      u2 = rep(0, sum(indep)),
      par = .copula_one_sided_par_deriv2(

        .copula_joe_pdf,
        p$u1[indep],
        p$u2[indep],
        par0 = 1,
        h = 1e-4
      ),
      stop("Unsupported Joe second derivative: ", deriv, call. = FALSE)
    )

    if (all(indep)) {
      out[!is.finite(out)] <- 0

      return(out)
    }
  } else {
    out <- numeric(length(p$theta))
  }


  dep <- !indep

  p_dep <- lapply(p, function(x) x[dep])

  out[dep] <- switch(deriv,
    u1 = {
      h <- pmin(1e-5, 0.25 * p_dep$u1, 0.25 * (1 - p_dep$u1))

      (

        .copula_joe_deriv(p_dep$u1 + h, p_dep$u2, p_dep$theta, deriv = "u1") -

          .copula_joe_deriv(p_dep$u1 - h, p_dep$u2, p_dep$theta, deriv = "u1")

      ) / (2 * h)
    },
    u2 = {
      h <- pmin(1e-5, 0.25 * p_dep$u2, 0.25 * (1 - p_dep$u2))

      (

        .copula_joe_deriv(p_dep$u1, p_dep$u2 + h, p_dep$theta, deriv = "u2") -

          .copula_joe_deriv(p_dep$u1, p_dep$u2 - h, p_dep$theta, deriv = "u2")

      ) / (2 * h)
    },
    par = {
      h <- pmin(1e-4, 0.25 * (p_dep$theta - 1))

      (

        .copula_joe_deriv(p_dep$u1, p_dep$u2, p_dep$theta + h, deriv = "par") -

          .copula_joe_deriv(p_dep$u1, p_dep$u2, p_dep$theta - h, deriv = "par")

      ) / (2 * h)
    },
    stop("Unsupported Joe second derivative: ", deriv, call. = FALSE)
  )


  out[!is.finite(out)] <- 0

  out
}
