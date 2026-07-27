.copula_t_deriv2 <- function(u1, u2, par, par2, deriv) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_gaussian_rho(par), .copula_t_df(par2))

  u1 <- vals[[1]]

  u2 <- vals[[2]]

  rho <- vals[[3]]

  df <- vals[[4]]

  out <- switch(deriv,
    u1 = {
      h <- pmin(1e-5, 0.25 * u1, 0.25 * (1 - u1))

      (

        .copula_t_deriv(u1 + h, u2, rho, df, deriv = "u1") -

          .copula_t_deriv(u1 - h, u2, rho, df, deriv = "u1")

      ) / (2 * h)
    },
    u2 = {
      h <- pmin(1e-5, 0.25 * u2, 0.25 * (1 - u2))

      (

        .copula_t_deriv(u1, u2 + h, rho, df, deriv = "u2") -

          .copula_t_deriv(u1, u2 - h, rho, df, deriv = "u2")

      ) / (2 * h)
    },
    par = {
      h <- pmin(1e-4, 0.25 * (1 - abs(rho)))

      (

        .copula_t_deriv(u1, u2, rho + h, df, deriv = "par") -

          .copula_t_deriv(u1, u2, rho - h, df, deriv = "par")

      ) / (2 * h)
    },
    par2 = {
      h <- .copula_t_df_step(df)

      (

        .copula_t_pdf(u1, u2, rho, df + h) -

          2 * .copula_t_pdf(u1, u2, rho, df) +

          .copula_t_pdf(u1, u2, rho, df - h)

      ) / (h^2)
    },
    par1par2 = {
      h <- .copula_t_df_step(df)

      (

        .copula_t_deriv(u1, u2, rho, df + h, deriv = "par") -

          .copula_t_deriv(u1, u2, rho, df - h, deriv = "par")

      ) / (2 * h)
    },
    stop("Unsupported t second derivative: ", deriv, call. = FALSE)
  )

  out[!is.finite(out)] <- 0

  out
}
