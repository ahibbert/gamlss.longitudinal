.copula_t_df <- function(par2) {

  pmax(as.numeric(par2), 2.000001)

}


.copula_t_df_step <- function(df, rel = 1e-3) {

  pmin(pmax(rel * abs(df), 1e-4), 0.25 * (df - 2))

}


.copula_t_pdf <- function(u1, u2, par, par2) {

  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_gaussian_rho(par), .copula_t_df(par2))

  u1 <- vals[[1]]

  u2 <- vals[[2]]

  rho <- vals[[3]]

  df <- vals[[4]]

  x <- stats::qt(u1, df = df)

  y <- stats::qt(u2, df = df)

  q <- (x^2 - 2 * rho * x * y + y^2) / (1 - rho^2)

  log_biv <- lgamma((df + 2) / 2) - lgamma(df / 2) -

    log(df * pi) - 0.5 * log1p(-rho^2) -

    (df + 2) / 2 * log1p(q / df)

  out <- exp(log_biv - stats::dt(x, df = df, log = TRUE) - stats::dt(y, df = df, log = TRUE))

  out[!is.finite(out)] <- 0

  out

}


.copula_t_hfunc1 <- function(u1, u2, par, par2) {

  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_gaussian_rho(par), .copula_t_df(par2))

  u1 <- vals[[1]]

  u2 <- vals[[2]]

  rho <- vals[[3]]

  df <- vals[[4]]

  x <- stats::qt(u1, df = df)

  y <- stats::qt(u2, df = df)

  scale <- sqrt((df + x^2) * (1 - rho^2) / (df + 1))

  .copula_clamp01(stats::pt((y - rho * x) / scale, df = df + 1))

}


.copula_t_cdf_one <- function(u1, u2, rho, df) {

  x <- stats::qt(u1, df = df)

  y <- stats::qt(u2, df = df)

  integrand <- function(z) {

    scale <- sqrt((df + z^2) * (1 - rho^2) / (df + 1))

    stats::dt(z, df = df) * stats::pt((y - rho * z) / scale, df = df + 1)

  }

  stats::integrate(integrand, lower = -Inf, upper = x, rel.tol = 1e-7)$value

}


.copula_t_cdf <- function(u1, u2, par, par2) {

  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_gaussian_rho(par), .copula_t_df(par2))

  out <- mapply(

    .copula_t_cdf_one,

    vals[[1]],

    vals[[2]],

    vals[[3]],

    vals[[4]],

    SIMPLIFY = TRUE,

    USE.NAMES = FALSE

  )

  .copula_clamp01(out)

}


.copula_t_deriv <- function(u1, u2, par, par2, deriv, log = FALSE) {

  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_gaussian_rho(par), .copula_t_df(par2))

  u1 <- vals[[1]]

  u2 <- vals[[2]]

  rho <- vals[[3]]

  df <- vals[[4]]

  x <- stats::qt(u1, df = df)

  y <- stats::qt(u2, df = df)

  rho_denom <- 1 - rho^2

  numerator <- x^2 - 2 * rho * x * y + y^2

  denom <- df * rho_denom + numerator

  q <- numerator / rho_denom

  log_biv <- lgamma((df + 2) / 2) - lgamma(df / 2) -

    log(df * pi) - 0.5 * log1p(-rho^2) -

    (df + 2) / 2 * log1p(q / df)

  log_density <- log_biv - stats::dt(x, df = df, log = TRUE) - stats::dt(y, df = df, log = TRUE)

  density <- exp(log_density)

  density[!is.finite(density)] <- 0


  score <- switch(

    deriv,

    u1 = {

      dlog_dx <- -(df + 2) * (x - rho * y) / denom +

        (df + 1) * x / (df + x^2)

      dlog_dx / stats::dt(x, df = df)

    },

    u2 = {

      dlog_dy <- -(df + 2) * (y - rho * x) / denom +

        (df + 1) * y / (df + y^2)

      dlog_dy / stats::dt(y, df = df)

    },

    par = {

      dq_drho <- (-2 * x * y * rho_denom + 2 * rho * numerator) / rho_denom^2

      rho / rho_denom - (df + 2) * dq_drho / (2 * (df + q))

    },

    par2 = {

      h <- .copula_t_df_step(df)

      dcd_df <- (

        .copula_t_pdf(u1, u2, rho, df + h) -

          .copula_t_pdf(u1, u2, rho, df - h)

      ) / (2 * h)

      dcd_df / density

    },

    stop("Unsupported t derivative: ", deriv, call. = FALSE)

  )


  out <- if (isTRUE(log)) score else density * score

  out[!is.finite(out)] <- 0

  out

}


.copula_t_deriv2 <- function(u1, u2, par, par2, deriv) {

  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_gaussian_rho(par), .copula_t_df(par2))

  u1 <- vals[[1]]

  u2 <- vals[[2]]

  rho <- vals[[3]]

  df <- vals[[4]]


  out <- switch(

    deriv,

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


