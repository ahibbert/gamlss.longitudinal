.copula_t_df <- function(par2) {
  pmax(as.numeric(par2), 2.000001)
}

.copula_t_df_step <- function(df, rel = 1e-3) {
  pmin(pmax(rel * abs(df), 1e-4), 0.25 * (df - 2))
}

.copula_t_gaussian_limit_df <- 1e6

.copula_t_components <- function(u1, u2, par, par2) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_gaussian_rho(par), .copula_t_df(par2))

  u1 <- vals[[1]]
  u2 <- vals[[2]]
  rho <- vals[[3]]
  df <- vals[[4]]

  gaussian_limit <- df >= .copula_t_gaussian_limit_df
  calculation_df <- df
  calculation_df[gaussian_limit] <- .copula_t_gaussian_limit_df
  x <- stats::qt(u1, df = calculation_df)
  y <- stats::qt(u2, df = calculation_df)
  if (any(gaussian_limit)) {
    x[gaussian_limit] <- stats::qnorm(u1[gaussian_limit])
    y[gaussian_limit] <- stats::qnorm(u2[gaussian_limit])
  }
  rho_denom <- 1 - rho^2
  numerator <- x^2 - 2 * rho * x * y + y^2
  denom <- df * rho_denom + numerator
  q <- numerator / rho_denom

  log_biv <- lgamma((calculation_df + 2) / 2) - lgamma(calculation_df / 2) -
    log(calculation_df * pi) - 0.5 * log1p(-rho^2) -
    (calculation_df + 2) / 2 * log1p(q / calculation_df)
  log_density <- log_biv -
    stats::dt(x, df = calculation_df, log = TRUE) -
    stats::dt(y, df = calculation_df, log = TRUE)
  if (any(gaussian_limit)) {
    gaussian_denom <- 1 - rho[gaussian_limit]^2
    log_density[gaussian_limit] <- (
      2 * rho[gaussian_limit] * x[gaussian_limit] * y[gaussian_limit] -
        rho[gaussian_limit]^2 *
          (x[gaussian_limit]^2 + y[gaussian_limit]^2)
    ) / (2 * gaussian_denom) - 0.5 * log(gaussian_denom)
  }
  density <- exp(log_density)
  density[!is.finite(density)] <- 0

  list(
    u1 = u1,
    u2 = u2,
    rho = rho,
    df = df,
    x = x,
    y = y,
    rho_denom = rho_denom,
    numerator = numerator,
    denom = denom,
    q = q,
    density = density,
    log_density = log_density
  )
}

.copula_t_pdf <- function(u1, u2, par, par2) {
  .copula_t_components(u1, u2, par, par2)$density
}

.copula_t_logpdf <- function(u1, u2, par, par2) {
  .copula_t_components(u1, u2, par, par2)$log_density
}

.copula_t_hfunc1 <- function(u1, u2, par, par2) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_gaussian_rho(par), .copula_t_df(par2))

  u1 <- vals[[1]]

  u2 <- vals[[2]]

  rho <- vals[[3]]

  df <- vals[[4]]

  gaussian_limit <- df >= .copula_t_gaussian_limit_df
  out <- numeric(length(df))
  if (any(gaussian_limit)) {
    out[gaussian_limit] <- .copula_gaussian_hfunc1(
      u1[gaussian_limit], u2[gaussian_limit], rho[gaussian_limit]
    )
  }
  ordinary <- !gaussian_limit
  if (any(ordinary)) {
    x <- stats::qt(u1[ordinary], df = df[ordinary])
    y <- stats::qt(u2[ordinary], df = df[ordinary])
    scale <- sqrt(
      (df[ordinary] + x^2) * (1 - rho[ordinary]^2) /
        (df[ordinary] + 1)
    )
    out[ordinary] <- stats::pt(
      (y - rho[ordinary] * x) / scale,
      df = df[ordinary] + 1
    )
  }
  .copula_clamp01(out)
}

.copula_t_cdf_one <- function(u1, u2, rho, df) {
  if (!all(is.finite(c(u1, u2, rho, df))) || df <= 2 || abs(rho) >= 1) {
    return(NA_real_)
  }
  if (df >= .copula_t_gaussian_limit_df) {
    return(.copula_gaussian_cdf_cached(u1, u2, rho))
  }
  x <- stats::qt(u1, df = df)

  y <- stats::qt(u2, df = df)
  if (!is.finite(x) || !is.finite(y)) return(NA_real_)

  integrand <- function(z) {
    scale <- sqrt((df + z^2) * (1 - rho^2) / (df + 1))

    stats::dt(z, df = df) * stats::pt((y - rho * z) / scale, df = df + 1)
  }

  tryCatch(
    stats::integrate(integrand, lower = -Inf, upper = x, rel.tol = 1e-7)$value,
    error = function(e) NA_real_
  )
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
