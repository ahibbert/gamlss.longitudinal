.copula_t_df <- function(par2) {
  pmax(as.numeric(par2), 2.000001)
}

.copula_t_df_step <- function(df, rel = 1e-3) {
  pmin(pmax(rel * abs(df), 1e-4), 0.25 * (df - 2))
}

.copula_t_components <- function(u1, u2, par, par2) {
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
    density = density
  )
}

.copula_t_pdf <- function(u1, u2, par, par2) {
  .copula_t_components(u1, u2, par, par2)$density
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
