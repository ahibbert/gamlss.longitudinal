.select_copula_fit_family <- function(u1, u2, family, t_df_grid, copula_time = NULL) {
  if (!is.null(copula_time)) {
    return(.select_copula_fit_family_by_time(u1, u2, family, t_df_grid, copula_time))
  }
  tau_start <- suppressWarnings(stats::cor(u1, u2, method = "kendall", use = "complete.obs"))
  if (!is.finite(tau_start)) tau_start <- 0

  fit <- switch(family,
    N = .select_copula_fit_one_par(u1, u2, family, lower = -0.95, upper = 0.95, start = .copula_tau_to_par("N", tau_start)),
    C = .select_copula_fit_one_par(u1, u2, family, lower = 1e-8, upper = 50, start = .copula_tau_to_par("C", pmax(tau_start, 0))),
    F = .select_copula_fit_frank(u1, u2, tau_start),
    G = .select_copula_fit_one_par(u1, u2, family, lower = 1 + 1e-8, upper = 50, start = .copula_tau_to_par("G", pmax(tau_start, 0))),
    J = .select_copula_fit_one_par(u1, u2, family, lower = 1 + 1e-8, upper = 50, start = .copula_tau_to_par("J", pmax(tau_start, 0))),
    t = .select_copula_fit_t(u1, u2, tau_start, t_df_grid),
    stop("Unsupported copula family: ", family, call. = FALSE)
  )

  k <- if (identical(family, "t")) 2 else 1
  data.frame(
    family = family,
    par = fit$par,
    par2 = fit$par2,
    tau = .copula_par_to_tau(family, fit$par, fit$par2),
    logLik = fit$logLik,
    AIC = -2 * fit$logLik + 2 * k,
    BIC = -2 * fit$logLik + log(length(u1)) * k,
    stringsAsFactors = FALSE
  )
}

.select_copula_fit_family_by_time <- function(u1, u2, family, t_df_grid, copula_time) {
  copula_time <- factor(copula_time, levels = unique(copula_time))
  if (length(copula_time) != length(u1)) {
    stop("'copula_time' must have one value per pseudo-observation pair.", call. = FALSE)
  }
  levels_time <- levels(copula_time)
  fits <- lapply(levels_time, function(level) {
    idx <- copula_time == level
    .select_copula_fit_family(u1[idx], u2[idx], family = family, t_df_grid = t_df_grid)
  })
  log_lik <- sum(vapply(fits, function(fit) fit$logLik[[1L]], numeric(1)))
  k_per_level <- if (identical(family, "t")) 2 else 1
  k_total <- length(levels_time) * k_per_level
  tau <- vapply(fits, function(fit) fit$tau[[1L]], numeric(1))
  data.frame(
    family = family,
    par = NA_real_,
    par2 = if (identical(family, "t")) NA_real_ else 0,
    tau = stats::weighted.mean(tau, w = as.numeric(table(copula_time)), na.rm = TRUE),
    logLik = log_lik,
    AIC = -2 * log_lik + 2 * k_total,
    BIC = -2 * log_lik + log(length(u1)) * k_total,
    n_copula_time_levels = length(levels_time),
    stringsAsFactors = FALSE
  )
}

.select_copula_fit_one_par <- function(u1, u2, family, lower, upper, start = NULL, par2 = 0) {
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  if (!is.null(start) && is.finite(start)) {
    start <- pmin(pmax(as.numeric(start), lower), upper)
  }
  obj <- function(par) -.select_copula_loglik(u1, u2, family = family, par = par, par2 = par2)
  opt <- stats::optimize(obj, interval = c(lower, upper))
  candidates <- c(opt$minimum, lower, upper, start)
  candidates <- unique(candidates[is.finite(candidates)])
  ll <- vapply(candidates, function(par) -obj(par), numeric(1))
  best <- which.max(ll)
  list(par = candidates[best], par2 = par2, logLik = ll[best])
}

.select_copula_fit_frank <- function(u1, u2, tau_start) {
  if (abs(tau_start) < 1e-4) {
    tau_start <- 0.05
  }
  start <- .copula_tau_to_par("F", tau_start)
  .select_copula_fit_one_par(u1, u2, "F", lower = -50, upper = 50, start = start)
}

.select_copula_fit_t <- function(u1, u2, tau_start, t_df_grid) {
  t_df_grid <- unique(as.numeric(t_df_grid))
  t_df_grid <- t_df_grid[is.finite(t_df_grid) & t_df_grid > 2]
  if (length(t_df_grid) < 1L) {
    stop("t_df_grid must contain at least one finite value greater than 2.", call. = FALSE)
  }
  fits <- lapply(t_df_grid, function(df) {
    .select_copula_fit_one_par(
      u1,
      u2,
      family = "t",
      lower = -0.95,
      upper = 0.95,
      start = .copula_tau_to_par("t", tau_start),
      par2 = df
    )
  })
  ll <- vapply(fits, `[[`, numeric(1), "logLik")
  fits[[which.max(ll)]]
}

.select_copula_loglik <- function(u1, u2, family, par, par2 = 0) {
  dens <- .copula_pdf(u1, u2, family = family, par = par, par2 = par2)
  dens <- pmax(dens, .Machine$double.xmin)
  sum(log(dens[is.finite(dens)]))
}
