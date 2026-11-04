.select_copula_fit_family <- function(
    u1,
    u2,
    family,
    t_df_grid,
    copula_time = NULL,
    upper1 = NULL,
    upper2 = NULL,
    lower1 = NULL,
    lower2 = NULL,
    log_mass1 = NULL,
    log_mass2 = NULL,
    log_survival1 = NULL,
    log_survival2 = NULL,
    bic_n = length(u1)) {
  if (!is.null(copula_time)) {
    return(.select_copula_fit_family_by_time(
      u1, u2, family, t_df_grid, copula_time,
      upper1 = upper1, upper2 = upper2, lower1 = lower1, lower2 = lower2,
      log_mass1 = log_mass1, log_mass2 = log_mass2,
      log_survival1 = log_survival1, log_survival2 = log_survival2,
      bic_n = bic_n
    ))
  }
  tau_start <- suppressWarnings(stats::cor(u1, u2, method = "kendall", use = "complete.obs"))
  if (!is.finite(tau_start)) tau_start <- 0

  fit <- switch(family,
    N = .select_copula_fit_one_par(u1, u2, family, lower = -0.95, upper = 0.95, start = .copula_tau_to_par("N", tau_start), upper1 = upper1, upper2 = upper2, lower1 = lower1, lower2 = lower2, log_mass1 = log_mass1, log_mass2 = log_mass2, log_survival1 = log_survival1, log_survival2 = log_survival2),
    C = .select_copula_fit_one_par(u1, u2, family, lower = 1e-8, upper = 50, start = .copula_tau_to_par("C", pmax(tau_start, 0)), upper1 = upper1, upper2 = upper2, lower1 = lower1, lower2 = lower2, log_mass1 = log_mass1, log_mass2 = log_mass2, log_survival1 = log_survival1, log_survival2 = log_survival2),
    F = .select_copula_fit_frank(u1, u2, tau_start, upper1 = upper1, upper2 = upper2, lower1 = lower1, lower2 = lower2, log_mass1 = log_mass1, log_mass2 = log_mass2, log_survival1 = log_survival1, log_survival2 = log_survival2),
    G = .select_copula_fit_one_par(u1, u2, family, lower = 1 + 1e-8, upper = 50, start = .copula_tau_to_par("G", pmax(tau_start, 0)), upper1 = upper1, upper2 = upper2, lower1 = lower1, lower2 = lower2, log_mass1 = log_mass1, log_mass2 = log_mass2, log_survival1 = log_survival1, log_survival2 = log_survival2),
    J = .select_copula_fit_one_par(u1, u2, family, lower = 1 + 1e-8, upper = 50, start = .copula_tau_to_par("J", pmax(tau_start, 0)), upper1 = upper1, upper2 = upper2, lower1 = lower1, lower2 = lower2, log_mass1 = log_mass1, log_mass2 = log_mass2, log_survival1 = log_survival1, log_survival2 = log_survival2),
    t = .select_copula_fit_t(u1, u2, tau_start, t_df_grid, upper1 = upper1, upper2 = upper2, lower1 = lower1, lower2 = lower2, log_mass1 = log_mass1, log_mass2 = log_mass2, log_survival1 = log_survival1, log_survival2 = log_survival2),
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
    BIC = -2 * fit$logLik + log(bic_n) * k,
    stringsAsFactors = FALSE
  )
}

.select_copula_fit_family_by_time <- function(
    u1,
    u2,
    family,
    t_df_grid,
    copula_time,
    upper1 = NULL,
    upper2 = NULL,
    lower1 = NULL,
    lower2 = NULL,
    log_mass1 = NULL,
    log_mass2 = NULL,
    log_survival1 = NULL,
    log_survival2 = NULL,
    bic_n = length(u1)) {
  copula_time <- factor(copula_time, levels = unique(copula_time))
  if (length(copula_time) != length(u1)) {
    stop("'copula_time' must have one value per pseudo-observation pair.", call. = FALSE)
  }
  levels_time <- levels(copula_time)
  fits <- lapply(levels_time, function(level) {
    idx <- copula_time == level
    .select_copula_fit_family(
      u1[idx], u2[idx], family = family, t_df_grid = t_df_grid,
      upper1 = if (is.null(upper1)) NULL else upper1[idx],
      upper2 = if (is.null(upper2)) NULL else upper2[idx],
      lower1 = if (is.null(lower1)) NULL else lower1[idx],
      lower2 = if (is.null(lower2)) NULL else lower2[idx],
      log_mass1 = if (is.null(log_mass1)) NULL else log_mass1[idx],
      log_mass2 = if (is.null(log_mass2)) NULL else log_mass2[idx],
      log_survival1 = if (is.null(log_survival1)) NULL else log_survival1[idx],
      log_survival2 = if (is.null(log_survival2)) NULL else log_survival2[idx]
    )
  })
  log_lik <- sum(vapply(fits, function(fit) fit$logLik[[1L]], numeric(1)))
  k_per_level <- if (identical(family, "t")) 2 else 1
  k_total <- length(levels_time) * k_per_level
  tau <- vapply(fits, function(fit) fit$tau[[1L]], numeric(1))
  out <- data.frame(
    family = family,
    par = NA_real_,
    par2 = if (identical(family, "t")) NA_real_ else 0,
    tau = stats::weighted.mean(tau, w = as.numeric(table(copula_time)), na.rm = TRUE),
    logLik = log_lik,
    AIC = -2 * log_lik + 2 * k_total,
    BIC = -2 * log_lik + log(bic_n) * k_total,
    n_copula_time_levels = length(levels_time),
    stringsAsFactors = FALSE
  )
  attr(out, "copula_time_fits") <- data.frame(
    copula_time = levels_time,
    family = family,
    par = vapply(fits, function(fit) fit$par[[1L]], numeric(1)),
    par2 = vapply(fits, function(fit) fit$par2[[1L]], numeric(1)),
    tau = tau,
    logLik = vapply(fits, function(fit) fit$logLik[[1L]], numeric(1)),
    stringsAsFactors = FALSE
  )
  out
}

.select_copula_fit_one_par <- function(
    u1,
    u2,
    family,
    lower,
    upper,
    start = NULL,
    par2 = 0,
    upper1 = NULL,
    upper2 = NULL,
    lower1 = NULL,
    lower2 = NULL,
    log_mass1 = NULL,
    log_mass2 = NULL,
    log_survival1 = NULL,
    log_survival2 = NULL) {
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  if (!is.null(start) && is.finite(start)) {
    start <- pmin(pmax(as.numeric(start), lower), upper)
  }
  obj <- function(par) -.select_copula_loglik(
    u1, u2, family = family, par = par, par2 = par2,
    upper1 = upper1, upper2 = upper2, lower1 = lower1, lower2 = lower2,
    log_mass1 = log_mass1, log_mass2 = log_mass2,
    log_survival1 = log_survival1, log_survival2 = log_survival2
  )
  opt <- stats::optimize(obj, interval = c(lower, upper))
  candidates <- c(opt$minimum, lower, upper, start)
  candidates <- unique(candidates[is.finite(candidates)])
  ll <- vapply(candidates, function(par) -obj(par), numeric(1))
  best <- which.max(ll)
  list(par = candidates[best], par2 = par2, logLik = ll[best])
}

.select_copula_fit_frank <- function(
    u1, u2, tau_start, upper1 = NULL, upper2 = NULL, lower1 = NULL, lower2 = NULL,
    log_mass1 = NULL, log_mass2 = NULL, log_survival1 = NULL, log_survival2 = NULL) {
  if (abs(tau_start) < 1e-4) {
    tau_start <- 0.05
  }
  start <- .copula_tau_to_par("F", tau_start)
  .select_copula_fit_one_par(
    u1, u2, "F", lower = -50, upper = 50, start = start,
    upper1 = upper1, upper2 = upper2, lower1 = lower1, lower2 = lower2,
    log_mass1 = log_mass1, log_mass2 = log_mass2,
    log_survival1 = log_survival1, log_survival2 = log_survival2
  )
}

.select_copula_fit_t <- function(
    u1, u2, tau_start, t_df_grid,
    upper1 = NULL, upper2 = NULL, lower1 = NULL, lower2 = NULL,
    log_mass1 = NULL, log_mass2 = NULL, log_survival1 = NULL, log_survival2 = NULL) {
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
      par2 = df,
      upper1 = upper1,
      upper2 = upper2,
      lower1 = lower1,
      lower2 = lower2,
      log_mass1 = log_mass1,
      log_mass2 = log_mass2,
      log_survival1 = log_survival1,
      log_survival2 = log_survival2
    )
  })
  ll <- vapply(fits, `[[`, numeric(1), "logLik")
  fits[[which.max(ll)]]
}

.select_copula_loglik <- function(
    u1,
    u2,
    family,
    par,
    par2 = 0,
    upper1 = NULL,
    upper2 = NULL,
    lower1 = NULL,
    lower2 = NULL,
    log_mass1 = NULL,
    log_mass2 = NULL,
    log_survival1 = NULL,
    log_survival2 = NULL) {
  bounds <- list(upper1, upper2, lower1, lower2)
  if (all(vapply(bounds, is.null, logical(1)))) {
    dens <- .copula_pdf(u1, u2, family = family, par = par, par2 = par2)
    dens <- pmax(dens, .Machine$double.xmin)
    return(sum(log(dens[is.finite(dens)])))
  }
  if (any(vapply(bounds, is.null, logical(1)))) {
    stop("All four discrete CDF bounds must be supplied together.", call. = FALSE)
  }
  rect <- .copula_rectangle_prob(
    upper1, upper2, lower1, lower2,
    family = family, par = par, par2 = par2,
    floor_probability = FALSE
  )
  if (is.null(log_mass1)) log_mass1 <- log(pmax(upper1 - lower1, .Machine$double.xmin))
  if (is.null(log_mass2)) log_mass2 <- log(pmax(upper2 - lower2, .Machine$double.xmin))
  contribution <- rep(-Inf, length(rect))
  valid_rect <- is.finite(rect) & rect > 0 & is.finite(log_mass1) & is.finite(log_mass2)
  contribution[valid_rect] <- log(rect[valid_rect]) -
    log_mass1[valid_rect] - log_mass2[valid_rect]
  unstable <- !is.finite(contribution) |
    (upper1 - lower1) <= 1e-10 | (upper2 - lower2) <= 1e-10
  if (any(unstable)) {
    par_values <- rep_len(par, length(upper1))
    par2_values <- rep_len(par2, length(upper1))
    stable <- .gl_discrete_copula_log_ratio_quadrature(
      family = family,
      par = par_values[unstable],
      par2 = par2_values[unstable],
      upper = cbind(upper1[unstable], upper2[unstable]),
      lower = cbind(lower1[unstable], lower2[unstable]),
      margin_log_d1 = log_mass1[unstable],
      margin_log_d2 = log_mass2[unstable],
      log_survival1 = if (is.null(log_survival1)) NULL else log_survival1[unstable],
      log_survival2 = if (is.null(log_survival2)) NULL else log_survival2[unstable]
    )
    contribution[unstable] <- stable
  }
  if (!all(is.finite(contribution))) return(-Inf)
  sum(contribution)
}
