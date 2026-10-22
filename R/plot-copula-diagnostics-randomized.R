.copula_v2_validate_randomized_controls <- function(randomized_residuals,
                                                     randomized_seed,
                                                     randomized_draws = 5L) {
  if (!is.logical(randomized_residuals) || length(randomized_residuals) != 1L ||
      is.na(randomized_residuals)) {
    stop("'randomized_residuals' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(randomized_seed) || length(randomized_seed) != 1L ||
      !is.finite(randomized_seed)) {
    stop("'randomized_seed' must be a single finite number.", call. = FALSE)
  }
  if (!is.numeric(randomized_draws) || length(randomized_draws) != 1L ||
      !is.finite(randomized_draws) || randomized_draws < 1 ||
      randomized_draws != as.integer(randomized_draws)) {
    stop("'randomized_draws' must be a positive integer.", call. = FALSE)
  }
  list(
    randomized_residuals = randomized_residuals,
    randomized_seed = as.integer(randomized_seed),
    randomized_draws = as.integer(randomized_draws)
  )
}

.copula_v2_latent_rectangle_draws <- function(pair_data, family_num,
                                               draws = 5L, seed = 1L,
                                               tolerance = 1e-10) {
  required <- c(
    "u1_lower", "u1_upper", "u2_lower", "u2_upper",
    "theta_pair", "zeta_pair"
  )
  if (!all(required %in% names(pair_data))) {
    stop("Discrete pair diagnostics do not contain fitted CDF rectangles.",
      call. = FALSE
    )
  }
  if (!nrow(pair_data)) return(pair_data)

  index <- rep(seq_len(nrow(pair_data)), each = draws)
  completed <- pair_data[index, , drop = FALSE]
  completed$observed_pair_id <- index
  completed$latent_draw <- rep(seq_len(draws), times = nrow(pair_data))

  a1 <- pmin(pmax(as.numeric(completed$u1_lower), 0), 1)
  b1 <- pmin(pmax(as.numeric(completed$u1_upper), 0), 1)
  a2 <- pmin(pmax(as.numeric(completed$u2_lower), 0), 1)
  b2 <- pmin(pmax(as.numeric(completed$u2_upper), 0), 1)
  par <- as.numeric(completed$theta_pair)
  par2 <- as.numeric(completed$zeta_pair)
  par2[!is.finite(par2)] <- 0

  bounds_ok <- is.finite(a1) & is.finite(b1) & is.finite(a2) &
    is.finite(b2) & is.finite(par) & b1 > a1 & b2 > a2
  if (!all(bounds_ok)) {
    stop("Discrete copula rectangles or fitted parameters are invalid.",
      call. = FALSE
    )
  }

  rectangle_mass <- .copula_rectangle_prob(
    u1 = b1, u2 = b2, l1 = a1, l2 = a2,
    family = family_num, par = par, par2 = par2,
    floor_probability = FALSE
  )
  if (any(!is.finite(rectangle_mass) | rectangle_mass <= 0)) {
    stop("A fitted discrete copula rectangle has non-positive probability.",
      call. = FALSE
    )
  }

  uniforms <- .gl_with_preserved_rng(seed, list(
    first = stats::runif(nrow(completed)),
    second = stats::runif(nrow(completed))
  ))

  # Under the fitted copula, the conditional CDF of U1 inside the observed
  # rectangle is the partial rectangle mass divided by its total mass.
  target_mass <- uniforms$first * rectangle_mass
  lower <- a1
  upper <- b1
  max_iterations <- ceiling(log2(1 / tolerance)) + 2L
  for (iteration in seq_len(max_iterations)) {
    active <- (upper - lower) > tolerance
    if (!any(active)) break
    middle <- (lower[active] + upper[active]) / 2
    partial_mass <- .copula_rectangle_prob(
      u1 = middle,
      u2 = b2[active],
      l1 = a1[active],
      l2 = a2[active],
      family = family_num,
      par = par[active],
      par2 = par2[active],
      floor_probability = FALSE
    )
    if (any(!is.finite(partial_mass))) {
      stop("Non-finite partial copula rectangle probability encountered.",
        call. = FALSE
      )
    }
    active_index <- which(active)
    move_lower <- partial_mass < target_mass[active]
    lower[active_index[move_lower]] <- middle[move_lower]
    upper[active_index[!move_lower]] <- middle[!move_lower]
  }
  u1 <- (lower + upper) / 2

  # Conditional on U1, invert the copula h-function inside the second-margin
  # interval. Together with the first inversion this is a draw from the fitted
  # copula conditional on the complete observed probability rectangle.
  h_lower <- .copula_hfunc1(
    u1, a2, family = family_num, par = par, par2 = par2
  )
  h_upper <- .copula_hfunc1(
    u1, b2, family = family_num, par = par, par2 = par2
  )
  if (any(!is.finite(h_lower) | !is.finite(h_upper) |
      h_upper + sqrt(.Machine$double.eps) < h_lower)) {
    stop("Invalid conditional copula bounds encountered during latent completion.",
      call. = FALSE
    )
  }
  h_target <- h_lower + uniforms$second * pmax(h_upper - h_lower, 0)
  u2 <- .sim_invert_hfunc1(
    u1, h_target, family = family_num, par = par, par2 = par2
  )

  completed$u1 <- pmin(pmax(.copula_clamp01(u1), a1), b1)
  completed$u2 <- pmin(pmax(.copula_clamp01(u2), a2), b2)
  attr(completed, "latent_completion") <- list(
    method = "copula_conditional_rectangle_completion",
    draws = as.integer(draws),
    seed = as.integer(seed),
    rng_state_preserved = TRUE
  )
  completed
}

.copula_v2_randomized_discrete_residuals <- function(fit_data, family_num,
                                                      seed = 1L) {
  required <- c(
    "subject", "time", "u_lower", "u_upper", "theta_fit", "zeta_fit"
  )
  if (!all(required %in% names(fit_data))) {
    stop("Discrete diagnostic data do not contain fitted CDF bounds.", call. = FALSE)
  }

  lower <- pmin(pmax(as.numeric(fit_data$u_lower), 0), 1)
  upper <- pmin(pmax(as.numeric(fit_data$u_upper), 0), 1)
  if (any(!is.finite(lower)) || any(!is.finite(upper)) ||
      any(upper + sqrt(.Machine$double.eps) < lower)) {
    stop("Discrete diagnostic CDF bounds are invalid.", call. = FALSE)
  }

  uniforms <- .gl_with_preserved_rng(seed, list(
    marginal = stats::runif(nrow(fit_data)),
    conditional = stats::runif(nrow(fit_data))
  ))
  randomized_fit <- fit_data
  randomized_fit$u <- .copula_v2_clamp01(
    lower + uniforms$marginal * pmax(upper - lower, 0)
  )

  pair_data <- .copula_v2_pair_data(randomized_fit, lags = 1L)
  denominator <- pair_data$u1_upper - pair_data$u1_lower
  upper_numerator <-
    .copula_v2_bicop_cdf(
      pair_data$u1_upper, pair_data$u2_upper,
      family_num, pair_data$theta_pair, pair_data$zeta_pair
    ) -
    .copula_v2_bicop_cdf(
      pair_data$u1_lower, pair_data$u2_upper,
      family_num, pair_data$theta_pair, pair_data$zeta_pair
    )
  lower_numerator <-
    .copula_v2_bicop_cdf(
      pair_data$u1_upper, pair_data$u2_lower,
      family_num, pair_data$theta_pair, pair_data$zeta_pair
    ) -
    .copula_v2_bicop_cdf(
      pair_data$u1_lower, pair_data$u2_lower,
      family_num, pair_data$theta_pair, pair_data$zeta_pair
    )

  conditional_lower <- lower_numerator / denominator
  conditional_upper <- upper_numerator / denominator
  invalid <- !is.finite(denominator) | denominator <= 0 |
    !is.finite(conditional_lower) | !is.finite(conditional_upper)
  conditional_lower <- pmin(pmax(conditional_lower, 0), 1)
  conditional_upper <- pmin(pmax(conditional_upper, conditional_lower), 1)
  conditional_lower[invalid] <- NA_real_
  conditional_upper[invalid] <- NA_real_

  fit_key <- paste(
    as.character(randomized_fit$subject), as.character(randomized_fit$time),
    sep = "::"
  )
  right_key <- paste(
    as.character(pair_data$subject), as.character(pair_data$time_right),
    sep = "::"
  )
  right_row <- match(right_key, fit_key)
  if (anyNA(right_row)) {
    stop("Could not align discrete conditional residuals to fitted rows.", call. = FALSE)
  }
  conditional_uniform <- uniforms$conditional[right_row]
  pair_data$conditional_lower <- conditional_lower
  pair_data$conditional_upper <- conditional_upper
  pair_data$rosenblatt <- .copula_v2_clamp01(
    conditional_lower + conditional_uniform *
      pmax(conditional_upper - conditional_lower, 0)
  )

  series <- randomized_fit[, c("subject", "time"), drop = FALSE]
  series$key <- fit_key
  series$rosenblatt <- NA_real_
  time_values <- randomized_fit$time
  time_levels <- if (is.factor(time_values)) {
    levels(time_values)[levels(time_values) %in% as.character(unique(time_values))]
  } else if (is.numeric(time_values) || is.integer(time_values)) {
    sort(unique(time_values))
  } else {
    sort(unique(as.character(time_values)))
  }
  time_lookup <- stats::setNames(seq_along(time_levels), as.character(time_levels))
  first_rows <- unlist(lapply(
    split(seq_len(nrow(series)), as.character(series$subject)),
    function(index) index[[which.min(unname(time_lookup[
      as.character(series$time[index])
    ]))]]
  ), use.names = FALSE)
  series$rosenblatt[first_rows] <- randomized_fit$u[first_rows]
  series$rosenblatt[right_row] <- pair_data$rosenblatt
  series$rosenblatt <- .copula_v2_clamp01(series$rosenblatt)
  series$z <- stats::qnorm(series$rosenblatt)
  series <- series[is.finite(series$z), c(
    "subject", "time", "rosenblatt", "z"
  ), drop = FALSE]

  series_key <- paste(
    as.character(series$subject), as.character(series$time), sep = "::"
  )
  left_key <- paste(
    as.character(pair_data$subject), as.character(pair_data$time_left),
    sep = "::"
  )
  pair_data$z_prev <- series$z[match(left_key, series_key)]
  pair_data$z_curr <- series$z[match(right_key, series_key)]
  attr(pair_data, "conditional_residual_pairs") <- TRUE

  list(
    fit_data = randomized_fit,
    rosenblatt = series,
    rosenblatt_pairs = pair_data,
    settings = list(
      method = "randomized_discrete_conditional_pit",
      seed = as.integer(seed),
      rng_state_preserved = TRUE
    )
  )
}

.copula_v2_discrete_unavailable <- function(fit_data, pair_data_plot,
                                            pair_data_uniform, type,
                                            dashboard_ncol, plot,
                                            randomized_seed) {
  unavailable <- function(title) {
    .copula_v2_message_plot(
      title = title,
      subtitle = "Requires randomized conditional residuals for a discrete margin",
      message = "Set randomized_residuals = TRUE"
    )
  }
  all_plots <- list(
    empirical_overlay = unavailable("Copula fit"),
    quartile_correlation = unavailable("Kendall's tau by fitted bin"),
    kendall_function = unavailable("Kendall calibration"),
    rosenblatt_by_time = unavailable("Rosenblatt scores by time"),
    rosenblatt_qq = unavailable("Rosenblatt Q-Q"),
    rosenblatt_lag = unavailable("Rosenblatt lag plot"),
    tail_cooccurrence = unavailable("Tail co-occurrence"),
    conditional_tail_exceedance = unavailable("Conditional tail exceedance"),
    residual_lag_correlation = unavailable("Residual correlation by lag")
  )
  dashboard_plots <- .copula_v2_select_dashboard_plots(all_plots, type)
  dashboard <- .copula_v2_arrange_dashboard(dashboard_plots, dashboard_ncol)
  if (isTRUE(plot)) print(dashboard)

  invisible(.copula_v2_plot_result(
    plots = all_plots,
    displayed_plots = dashboard_plots,
    type = type,
    dashboard = dashboard,
    fit_data = fit_data,
    pair_data = pair_data_plot,
    pair_data_uniform = pair_data_uniform,
    rosenblatt = data.frame(),
    rosenblatt_pairs = data.frame(),
    quartile_summary = data.frame(),
    kendall_summary = data.frame(),
    tail_summary = data.frame(),
    conditional_tail_summary = data.frame(),
    residual_lag_summary = data.frame(),
    randomized_settings = list(
      method = "unavailable_for_discrete_margin",
      enabled = FALSE,
      seed = as.integer(randomized_seed)
    )
  ))
}
