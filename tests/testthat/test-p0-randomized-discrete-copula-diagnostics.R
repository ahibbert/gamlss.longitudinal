test_that("randomized discrete diagnostic controls reject invalid inputs", {
  expect_error(
    .copula_v2_validate_randomized_controls(NA, 1),
    "randomized_residuals"
  )
  expect_error(
    .copula_v2_validate_randomized_controls(TRUE, Inf),
    "randomized_seed"
  )
  expect_error(
    .copula_v2_validate_randomized_controls(TRUE, 1, 0),
    "randomized_draws"
  )
  expect_error(
    .copula_v2_validate_randomized_controls(TRUE, 1, 1.5),
    "randomized_draws"
  )
})

test_that("joint latent completions stay inside and follow fitted rectangles", {
  pair_data <- data.frame(
    subject = 1,
    time_left = "1",
    time_right = "2",
    time_pair = "T1 vs T2",
    lag = 1L,
    u1 = 0.4,
    u2 = 0.8,
    u1_lower = 0.1,
    u1_upper = 0.4,
    u2_lower = 0.55,
    u2_upper = 0.8,
    theta_pair = 0.6,
    zeta_pair = NA_real_,
    tau_fit = 2 / pi * asin(0.6)
  )

  set.seed(431)
  state <- .Random.seed
  completed <- .copula_v2_latent_rectangle_draws(
    pair_data, "N", draws = 1000L, seed = 72L
  )
  expect_identical(.Random.seed, state)
  expect_equal(nrow(completed), 1000L)
  expect_true(all(completed$u1 >= completed$u1_lower &
                    completed$u1 <= completed$u1_upper))
  expect_true(all(completed$u2 >= completed$u2_lower &
                    completed$u2 <= completed$u2_upper))

  mass <- .copula_rectangle_prob(
    u1 = completed$u1_upper,
    u2 = completed$u2_upper,
    l1 = completed$u1_lower,
    l2 = completed$u2_lower,
    family = "N",
    par = completed$theta_pair,
    floor_probability = FALSE
  )
  first_uniform <- .copula_rectangle_prob(
    u1 = completed$u1,
    u2 = completed$u2_upper,
    l1 = completed$u1_lower,
    l2 = completed$u2_lower,
    family = "N",
    par = completed$theta_pair,
    floor_probability = FALSE
  ) / mass
  h_lower <- .copula_hfunc1(
    completed$u1, completed$u2_lower, "N", completed$theta_pair
  )
  h_upper <- .copula_hfunc1(
    completed$u1, completed$u2_upper, "N", completed$theta_pair
  )
  second_uniform <- (
    .copula_hfunc1(completed$u1, completed$u2, "N", completed$theta_pair) -
      h_lower
  ) / (h_upper - h_lower)
  expect_equal(mean(first_uniform), 0.5, tolerance = 0.03)
  expect_equal(mean(second_uniform), 0.5, tolerance = 0.03)
  expect_lt(abs(stats::cor(first_uniform, second_uniform)), 0.08)
  expect_identical(
    attr(completed, "latent_completion")$method,
    "copula_conditional_rectangle_completion"
  )

  repeated <- .copula_v2_latent_rectangle_draws(
    pair_data, "N", draws = 1000L, seed = 72L
  )
  expect_equal(completed, repeated)
})

test_that("fast empirical copula values match the direct definition", {
  x <- c(0.2, 0.2, 0.7, 0.5, NA)
  y <- c(0.4, 0.1, 0.8, 0.4, 0.2)
  direct <- vapply(seq_len(4L), function(i) {
    mean(x[1:4] <= x[i] & y[1:4] <= y[i])
  }, numeric(1))
  expect_equal(.copula_v2_empirical_cdf_2d(x, y), c(direct, NA_real_))
})

test_that("latent rectangle completion supports every registered copula", {
  pair_data <- data.frame(
    subject = 1,
    time_left = "1",
    time_right = "2",
    time_pair = "T1 vs T2",
    lag = 1L,
    u1 = 0.4,
    u2 = 0.8,
    u1_lower = 0.1,
    u1_upper = 0.4,
    u2_lower = 0.55,
    u2_upper = 0.8,
    theta_pair = NA_real_,
    zeta_pair = NA_real_,
    tau_fit = NA_real_
  )
  cases <- list(
    N = c(0.4, NA), C = c(1.2, NA), F = c(2, NA),
    G = c(1.5, NA), J = c(1.5, NA), t = c(0.4, 5)
  )

  for (family in names(cases)) {
    pair_data$theta_pair <- cases[[family]][[1L]]
    pair_data$zeta_pair <- cases[[family]][[2L]]
    completed <- .copula_v2_latent_rectangle_draws(
      pair_data, family, draws = 8L, seed = 4L
    )
    expect_equal(nrow(completed), 8L, info = family)
    expect_true(all(completed$u1 >= 0.1 & completed$u1 <= 0.4),
      info = family
    )
    expect_true(all(completed$u2 >= 0.55 & completed$u2 <= 0.8),
      info = family
    )
  }
})

test_that("randomized conditional PIT uses exact discrete conditional bounds", {
  fit_data <- data.frame(
    subject = c(1, 1),
    time = c(1, 2),
    response = c(0, 1),
    u = c(0.4, 0.7),
    u_lower = c(0, 0.4),
    u_upper = c(0.4, 0.7),
    theta_fit = c(0, NA),
    zeta_fit = c(NA, NA),
    tau_fit = c(0, NA)
  )
  out <- .copula_v2_randomized_discrete_residuals(fit_data, "N", seed = 9L)
  pair <- out$rosenblatt_pairs

  expect_equal(pair$conditional_lower, 0.4, tolerance = 1e-10)
  expect_equal(pair$conditional_upper, 0.7, tolerance = 1e-10)
  expect_gte(pair$rosenblatt, pair$conditional_lower)
  expect_lte(pair$rosenblatt, pair$conditional_upper)
  expect_identical(out$settings$method, "randomized_discrete_conditional_pit")

  fit_data$theta_fit[[1L]] <- 0.5
  dependent <- .copula_v2_randomized_discrete_residuals(
    fit_data, "N", seed = 9L
  )$rosenblatt_pairs
  rectangle_mass <- .copula_rectangle_prob(
    u1 = 0.4, u2 = 0.7, l1 = 0, l2 = 0.4,
    family = "N", par = 0.5, floor_probability = FALSE
  )
  expect_equal(
    (dependent$conditional_upper - dependent$conditional_lower) * 0.4,
    rectangle_mass,
    tolerance = 1e-10
  )
})

test_that("discrete copula diagnostics are opt-in, reproducible, and keep old views", {
  dat <- simulate_longitudinal_dataset(
    n = 14L,
    times = 0:3,
    margin_dist = gamlss.dist::PO(),
    copula_dist = "N",
    margin_params = list(mu = 3),
    copula_params = list(theta = 0.35),
    seed = 2401
  )
  fit <- suppressWarnings(gamlss_longitudinal(
    dataset = dat,
    margin_dist = gamlss.dist::PO(),
    copula_dist = "N",
    time_var = "time",
    subject_var = "subject",
    mu.formula = "response ~ 1",
    theta.formula = "~ 1",
    method = "RS",
    max_outer_iter = 3L,
    max_inner_iter = 4L,
    compute_vcov = FALSE,
    verbose = 0
  ))

  unavailable <- plot_copula_diagnostics(
    fit,
    type = "all",
    randomized_residuals = FALSE,
    plot = FALSE
  )
  expect_identical(unavailable$randomized_settings$enabled, FALSE)
  expect_equal(nrow(unavailable$rosenblatt), 0L)

  set.seed(991)
  seed_before <- .Random.seed
  randomized <- plot_copula_diagnostics(
    fit,
    type = "all",
    ncol = 3,
    randomized_residuals = TRUE,
    randomized_seed = 17L,
    plot = FALSE
  )
  expect_identical(.Random.seed, seed_before)
  expect_identical(length(randomized$displayed_plots), 9L)
  expect_identical(names(randomized$displayed_plots), c(
    "empirical_overlay", "quartile_correlation", "kendall_function",
    "rosenblatt_by_time", "rosenblatt_qq", "rosenblatt_lag",
    "tail_cooccurrence", "conditional_tail_exceedance",
    "residual_lag_correlation"
  ))
  expect_true(nrow(randomized$rosenblatt) > 0L)
  expect_true(all(randomized$rosenblatt$rosenblatt > 0 &
                    randomized$rosenblatt$rosenblatt < 1))
  expect_identical(
    randomized$randomized_settings$method,
    "randomized_discrete_conditional_pit"
  )
  expect_identical(
    randomized$randomized_settings$latent_pair_method,
    "copula_conditional_rectangle_completion"
  )
  expect_identical(randomized$randomized_settings$latent_draws, 5L)
  expect_true(all(randomized$pair_data_uniform$u1 >=
                    randomized$pair_data_uniform$u1_lower))
  expect_true(all(randomized$pair_data_uniform$u1 <=
                    randomized$pair_data_uniform$u1_upper))
  expect_true(all(randomized$pair_data_uniform$u2 >=
                    randomized$pair_data_uniform$u2_lower))
  expect_true(all(randomized$pair_data_uniform$u2 <=
                    randomized$pair_data_uniform$u2_upper))
  expect_match(
    randomized$plots$quartile_correlation$labels$title,
    "Latent Kendall"
  )
  expect_match(
    randomized$plots$tail_cooccurrence$labels$title,
    "Latent tail"
  )

  default_plot <- plot_copula_diagnostics(
    fit,
    type = "all",
    ncol = 3,
    randomized_seed = 17L,
    plot = FALSE
  )
  expect_equal(default_plot$rosenblatt, randomized$rosenblatt)
  expect_identical(
    default_plot$randomized_settings$method,
    "randomized_discrete_conditional_pit"
  )

  repeated <- plot_copula_diagnostics(
    fit,
    type = "all",
    ncol = 3,
    randomized_residuals = TRUE,
    randomized_seed = 17L,
    plot = FALSE
  )
  expect_equal(randomized$rosenblatt, repeated$rosenblatt)
  expect_equal(randomized$pair_data_uniform, repeated$pair_data_uniform)

  changed <- plot_copula_diagnostics(
    fit,
    randomized_residuals = TRUE,
    randomized_seed = 18L,
    plot = FALSE
  )
  expect_false(isTRUE(all.equal(
    randomized$rosenblatt$rosenblatt,
    changed$rosenblatt$rosenblatt
  )))
  expect_false(isTRUE(all.equal(
    randomized$pair_data_uniform[, c("u1", "u2")],
    changed$pair_data_uniform[, c("u1", "u2")]
  )))

  check_off <- check_model(fit)
  copula_off <- check_off$checks[
    check_off$checks$area == "Copula fit", , drop = FALSE
  ]
  expect_identical(copula_off$status, "unavailable")
  expect_true(all(check_off$residual_dependence$scope_status ==
                    "unavailable_for_discrete_margin"))

  check_on <- check_model(
    fit,
    randomized_residuals = TRUE,
    pit_seed = 17L
  )
  expect_true(all(check_on$residual_dependence$scope_status == "available"))
  expect_identical(
    check_on$conditional_residual_method$method,
    "randomized_discrete_conditional_pit"
  )
  expect_identical(check_on$conditional_residual_method$seed, 17L)
})
