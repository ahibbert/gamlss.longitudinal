test_that("simulate_longitudinal_dataset returns balanced long data", {
  dat <- simulate_longitudinal_dataset(
    n = 6,
    times = c("t1", "t2", "t3"),
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "N",
    copula_params = list(theta = 0.2),
    seed = 123
  )

  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 18)
  expect_true(all(c("subject", "time", "response", "u", "true_mu", "true_sigma", "true_theta") %in% names(dat)))
  expect_equal(levels(dat$subject), as.character(seq_len(6)))
  expect_equal(unique(dat$time), c("t1", "t2", "t3"))
  expect_true(all(is.finite(dat$response)))
  expect_true(all(dat$u > 0 & dat$u < 1))
  expect_equal(dat$true_theta[dat$time == "t1"], rep(NA_real_, 6))
  expect_equal(dat$true_theta[dat$time != "t1"], rep(0.2, 12))
})

test_that("simulation seed helper restores caller RNG state", {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit({
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  set.seed(999)
  before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  restore <- gamlss.longitudinal:::.sim_set_seed_restore(123)
  expect_false(identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), before))
  expect_invisible(restore())
  expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), before)

  rm(".Random.seed", envir = .GlobalEnv)
  restore <- gamlss.longitudinal:::.sim_set_seed_restore(123)
  expect_true(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_invisible(restore())
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("simulation uniform bounds helper validates and clamps uniforms", {
  u <- c(0, 0.1, 0.9, 1)

  expect_equal(gamlss.longitudinal:::.sim_apply_u_bounds(u, NULL), u)
  expect_equal(
    gamlss.longitudinal:::.sim_apply_u_bounds(u, c(0.05, 0.95)),
    c(0.05, 0.1, 0.9, 0.95)
  )
  expect_error(
    gamlss.longitudinal:::.sim_validate_u_bounds(c(0.5, 0.5)),
    "finite increasing",
    fixed = TRUE
  )
  expect_error(
    gamlss.longitudinal:::.sim_validate_u_bounds(c(-0.1, 0.9)),
    "inside [0, 1]",
    fixed = TRUE
  )
})

test_that("simulate_longitudinal_dataset supports function and time-varying parameters", {
  dat <- simulate_longitudinal_dataset(
    n = 5,
    times = 1:4,
    margin_dist = gamlss.dist::GA(mu.link = "log", sigma.link = "log"),
    margin_params = list(
      mu = function(data) exp(0.2 + 0.1 * data$.sim_time_index),
      sigma = c(0.4, 0.5, 0.6, 0.7)
    ),
    copula_dist = "C",
    copula_params = list(tau = c(0.1, 0.2, 0.3)),
    seed = 456
  )

  expect_equal(nrow(dat), 20)
  expect_true(all(dat$response > 0))
  expect_equal(dat$true_sigma, rep(c(0.4, 0.5, 0.6, 0.7), times = 5))
  expect_equal(
    unique(dat$true_theta[dat$time == 2]),
    .copula_tau_to_par("C", 0.1),
    tolerance = 1e-12
  )
  expect_equal(
    unique(dat$true_theta[dat$time == 4]),
    .copula_tau_to_par("C", 0.3),
    tolerance = 1e-12
  )
})

test_that("simulation smooth and factor helper functions are readable building blocks", {
  x <- c(2, 4, 6)
  expect_equal(sim_rescale01(x), c(0, 0.5, 1))
  expect_equal(sim_rescale01(c(3, 3, 3)), c(0, 0, 0))

  expect_equal(
    sim_factor_effect(
      factor(c("control", "active", "active"), levels = c("control", "active")),
      c(control = 0, active = 0.4)
    ),
    c(0, 0.4, 0.4)
  )
  expect_error(
    sim_factor_effect(c("control", "missing"), c(control = 0)),
    "Missing factor effect",
    fixed = TRUE
  )

  smooth_values <- cbind(
    linear = sim_smooth_linear(x = c(0, 0.5, 1)),
    sin = sim_smooth_sin(x = c(0, 0.25, 0.5)),
    bump = sim_smooth_bump(x = c(0, 0.5, 1)),
    sigmoid = sim_smooth_sigmoid(x = c(0, 0.5, 1)),
    u = sim_smooth_u(x = c(0, 0.5, 1)),
    wiggle = sim_smooth_wiggle(x = c(0, 0.5, 1))
  )
  expect_true(all(is.finite(smooth_values)))
})

test_that("simulate_longitudinal_covariates expands subject and observation covariates", {
  dat <- simulate_longitudinal_dataset(
    n = 4,
    times = 1:3,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "N",
    copula_params = list(theta = 0.1),
    covariates = function(base) {
      simulate_longitudinal_covariates(
        base,
        subject = list(group = function(subject_data) factor(rep(c("A", "B"), length.out = nrow(subject_data)))),
        observation = list(time_scaled = function(long_data) sim_rescale01(long_data$time))
      )
    },
    seed = 789
  )

  expect_true(all(c("group", "time_scaled") %in% names(dat)))
  expect_equal(length(unique(dat$group[dat$subject == "1"])), 1)
  expect_equal(unique(dat$time_scaled), c(0, 0.5, 1))
})

test_that("fitted-model newdata simulation builds a deterministic time grid", {
  nd <- data.frame(
    subject = c("a", "a", "b", "b"),
    time = c(2, 1, 2, 1),
    response = NA_real_
  )
  object <- list(response_margin = 1:3)

  grid <- .gl_simulation_time_grid(nd, object)

  expect_equal(grid$time_for_grid, nd$time)
  expect_equal(grid$time_levels, 1:3)
  expect_equal(grid$time_idx, c(2, 1, 2, 1))
})

test_that("fitted-model newdata simulation respects factor time order", {
  nd <- data.frame(
    subject = c("a", "a", "b", "b"),
    time = c(1, 2, 1, 2),
    time_covariate = factor(c("followup", "baseline", "followup", "baseline"),
      levels = c("baseline", "followup", "later")
    ),
    response = NA_real_
  )
  object <- list(response_margin = factor(character(), levels = c("baseline", "followup", "later")))

  grid <- .gl_simulation_time_grid(nd, object)

  expect_equal(grid$time_for_grid, nd$time_covariate)
  expect_equal(grid$time_levels, c("baseline", "followup", "later"))
  expect_equal(grid$time_idx, c(2, 1, 2, 1))
})

test_that("fitted-model newdata simulation rejects duplicate subject-time rows", {
  nd <- data.frame(
    subject = c("a", "a"),
    time = c(1, 1),
    response = NA_real_
  )
  object <- list(response_margin = 1:2)

  expect_error(
    .gl_simulation_time_grid(nd, object),
    "at most one row per subject/time",
    fixed = TRUE
  )
})

test_that("fitted-model newdata simulation aligns dependence parameters", {
  nd_eval <- data.frame(
    .gl_sim_time_idx = c(1L, 2L, 3L, 1L, 2L, 3L)
  )
  time_levels <- 1:3

  expect_equal(
    .gl_simulation_align_dependence_parameter(numeric(0), 6L, nd_eval, time_levels),
    rep(NA_real_, 6)
  )
  expect_equal(
    .gl_simulation_align_dependence_parameter(11:16, 6L, nd_eval, time_levels),
    as.numeric(11:16)
  )
  expect_equal(
    .gl_simulation_align_dependence_parameter(c(0.2, 0.4, 0.6, 0.8), 6L, nd_eval, time_levels),
    c(0.2, 0.4, NA, 0.6, 0.8, NA)
  )
  expect_equal(
    .gl_simulation_align_dependence_parameter(c(3, 4), 6L, nd_eval, time_levels),
    c(3, 4, 3, 4, 3, 4)
  )
})

test_that("simulate_longitudinal_dataset runs across native copula families", {
  cases <- list(
    N = list(theta = 0.25),
    C = list(theta = 1.2),
    F = list(theta = 2.5),
    G = list(theta = 1.4),
    J = list(theta = 1.6),
    t = list(theta = 0.25, zeta = 5)
  )

  for (family in names(cases)) {
    dat <- simulate_longitudinal_dataset(
      n = 4,
      times = 1:3,
      margin_dist = gamlss.dist::NO(),
      margin_params = list(mu = 0, sigma = 1),
      copula_dist = family,
      copula_params = cases[[family]],
      seed = 10
    )
    expect_equal(nrow(dat), 12)
    expect_true(all(is.finite(dat$response)))
    expect_true(all(dat$u > 0 & dat$u < 1))
  }
})

test_that("simulated fixed Normal Gaussian parameters recover under separate RS", {
  suppressPackageStartupMessages(library(gamlss.dist))

  true_mu <- 1.25
  true_sigma <- 0.7
  true_theta <- 0.35

  dat <- simulate_longitudinal_dataset(
    n = 600,
    times = 1:4,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    margin_params = list(mu = true_mu, sigma = true_sigma),
    copula_params = list(theta = true_theta),
    seed = 101
  )

  fit <- gamlss_longitudinal(
    dataset = dat,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    time_var = "time",
    subject_var = "subject",
    mu.formula = response ~ 1,
    sigma.formula = ~ 1,
    theta.formula = ~ 1,
    include_dlcopdpar = FALSE,
    method = "RS",
    warm_start_joint = FALSE,
    max_outer_iter = 80,
    max_inner_iter = 80,
    outer_stop_crit = 1e-5,
    inner_stop_crit = 1e-5,
    verbose = 0,
    compute_vcov = FALSE
  )

  theta_link <- gamlss.longitudinal:::get_copula_dist("N")$copula_link$theta.linkinv
  estimate <- c(
    mu = unname(fit$par["mu.intercept"]),
    sigma = exp(unname(fit$par["sigma.intercept"])),
    theta = theta_link(unname(fit$par["theta.intercept"]))
  )
  truth <- c(mu = true_mu, sigma = true_sigma, theta = true_theta)

  expect_equal(estimate["mu"], truth["mu"], tolerance = 0.03)
  expect_equal(estimate["sigma"], truth["sigma"], tolerance = 0.03)
  expect_equal(estimate["theta"], truth["theta"], tolerance = 0.06)
  expect_true(all(is.finite(fit$calc_lik_out_end$log_lik[c("marginal", "copula", "joint")])))
})
