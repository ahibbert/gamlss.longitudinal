test_that("prediction distinguishes response mean, median, and mu for non-Normal margins", {
  dat <- simulate_longitudinal_dataset(
    n = 18,
    times = 0:2,
    margin_dist = gamlss.dist::LOGNO(),
    copula_dist = "N",
    margin_params = list(mu = 1, sigma = 0.35),
    copula_params = list(theta = 0.2),
    seed = 501
  )

  fit <- gamlss_longitudinal(
    dataset = dat,
    margin_dist = gamlss.dist::LOGNO(),
    copula_dist = "N",
    time_var = "time",
    subject_var = "subject",
    mu.formula = response ~ 1,
    sigma.formula = ~ 1,
    theta.formula = ~ 1,
    max_outer_iter = 2,
    max_inner_iter = 2,
    outer_stop_crit = 1,
    inner_stop_crit = 1,
    compute_vcov = FALSE,
    verbose = 0
  )

  p_mu <- predict(fit, type = "mu")
  p_response <- predict(fit, type = "response")
  p_mean <- predict(fit, type = "mean")
  p_median <- predict(fit, type = "median")

  expect_equal(p_response, p_mu)
  expect_true(all(is.finite(p_mean)))
  expect_true(all(is.finite(p_median)))
  expect_true(mean(abs(p_mean - p_mu)) > 0.5)
  expect_true(mean(abs(p_median - exp(p_mu))) < 1e-6)
})

test_that("benchmark scoring uses response means for gamlss.longitudinal fits", {
  dat <- simulate_longitudinal_dataset(
    n = 16,
    times = 0:2,
    margin_dist = gamlss.dist::LOGNO(),
    copula_dist = "N",
    margin_params = list(mu = 1, sigma = 0.25),
    copula_params = list(theta = 0.2),
    seed = 502
  )

  fit <- gamlss_longitudinal(
    dataset = dat,
    margin_dist = gamlss.dist::LOGNO(),
    copula_dist = "N",
    time_var = "time",
    subject_var = "subject",
    mu.formula = response ~ 1,
    sigma.formula = ~ 1,
    theta.formula = ~ 1,
    max_outer_iter = 2,
    max_inner_iter = 2,
    outer_stop_crit = 1,
    inner_stop_crit = 1,
    compute_vcov = FALSE,
    verbose = 0
  )

  bench <- benchmark_standard_models(
    data = dat,
    formula = response ~ 1,
    subject_var = "subject",
    family = "gaussian",
    comparators = character(0),
    fit = fit
  )

  pred_mean <- predict(fit, newdata = dat, type = "mean")
  expected_rmse <- sqrt(mean((pred_mean - dat$response)^2))
  expect_equal(bench$results$benchmark_rmse, expected_rmse, tolerance = 1e-10)
})

test_that("verbose = 0 suppresses fitter iteration and model-fit banners", {
  dat <- simulate_longitudinal_dataset(
    n = 8,
    times = 0:2,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    seed = 503
  )

  output <- capture.output(
    suppressWarnings(gamlss_longitudinal(
      dataset = dat,
      margin_dist = gamlss.dist::NO(),
      copula_dist = "N",
      time_var = "time",
      subject_var = "subject",
      mu.formula = response ~ 1,
      sigma.formula = ~ 1,
      theta.formula = ~ 1,
      max_outer_iter = 2,
      max_inner_iter = 2,
      outer_stop_crit = 1,
      inner_stop_crit = 1,
      compute_vcov = FALSE,
      warm_start_joint = FALSE,
      verbose = 0
    ))
  )

  expect_false(any(grepl("OUTER ITERATION|MODEL FIT", output)))
})

test_that("prediction evaluation value helper preserves q/y policy", {
  diag_data <- list(
    response = c(1, 2, NA),
    subject = c(1, 1, 2)
  )

  expect_equal(
    gamlss.longitudinal:::.gl_prediction_eval_values(5, "q", "cdf", diag_data, require_response = FALSE),
    c(5, 5, 5)
  )
  expect_equal(
    gamlss.longitudinal:::.gl_prediction_eval_values(c(1, 2, 3), "q", "cdf", diag_data, require_response = FALSE),
    c(1, 2, 3)
  )
  expect_equal(
    gamlss.longitudinal:::.gl_prediction_eval_values(NULL, "q", "cdf", diag_data, require_response = TRUE),
    c(1, 2, NA)
  )
  expect_error(
    gamlss.longitudinal:::.gl_prediction_eval_values(NULL, "q", "cdf", diag_data, require_response = FALSE),
    "'q' is required for type = 'cdf' when no response is available.",
    fixed = TRUE
  )
  expect_error(
    gamlss.longitudinal:::.gl_prediction_eval_values(c(1, 2), "y", "density", diag_data, require_response = FALSE),
    "'y' must be length 1 or match the number of prediction rows.",
    fixed = TRUE
  )
})

test_that("simulation covariate callback may return input data plus new columns", {
  dat <- simulate_longitudinal_dataset(
    n = 6,
    times = 0:2,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    covariates = function(d) {
      d$group <- factor(rep(c("A", "B"), length.out = nrow(d)))
      d
    },
    margin_params = list(mu = function(d) ifelse(d$group == "B", 1, 0), sigma = 1),
    seed = 504
  )

  expect_true("group" %in% names(dat))
  expect_false(any(startsWith(names(dat), ".sim_")))
})
