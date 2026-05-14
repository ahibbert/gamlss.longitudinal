test_that("T101 required-input validation errors are explicit", {
  dat <- make_fixture_factor_time_interaction(n_subject = 10L)

  expect_error(
    fit_fixture_model(dat, time_var = "not_a_col"),
    "time_var='not_a_col' not found",
    fixed = TRUE
  )

  expect_error(
    fit_fixture_model(dat, subject_var = "not_a_col"),
    "subject_var='not_a_col' not found",
    fixed = TRUE
  )

  expect_error(
    fit_fixture_model(dat, mu_formula = "response_not_present ~ time_raw + gender"),
    "response variable 'response_not_present' not found",
    fixed = TRUE
  )
})

test_that("T102 duplicate subject-time combinations are rejected", {
  dat <- make_fixture_with_duplicate_subject_time()

  expect_error(
    fit_fixture_model(dat),
    "Duplicate subject/time combinations found",
    fixed = TRUE
  )
})

test_that("T103 full-grid expansion restores structural missing rows", {
  dat_full <- make_fixture_factor_time_interaction(n_subject = 16L)
  dat <- make_fixture_with_structural_missing_rows()
  expected_full_n <- length(unique(dat$id)) * length(unique(dat$time_raw))

  # First fit on complete grid to obtain stable starting values.
  fit_full <- fit_fixture_model(dat_full, include_dlcopdpar = TRUE)

  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = TRUE,
    start_from = fit_full$par
  )

  expect_equal(length(fit$response), expected_full_n)
  expect_true(any(is.na(fit$response)))
})

test_that("T104 missingness hard stops trigger for margin and pair failures", {
  dat_margin_missing <- make_fixture_missingness_margin_all_missing()
  expect_error(
    fit_fixture_model(dat_margin_missing),
    "100% missing response values detected",
    fixed = TRUE
  )

  dat_pair_missing <- make_fixture_missingness_zero_complete_pairs()
  expect_error(
    fit_fixture_model(dat_pair_missing),
    "100% missing complete copula pairs detected",
    fixed = TRUE
  )
})

test_that("T105 formula normalization accepts compact string formulas", {
  dat <- make_fixture_factor_time_interaction(n_subject = 12L)

  fit <- fit_fixture_model(
    dat,
    mu_formula = "y ~ time_raw * gender + age",
    sigma_formula = "time_raw + gender",
    nu_formula = "1",
    tau_formula = "1",
    theta_formula = "time_raw",
    include_dlcopdpar = TRUE
  )

  expect_s3_class(fit, "gamlss.longitudinal")
  expect_true(all(c("mu", "sigma", "theta") %in% names(fit$model_matrix$x)))
})

test_that("T106 smooth and fixed terms coexist in model setup", {
  dat <- make_fixture_factor_time_interaction(n_subject = 20L)

  fit <- fit_fixture_model(
    dat,
    mu_formula = "y ~ time_raw * gender + s(age, bs='ps')",
    sigma_formula = "~ time_raw + gender",
    include_dlcopdpar = FALSE
  )

  expect_true(ncol(fit$model_matrix$x$mu) > 0)
  expect_true(length(fit$model_matrix$s$mu) >= 1)
  s_name <- names(fit$model_matrix$s$mu)[1]
  expect_true(nrow(fit$model_matrix$s$mu[[s_name]]) == length(fit$response))
})

test_that("T107 starting values fallback warning is emitted for non-finite tau", {
  suppressPackageStartupMessages({
    library(gamlss.dist)
    library(VineCopula)
  })

  dat <- make_fixture_factor_time_interaction(n_subject = 14L)
  dat$y[] <- 1

  warn_out <- capture_warnings(
    gamlss.longitudinal::get_starting_values(
      copula_dist = "N",
      margin_dist = gamlss.dist::NO(),
      dataset = data.frame(time = as.integer(dat$time_raw), response = dat$y),
      eta_transform = FALSE
    )
  )

  expect_true(any(grepl("Non-finite Kendall tau", warn_out$warnings, fixed = TRUE)))
})

test_that("T108 intercept-only formulas are accepted across parameters", {
  dat <- make_fixture_factor_time_interaction(n_subject = 12L)

  fit <- fit_fixture_model(
    dat,
    mu_formula = "y ~ 1",
    sigma_formula = "~ 1",
    nu_formula = "~ 1",
    tau_formula = "~ 1",
    theta_formula = "~ 1",
    include_dlcopdpar = TRUE
  )

  expect_s3_class(fit, "gamlss.longitudinal")
  # For NO() margins the active parameters are mu/sigma plus copula theta.
  expected_pars <- c("mu", "sigma", "theta")
  expect_true(all(expected_pars %in% names(fit$model_matrix$x)))
  expect_true(all(vapply(fit$model_matrix$x[expected_pars], ncol, integer(1)) >= 1L))
})

test_that("T109 CG optimizer runs with fixed and smooth terms", {
  dat <- make_fixture_factor_time_interaction(n_subject = 14L)

  fit <- fit_fixture_model(
    dat,
    method = "CG",
    include_dlcopdpar = TRUE,
    mu_formula = "y ~ time_raw * gender + s(age, bs='ps')",
    sigma_formula = "~ time_raw + gender",
    max_outer_iter = 2,
    outer_stop_crit = 1
  )

  expect_s3_class(fit, "gamlss.longitudinal")
  expect_identical(fit$optim_method, "CG")
  expect_true(is.finite(as.numeric(fit$calc_lik_out_end$log_lik["joint"])))
})
