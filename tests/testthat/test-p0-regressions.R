test_that("T001 factor time is preserved as model covariate", {
  dat <- make_fixture_factor_time_interaction()

  out <- capture_warnings(fit_fixture_model(dat, include_dlcopdpar = TRUE))
  fit <- out$value

  expect_s3_class(fit, "gamlss.longitudinal")

  warn_txt <- paste(out$warnings, collapse = "\n")
  expect_false(grepl("not meaningful for factors", warn_txt, fixed = TRUE))

  mu_cols <- colnames(fit$model_matrix$x$mu)
  expect_true(any(grepl("^time_covariate", mu_cols)))
})

test_that("T002 interaction expansion contains base and interaction columns", {
  dat <- make_fixture_factor_time_interaction()
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  mu_cols <- colnames(fit$model_matrix$x$mu)

  expect_true(any(grepl("^time_covariate", mu_cols)))
  expect_true(any(grepl("^gender", mu_cols)))
  expect_true(any(grepl(":", mu_cols, fixed = TRUE)))

  par_names <- names(fit$par)
  expect_true(any(grepl("mu\\..*:", par_names)))
})

test_that("T003 include_dlcopdpar TRUE path yields finite likelihood outputs", {
  dat <- make_fixture_factor_time_interaction()
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  ll <- fit$calc_lik_out_end$log_lik
  expect_true(all(is.finite(ll[c("marginal", "copula", "joint")])))
})

test_that("T004 dlcopdpar TRUE/FALSE parity smoke test", {
  dat <- make_fixture_factor_time_interaction()

  fit_false <- fit_fixture_model(dat, include_dlcopdpar = FALSE)
  fit_true <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  expect_true(is.finite(fit_false$calc_lik_out_end$log_lik["joint"]))
  expect_true(is.finite(fit_true$calc_lik_out_end$log_lik["joint"]))

  expect_equal(length(fit_false$par), length(fit_true$par))
  expect_setequal(names(fit_false$par), names(fit_true$par))
})

test_that("T005 baseline fit fingerprint stays stable with use_backtracking FALSE", {
  dat <- make_fixture_factor_time_interaction(n_subject = 24L)

  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = FALSE,
    use_backtracking = FALSE,
    max_outer_iter = 3L,
    max_inner_iter = 3L,
    outer_stop_crit = 0.5,
    inner_stop_crit = 0.5,
    verbose = 0
  )

  expected_par <- c(
    "theta.intercept" = 1.5641422709077644,
    "theta.time_covariate.L" = 0.00036952677386012129,
    "sigma.intercept" = -1.4662507430579528,
    "sigma.time_covariate.L" = 0.20867796595747662,
    "sigma.time_covariate.Q" = 0.077153524488408259,
    "sigma.genderM" = 0.3370667005163685,
    "mu.intercept" = 2.4068416883317543,
    "mu.time_covariate.L" = 0.30616664489726197,
    "mu.time_covariate.Q" = 0.13587596959788967,
    "mu.genderM" = 0.8220603720064974,
    "mu.age" = 0.0011928517367271545,
    "mu.time_covariate.L:genderM" = 0.25179058880563637,
    "mu.time_covariate.Q:genderM" = -0.38387806349188175
  )

  expected_loglik <- c(
    "marginal" = -5.6475829180525805,
    "copula" = -164.58474063030087,
    "joint" = -170.23232354835346
  )

  expect_identical(names(fit$par), names(expected_par))
  expect_equal(unname(fit$par), unname(expected_par), tolerance = 1e-9)

  ll <- fit$calc_lik_out_end$log_lik[c("marginal", "copula", "joint")]
  expect_equal(unname(ll), unname(expected_loglik), tolerance = 1e-9)
})

test_that("T006 baseline fit fingerprint stays stable with use_backtracking TRUE", {
  dat <- make_fixture_factor_time_interaction(n_subject = 24L)

  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = FALSE,
    use_backtracking = TRUE,
    max_outer_iter = 3L,
    max_inner_iter = 3L,
    outer_stop_crit = 0.5,
    inner_stop_crit = 0.5,
    verbose = 0
  )

  expected_par <- c(
    "theta.intercept" = 1.406839689988608,
    "theta.time_covariate.L" = 0.025319210421813636,
    "sigma.intercept" = -0.5745834057904562,
    "sigma.time_covariate.L" = 0.14056958293200163,
    "sigma.time_covariate.Q" = 0.20124015862764695,
    "sigma.genderM" = 0.14987686912288116,
    "mu.intercept" = 2.936366742574324,
    "mu.time_covariate.L" = 0,
    "mu.time_covariate.Q" = 0,
    "mu.genderM" = 0,
    "mu.age" = 0,
    "mu.time_covariate.L:genderM" = 0,
    "mu.time_covariate.Q:genderM" = 0
  )

  expected_loglik <- c(
    "marginal" = -66.50547493356043,
    "copula" = 9.98291714810978,
    "joint" = -56.52255778545065
  )

  expect_identical(names(fit$par), names(expected_par))
  expect_equal(unname(fit$par), unname(expected_par), tolerance = 1e-6)

  ll <- fit$calc_lik_out_end$log_lik[c("marginal", "copula", "joint")]
  expect_equal(unname(ll), unname(expected_loglik), tolerance = 1e-6)
})
