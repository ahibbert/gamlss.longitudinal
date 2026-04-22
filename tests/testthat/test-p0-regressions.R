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
    "theta.intercept" = 1.5638809760,
    "theta.time_covariatet2" = 0.0005225898,
    "sigma.intercept" = -1.5823105534,
    "sigma.time_covariatet2" = 0.0530642214,
    "sigma.time_covariatet3" = 0.2951152096,
    "sigma.genderM" = 0.3370667005,
    "mu.intercept" = 2.2458203099,
    "mu.time_covariatet2" = 0.0500791139,
    "mu.time_covariatet3" = 0.4329850216,
    "mu.genderM" = 0.4872999761,
    "mu.age" = 0.0011928517,
    "mu.time_covariatet2:genderM" = 0.6481955223,
    "mu.time_covariatet3:genderM" = 0.3560856656
  )

  expected_loglik <- c(
    "marginal" = -5.647582918,
    "copula" = -164.584740630,
    "joint" = -170.232323548
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
    "theta.intercept" = 1.38893630,
    "theta.time_covariatet2" = 0.03580677,
    "sigma.intercept" = -0.59182516,
    "sigma.time_covariatet2" = -0.14707015,
    "sigma.time_covariatet3" = 0.19879541,
    "sigma.genderM" = 0.14987687,
    "mu.intercept" = 2.93636674,
    "mu.time_covariatet2" = 0.00000000,
    "mu.time_covariatet3" = 0.00000000,
    "mu.genderM" = 0.00000000,
    "mu.age" = 0.00000000,
    "mu.time_covariatet2:genderM" = 0.00000000,
    "mu.time_covariatet3:genderM" = 0.00000000
  )

  expected_loglik <- c(
    "marginal" = -66.505475,
    "copula" = 9.982917,
    "joint" = -56.522558
  )

  expect_identical(names(fit$par), names(expected_par))
  expect_equal(unname(fit$par), unname(expected_par), tolerance = 1e-6)

  ll <- fit$calc_lik_out_end$log_lik[c("marginal", "copula", "joint")]
  expect_equal(unname(ll), unname(expected_loglik), tolerance = 1e-6)
})

test_that("T007 copula link functions cover Frank Joe and Gumbel", {
  frank <- get_copula_dist("Frank")
  joe <- get_copula_dist("Joe")
  gumbel <- get_copula_dist("Gumbel")

  expect_equal(frank$parameters, "theta")
  expect_equal(frank$copula_link$theta.linkfun(2.5), 2.5)
  expect_equal(frank$copula_link$theta.linkinv(-1.75), -1.75)

  expect_equal(joe$parameters, "theta")
  expect_equal(joe$copula_link$theta.linkinv(log(4)), 5)
  expect_equal(gumbel$copula_link$theta.linkinv(log(2)), 3)

  expect_equal(gumbel$parameters, "theta")
  expect_equal(joe$copula_link$theta.linkfun(3), log(2))
  expect_equal(gumbel$copula_link$theta.linkfun(3), log(2))
})
