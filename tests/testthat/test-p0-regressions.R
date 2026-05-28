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
  expect_true(isTRUE(fit_true$warm_start_joint$used))

  expect_equal(length(fit_false$par), length(fit_true$par))
  expect_false(anyNA(names(fit_true$par)))
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
    "theta.intercept" = 1.54578346770453,
    "theta.time_covariatet2" = 0.0197054581248683,
    "sigma.intercept" = -1.75640640984722,
    "sigma.time_covariatet2" = 0.0480181791101031,
    "sigma.time_covariatet3" = -0.129516565774706,
    "sigma.genderM" = -0.0395172530114576,
    "mu.intercept" = 2.07667100817415,
    "mu.time_covariatet2" = 0.0677478953258934,
    "mu.time_covariatet3" = 0.402764578386612,
    "mu.genderM" = 0.666986768172987,
    "mu.age" = 0.00474355030623712,
    "mu.time_covariatet2:genderM" = 0.353837664435666,
    "mu.time_covariatet3:genderM" = 0.510347291657982
  )

  expected_loglik <- c(
    "marginal" = 18.9898656945639,
    "copula" = -210.533686658939,
    "joint" = -191.543820964375
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
    "theta.intercept" = 1.31117942359909,
    "theta.time_covariatet2" = -0.121173989509146,
    "sigma.intercept" = -0.585801217424111,
    "sigma.time_covariatet2" = -0.147767622745633,
    "sigma.time_covariatet3" = 0.218013053439918,
    "sigma.genderM" = 0.147337617843285,
    "mu.intercept" = 2.93636674301033,
    "mu.time_covariatet2" = -4.97969501316779e-17,
    "mu.time_covariatet3" = -6.51741497346979e-16,
    "mu.genderM" = 6.55850753674967e-16,
    "mu.age" = -3.89801216142305e-18,
    "mu.time_covariatet2:genderM" = 2.76394519732804e-16,
    "mu.time_covariatet3:genderM" = 8.20947302630695e-16
  )

  expected_loglik <- c(
    "marginal" = -66.4886588595176,
    "copula" = 16.0021762532036,
    "joint" = -50.486482606314
  )

  expect_identical(names(fit$par), names(expected_par))
  expect_equal(unname(fit$par), unname(expected_par), tolerance = 1e-6)

  ll <- fit$calc_lik_out_end$log_lik[c("marginal", "copula", "joint")]
  expect_equal(unname(ll), unname(expected_loglik), tolerance = 1e-6)
})

test_that("T007 copula link functions cover Frank Joe and Gumbel", {
  frank <- gamlss.longitudinal:::get_copula_dist("F")
  joe <- gamlss.longitudinal:::get_copula_dist("J")
  gumbel <- gamlss.longitudinal:::get_copula_dist("G")

  expect_equal(frank$parameters, "theta")
  expect_equal(frank$copula_link$theta.linkfun(2.5), 2.5)
  expect_equal(frank$copula_link$theta.linkinv(-1.75), -1.75)

  expect_equal(joe$parameters, "theta")
  expect_equal(joe$copula_link$theta.linkinv(log(4)), 5)
  expect_equal(gumbel$copula_link$theta.linkinv(log(2)), 3)
  expect_equal(gumbel$copula_link$theta.linkinv(log(100)), 17)

  expect_equal(gumbel$parameters, "theta")
  expect_equal(joe$copula_link$theta.linkfun(3), log(2))
  expect_equal(gumbel$copula_link$theta.linkfun(3), log(2))
  expect_equal(gumbel$copula_link$theta.linkfun(100), log(16))
})

test_that("T008 t-copula starting values include theta and zeta", {
  dat <- make_fixture_factor_time_interaction(n_subject = 16L)
  start_dat <- data.frame(
    subject = dat$id,
    time = dat$time_raw,
    response = dat$y
  )

  start_vals <- gamlss.longitudinal:::get_starting_values(
    copula_dist = "t",
    margin_dist = gamlss.dist::NO(),
    dataset = start_dat,
    eta_transform = FALSE
  )

  expect_true(all(c("theta", "zeta") %in% names(start_vals)))
  expect_true(is.finite(unname(start_vals["theta"])))
  expect_true(is.finite(unname(start_vals["zeta"])))
  expect_gt(unname(start_vals["zeta"]), 2)
})

test_that("T009 CG optimizer produces finite improving fixture fit", {
  dat <- make_fixture_factor_time_interaction(n_subject = 24L)

  fit_rs <- fit_fixture_model(
    dat,
    method = "RS",
    include_dlcopdpar = TRUE,
    max_outer_iter = 4L,
    max_inner_iter = 4L,
    outer_stop_crit = 0.2,
    inner_stop_crit = 0.2,
    use_backtracking = TRUE,
    verbose = 0
  )

  fit_cg <- fit_fixture_model(
    dat,
    method = "CG",
    include_dlcopdpar = TRUE,
    max_outer_iter = 8L,
    outer_stop_crit = 0.05,
    use_backtracking = TRUE,
    verbose = 0
  )

  ll_cg <- as.numeric(fit_cg$calc_lik_out_end$log_lik["joint"])
  ll_cg_start <- as.numeric(fit_cg$calc_lik_out$log_lik["joint"])

  expect_true(all(is.finite(fit_rs$calc_lik_out_end$log_lik[c("marginal", "copula", "joint")])))
  expect_true(is.finite(ll_cg))
  expect_true(is.finite(ll_cg_start))
  expect_gte(ll_cg, ll_cg_start - 1e-8)
  expect_true(is.list(fit_cg$convergence))
  expect_false(isTRUE(fit_cg$convergence$raw_loglik_drop))
})
