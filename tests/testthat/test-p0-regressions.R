test_that("T000 underscore constructor preserves dotted S3 surface", {
  suppressPackageStartupMessages(library(gamlss.dist))

  dat <- make_fixture_factor_time_interaction(n_subject = 12L)
  fit <- capture_warnings(
    gamlss.longitudinal::gamlss_longitudinal(
      dataset = dat,
      margin_dist = gamlss.dist::NO(),
      copula_dist = "N",
      time_var = "time_raw",
      subject_var = "id",
      mu.formula = y ~ time_raw * gender + age,
      sigma.formula = ~ time_raw + gender,
      theta.formula = ~ time_raw,
      include_dlcopdpar = FALSE,
      method = "RS",
      warm_start_joint = FALSE,
      max_outer_iter = 2,
      max_inner_iter = 2,
      outer_stop_crit = 1,
      inner_stop_crit = 1,
      compute_vcov = FALSE,
      verbose = 0
    )
  )$value

  expect_s3_class(fit, "gamlss.longitudinal")
  expect_false(inherits(fit, "gamlss_longitudinal"))
  expect_identical(stats::coef(fit), fit$par)

  s <- summary(fit, include_vcov = FALSE)
  expect_s3_class(s, "summary.gamlss.longitudinal")

  pred <- stats::predict(fit)
  expect_type(pred, "double")
  expect_equal(length(pred), length(fit$response))
})

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

test_that("T004b warm-start convergence warnings do not escape final converged fit", {
  dat <- make_fixture_factor_time_interaction(n_subject = 12L)
  escaped_warnings <- character(0)
  fit <- NULL

  capture.output({
    fit <- withCallingHandlers(
      gamlss.longitudinal::gamlss_longitudinal(
        dataset = dat,
        margin_dist = gamlss.dist::NO(),
        copula_dist = "N",
        time_var = "time_raw",
        subject_var = "id",
        mu.formula = y ~ time_raw * gender + age,
        sigma.formula = ~ time_raw + gender,
        theta.formula = ~ time_raw,
        include_dlcopdpar = TRUE,
        warm_start_joint = TRUE,
        warm_start_joint_iter = 2L,
        max_outer_iter = 10L,
        max_inner_iter = 2L,
        outer_stop_crit = 1,
        inner_stop_crit = 1,
        compute_vcov = FALSE,
        verbose = 0
      ),
      warning = function(w) {
        escaped_warnings <<- c(escaped_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  })

  expect_true(isTRUE(fit$convergence$converged))
  expect_false(any(grepl(
    "Model returned without satisfying the optimizer convergence contract",
    escaped_warnings,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "stop reason: max_iterations",
    fit$warm_start_joint$captured_warnings,
    fixed = TRUE
  )))
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
    "theta.intercept" = 1.52432081401975,
    "theta.time_covariatet2" = -0.137026522361319,
    "sigma.intercept" = -2.0701396400213,
    "sigma.time_covariatet2" = -0.0395820535922465,
    "sigma.time_covariatet3" = 0.0339199214749465,
    "sigma.genderM" = 0.197718129854906,
    "mu.intercept" = 1.76625280226093,
    "mu.time_covariatet2" = 0.264278685913664,
    "mu.time_covariatet3" = 0.501336192479701,
    "mu.genderM" = 0.713406173253638,
    "mu.age" = 0.00989723440416187,
    "mu.time_covariatet2:genderM" = 0.181226838909403,
    "mu.time_covariatet3:genderM" = 0.597042016726706
  )

  expected_loglik <- c(
    "marginal" = 40.7567885070402,
    "copula" = -193.140969909835,
    "joint" = -152.384181402795
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
    "theta.intercept" = 0.779286043708592,
    "theta.time_covariatet2" = 0.141012935738407,
    "sigma.intercept" = -1.16217675162619,
    "sigma.time_covariatet2" = -0.150525886015387,
    "sigma.time_covariatet3" = 0.224941716928655,
    "sigma.genderM" = 0.141584174367881,
    "mu.intercept" = 2.26606264132697,
    "mu.time_covariatet2" = 0.109926147535223,
    "mu.time_covariatet3" = 0.427381728677277,
    "mu.genderM" = 0.819973822728021,
    "mu.age" = 0.00244574299863918,
    "mu.time_covariatet2:genderM" = 0.118681177183672,
    "mu.time_covariatet3:genderM" = 0.182488704562899
  )

  expected_loglik <- c(
    "marginal" = -5.2728691724679,
    "copula" = 3.79673322608529,
    "joint" = -1.47613594638261
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
