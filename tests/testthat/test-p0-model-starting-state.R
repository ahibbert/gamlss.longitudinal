test_that("initial parameter state builds fixed starts from family starts", {
  prepared <- gamlss.longitudinal:::.gl_prepare_fit_data(
    dataset = data.frame(
      id = c("a", "a", "b", "b"),
      visit = c(1, 2, 1, 2),
      y = c(1.0, 2.0, 1.5, 2.5),
      x = c(0, 1, 0, 1)
    ),
    time_var = "visit",
    subject_var = "id",
    mu.formula = y ~ x,
    sigma.formula = ~ 1,
    nu.formula = ~ 1,
    tau.formula = ~ 1,
    theta.formula = ~ 1,
    zeta.formula = ~ 1,
    verbose = 0
  )
  bundle <- gamlss.longitudinal:::.gl_build_model_matrix_bundle(
    formulas_int = prepared$formulas_int,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    dataset = prepared$dataset
  )

  state <- gamlss.longitudinal:::.gl_build_initial_parameter_state(
    start_from = NA,
    mm = bundle$mm,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    dataset = prepared$dataset
  )

  expect_true(any(grepl("^mu[.]", names(state$par_cov))))
  expect_true(any(grepl("^sigma[.]", names(state$par_cov))))
  expect_true(any(grepl("^theta[.]", names(state$par_cov))))
  expect_true(any(grepl("^mu[.].*x", names(state$par_cov))))
  expect_equal(unname(state$par_cov[grep("^mu[.].*x", names(state$par_cov))[1]]), 0)
  expect_named(state$par_s, names(bundle$mm$x))
  expect_named(state$df_s, names(bundle$mm$x))
  expect_named(state$lambda_s, names(bundle$mm$x))
})

test_that("starting-value moment helpers handle degenerate and finite inputs", {
  expect_equal(gamlss.longitudinal:::.starting_moment_skewness(c(1, 1, 1)), 0)
  expect_equal(gamlss.longitudinal:::.starting_moment_skewness(c(1, NA, Inf)), 0)
  expect_equal(gamlss.longitudinal:::.starting_moment_kurtosis(c(1, 1, 1, 1)), 3)
  expect_true(is.finite(gamlss.longitudinal:::.starting_moment_skewness(c(1, 2, 5, 10))))
  expect_true(is.finite(gamlss.longitudinal:::.starting_moment_kurtosis(c(1, 2, 5, 10))))
})

test_that("starting margin helper returns family-parameter starts", {
  dat <- data.frame(response = c(1, 2, 3, 4), time = c(1, 2, 1, 2), subject = c("a", "a", "b", "b"))

  no_start <- gamlss.longitudinal:::.starting_margin_parameter_values(
    margin_dist = gamlss.dist::NO(),
    finite_response = dat$response,
    dataset = dat
  )
  po_start <- gamlss.longitudinal:::.starting_margin_parameter_values(
    margin_dist = gamlss.dist::PO(),
    finite_response = dat$response,
    dataset = dat
  )

  expect_equal(names(no_start$margin_par), c("mu", "sigma"))
  expect_equal(unname(no_start$margin_par), c(mean(dat$response), stats::sd(dat$response)))
  expect_false(no_start$margin_par_already_eta)
  expect_equal(names(po_start$margin_par), "mu")
  expect_equal(unname(po_start$margin_par), mean(dat$response))
})

test_that("initial parameter state preserves supplied starts and warm smooths", {
  B <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3)
  attr(B, "penalty") <- diag(2)
  mm <- list(
    x = list(mu = cbind(`(Intercept)` = 1, x = c(0, 1, 2))),
    s = list(mu = list(`s(x)` = B))
  )
  start_from <- c(`mu.(Intercept)` = 1, mu.x = 2)
  warm_start_par_s <- list(mu = list(`s(x)` = c(`mu.s(x).1` = 3, `mu.s(x).2` = 4)))

  state <- gamlss.longitudinal:::.gl_build_initial_parameter_state(
    start_from = start_from,
    warm_start_par_s = warm_start_par_s,
    mm = mm,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    dataset = data.frame(response = c(1, 2, 3), time = c(1, 2, 3), subject = c("a", "a", "a")),
    lambda_start = NA
  )

  expect_equal(state$par_cov, start_from)
  expect_equal(unname(state$par_s$mu$`s(x)`), c(3, 4))
  expect_equal(names(state$par_s$mu$`s(x)`), c("mu.s(x).1", "mu.s(x).2"))
  expect_equal(unname(state$df_s$mu$`s(x)`), 0)
  expect_equal(unname(state$lambda_s$mu$`s(x)`), sum(diag(t(B) %*% B)) / 2)
})
