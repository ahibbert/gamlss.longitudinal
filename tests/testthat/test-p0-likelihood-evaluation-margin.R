test_that("likelihood margin evaluation handles continuous densities and missing responses", {
  skip_if_not_installed("gamlss.dist")
  suppressPackageStartupMessages(library(gamlss.dist))

  response <- c(0, NA_real_, 1)
  margin_dist <- gamlss.dist::NO()
  margin_eval_cache <- gamlss.longitudinal:::.build_margin_eval_cache(margin_dist, calc_d2 = FALSE)
  out <- gamlss.longitudinal:::.gl_likelihood_evaluate_margins(
    eta_inv = list(mu = c(0, 0, 0), sigma = c(1, 1, 1)),
    mm = list(mu = matrix(1, 3, 1), sigma = matrix(1, 3, 1)),
    margin_dist = margin_dist,
    margin_eval_cache = margin_eval_cache,
    response = response,
    discrete_margin = FALSE,
    calc_margin_deriv = TRUE
  )

  expect_identical(out$likelihood_type, "continuous_density")
  expect_null(out$margin_p_lower)
  expect_equal(out$margin_p[c(1, 3)], stats::pnorm(c(0, 1)), tolerance = 1e-12)
  expect_equal(out$margin_d[c(1, 3)], stats::dnorm(c(0, 1)), tolerance = 1e-12)
  expect_true(is.na(out$margin_p[2]))
  expect_true(is.na(out$margin_d[2]))
  expect_true(length(out$margin_deriv) > 0)
})

test_that("likelihood margin evaluation can restrict marginal derivatives", {
  skip_if_not_installed("gamlss.dist")
  suppressPackageStartupMessages(library(gamlss.dist))

  response <- c(-1, 0, 1)
  margin_dist <- gamlss.dist::NO()
  eta_inv <- list(mu = c(0, 0.1, 0.2), sigma = c(1, 1.1, 1.2))
  mm <- list(mu = matrix(1, 3, 1), sigma = matrix(1, 3, 1))
  margin_eval_cache <- gamlss.longitudinal:::.build_margin_eval_cache(margin_dist, calc_d2 = FALSE)

  full <- gamlss.longitudinal:::.gl_likelihood_evaluate_margins(
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = margin_dist,
    margin_eval_cache = margin_eval_cache,
    response = response,
    discrete_margin = FALSE,
    calc_margin_deriv = TRUE
  )
  subset <- gamlss.longitudinal:::.gl_likelihood_evaluate_margins(
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = margin_dist,
    margin_eval_cache = margin_eval_cache,
    response = response,
    discrete_margin = FALSE,
    calc_margin_deriv = TRUE,
    margin_deriv_names = "dldm"
  )
  none <- gamlss.longitudinal:::.gl_likelihood_evaluate_margins(
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = margin_dist,
    margin_eval_cache = margin_eval_cache,
    response = response,
    discrete_margin = FALSE,
    calc_margin_deriv = TRUE,
    margin_deriv_names = character(0)
  )

  expect_identical(names(subset$margin_deriv), "dldm")
  expect_equal(subset$margin_deriv$dldm, full$margin_deriv$dldm)
  expect_equal(subset$margin_d, full$margin_d)
  expect_equal(subset$margin_p, full$margin_p)
  expect_identical(none$margin_deriv, list())
})

test_that("likelihood margin evaluation builds lower CDFs for discrete rectangles", {
  skip_if_not_installed("gamlss.dist")
  suppressPackageStartupMessages(library(gamlss.dist))

  response <- c(0, 3)
  margin_dist <- gamlss.dist::DEL()
  eta_inv <- list(
    mu = c(3, 3),
    sigma = c(0.8, 0.8),
    nu = c(0.5, 0.5)
  )
  margin_eval_cache <- gamlss.longitudinal:::.build_margin_eval_cache(margin_dist, calc_d2 = FALSE)
  out <- gamlss.longitudinal:::.gl_likelihood_evaluate_margins(
    eta_inv = eta_inv,
    mm = list(
      mu = matrix(1, 2, 1),
      sigma = matrix(1, 2, 1),
      nu = matrix(1, 2, 1)
    ),
    margin_dist = margin_dist,
    margin_eval_cache = margin_eval_cache,
    response = response,
    discrete_margin = TRUE,
    calc_margin_deriv = TRUE
  )

  expect_identical(out$likelihood_type, "discrete_rectangle")
  expect_equal(out$margin_p_lower[1], 0)
  expect_equal(
    out$margin_p_lower[2],
    gamlss.dist::pDEL(2, mu = 3, sigma = 0.8, nu = 0.5),
    tolerance = 1e-12
  )
  expect_equal(
    out$margin_d[1],
    gamlss.dist::dDEL(response[1], mu = 3, sigma = 0.8, nu = 0.5),
    tolerance = 1e-12
  )
  expect_equal(
    out$margin_d[2],
    gamlss.dist::dDEL(response[2], mu = 3, sigma = 0.8, nu = 0.5),
    tolerance = 1e-12
  )
})
