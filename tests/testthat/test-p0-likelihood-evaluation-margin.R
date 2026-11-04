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

test_that("non-finite family scores fall back to numerical log-density derivatives", {
  skip_if_not_installed("gamlss.dist")

  margin_dist <- gamlss.dist::NO()
  margin_dist$dldm <- function(y, mu, sigma) rep(NA_real_, length(y))
  response <- c(-1, 0.5, 2)
  mu <- c(0, 0.1, 0.2)
  sigma <- c(1, 1.2, 0.8)
  out <- gamlss.longitudinal:::.gl_likelihood_evaluate_margins(
    eta_inv = list(mu = mu, sigma = sigma),
    mm = list(mu = matrix(1, 3, 1), sigma = matrix(1, 3, 1)),
    margin_dist = margin_dist,
    margin_eval_cache = gamlss.longitudinal:::.build_margin_eval_cache(
      margin_dist, calc_d2 = FALSE
    ),
    response = response,
    discrete_margin = FALSE,
    calc_margin_deriv = TRUE,
    margin_deriv_names = "dldm"
  )

  expect_equal(
    out$margin_deriv$dldm,
    (response - mu) / sigma^2,
    tolerance = 1e-8
  )
})

test_that("known non-vector-safe GAMLSS kernels are evaluated row by row", {
  skip_if_not_installed("gamlss.dist")

  cases <- list(
    ZANBI = list(mu = c(1, 2, 3), sigma = c(0.5, 0.8, 1.1), nu = c(0.2, 0.3, 0.4)),
    ZAP = list(mu = c(1, 2, 3), sigma = c(0.1, 0.2, 0.3)),
    ZAPIG = list(mu = c(1, 2, 3), sigma = c(0.5, 0.8, 1.1), nu = c(0.2, 0.3, 0.4)),
    ZABI = list(mu = c(0.2, 0.4, 0.6), sigma = c(0.1, 0.2, 0.3))
  )
  response <- c(0, 1, 3)

  for (family in names(cases)) {
    object <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())
    parameters <- cases[[family]]
    fixed <- if (identical(family, "ZABI")) list(bd = c(10, 10, 10)) else list()
    attr(object, "fixed_family_args") <- fixed
    mm <- lapply(parameters, function(value) matrix(1, length(value), 1L))
    out <- gamlss.longitudinal:::.gl_likelihood_evaluate_margins(
      eta_inv = parameters,
      mm = mm,
      margin_dist = object,
      margin_eval_cache = gamlss.longitudinal:::.build_margin_eval_cache(
        object, calc_d2 = FALSE
      ),
      response = response,
      discrete_margin = TRUE,
      calc_margin_deriv = FALSE
    )
    dfun <- get(paste0("d", family), envir = asNamespace("gamlss.dist"))
    pfun <- get(paste0("p", family), envir = asNamespace("gamlss.dist"))
    scalar_density <- vapply(seq_along(response), function(ii) {
      args <- c(list(x = response[[ii]]), lapply(parameters, `[`, ii), lapply(fixed, `[`, ii))
      do.call(dfun, args)
    }, numeric(1))
    scalar_probability <- vapply(seq_along(response), function(ii) {
      args <- c(list(q = response[[ii]]), lapply(parameters, `[`, ii), lapply(fixed, `[`, ii))
      do.call(pfun, args)
    }, numeric(1))

    expect_equal(out$margin_d, scalar_density, tolerance = 1e-12, info = family)
    expect_equal(out$margin_p, scalar_probability, tolerance = 1e-12, info = family)
  }
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
  expect_equal(
    out$margin_log_survival,
    gamlss.dist::pDEL(
      response, mu = 3, sigma = 0.8, nu = 0.5,
      lower.tail = FALSE, log.p = TRUE
    ),
    tolerance = 1e-12
  )
  expect_equal(
    out$margin_log_survival_lower,
    gamlss.dist::pDEL(
      response - 1, mu = 3, sigma = 0.8, nu = 0.5,
      lower.tail = FALSE, log.p = TRUE
    ),
    tolerance = 1e-12
  )
})

test_that("discrete CDF companions retain directly evaluated extreme tails", {
  upper <- c(0.2, 1, NA_real_)
  mass <- c(0.1, 1e-30, NA_real_)
  log_mass <- c(log(0.1), log(1e-30), NA_real_)
  direct <- c(NA_real_, log(2e-30), NA_real_)

  out <- gamlss.longitudinal:::.gl_discrete_cdf_companions(
    upper, mass, log_mass, direct
  )

  expect_equal(out$lower[1:2], c(0.1, 1 - 1e-30))
  expect_equal(out$log_survival[1], log(0.8))
  expect_equal(out$log_survival[2], log(2e-30))
  expect_equal(out$log_survival_lower[1], log(0.9))
  expect_equal(out$log_survival_lower[2], log(3e-30), tolerance = 1e-12)
  expect_true(all(vapply(out, function(x) is.na(x[3]), logical(1))))
})
