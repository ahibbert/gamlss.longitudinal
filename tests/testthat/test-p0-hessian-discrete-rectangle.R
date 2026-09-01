test_that("discrete endpoint evaluation supports count, zero-inflated, and binary margins", {
  skip_if_not_installed("gamlss.dist")

  cases <- list(
    DEL = list(
      family = gamlss.dist::DEL(),
      eta = list(mu = c(3, 4), sigma = c(0.7, 0.8), nu = c(0.45, 0.5)),
      response = c(0, 2)
    ),
    NBI = list(
      family = gamlss.dist::NBI(),
      eta = list(mu = c(3, 4), sigma = c(0.6, 0.7)),
      response = c(0, 2)
    ),
    ZIP = list(
      family = gamlss.dist::ZIP(),
      eta = list(mu = c(3, 4), sigma = c(0.25, 0.35)),
      response = c(0, 2)
    ),
    ZINBI = list(
      family = gamlss.dist::ZINBI(),
      eta = list(mu = c(3, 4), sigma = c(0.6, 0.7), nu = c(0.25, 0.35)),
      response = c(0, 2)
    ),
    BI = list(
      family = gamlss.dist::BI(mu.link = "logit"),
      eta = list(mu = c(0.35, 0.65)),
      response = c(0, 1)
    )
  )

  for (case in cases) {
    pfun <- get(paste0("p", case$family$family[1]), envir = asNamespace("gamlss.dist"))
    dfun <- get(paste0("d", case$family$family[1]), envir = asNamespace("gamlss.dist"))

    zero <- gamlss.longitudinal:::.discrete_hessian_endpoint_eval(
      1L, case$eta, case$family, case$response, pfun, dfun
    )
    upper <- gamlss.longitudinal:::.discrete_hessian_endpoint_eval(
      2L, case$eta, case$family, case$response, pfun, dfun
    )

    expect_equal(zero$lower, 0)
    expect_true(all(is.finite(unlist(zero))))
    expect_true(all(is.finite(unlist(upper))))
    expect_lte(zero$lower, zero$upper)
    expect_lte(upper$lower, upper$upper)
  }
})

test_that("bounded-binomial families beyond Bernoulli are not silently marked supported", {
  skip_if_not_installed("gamlss.dist")

  expect_false(gamlss.longitudinal:::.is_discrete_margin(gamlss.dist::BI()))
  expect_false(gamlss.longitudinal:::.is_discrete_margin(gamlss.dist::BB()))
  expect_false(gamlss.longitudinal:::.is_discrete_margin(gamlss.dist::ZIBI()))
  expect_false(gamlss.longitudinal:::.is_discrete_margin(gamlss.dist::ZABB()))
})

test_that("discrete rectangle Hessian helper returns assembler-compatible shapes", {
  skip_if_not_installed("gamlss.dist")

  response <- c(0, 1, 1, 0, 0, 1)
  time <- rep(1:3, times = 2L)
  subject <- rep(1:2, each = 3L)
  n <- length(response)
  pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(response, time, subject)

  mm_x <- list(
    mu = matrix(1, n, 1, dimnames = list(NULL, "(Intercept)")),
    theta = matrix(1, length(pair_cache$row_id1), 1, dimnames = list(NULL, "(Intercept)"))
  )
  eta_inv <- list(
    mu = rep(0.45, n),
    theta = rep(0.25, length(pair_cache$row_id1))
  )
  eta_dr <- list(mu = rep(0.25, n), theta = rep(0.94, length(pair_cache$row_id1)))
  eta_d2 <- list(mu = rep(0, n), theta = rep(0, length(pair_cache$row_id1)))

  calc_lik <- gamlss.longitudinal:::calc_likelihood_minimal(
    eta_inv,
    mm = mm_x,
    margin_dist = gamlss.dist::BI(mu.link = "logit"),
    copula_dist = "N",
    response = response,
    response_margin = time,
    response_subject = subject,
    pair_cache = pair_cache
  )

  out <- gamlss.longitudinal:::.calc_discrete_rectangle_hessian_contributions(
    eta_inv = eta_inv,
    eta_dr = eta_dr,
    eta_d2 = eta_d2,
    pair_cache = pair_cache,
    margin_dist = gamlss.dist::BI(mu.link = "logit"),
    copula_dist = "N",
    response = response,
    calc_lik = calc_lik,
    mm = list(x = mm_x),
    h = 1e-4
  )

  expect_named(out$cop_d1l_margin, "mu")
  expect_named(out$cop_d2l_margin$mu, "mu")
  expect_length(out$cop_d1l_margin$mu, n)
  expect_length(out$cop_d2l_margin$mu$mu, n)
  expect_length(out$cop_d2l_theta, length(pair_cache$row_id1))
  expect_length(out$cross_pair_contribs$mu$mu, length(pair_cache$row_id1))
  expect_true(all(is.finite(unlist(out$cop_d1l_margin))))
  expect_true(all(is.finite(unlist(out$cop_d2l_margin))))
  expect_true(all(is.finite(out$cop_d2l_theta)))
})

test_that("analytical discrete vcov succeeds for a registered DEL fixture and BI fails preflight", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  set.seed(7301)
  dat_del <- expand.grid(id = seq_len(8L), time = seq_len(2L), KEEP.OUT.ATTRS = FALSE)
  dat_del <- dat_del[order(dat_del$id, dat_del$time), ]
  dat_del$x <- stats::rnorm(nrow(dat_del))
  dat_del$y <- gamlss.dist::rDEL(nrow(dat_del), mu = exp(1 + 0.1 * dat_del$x), sigma = 0.7, nu = 0.45)

  fit_del <- suppressWarnings(gamlss.longitudinal::gamlss_longitudinal(
    dataset = dat_del,
    margin_dist = gamlss.dist::DEL(),
    copula_dist = "N",
    time_var = "time",
    subject_var = "id",
    mu.formula = "y ~ x",
    sigma.formula = "~ 1",
    nu.formula = "~ 1",
    theta.formula = "~ 1",
    zeta.formula = "~ 1",
    include_dlcopdpar = TRUE,
    warm_start_joint = FALSE,
    method = "RS",
    max_outer_iter = 2L,
    max_inner_iter = 2L,
    outer_stop_crit = 1,
    inner_stop_crit = 1,
    compute_vcov = FALSE,
    verbose = 0
  ))
  vc_del <- vcov(fit_del, method = "analytical_only", progress = FALSE)
  expect_equal(vc_del$method, "analytical")
  expect_true(all(is.finite(vc_del$se$overall)))

  set.seed(7302)
  dat_bi <- expand.grid(id = seq_len(16L), time = seq_len(3L), KEEP.OUT.ATTRS = FALSE)
  dat_bi <- dat_bi[order(dat_bi$id, dat_bi$time), ]
  dat_bi$x <- stats::rnorm(nrow(dat_bi))
  dat_bi$y <- stats::rbinom(nrow(dat_bi), size = 1, prob = stats::plogis(-0.1 + 0.4 * dat_bi$x))

  expect_error(
    gamlss.longitudinal::gamlss_longitudinal(
      dataset = dat_bi,
      margin_dist = gamlss.dist::BI(mu.link = "logit"),
      copula_dist = "N",
      time_var = "time",
      subject_var = "id",
      mu.formula = "y ~ x"
    ),
    class = "gamlss_longitudinal_unsupported_margin_error"
  )
})

test_that("fit-time analytical discrete vcov is cached and reused", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  set.seed(7303)
  dat <- expand.grid(id = seq_len(12L), time = seq_len(2L), KEEP.OUT.ATTRS = FALSE)
  dat <- dat[order(dat$id, dat$time), ]
  dat$x <- stats::rnorm(nrow(dat))
  dat$y <- gamlss.dist::rNBI(nrow(dat), mu = exp(1 + 0.1 * dat$x), sigma = 0.7)

  fit <- suppressWarnings(gamlss.longitudinal::gamlss_longitudinal(
    dataset = dat,
    margin_dist = gamlss.dist::NBI(),
    copula_dist = "N",
    time_var = "time",
    subject_var = "id",
    mu.formula = "y ~ x",
    sigma.formula = "~ 1",
    theta.formula = "~ 1",
    zeta.formula = "~ 1",
    include_dlcopdpar = TRUE,
    warm_start_joint = FALSE,
    method = "RS",
    max_outer_iter = 2L,
    max_inner_iter = 2L,
    outer_stop_crit = 1,
    inner_stop_crit = 1,
    compute_vcov = TRUE,
    vcov_method = "analytical",
    verbose = 0
  ))

  expect_true(fit$vcov_meta$precomputed)
  expect_equal(fit$vcov_meta$method_used, "analytical")
  resolved <- gamlss.longitudinal:::.resolve_vcov(fit, extra_args = list(method = "analytical"))
  expect_identical(resolved, fit$vcov)
})
