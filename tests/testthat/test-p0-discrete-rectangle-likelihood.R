test_that("DEL Clayton likelihood uses exact rectangle probabilities", {
  skip_if_not_installed("gamlss.dist")
  suppressPackageStartupMessages(library(gamlss.dist))

  set.seed(11)
  n_subject <- 500L
  true_theta <- 2
  mu <- 3
  sigma <- 0.8
  nu <- 0.5

  u1 <- stats::runif(n_subject)
  target <- stats::runif(n_subject)
  u2 <- vapply(seq_len(n_subject), function(i) {
    objective <- function(u) {
      gamlss.longitudinal:::.copula_hfunc1(u1[i], u, family = "C", par = true_theta) - target[i]
    }
    stats::uniroot(objective, c(1e-8, 1 - 1e-8), tol = 1e-10)$root
  }, numeric(1))

  y1 <- gamlss.dist::qDEL(u1, mu = mu, sigma = sigma, nu = nu)
  y2 <- gamlss.dist::qDEL(u2, mu = mu, sigma = sigma, nu = nu)
  response <- as.numeric(rbind(y1, y2))
  time <- rep(1:2, times = n_subject)
  subject <- rep(seq_len(n_subject), each = 2L)

  mm <- list(
    mu = data.frame(intercept = rep(1, length(response))),
    sigma = data.frame(intercept = rep(1, length(response))),
    nu = data.frame(intercept = rep(1, length(response))),
    theta = data.frame(intercept = rep(1, n_subject))
  )
  pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(response, time, subject)

  loglik_at <- function(theta) {
    eta_inv <- list(
      mu = rep(mu, length(response)),
      sigma = rep(sigma, length(response)),
      nu = rep(nu, length(response)),
      theta = rep(theta, n_subject)
    )
    gamlss.longitudinal:::calc_likelihood_minimal(
      eta_inv,
      mm = mm,
      margin_dist = gamlss.dist::DEL(),
      copula_dist = "C",
      response = response,
      response_margin = time,
      response_subject = subject,
      pair_cache = pair_cache
    )$log_lik["joint"]
  }

  theta_hat <- stats::optimize(function(th) -loglik_at(th), c(1e-4, 8))$minimum
  expect_equal(theta_hat, true_theta, tolerance = 0.35)

  eta_inv <- list(
    mu = rep(mu, length(response)),
    sigma = rep(sigma, length(response)),
    nu = rep(nu, length(response)),
    theta = rep(true_theta, n_subject)
  )
  lik <- gamlss.longitudinal:::calc_likelihood_minimal(
    eta_inv,
    mm = mm,
    margin_dist = gamlss.dist::DEL(),
    copula_dist = "C",
    response = response,
    response_margin = time,
    response_subject = subject,
    pair_cache = pair_cache
  )

  upper1 <- gamlss.dist::pDEL(y1, mu = mu, sigma = sigma, nu = nu)
  upper2 <- gamlss.dist::pDEL(y2, mu = mu, sigma = sigma, nu = nu)
  lower1 <- gamlss.dist::pDEL(y1 - 1, mu = mu, sigma = sigma, nu = nu)
  lower2 <- gamlss.dist::pDEL(y2 - 1, mu = mu, sigma = sigma, nu = nu)
  rect <- gamlss.longitudinal:::.copula_cdf(upper1, upper2, "C", true_theta) -
    gamlss.longitudinal:::.copula_cdf(lower1, upper2, "C", true_theta) -
    gamlss.longitudinal:::.copula_cdf(upper1, lower2, "C", true_theta) +
    gamlss.longitudinal:::.copula_cdf(lower1, lower2, "C", true_theta)
  pmf_ratio <- rect / (lik$margin_d[lik$copula_row_id1] * lik$margin_d[lik$copula_row_id2])

  expect_equal(lik$copula_d, pmf_ratio, tolerance = 1e-10)
  expect_equal(as.numeric(lik$log_lik["joint"]), sum(log(pmax(rect, 1e-300))), tolerance = 1e-8)
})

test_that("discrete rectangle analytical scores match finite differences", {
  skip_if_not_installed("gamlss.dist")
  suppressPackageStartupMessages(library(gamlss.dist))

  set.seed(12)
  n_subject <- 5L
  n_time <- 3L
  response <- gamlss.dist::rDEL(n_subject * n_time, mu = 3, sigma = 0.75, nu = 0.45)
  time <- rep(seq_len(n_time), times = n_subject)
  subject <- rep(seq_len(n_subject), each = n_time)

  mm <- list(
    mu = data.frame(intercept = rep(1, length(response))),
    sigma = data.frame(intercept = rep(1, length(response))),
    nu = data.frame(intercept = rep(1, length(response))),
    theta = data.frame(intercept = rep(1, n_subject * (n_time - 1L)))
  )
  eta_inv <- list(
    mu = rep(3, length(response)),
    sigma = rep(0.75, length(response)),
    nu = rep(0.45, length(response)),
    theta = rep(1.8, n_subject * (n_time - 1L))
  )
  pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(response, time, subject)
  lik <- gamlss.longitudinal:::calc_likelihood_minimal(
    eta_inv,
    mm = mm,
    margin_dist = gamlss.dist::DEL(),
    copula_dist = "C",
    response = response,
    response_margin = time,
    response_subject = subject,
    pair_cache = pair_cache
  )

  analytical <- gamlss.longitudinal:::.calc_discrete_rectangle_scores(
    eta_inv, mm, gamlss.dist::DEL(), "C", response, time, subject,
    pair_cache = pair_cache, calc_lik = lik, method = "analytical", h = 1e-4
  )
  finite <- gamlss.longitudinal:::.calc_discrete_rectangle_scores(
    eta_inv, mm, gamlss.dist::DEL(), "C", response, time, subject,
    pair_cache = pair_cache, calc_lik = lik, method = "finite", h = 1e-4
  )

  for (nm in names(finite)) {
    expect_equal(analytical[[nm]], finite[[nm]], tolerance = 1e-2)
  }
})

test_that("RS discrete rectangle score modes both fit finite tiny DEL models", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  suppressPackageStartupMessages({
    library(gamlss)
    library(gamlss.dist)
  })

  set.seed(13)
  dat <- expand.grid(id = seq_len(6L), time = seq_len(2L), KEEP.OUT.ATTRS = FALSE)
  dat <- dat[order(dat$id, dat$time), ]
  dat$x <- stats::rnorm(nrow(dat))
  dat$y <- gamlss.dist::rDEL(nrow(dat), mu = exp(1 + 0.1 * dat$x), sigma = 0.7, nu = 0.45)

  fit_once <- function(score_method) {
    suppressWarnings(gamlss.longitudinal::gamlss.longitudinal(
      dataset = dat,
      margin_dist = gamlss.dist::DEL(),
      copula_dist = "C",
      time_var = "time",
      subject_var = "id",
      mu.formula = "y ~ x",
      sigma.formula = "~ 1",
      nu.formula = "~ 1",
      theta.formula = "~ 1",
      zeta.formula = "~ 1",
      include_dlcopdpar = TRUE,
      warm_start_joint = FALSE,
      discrete_score_method = score_method,
      method = "RS",
      max_outer_iter = 2L,
      max_inner_iter = 2L,
      outer_stop_crit = 1,
      inner_stop_crit = 1,
      compute_vcov = FALSE,
      verbose = 0
    ))
  }

  fit_analytical <- fit_once("analytical")
  fit_finite <- fit_once("finite")

  expect_identical(fit_analytical$calc_lik_out_end$likelihood_type, "discrete_rectangle")
  expect_identical(fit_finite$calc_lik_out_end$likelihood_type, "discrete_rectangle")
  expect_true(all(is.finite(fit_analytical$calc_lik_out_end$log_lik)))
  expect_true(all(is.finite(fit_finite$calc_lik_out_end$log_lik)))
})
