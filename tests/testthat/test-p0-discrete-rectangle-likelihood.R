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

test_that("Gaussian independence keeps extreme discrete likelihoods finite", {
  skip_if_not_installed("gamlss.dist")
  suppressPackageStartupMessages(library(gamlss.dist))

  response <- c(50, 50)
  time <- 1:2
  subject <- c(1, 1)
  mm <- list(
    mu = data.frame(intercept = c(1, 1)),
    theta = data.frame(intercept = 1)
  )
  pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(response, time, subject)

  lik <- gamlss.longitudinal:::calc_likelihood_minimal(
    eta_inv = list(mu = c(1, 1), theta = 0),
    mm = mm,
    margin_dist = gamlss.dist::PO(),
    copula_dist = "N",
    response = response,
    response_margin = time,
    response_subject = subject,
    pair_cache = pair_cache
  )

  oracle <- sum(stats::dpois(response, lambda = 1, log = TRUE))
  expect_true(is.finite(lik$log_lik[["joint"]]))
  expect_equal(unname(lik$log_lik[["joint"]]), oracle, tolerance = 1e-12)
  expect_equal(lik$copula_log_d, 0)
  expect_equal(lik$copula_d, 1)
  expect_identical(unname(lik$contribution_counts[["marginal_included"]]), 2L)
  expect_identical(unname(lik$contribution_counts[["pair_included"]]), 1L)
  expect_true(lik$valid)
})

test_that("dependent Gaussian tail rectangles match a direct bivariate-normal oracle", {
  skip_if_not_installed("gamlss.dist")
  skip_if_not_installed("mvtnorm")
  suppressPackageStartupMessages(library(gamlss.dist))

  response <- c(15, 15)
  rho <- 0.4
  pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(
    response, 1:2, c(1, 1)
  )
  lik <- gamlss.longitudinal:::calc_likelihood_minimal(
    eta_inv = list(mu = c(1, 1), theta = rho),
    mm = list(
      mu = data.frame(intercept = c(1, 1)),
      theta = data.frame(intercept = 1)
    ),
    margin_dist = gamlss.dist::PO(),
    copula_dist = "N",
    response = response,
    response_margin = 1:2,
    response_subject = c(1, 1),
    pair_cache = pair_cache
  )

  log_survival_upper <- stats::ppois(15, 1, lower.tail = FALSE, log.p = TRUE)
  log_survival_lower <- stats::ppois(14, 1, lower.tail = FALSE, log.p = TRUE)
  z_upper <- stats::qnorm(log_survival_upper, lower.tail = FALSE, log.p = TRUE)
  z_lower <- stats::qnorm(log_survival_lower, lower.tail = FALSE, log.p = TRUE)
  oracle_mass <- as.numeric(mvtnorm::pmvnorm(
    lower = c(z_lower, z_lower),
    upper = c(z_upper, z_upper),
    corr = matrix(c(1, rho, rho, 1), 2, 2)
  ))

  expect_true(is.finite(lik$log_lik[["joint"]]))
  expect_equal(unname(lik$log_lik[["joint"]]), log(oracle_mass), tolerance = 1e-5)
  expect_true(lik$valid)
})

test_that("vectorized stable rectangle quadrature matches a scalar reference", {
  scalar_reference <- function(family, par, par2, upper, lower,
                               log_mass1, log_mass2,
                               log_survival1 = NULL, log_survival2 = NULL) {
    nodes <- (c(
      -0.9602898564975363, -0.7966664774136267,
      -0.5255324099163290, -0.1834346424956498,
       0.1834346424956498,  0.5255324099163290,
       0.7966664774136267,  0.9602898564975363
    ) + 1) / 2
    weights <- c(
      0.1012285362903763, 0.2223810344533745,
      0.3137066458778873, 0.3626837833783620,
      0.3626837833783620, 0.3137066458778873,
      0.2223810344533745, 0.1012285362903763
    ) / 2
    grid <- expand.grid(i = seq_along(nodes), j = seq_along(nodes))
    log_weights <- log(weights[grid$i]) + log(weights[grid$j])
    logspace_add <- function(a, b) {
      maximum <- pmax(a, b)
      maximum + log(exp(a - maximum) + exp(b - maximum))
    }
    node_probability <- function(k, node, log_survival, log_mass, lo, hi) {
      if (!is.null(log_survival) && is.finite(log_survival[[k]]) &&
          is.finite(log_mass[[k]])) {
        log_s <- logspace_add(log_survival[[k]], log1p(-node) + log_mass[[k]])
        return(list(u = -expm1(log_s), log_survival = log_s))
      }
      list(u = lo + node * (hi - lo), log_survival = NA_real_)
    }
    vapply(seq_along(par), function(k) {
      log_density <- vapply(seq_len(nrow(grid)), function(g) {
        p1 <- node_probability(
          k, nodes[grid$i[g]], log_survival1, log_mass1,
          lower[k, 1], upper[k, 1]
        )
        p2 <- node_probability(
          k, nodes[grid$j[g]], log_survival2, log_mass2,
          lower[k, 2], upper[k, 2]
        )
        if (identical(family, "N") && is.finite(p1$log_survival) &&
            is.finite(p2$log_survival)) {
          z1 <- stats::qnorm(p1$log_survival, lower.tail = FALSE, log.p = TRUE)
          z2 <- stats::qnorm(p2$log_survival, lower.tail = FALSE, log.p = TRUE)
          rho <- .copula_gaussian_rho(par[[k]])
          denominator <- 1 - rho^2
          return((2 * rho * z1 * z2 - rho^2 * (z1^2 + z2^2)) /
            (2 * denominator) - 0.5 * log(denominator))
        }
        .copula_logpdf(p1$u, p2$u, family, par[[k]], par2[[k]])
      }, numeric(1))
      maximum <- max(log_weights + log_density)
      maximum + log(sum(exp(log_weights + log_density - maximum)))
    }, numeric(1))
  }

  lower <- cbind(c(0.2, 0.9, 0.999999), c(0.3, 0.85, 0.99999))
  upper <- cbind(c(0.4, 0.95, 0.9999999), c(0.55, 0.93, 0.999999))
  log_mass1 <- log(upper[, 1] - lower[, 1])
  log_mass2 <- log(upper[, 2] - lower[, 2])
  log_survival1 <- log1p(-upper[, 1])
  log_survival2 <- log1p(-upper[, 2])

  for (family in c("N", "C")) {
    par <- if (identical(family, "N")) c(0.2, 0.6, 0.8) else c(0.5, 1.5, 3)
    expected <- scalar_reference(
      family, par, rep(0, 3), upper, lower,
      log_mass1, log_mass2,
      if (identical(family, "N")) log_survival1 else NULL,
      if (identical(family, "N")) log_survival2 else NULL
    )
    actual <- .gl_discrete_copula_log_ratio_quadrature(
      family, par, rep(0, 3), upper, lower,
      log_mass1, log_mass2,
      if (identical(family, "N")) log_survival1 else NULL,
      if (identical(family, "N")) log_survival2 else NULL
    )
    expect_equal(actual, expected, tolerance = 1e-12)
  }
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
    suppressWarnings(gamlss.longitudinal::gamlss_longitudinal(
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

test_that("discrete likelihood evaluates observed rows when responses contain gaps", {
  response <- c(2, NA, 4, 3, 1, 2, NA, 5)
  time <- rep(1:4, 2)
  subject <- rep(1:2, each = 4)
  eta_inv <- list(
    mu = rep(3, length(response)),
    sigma = rep(0.7, length(response)),
    nu = rep(0.4, length(response)),
    theta = rep(0.25, length(response))
  )
  mm <- list(mu = NULL, sigma = NULL, nu = NULL, theta = NULL)

  out <- gamlss.longitudinal:::calc_likelihood_minimal(
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = gamlss.dist::DEL(),
    copula_dist = "N",
    response = response,
    response_margin = time,
    response_subject = subject
  )

  expect_true(isTRUE(out$valid))
  expect_true(all(is.finite(out$margin_log_d[!is.na(response)])))
  expect_true(all(is.na(out$margin_log_d[is.na(response)])))
  expect_true(is.finite(unname(out$log_lik["joint"])))
  expect_identical(out$contribution_counts[["marginal_included"]], 6L)
  expect_identical(out$contribution_counts[["pair_included"]], 2L)
})
