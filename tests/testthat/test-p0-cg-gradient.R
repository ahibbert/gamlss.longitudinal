test_that("CG analytical gradient matches finite differences on a tiny model", {
  skip_if_not_installed("gamlss.dist")
  suppressPackageStartupMessages({
    library(gamlss.dist)
  })

  set.seed(1)
  n_subject <- 8L
  n_time <- 3L
  dat <- expand.grid(id = seq_len(n_subject), time = seq_len(n_time))
  dat <- dat[order(dat$id, dat$time), ]
  dat$x <- stats::rnorm(nrow(dat))
  dat$response <- 1 + 0.2 * dat$x + 0.1 * dat$time + stats::rnorm(nrow(dat), sd = 0.4)

  fit <- withCallingHandlers(
    gamlss_longitudinal(
      dat,
      gamlss.dist::NO(),
      "N",
      time_var = "time",
      subject_var = "id",
      mu.formula = "response ~ x + time",
      sigma.formula = "~ 1",
      theta.formula = "~ 1",
      zeta.formula = "~ 1",
      method = "CG",
      max_outer_iter = 2,
      compute_vcov = FALSE,
      include_dlcopdpar = TRUE,
      cg_gradient_method = "forward",
      verbose = 0
    ),
    warning = function(w) {
      if (grepl("Model stopped at max_outer_iter", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )

  mm <- fit$model_matrix
  beta <- fit$par
  copula_link <- gamlss.longitudinal:::get_copula_dist(fit$copula_dist)$copula_link
  eta <- gamlss.longitudinal:::calc_eta(beta, mm, fit$margin_dist, copula_link, fit$par_s)
  pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(fit$response, fit$response_margin, fit$response_subject)
  lik <- gamlss.longitudinal:::calc_likelihood_minimal(
    eta$eta_inv,
    mm = mm$x,
    margin_dist = fit$margin_dist,
    copula_dist = fit$copula_dist,
    calc_d2 = FALSE,
    response = fit$response,
    response_margin = fit$response_margin,
    response_subject = fit$response_subject,
    pair_cache = pair_cache
  )

  analytical <- gamlss.longitudinal:::.cg_analytical_gradient(
    beta,
    mm,
    eta,
    lik,
    fit$margin_dist,
    fit$copula_dist,
    TRUE,
    fit$response,
    fit$response_margin,
    fit$response_subject
  )

  finite_diff <- beta * 0
  for(i in seq_along(beta)) {
    h <- 1e-6 * max(1, abs(beta[i]))
    beta_plus <- beta_minus <- beta
    beta_plus[i] <- beta_plus[i] + h
    beta_minus[i] <- beta_minus[i] - h

    eta_plus <- gamlss.longitudinal:::calc_eta(beta_plus, mm, fit$margin_dist, copula_link, fit$par_s)
    eta_minus <- gamlss.longitudinal:::calc_eta(beta_minus, mm, fit$margin_dist, copula_link, fit$par_s)

    ll_plus <- gamlss.longitudinal:::calc_likelihood_minimal(
      eta_plus$eta_inv,
      mm = mm$x,
      margin_dist = fit$margin_dist,
      copula_dist = fit$copula_dist,
      calc_d2 = FALSE,
      response = fit$response,
      response_margin = fit$response_margin,
      response_subject = fit$response_subject,
      pair_cache = pair_cache
    )$log_lik["joint"]
    ll_minus <- gamlss.longitudinal:::calc_likelihood_minimal(
      eta_minus$eta_inv,
      mm = mm$x,
      margin_dist = fit$margin_dist,
      copula_dist = fit$copula_dist,
      calc_d2 = FALSE,
      response = fit$response,
      response_margin = fit$response_margin,
      response_subject = fit$response_subject,
      pair_cache = pair_cache
    )$log_lik["joint"]

    finite_diff[i] <- (ll_plus - ll_minus) / (2 * h)
  }

  expect_equal(unname(analytical), unname(finite_diff), tolerance = 1e-4)
})

test_that("copula parameter gradient maps pair derivatives onto beta coefficients", {
  mm_cg <- list(
    x = list(
      theta = cbind(theta = rep(1, 4), x = c(0, 1, 0, 1))
    )
  )
  eta <- list(theta = rep(0, 4))
  eta_dr <- list(theta = rep(2, 4))
  calc_lik <- list(
    copula_row_id1 = c(1L, 2L, 2L, 4L),
    copula_theta_index_map = seq_len(4)
  )

  grad <- gamlss.longitudinal:::.cg_copula_parameter_gradient(
    "theta",
    derivative = c(1, 2, 3, 4),
    eta = eta,
    eta_dr = eta_dr,
    mm_cg = mm_cg,
    calc_lik = calc_lik,
    response = rep(0, 4)
  )

  expect_equal(grad, c(theta.theta = 20, theta.x = 18))
})

test_that("CG margin natural scores use margin derivatives and zero missing pieces", {
  margin_par <- c("mu", "sigma", "nu")
  eta <- list(
    mu = rep(0, 3),
    sigma = rep(0, 3),
    nu = rep(0, 3)
  )
  calc_lik <- list(
    margin_deriv = list(
      dldm = c(1, 2, 3),
      dldd.extra = c(4, 5, 6)
    )
  )

  scores <- gamlss.longitudinal:::.cg_margin_natural_scores(
    margin_par = margin_par,
    eta = eta,
    eta_inv = list(mu = rep(0, 3)),
    mm_cg = list(x = list()),
    calc_lik = calc_lik,
    margin_dist = NULL,
    copula_derivatives = list(),
    include_dlcopdpar = FALSE,
    response = rep(0, 3),
    response_margin = rep(1, 3),
    response_subject = seq_len(3)
  )

  expect_equal(scores$mu, c(1, 2, 3))
  expect_equal(scores$sigma, c(4, 5, 6))
  expect_equal(scores$nu, c(0, 0, 0))
})
