copula_zero_cell <- function(u, v, family, theta, df) {
  switch(family,
    N = pbivnorm::pbivnorm(qnorm(u), qnorm(v), theta),
    C = (u^(-theta) + v^(-theta) - 1)^(-1 / theta),
    F = -log1p(expm1(-theta * u) * expm1(-theta * v) / expm1(-theta)) / theta,
    G = exp(-((-log(u))^theta + (-log(v))^theta)^(1 / theta)),
    J = 1 - ((1 - u)^theta + (1 - v)^theta -
      ((1 - u) * (1 - v))^theta)^(1 / theta),
    t = integrate(function(x) {
      scale <- sqrt((df + x^2) * (1 - theta^2) / (df + 1))
      dt(x, df) * pt((qt(v, df) - theta * x) / scale, df + 1)
    }, -Inf, qt(u, df), rel.tol = 1e-10)$value
  )
}

binary_markov_cells <- function(zero, family, theta, df = c(0, 0)) {
  rectangle <- function(edge) {
    u <- zero[edge]
    v <- zero[edge + 1L]
    q00 <- copula_zero_cell(u, v, family, theta[edge], df[edge])
    matrix(c(q00, v - q00, u - q00, 1 - u - v + q00), nrow = 2L)
  }
  q12 <- rectangle(1L)
  q23 <- rectangle(2L)
  cell <- as.matrix(expand.grid(1:2, 1:2, 1:2))
  q12[cell[, 1:2]] * q23[cell[, 2:3]] / c(zero[2], 1 - zero[2])[cell[, 2]]
}

test_that("discrete dataset simulation agrees with all three-visit rectangle cells", {
  zero <- c(.43, .57, .39)
  par <- list(N = c(.76, .71), C = c(2.9, 2.2), F = c(6, 5),
              G = c(2.6, 2.1), J = c(2.8, 2.2), t = c(.73, .68))
  n <- 16000L
  for (family in names(par)) {
    df <- if (family == "t") c(4, 7) else c(0, 0)
    oracle <- binary_markov_cells(zero, family, par[[family]], df)
    expect_equal(sum(oracle), 1, tolerance = 1e-10, info = family)
    expect_true(all(oracle > 0), info = family)
    dat <- simulate_longitudinal_dataset(
      n = n, times = c(2, 4, 7), seed = 18021L,
      margin_dist = gamlss.dist::BI(), copula_dist = family,
      margin_params = list(mu = 1 - zero, bd = 1),
      copula_params = list(theta = par[[family]], zeta = df),
      include_truth = FALSE
    )
    y <- matrix(dat$response, ncol = 3L, byrow = TRUE)
    observed <- tabulate(1 + y[, 1L] + 2 * y[, 2L] + 4 * y[, 3L], 8L) / n
    se <- sqrt(oracle * (1 - oracle) / n)
    expect_true(all(abs(observed - oracle) < 5 * se), info = family)
  }
})

discrete_simulation_object <- function(margin = gamlss.dist::PO(), rho = .7,
                                       mu = log(2), response = c(0, 1, 0)) {
  dat <- data.frame(subject = 1L, time = seq_along(response),
                    time_covariate = seq_along(response), response = response)
  formulas <- list(mu = response ~ 1, sigma = ~1, nu = ~1, tau = ~1,
                   theta = ~1, zeta = ~1)
  link <- get_copula_dist("N")$copula_link
  mm <- create_model_matrices(
    mu.formula = formulas$mu, sigma.formula = formulas$sigma,
    nu.formula = formulas$nu, tau.formula = formulas$tau,
    theta.formula = formulas$theta, zeta.formula = formulas$zeta,
    margin.family = margin, copula.family = "N", copula.link = link,
    dataset = dat, preserve_factor_levels = TRUE
  )
  pars <- c(mu.intercept = margin$mu.linkfun(mu), sigma.intercept = 0,
            theta.intercept = link$theta.linkfun(rho))
  structure(list(
    margin_dist = margin, copula_dist = "N", formulas_int = formulas,
    var_map = c(y = "response", id = "subject", visit = "time"),
    dataset = dat, response = dat$response, response_margin = dat$time,
    response_subject = dat$subject, model_matrix = mm, par = pars,
    par_s = list(mu = NA, sigma = NA, nu = NA, tau = NA, theta = NA, zeta = NA),
    convergence = list(converged = TRUE)
  ), class = "gamlss.longitudinal")
}

test_that("fitted discrete trajectories match the observed-count law in row order", {
  n <- 120000L
  fit <- discrete_simulation_object()
  q00 <- 1/4 + asin(.7)/(2*pi)
  p000 <- 2 * q00^2
  old <- 1/8 + (2 * asin(.7) + asin(.7^2))/(4*pi)
  sim <- as.matrix(simulate(fit, nsim = n, seed = 18011L))
  p_hat <- mean(colSums(sim == 0) == 3L)
  expect_lt(abs(mean(colSums(sim[1:2, ] == 0) == 2L) - q00),
            5 * sqrt(q00 * (1 - q00) / n))
  expect_lt(abs(p_hat - p000), 5 * sqrt(p000 * (1 - p000) / n))
  expect_gt(abs(p_hat - old), 5 * sqrt(old * (1 - old) / n))

  panel <- data.frame(id = 1L, visit = 1:3)
  order <- c(3L, 1L, 2L)
  set.seed(198)
  old_seed <- .Random.seed
  ordered <- simulate(fit, nsim = 30L, seed = 710, newdata = panel)
  shuffled <- simulate(fit, nsim = 30L, seed = 710, newdata = panel[order, ])
  expect_identical(as.matrix(ordered)[order, ], as.matrix(shuffled))
  expect_identical(.Random.seed, old_seed)

  independent <- simulate(discrete_simulation_object(rho = 0), nsim = 8000L,
                          seed = 847, newdata = panel)
  expect_lt(abs(mean(colSums(independent == 0) == 3L) - .125), .02)
})

test_that("discrete conditional forecasts use the latest count at every step", {
  n <- 14000L
  fit <- discrete_simulation_object(mu = 1, response = c(1, 0, NA, NA, NA))
  history <- data.frame(id = 1L, visit = 1:2, y = c(1, 0))
  future <- data.frame(id = 1L, visit = 3:4)
  prob_zero <- exp(-1)
  r00 <- copula_zero_cell(prob_zero, prob_zero, "N", .7, 0)
  step <- r00 / prob_zero
  sample <- as.matrix(simulate(fit, n, seed = 21811L, newdata = future,
                               condition_on = history))
  expect_lt(abs(mean(sample[1L, ] == 0) - step), 4 * sqrt(step * (1 - step) / n))
  expect_lt(abs(mean(colSums(sample == 0) == 2L) - step^2),
            4 * sqrt(step^2 * (1 - step^2) / n))
  history$y[1L] <- 6L
  expect_identical(
    as.matrix(simulate(fit, 20L, seed = 9L, newdata = future, condition_on = history)),
    as.matrix(simulate(fit, 20L, seed = 9L, newdata = future, condition_on = "fitted"))
  )
  gap <- suppressWarnings(simulate(fit, 6000L, seed = 5L, newdata =
    data.frame(id = 1L, visit = c(3L, 5L)), condition_on = history))
  expect_equal(attr(gap, "conditioning")$starts$start_type,
               c("conditioned_history", "unconditional_gap"))
  expect_lt(abs(mean(as.matrix(gap)[2, ] == 0) - prob_zero), .03)
  fitted <- simulate(discrete_simulation_object(response = c(0, NA, 2)),
                     5L, seed = 8L)
  expect_true(all(is.na(fitted[2, ])))
})

test_that("discrete truth and redraw use row-specific family arguments", {
  evaluations <- 0L
  margins <- list(
    mu = function(data) {
      evaluations <<- evaluations + 1L
      c(.27, .44, .77)
    },
    bd = c(1L, 3L, 2L)
  )
  simulate_bi <- function(bounds = NULL) simulate_longitudinal_dataset(
    n = 30L, margin_dist = gamlss.dist::BI(), copula_dist = "C",
    margin_params = margins, copula_params = list(theta = c(1.3, 2.2)),
    u_bounds = bounds, seed = 1821L
  )
  raw <- simulate_bi()
  expect_identical(evaluations, 1L)
  bounded <- simulate_bi(c(.2, .8))
  expect_identical(evaluations, 2L)
  expect_identical(bounded$u, pmin(pmax(raw$u, .2), .8))
  expect_equal(bounded$response, qbinom(bounded$u, bounded$true_bd, bounded$true_mu))
  expect_true(all(raw$response >= 0 & raw$response <= raw$true_bd))

  margin <- .gl_attach_margin_args_spec(gamlss.dist::BI(), list(bd = "trials"))
  fit <- discrete_simulation_object(margin = margin, mu = .5)
  fit$dataset$trials <- 1:3
  fit$margin_dist <- .gl_bind_margin_args(margin, fit$dataset)
  panel <- data.frame(id = 1L, visit = 1:3, trials = 1:3)
  sim <- as.matrix(simulate(fit, 120L, seed = 5L, newdata = panel))
  expect_true(all(sim >= 0 & sim <= 1:3))
  history <- transform(panel[1L, ], y = 1L)
  forecast <- as.matrix(simulate(fit, 100L, seed = 5L, newdata = panel[2:3, ],
                                 condition_on = history))
  expect_true(all(forecast >= 0 & forecast <= 2:3))
  history$y <- 2L
  expect_error(simulate(fit, 1L, seed = 5L, newdata = panel[2:3, ],
                        condition_on = history), class = "gamlss_longitudinal_response_domain_error")
  expect_error(.gl_simulation_redraw_discrete(c(0, 1), c(.2, 1), "BI", c(5, 6)),
               "interior BI CDF interval.*row\\(s\\) 6")
})

test_that("parametric bootstrap passes simulated counts to its refit path", {
  fit <- discrete_simulation_object(response = c(0, NA, 1))
  calls <- list()
  testthat::local_mocked_bindings(
    .gl_bootstrap_refit = function(object, response_b, fit_args) {
      calls[[length(calls) + 1L]] <<- response_b
      object
    }
  )
  expected <- simulate(fit, 3L, seed = 519)
  ans <- bootstrap_inference(fit, R = 3L, seed = 519, type = "parametric")
  expect_identical(calls, unname(as.list(expected)))
  expect_equal(ans$successful_replicates, 3L)
})

test_that("rare discrete conditioning does not clip away outcome categories", {
  rho <- .9
  fit <- discrete_simulation_object(rho = rho, mu = 20)
  history <- data.frame(id = 1L, visit = 1L, y = 0L)
  future <- data.frame(id = 1L, visit = 2L)
  p0 <- exp(-20)
  expected <- integrate(function(v) {
    pnorm((qnorm(p0) - rho * qnorm(p0 * v)) / sqrt(1 - rho^2))
  }, 0, 1, rel.tol = 1e-9)$value
  n <- 3500L
  forecast <- as.matrix(simulate(fit, n, seed = 661L, newdata = future,
                                 condition_on = history))
  expect_gt(expected, .01)
  expect_lt(abs(mean(forecast == 0) - expected), 5 * sqrt(expected * (1 - expected) / n))

  u <- mean(ppois(c(14, 15), lambda = 1))
  next_u <- .sim_invert_hfunc1(u, .5, "N", .9, discrete = TRUE)
  expect_equal(qpois(next_u, lambda = 1), 13)
})

test_that("native Gumbel and Joe limits resolve an unclipped lower-tail inverse", {
  skip_if(!identical(.copula_backend(), "native"))
  u <- .2
  theta <- 2
  v <- 1e-10
  # Analytic derivative of each copula CDF with respect to its first argument.
  x <- -log(u)
  w <- sqrt(x^2 + (-log(v))^2)
  gumbel <- exp(-w) * x / (u * w)
  joe <- (1 - u) * (2 * v - v^2) /
    sqrt((1 - u)^2 + (1 - v)^2 * (2 * u - u^2))
  g <- .sim_invert_hfunc1(u, gumbel, "G", theta, discrete = TRUE)
  j <- .sim_invert_hfunc1(u, joe, "J", theta, discrete = TRUE)
  expect_lt(abs(g / v - 1), 1e-7)
  expect_lt(abs(j / v - 1), 1e-7)
})
