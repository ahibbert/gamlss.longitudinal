make_conditional_simulation_setup <- function(
    n,
    family = "NO",
    params,
    history_u_lower,
    history_u_upper = history_u_lower,
    theta = 0.7,
    zeta = 0) {
  keys <- as.character(seq_len(n))
  list(
    future_diag = list(family = family, params = params),
    future_fit = data.frame(
      .time_idx = rep(2L, n),
      .future_row = seq_len(n),
      .subject_key = keys,
      theta_fit = rep(NA_real_, n),
      zeta_fit = rep(NA_real_, n)
    ),
    history = data.frame(
      subject = keys,
      time = rep("1", n),
      time_idx = rep(1L, n),
      response = rep(NA_real_, n),
      cdf_lower = rep(history_u_lower, n),
      cdf_upper = rep(history_u_upper, n),
      theta_fit = rep(theta, n),
      zeta_fit = rep(zeta, n)
    ),
    plan = list(
      groups = stats::setNames(as.list(seq_len(n)), keys),
      history_row = stats::setNames(seq_len(n), keys),
      starts = data.frame(
        subject = keys,
        future_row = seq_len(n),
        time = rep("2", n),
        previous_time = rep("1", n),
        start_type = rep("conditioned_history", n),
        reason = rep("adjacent_to_history", n)
      )
    ),
    history_source = "test"
  )
}

test_that("conditional simulation keeps the unconditional default and restores RNG", {
  dat <- make_fixture_numeric_time(n_subject = 10L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE, theta_formula = "~ time_raw")

  expect_equal(
    simulate(fit, nsim = 2L, seed = 922L),
    simulate(fit, nsim = 2L, seed = 922L, condition_on = NULL)
  )

  set.seed(823L)
  before <- .Random.seed
  invisible(simulate(fit, nsim = 1L, seed = 923L))
  expect_identical(.Random.seed, before)
})

test_that("conditional forward simulation records adjacent, new, and gap starts", {
  dat <- make_fixture_numeric_time(n_subject = 12L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE, theta_formula = "~ time_raw")

  template <- dat[match(c(1L, 2L), dat$id), , drop = FALSE]
  future <- rbind(
    transform(template[1L, , drop = FALSE], time_raw = 4L),
    transform(template[1L, , drop = FALSE], time_raw = 6L),
    transform(template[2L, , drop = FALSE], time_raw = 5L),
    transform(template[1L, , drop = FALSE], time_raw = 4L),
    transform(template[1L, , drop = FALSE], time_raw = 5L)
  )
  future$id[4:5] <- 1001L
  future$y <- NA_real_

  expect_warning(
    simulate(fit, nsim = 2L, seed = 924L, newdata = future, condition_on = "fitted"),
    class = "gamlss.longitudinal_conditional_gap_warning"
  )
  out <- suppressWarnings(
    simulate(fit, nsim = 2L, seed = 924L, newdata = future, condition_on = "fitted")
  )
  metadata <- attr(out, "conditioning")

  expect_equal(nrow(out), nrow(future))
  expect_s3_class(metadata, "gamlss_longitudinal_conditioning_metadata")
  expect_equal(unname(metadata$counts[["conditioned_starts"]]), 1L)
  expect_equal(unname(metadata$counts[["unconditional_new_subject_starts"]]), 1L)
  expect_equal(unname(metadata$counts[["segmented_gap_starts"]]), 2L)
  expect_true(any(metadata$starts$subject == "1" & metadata$starts$start_type == "unconditional_gap"))
})

test_that("custom histories and known factor levels are supported", {
  dat <- make_fixture_factor_time_interaction(n_subject = 12L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)
  lev <- levels(dat$time_raw)
  history <- dat[dat$id %in% 1:3 & dat$time_raw == lev[1L], , drop = FALSE]
  future <- history
  future$time_raw <- factor(lev[2L], levels = lev, ordered = TRUE)
  future$y <- NA_real_

  out <- simulate(fit, nsim = 2L, seed = 925L, newdata = future, condition_on = history)
  expect_equal(nrow(out), nrow(future))
  expect_equal(attr(out, "conditioning")$history_source, "data.frame")

  bad <- future
  bad$gender <- factor("unseen", levels = c(levels(dat$gender), "unseen"))
  expect_error(
    simulate(fit, nsim = 1L, newdata = bad, condition_on = history),
    "contains level(s) not seen during fitting",
    fixed = TRUE
  )
})

test_that("continuous conditional simulation agrees with the Gaussian h-function", {
  n <- 2500L
  u1 <- 0.8
  rho <- 0.65
  setup <- make_conditional_simulation_setup(
    n = n,
    params = list(mu = rep(0, n), sigma = rep(1, n)),
    history_u_lower = u1,
    theta = rho
  )
  actual <- gamlss.longitudinal:::.gl_simulate_conditional_copula_matrix(
    object = list(copula_dist = "N"), setup = setup, nsim = 1L
  )$values[, 1L]

  expect_equal(mean(actual), rho * stats::qnorm(u1), tolerance = 0.06)
})

test_that("discrete conditional simulation agrees with exact copula rectangles", {
  n <- 4000L
  mu <- 2.1
  y1 <- 2L
  y2 <- 3L
  rho <- 0.55
  lower1 <- gamlss.dist::pPO(y1 - 1L, mu = mu)
  upper1 <- gamlss.dist::pPO(y1, mu = mu)
  setup <- make_conditional_simulation_setup(
    n = n,
    family = "PO",
    params = list(mu = rep(mu, n)),
    history_u_lower = lower1,
    history_u_upper = upper1,
    theta = rho
  )
  actual <- gamlss.longitudinal:::.gl_simulate_conditional_copula_matrix(
    object = list(copula_dist = "N"), setup = setup, nsim = 1L
  )$values[, 1L]
  lower2 <- gamlss.dist::pPO(y2 - 1L, mu = mu)
  upper2 <- gamlss.dist::pPO(y2, mu = mu)
  expected <- gamlss.longitudinal:::.copula_rectangle_prob(
    upper1, upper2, lower1, lower2, family = "N", par = rho,
    floor_probability = FALSE
  ) / (upper1 - lower1)

  expect_equal(mean(actual == y2), expected, tolerance = 0.035)
})

test_that("all supported copulas can start an adjacent conditional segment", {
  n <- 8L
  for (copula in c("N", "C", "F", "G", "J", "t")) {
    spec <- get_copula_dist(copula)
    theta <- .copula_tau_to_par(copula, 0.25)
    setup <- make_conditional_simulation_setup(
      n = n,
      params = list(mu = rep(0, n), sigma = rep(1, n)),
      history_u_lower = 0.6,
      theta = theta,
      zeta = if (identical(copula, "t")) 6 else 0
    )
    out <- gamlss.longitudinal:::.gl_simulate_conditional_copula_matrix(
      object = list(copula_dist = copula), setup = setup, nsim = 1L
    )$values
    expect_true(all(is.finite(out)), info = spec$name)
  }
})
