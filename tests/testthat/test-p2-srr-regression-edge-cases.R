test_that("srr regression edge cases cover noiseless predictors and response", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  dat$age_copy <- dat$age
  expect_warning(
    fit_fixture_model(
      dat,
      include_dlcopdpar = FALSE,
      mu_formula = "y ~ age + age_copy",
      sigma_formula = "~ 1",
      theta_formula = "~ 1",
      max_outer_iter = 2
    ),
    "rank deficient",
    fixed = TRUE
  )

  dat_constant <- make_fixture_factor_time_interaction(n_subject = 10L)
  dat_constant$y <- 2
  warn_out <- capture_warnings(
    gamlss.longitudinal:::get_starting_values(
      copula_dist = "N",
      margin_dist = gamlss.dist::NO(),
      dataset = data.frame(time = as.integer(dat_constant$time_raw), response = dat_constant$y),
      eta_transform = FALSE
    )
  )
  expect_true(any(grepl("Non-finite Kendall tau", warn_out$warnings, fixed = TRUE)))
})

test_that("srr accessor contracts and row/case metadata are stable", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  dat <- make_fixture_with_structural_missing_rows()
  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = FALSE,
    max_outer_iter = 2,
    theta_formula = "~ 1"
  )

  expect_true(all(c("converged", "hit_outer_limit") %in% names(fit$convergence)))
  expect_s3_class(stats::formula(fit), "formula")
  expect_s3_class(stats::terms(fit), "terms")
  expect_equal(stats::nobs(fit, type = "submitted"), nrow(dat))
  expect_equal(stats::nobs(fit, type = "expanded"), length(unique(dat$id)) * length(unique(dat$time_raw)))
  expect_equal(nrow(stats::model.frame(fit, type = "observed")), stats::nobs(fit))
  expect_equal(rownames(stats::model.frame(fit, type = "expanded")), as.character(seq_len(stats::nobs(fit, type = "expanded"))))
  expect_length(stats::fitted(fit), stats::nobs(fit, type = "expanded"))
  expect_length(stats::residuals(fit), stats::nobs(fit))
})

test_that("srr newdata policy rejects unseen factor levels and supports future panels", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = FALSE,
    max_outer_iter = 2,
    theta_formula = "~ 1"
  )

  nd_bad <- dat[seq_len(3), , drop = FALSE]
  nd_bad$gender <- factor("X", levels = c(levels(dat$gender), "X"))
  expect_error(
    stats::predict(fit, newdata = nd_bad),
    "contains level(s) not seen during fitting",
    fixed = TRUE
  )

  nd_future <- dat[seq_len(6), , drop = FALSE]
  nd_future$id <- nd_future$id + 1000L
  nd_future$y <- NA_real_
  pred <- stats::predict(fit, newdata = nd_future, type = "quantile", probs = c(0.1, 0.9))
  expect_s3_class(pred, "data.frame")
  expect_equal(nrow(pred), nrow(nd_future))
  expect_true(all(c("q01", "q09") %in% names(pred)))
  expect_true(all(is.finite(pred$q01)))
  expect_true(all(is.finite(pred$q09)))
  expect_true(all(pred$q01 <= pred$q09))

  pred_interval <- stats::predict(
    fit,
    newdata = nd_future,
    type = "mean",
    se.fit = TRUE,
    interval = "confidence"
  )
  expect_s3_class(pred_interval, "data.frame")
  expect_equal(nrow(pred_interval), nrow(nd_future))
  expect_true(all(c("fit", "se.fit", "conf.low", "conf.high") %in% names(pred_interval)))
})
