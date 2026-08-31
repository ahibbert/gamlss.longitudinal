test_that("benchmark predictive distribution returns Gaussian intervals and PIT values", {
  y <- c(1, 2, 3, 4)
  fitted <- c(1.1, 1.9, 3.1, 3.9)

  pred <- .benchmark_predictive_distribution(
    y = y,
    fitted = fitted,
    family = stats::gaussian(),
    p = 0.9,
    interval_level = 0.8
  )

  expect_named(pred, c("q_p", "lower", "upper", "pit", "density"))
  expect_length(pred$q_p, length(y))
  expect_true(all(is.finite(pred$q_p)))
  expect_true(all(pred$lower < pred$upper))
  expect_true(all(pred$pit >= 0 & pred$pit <= 1))
  expect_true(all(pred$density > 0))
})

test_that("benchmark predictive distribution returns Poisson count summaries", {
  y <- c(0, 1, 2)
  fitted <- c(0.5, 1.5, 3)

  pred <- .benchmark_predictive_distribution(
    y = y,
    fitted = fitted,
    family = stats::poisson(),
    p = 0.9,
    interval_level = 0.8
  )

  expect_equal(pred$q_p, stats::qpois(0.9, lambda = fitted))
  expect_equal(pred$lower, stats::qpois(0.1, lambda = fitted))
  expect_equal(pred$upper, stats::qpois(0.9, lambda = fitted))
  expect_equal(
    pred$pit,
    stats::ppois(y, lambda = fitted) - 0.5 * stats::dpois(y, lambda = fitted)
  )
  expect_equal(pred$density, stats::dpois(y, lambda = fitted))
})

test_that("benchmark predictive distribution returns empty metrics for unsupported data", {
  pred <- .benchmark_predictive_distribution(
    y = c(NA_real_, Inf),
    fitted = c(1, 2),
    family = stats::gaussian()
  )

  expect_true(all(is.na(unlist(pred, use.names = FALSE))))
})
