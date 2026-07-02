make_prediction_interval_object <- function() {
  X <- matrix(1, nrow = 3, ncol = 1)
  colnames(X) <- "mu.intercept"
  X_sigma <- matrix(1, nrow = 3, ncol = 1)
  colnames(X_sigma) <- "sigma.intercept"

  structure(
    list(
      response = c(1, 2, NA),
      response_margin = c(1, 2, 3),
      response_subject = c("a", "a", "b"),
      margin_dist = gamlss.dist::NO(),
      copula_dist = "N",
      par = c(mu.intercept = 10, sigma.intercept = 0),
      par_s = list(mu = NA, sigma = NA),
      model_matrix = list(
        x = list(mu = X, sigma = X_sigma),
        s = list(mu = list(), sigma = list())
      )
    ),
    class = "gamlss.longitudinal"
  )
}

test_that("prediction interval helper returns standard-error frame without intervals", {
  object <- make_prediction_interval_object()

  out <- gamlss.longitudinal:::.gl_prediction_interval_frame(
    object = object,
    newdata = NULL,
    fit_values = c(10, 10, 10),
    interval = "none",
    level = 0.95,
    vcov_method = "analytical",
    se_fn = function(object, newdata, method, ...) c(0.1, 0.2, 0.3)
  )

  expect_equal(names(out), c("subject", "time", "response", "fit", "se.fit"))
  expect_equal(out$fit, c(10, 10, 10))
  expect_equal(out$se.fit, c(0.1, 0.2, 0.3))
  expect_equal(out$response, c(1, 2, NA))
})

test_that("prediction interval helper adds confidence limits when requested", {
  object <- make_prediction_interval_object()

  out <- gamlss.longitudinal:::.gl_prediction_interval_frame(
    object = object,
    newdata = NULL,
    fit_values = c(10, 10, 10),
    interval = "confidence",
    level = 0.8,
    vcov_method = "numderiv",
    se_fn = function(object, newdata, method, ...) {
      expect_equal(method, "numderiv")
      c(1, 2, 3)
    }
  )

  z <- stats::qnorm(0.9)
  expect_equal(names(out), c("subject", "time", "response", "fit", "se.fit", "conf.low", "conf.high"))
  expect_equal(out$conf.low, c(10, 10, 10) - z * c(1, 2, 3))
  expect_equal(out$conf.high, c(10, 10, 10) + z * c(1, 2, 3))
})

test_that("distributional prediction helpers build quantile and cdf frames", {
  object <- make_prediction_interval_object()
  diag_data <- gamlss.longitudinal:::.gl_fitted_distribution(object, require_response = FALSE)
  pred <- gamlss.longitudinal:::.gl_prediction_frame(object, require_response = FALSE)

  quant <- gamlss.longitudinal:::.gl_prediction_quantile_frame(
    pred = pred,
    params = diag_data$params,
    family = diag_data$family,
    probs = c(0.25, 0.5)
  )
  cdf <- gamlss.longitudinal:::.gl_prediction_cdf_frame(
    pred = pred,
    diag_data = diag_data,
    q = c(10, 11, 12),
    require_response = FALSE
  )

  expect_equal(names(quant), c("subject", "time", "response", "q025", "q05"))
  expect_equal(quant$q05, c(10, 10, 10), tolerance = 1e-12)
  expect_equal(names(cdf), c("subject", "time", "response", "q", "cdf"))
  expect_equal(cdf$q, c(10, 11, 12))
  expect_true(all(cdf$cdf >= 0 & cdf$cdf <= 1))
})

test_that("distributional prediction helpers build probability and density frames", {
  object <- make_prediction_interval_object()
  diag_data <- gamlss.longitudinal:::.gl_fitted_distribution(object, require_response = FALSE)
  pred <- gamlss.longitudinal:::.gl_prediction_frame(object, require_response = FALSE)

  below <- gamlss.longitudinal:::.gl_prediction_probability_frame(
    pred = pred,
    diag_data = diag_data,
    q = 10,
    direction = "below",
    require_response = FALSE
  )
  above <- gamlss.longitudinal:::.gl_prediction_probability_frame(
    pred = pred,
    diag_data = diag_data,
    q = 10,
    direction = "above",
    require_response = FALSE
  )
  density <- gamlss.longitudinal:::.gl_prediction_density_frame(
    pred = pred,
    diag_data = diag_data,
    y = c(10, 11, 12),
    require_response = FALSE
  )

  expect_equal(names(below), c("subject", "time", "response", "q", "direction", "probability"))
  expect_equal(below$direction, rep("below", 3))
  expect_equal(above$probability, 1 - below$probability)
  expect_equal(names(density), c("subject", "time", "response", "y", "density"))
  expect_equal(density$y, c(10, 11, 12))
  expect_true(all(density$density >= 0))
  expect_error(
    gamlss.longitudinal:::.gl_prediction_probability_frame(pred, diag_data, q = NULL, direction = "below", require_response = FALSE),
    "'q' is required for type = 'probability'.",
    fixed = TRUE
  )
})
