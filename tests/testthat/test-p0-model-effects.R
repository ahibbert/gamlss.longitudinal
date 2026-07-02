test_that("marginal-effects default values follow column type policy", {
  expect_equal(
    .gl_effect_counterfactual_values(factor(c("b", "a"), levels = c("a", "b", "c"))),
    c("a", "b", "c")
  )
  expect_equal(
    .gl_effect_counterfactual_values(c("z", "a", "z")),
    c("a", "z")
  )
  expect_equal(
    .gl_effect_counterfactual_values(1:4),
    as.numeric(stats::quantile(1:4, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE))
  )
})

test_that("marginal-effects calibration rows add missing factor levels", {
  nd <- data.frame(
    arm = factor(c("active", "active"), levels = c("control", "active")),
    visit = ordered(c("week1", "week1"), levels = c("week1", "week2")),
    x = c(1, 2)
  )
  template <- data.frame(
    arm = factor("active", levels = c("control", "active")),
    visit = ordered("week1", levels = c("week1", "week2")),
    x = 1
  )

  out <- .gl_effect_add_factor_calibration_rows(nd, template)

  expect_equal(nrow(out), 4)
  expect_true(all(levels(out$arm) == c("control", "active")))
  expect_true(is.ordered(out$visit))
  expect_true("control" %in% as.character(out$arm))
  expect_true("week2" %in% as.character(out$visit))
})

test_that("marginal-effects calibration rows are unchanged when levels are present", {
  nd <- data.frame(group = factor(c("a", "b"), levels = c("a", "b")))

  out <- .gl_effect_add_factor_calibration_rows(nd, nd)

  expect_identical(out, nd)
})

test_that("marginal-effects row helper averages requested parameter predictions", {
  nd <- data.frame(
    arm = factor(c("control", "active"), levels = c("control", "active")),
    x = c(1, 2)
  )
  captured <- list()

  out <- gamlss.longitudinal:::.gl_effect_counterfactual_row(
    object = structure(list(), class = "gamlss.longitudinal"),
    newdata = nd,
    variable = "arm",
    value = "active",
    parameter = "sigma",
    se.fit = FALSE,
    vcov_method = "analytical",
    predict_fn = function(object, newdata, type, ...) {
      captured$newdata <<- newdata
      captured$type <<- type
      data.frame(subject = 1:2, time = 1:2, response = c(NA, NA), mu = c(1, 2), sigma = c(3, 5))
    }
  )

  expect_equal(captured$type, "parameters")
  expect_equal(as.character(captured$newdata$arm[seq_len(nrow(nd))]), c("active", "active"))
  expect_equal(out$variable, "arm")
  expect_equal(out$value, "active")
  expect_equal(out$parameter, "sigma")
  expect_equal(out$estimate, 4)
  expect_true(is.na(out$std_error))
})

test_that("marginal-effects row helper summarizes mu standard errors", {
  nd <- data.frame(x = c(1, 2, 3))

  out <- gamlss.longitudinal:::.gl_effect_counterfactual_row(
    object = structure(list(), class = "gamlss.longitudinal"),
    newdata = nd,
    variable = "x",
    value = 10,
    parameter = "mu",
    se.fit = TRUE,
    vcov_method = "numderiv",
    predict_fn = function(object, newdata, type, se.fit, interval, vcov_method, ...) {
      expect_equal(type, "response")
      expect_true(se.fit)
      expect_equal(interval, "none")
      expect_equal(vcov_method, "numderiv")
      data.frame(
        subject = seq_len(nrow(newdata)),
        time = seq_len(nrow(newdata)),
        response = NA_real_,
        fit = c(2, 4, 6),
        se.fit = c(0.3, 0.4, NA_real_)
      )
    }
  )

  expect_equal(out$value, "10")
  expect_equal(out$estimate, 4)
  expect_equal(out$std_error, sqrt(0.3^2 + 0.4^2) / 2)
})

test_that("marginal-effects row helper validates requested parameter", {
  nd <- data.frame(x = c(1, 2))

  expect_error(
    gamlss.longitudinal:::.gl_effect_counterfactual_row(
      object = structure(list(), class = "gamlss.longitudinal"),
      newdata = nd,
      variable = "x",
      value = 10,
      parameter = "sigma",
      se.fit = FALSE,
      vcov_method = "analytical",
      predict_fn = function(object, newdata, type, ...) {
        data.frame(mu = c(1, 2))
      }
    ),
    "Parameter 'sigma' is not available",
    fixed = TRUE
  )
})

test_that("marginal-effects result helper adds reference contrasts", {
  rows <- list(
    data.frame(variable = "arm", value = "control", parameter = "mu", estimate = 2, std_error = NA_real_),
    data.frame(variable = "arm", value = "active", parameter = "mu", estimate = 5, std_error = NA_real_)
  )

  out <- gamlss.longitudinal:::.gl_finalize_marginal_effects(
    rows = rows,
    reference = "control",
    se.fit = FALSE,
    level = 0.95
  )

  expect_equal(out$reference, c("control", "control"))
  expect_equal(out$contrast, c(0, 3))
  expect_false(any(c("conf.low", "conf.high") %in% names(out)))
})

test_that("marginal-effects result helper adds confidence intervals when available", {
  rows <- list(
    data.frame(variable = "x", value = "1", parameter = "mu", estimate = 2, std_error = 0.5),
    data.frame(variable = "x", value = "2", parameter = "mu", estimate = 5, std_error = 1)
  )

  out <- gamlss.longitudinal:::.gl_finalize_marginal_effects(
    rows = rows,
    reference = "1",
    se.fit = TRUE,
    level = 0.8
  )

  z <- stats::qnorm(0.9)
  expect_equal(out$contrast, c(0, 3))
  expect_equal(out$conf.low, out$estimate - z * out$std_error)
  expect_equal(out$conf.high, out$estimate + z * out$std_error)
})

test_that("marginal-effects result helper validates reference values", {
  rows <- list(
    data.frame(variable = "x", value = "1", parameter = "mu", estimate = 2, std_error = NA_real_)
  )

  expect_error(
    gamlss.longitudinal:::.gl_finalize_marginal_effects(
      rows = rows,
      reference = "missing",
      se.fit = FALSE,
      level = 0.95
    ),
    "'reference' must be one of the counterfactual values.",
    fixed = TRUE
  )
})
