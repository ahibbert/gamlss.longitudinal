test_that("plot method quantile column helpers preserve naming contract", {
  expect_equal(
    gamlss.longitudinal:::.plot_method_quantile_col_name(0.1),
    "q1"
  )
  expect_equal(
    gamlss.longitudinal:::.plot_method_quantile_col_name(0.95),
    "q95"
  )

  q_cols <- gamlss.longitudinal:::.plot_method_quantile_columns(c(0.1, 0.9, 0.5))

  expect_equal(q_cols$q_low, "q1")
  expect_equal(q_cols$q_high, "q9")
  expect_equal(q_cols$q_mid, "q5")
})

test_that("plot method empty panel helper returns a ggplot", {
  p <- gamlss.longitudinal:::.plot_method_empty_plot(
    "Newdata Forecast Quantiles",
    "Provide 'newdata' or 'data' for this panel"
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot method newdata helper adds response fallback columns", {
  x <- list(var_map = c(id = "subject", visit = "time", y = "response"))
  data <- data.frame(id = 1:3, visit = 1:3)

  nd <- gamlss.longitudinal:::.plot_method_newdata_quantile_input(
    x,
    data = data,
    newdata = NULL,
    newdata_n = 2
  )

  expect_equal(nrow(nd), 2L)
  expect_true("y" %in% names(nd))
  expect_true(all(is.na(nd$y)))
})

test_that("plot method newdata helper respects explicit newdata", {
  x <- list(var_map = NULL)
  newdata <- data.frame(subject = 1:2, time = 1:2)

  nd <- gamlss.longitudinal:::.plot_method_newdata_quantile_input(
    x,
    data = data.frame(subject = 1:5, time = 1:5),
    newdata = newdata,
    newdata_n = 1
  )

  expect_identical(nd, newdata)
})
