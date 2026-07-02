test_that("select-margin input helper normalizes vector and data-frame responses", {
  counts <- .select_margin_inputs(response = c(0, 1, 2, NA_real_))
  expect_equal(counts$response, c(0, 1, 2))
  expect_equal(counts$response_all, c(0, 1, 2, NA_real_))
  expect_equal(counts$type, "counts")

  positive <- .select_margin_inputs(response = c(0.5, 1.2, 3.4))
  expect_equal(positive$type, "realplus")

  signed <- .select_margin_inputs(response = c(-1, 0, 1))
  expect_equal(signed$type, "realAll")

  dat <- data.frame(y = c(1, 2, 3), time = c("a", "b", "c"))
  from_data <- .select_margin_inputs(
    response = dat,
    response_var = "y",
    time_intercepts = TRUE,
    time_var = "time"
  )
  expect_equal(from_data$response, c(1, 2, 3))
  expect_true(from_data$time_intercepts)
  expect_equal(from_data$time_var, "time")
  expect_s3_class(from_data$data, "data.frame")
})

test_that("select-margin input helper preserves reviewer-facing errors", {
  expect_error(
    .select_margin_inputs(),
    "Supply either 'response' or both 'data' and 'response_var'.",
    fixed = TRUE
  )
  expect_error(
    .select_margin_inputs(response = c(1, NA_real_, NA_real_)),
    "Need at least three finite response values",
    fixed = TRUE
  )
  expect_error(
    .select_margin_inputs(data = data.frame(y = 1:3), response_var = "missing"),
    "response_var='missing' not found",
    fixed = TRUE
  )
  expect_error(
    .select_margin_inputs(
      data = data.frame(y = 1:3),
      response_var = "y",
      time_intercepts = TRUE
    ),
    "time_var",
    fixed = TRUE
  )
  expect_error(
    .select_margin_inputs(
      data = data.frame(y = 1:3, time = 1:3),
      response_var = "y",
      time_intercepts = TRUE,
      time_var = "missing_time"
    ),
    "time_var='missing_time' not found",
    fixed = TRUE
  )
})
