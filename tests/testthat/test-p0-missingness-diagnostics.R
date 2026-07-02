test_that("missingness helper validation normalizes inputs and bounds", {
  out <- .missingness_validate_args(
    data = list(y = c(1, NA), x = c(1, 2)),
    response_var = "y",
    time_var = NULL,
    subject_var = NULL,
    alpha = "0.1",
    max_factor_levels = 5
  )

  expect_s3_class(out$data, "data.frame")
  expect_equal(out$alpha, 0.1)
  expect_equal(out$max_factor_levels, 5L)

  expect_error(
    .missingness_validate_args(data.frame(y = 1), "missing", NULL, NULL, 0.1, 5),
    "response_var='missing' not found",
    fixed = TRUE
  )
  expect_error(
    .missingness_validate_args(data.frame(y = 1), "y", NULL, NULL, 1, 5),
    "between 0 and 1",
    fixed = TRUE
  )
})

test_that("missingness predictor resolver excludes response subject and truth columns", {
  data <- data.frame(
    y = c(1, NA),
    id = 1:2,
    visit = 1:2,
    true_mu = c(0, 0),
    x = c(1, 2)
  )

  predictors <- .missingness_resolve_predictors(
    data = data,
    response_var = "y",
    predictors = NULL,
    time_var = "visit",
    subject_var = "id"
  )

  expect_equal(predictors, c("visit", "x"))
})

test_that("missingness early result returns reviewer assessments", {
  response_summary <- data.frame(n = 3L, n_missing = 0L, prop_missing = 0)
  predictor_summary <- data.frame(
    predictor = "x",
    used = TRUE,
    stringsAsFactors = FALSE
  )

  out <- .missingness_early_result(response_summary, predictor_summary, alpha = 0.05)

  expect_s3_class(out, "gamlss_longitudinal_missingness_check")
  expect_equal(out$assessment, "no_missing_responses")
})

test_that("missingness assessment labels flagged and unflagged models", {
  flagged <- data.frame(term = "x", p_value = 0.01)
  unflagged <- data.frame(term = "x", p_value = 0.5)

  expect_equal(
    .missingness_assessment(flagged, alpha = 0.05)$assessment,
    "covariate_related_missingness"
  )
  expect_equal(
    .missingness_assessment(unflagged, alpha = 0.05)$assessment,
    "no_detected_covariate_association"
  )
})
