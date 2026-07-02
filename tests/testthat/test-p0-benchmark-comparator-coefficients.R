test_that("benchmark coefficient table helpers normalize and reshape estimates", {
  mat <- matrix(
    c(1.0, 0.2, 2.0, 0.4),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("(Intercept)", "time_covariate"), c("Estimate", "Std. Error"))
  )

  tbl <- .benchmark_coef_matrix_table(mat, method = "glm", parameter = "mu", level = 0.95)
  tbl$term <- .benchmark_normalize_coef_term(tbl$term, time_var = "visit")
  wide <- .benchmark_wide_values(tbl, "estimate")

  expect_equal(tbl$term, c("intercept", "visit"))
  expect_equal(tbl$estimate, c(1, 2))
  expect_equal(tbl$std_error, c(0.2, 0.4))
  expect_equal(wide$term, c("intercept", "visit"))
  expect_equal(wide$glm, c(1, 2))
})

test_that("benchmark coefficient comparison computes reference differences", {
  fit_a <- stats::lm(mpg ~ wt, data = mtcars)
  fit_b <- stats::lm(mpg ~ wt + hp, data = mtcars)
  fits <- list(primary = fit_a, richer = fit_b)
  results <- data.frame(
    method = c("primary", "richer"),
    success = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  cmp <- .benchmark_coefficient_comparison(
    fits = fits,
    results = results,
    primary_method = "primary",
    parameter = "mu"
  )

  expect_equal(cmp$reference_method, "primary")
  expect_true(all(c("long", "estimates", "uncertainty") %in% names(cmp)))
  expect_true("estimate_diff_vs_reference" %in% names(cmp$uncertainty))
  primary_intercept <- cmp$uncertainty[
    cmp$uncertainty$method == "primary" & cmp$uncertainty$term == "intercept",
    ,
    drop = FALSE
  ]
  expect_equal(primary_intercept$estimate_diff_vs_reference, 0, tolerance = 1e-12)
  expect_true(isTRUE(primary_intercept$ci_overlap_with_reference))
})
