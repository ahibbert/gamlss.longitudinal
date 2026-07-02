test_that("model-check status helper maps statuses to severities", {
  expect_equal(gamlss.longitudinal:::.gl_check_status_severity("PASS"), "ok")
  expect_equal(gamlss.longitudinal:::.gl_check_status_severity("REVIEW"), "review")
  expect_equal(gamlss.longitudinal:::.gl_check_status_severity("FAIL"), "concern")
  expect_equal(gamlss.longitudinal:::.gl_check_status_severity("UNKNOWN"), "ok")
})

test_that("model-check dependence cutoff validator preserves threshold contract", {
  expect_equal(gamlss.longitudinal:::.gl_validate_dependence_cor_cutoff("0.25"), 0.25)
  expect_equal(gamlss.longitudinal:::.gl_validate_dependence_cor_cutoff(0.999), 0.999)

  expect_error(
    gamlss.longitudinal:::.gl_validate_dependence_cor_cutoff(0),
    "between 0 and 1",
    fixed = TRUE
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_dependence_cor_cutoff(1),
    "between 0 and 1",
    fixed = TRUE
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_dependence_cor_cutoff(c(0.2, 0.3)),
    "between 0 and 1",
    fixed = TRUE
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_dependence_cor_cutoff(NA_real_),
    "between 0 and 1",
    fixed = TRUE
  )
})

test_that("model-check calibration helpers summarize PIT and tail ratios", {
  pit <- c(0.01, 0.05, 0.5, 0.95, 0.99)

  pit_stats <- gamlss.longitudinal:::.gl_pit_calibration_stats(pit)
  tail <- gamlss.longitudinal:::.gl_tail_calibration_stats(pit)

  expect_s3_class(pit_stats, "data.frame")
  expect_equal(nrow(pit_stats), 1L)
  expect_named(pit_stats, c("n", "mean", "sd", "expected_sd", "ks_p_value"))
  expect_equal(pit_stats$n, length(pit))
  expect_equal(pit_stats$expected_sd, sqrt(1 / 12))
  expect_true(is.finite(pit_stats$ks_p_value) || is.na(pit_stats$ks_p_value))

  expect_s3_class(tail$tail_summary, "data.frame")
  expect_s3_class(tail$tail_stats, "data.frame")
  expect_named(tail$tail_stats, "tail_ratio_max")
  expect_equal(
    tail$tail_stats$tail_ratio_max,
    attr(tail$tail_summary, "tail_ratio_max")
  )
})

test_that("model-check row helper preserves reviewer table schema", {
  out <- gamlss.longitudinal:::.gl_check_row(
    area = "Convergence",
    quantity_checked = "object$convergence$converged",
    value = "not TRUE",
    threshold_condition = "Not TRUE",
    default = "n/a",
    status = "FAIL",
    message = "Convergence was not confirmed.",
    action = "Refit."
  )

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  expect_named(
    out,
    c(
      "area",
      "quantity_checked",
      "value",
      "threshold_condition",
      "default",
      "status",
      "severity",
      "message",
      "action"
    )
  )
  expect_equal(out$severity, "concern")
})

test_that("model-check table applies threshold statuses consistently", {
  checks <- gamlss.longitudinal:::.gl_check_table(
    summary_obj = list(convergence = list(converged = FALSE)),
    scores = list(),
    pit_stats = list(ks_p_value = 0.01),
    tail_stats = list(tail_ratio_max = 2.5),
    lag1_cor = 0.4,
    dependence_cor_cutoff = 0.25,
    vcov_method = "numderiv"
  )

  expect_equal(nrow(checks), 5L)
  expect_true(all(c("FAIL", "REVIEW") %in% checks$status))
  expect_equal(
    checks$status[match(c("Convergence", "Variance calculation"), checks$area)],
    c("FAIL", "REVIEW")
  )
  expect_equal(
    checks$severity[match(c("Convergence", "Variance calculation"), checks$area)],
    c("concern", "review")
  )
})

test_that("basic check result prioritizes failures then review statuses", {
  display <- gamlss.longitudinal:::.gl_basic_checks(data.frame(
    area = "A",
    status = "PASS",
    value = "ok",
    threshold_condition = "none",
    message = "message",
    action = "hidden"
  ))

  expect_named(display, c("area", "status", "value", "threshold_condition", "message"))

  checks <- data.frame(status = c("PASS", "REVIEW", "FAIL"))
  expect_equal(gamlss.longitudinal:::.gl_basic_checks_result(checks), "failed")

  checks <- data.frame(status = c("PASS", "REVIEW"))
  expect_equal(gamlss.longitudinal:::.gl_basic_checks_result(checks), "review")

  checks <- data.frame(status = c("PASS", "PASS"))
  expect_equal(gamlss.longitudinal:::.gl_basic_checks_result(checks), "passed")
})
