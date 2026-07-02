test_that("longitudinal coverage fitter selects the first finite fit attempt", {
  finite_fit <- structure(
    list(calc_lik_out_end = list(log_lik = c(marginal = -1, copula = -2, joint = -3))),
    class = "gamlss.longitudinal"
  )
  nonfinite_fit <- structure(
    list(calc_lik_out_end = list(log_lik = c(marginal = -1, copula = NA_real_, joint = -3))),
    class = "gamlss.longitudinal"
  )
  failed_fit <- structure(
    list(message = "failed"),
    class = c("simpleError", "error", "condition")
  )
  attempts <- list(
    list(
      value = nonfinite_fit,
      elapsed = 1,
      warnings = c("first warning"),
      attempt_label = "nonfinite"
    ),
    list(
      value = finite_fit,
      elapsed = 2,
      warnings = c("second warning", "first warning"),
      attempt_label = "finite"
    ),
    list(
      value = failed_fit,
      elapsed = 4,
      warnings = character(),
      attempt_label = "failed"
    )
  )

  expect_equal(.coverage_longitudinal_successful_attempts(attempts), c(FALSE, TRUE, FALSE))
  expect_equal(.coverage_longitudinal_chosen_index(attempts), 2L)
  expect_equal(.coverage_longitudinal_elapsed(attempts, 2L), 3)
  expect_equal(
    .coverage_longitudinal_attempt_warnings(attempts, 2L),
    "first warning | second warning"
  )
  expect_equal(
    .coverage_longitudinal_attempt_trace(attempts, 2L),
    "nonfinite > finite"
  )
})

test_that("longitudinal coverage fitter falls back to the final attempt", {
  failed_fit <- structure(
    list(message = "failed"),
    class = c("simpleError", "error", "condition")
  )
  attempts <- list(
    list(value = failed_fit, elapsed = 1, warnings = "a", attempt_label = "first"),
    list(value = failed_fit, elapsed = 2, warnings = "b", attempt_label = "second")
  )

  expect_equal(.coverage_longitudinal_successful_attempts(attempts), c(FALSE, FALSE))
  expect_equal(.coverage_longitudinal_chosen_index(attempts), 2L)
  expect_equal(.coverage_longitudinal_elapsed(attempts, 2L), 3)
  expect_equal(.coverage_longitudinal_attempt_trace(attempts, 2L), "first > second")
})

test_that("longitudinal coverage result row records attempt and benchmark metadata", {
  fit <- structure(
    list(convergence = list(converged = TRUE)),
    class = "gamlss.longitudinal"
  )
  attempts <- list(
    list(value = fit, elapsed = 1, warnings = c("first warning"), attempt_label = "default"),
    list(value = fit, elapsed = 2, warnings = "unused", attempt_label = "cold_start")
  )
  truth_metrics <- c(
    benchmark_mae = 0.1,
    benchmark_rmse = 0.2,
    benchmark_mean_bias = 0.3,
    benchmark_mean_mae = 0.4,
    benchmark_mean_rmse = 0.5,
    benchmark_q90_mae = 0.6,
    benchmark_neg_log_score = 0.7,
    benchmark_upper_tail_error_90 = 0.8,
    benchmark_interval_coverage_95 = 0.9,
    benchmark_interval_width_95 = 1.0,
    benchmark_pit_ks_p_value = 0.11,
    benchmark_pit_mean_abs_error = 0.12,
    benchmark_tail_error_lower_05 = 0.13,
    benchmark_tail_error_upper_05 = 0.14
  )

  out <- .coverage_longitudinal_result_row(
    method = "cg",
    success = TRUE,
    is_fit = TRUE,
    fit = fit,
    chosen = attempts[[1L]],
    chosen_idx = 1L,
    attempts = attempts,
    elapsed = 1,
    max_elapsed_sec = 20,
    loglik = c(marginal = -1, copula = -2, joint = -3),
    abs_error = 0.01,
    rel_error = 0.02,
    theta_tau = 0.25,
    truth = c(copula_tau = 0.3),
    smooth_metrics = c(
      smooth_eta_rmse = 0.4,
      smooth_eta_max_abs_error = 0.5,
      smooth_eta_n = 6
    ),
    truth_metrics = truth_metrics
  )

  expect_equal(out$method, "cg")
  expect_true(out$success)
  expect_true(out$converged)
  expect_equal(out$warnings, "first warning")
  expect_equal(out$fit_attempt, "default")
  expect_equal(out$fit_attempt_count, 1L)
  expect_equal(out$fit_attempt_trace, "default")
  expect_equal(out$benchmark_estimator, "gamlss.longitudinal::cg")
  expect_equal(out$benchmark_rmse, 0.2)
  expect_equal(out$smooth_eta_n, 6)
  expect_equal(out$true_copula_tau, 0.3)
})
