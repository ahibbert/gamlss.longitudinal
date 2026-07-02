test_that("copula diagnostics Rosenblatt QQ data is built from finite scores", {
  rosenblatt_df <- data.frame(
    subject = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    z = c(0.5, NA, -1, Inf)
  )

  qq <- .copula_v2_rosenblatt_qq_data(rosenblatt_df)

  expect_equal(nrow(qq), 2L)
  expect_equal(qq$observed, c(-1, 0.5))
  expect_true(all(is.finite(qq$theoretical)))
})

test_that("copula diagnostics Rosenblatt QQ data handles empty scores", {
  expect_equal(
    nrow(.copula_v2_rosenblatt_qq_data(data.frame(z = c(NA_real_, Inf)))),
    0L
  )
  expect_equal(
    nrow(.copula_v2_rosenblatt_qq_data(data.frame())),
    0L
  )
})

test_that("copula diagnostics tail summaries are converted to long plotting data", {
  tail_df <- data.frame(
    threshold = c(0.05, 0.05),
    tail = c("Lower", "Upper"),
    empirical = c(0.01, 0.02),
    fitted = c(0.015, 0.025)
  )

  long <- .copula_v2_tail_long_data(tail_df)

  expect_equal(nrow(long), 4L)
  expect_setequal(long$source, c("Empirical", "Fitted"))
  expect_equal(long$probability[long$source == "Empirical"], tail_df$empirical)
  expect_equal(long$probability[long$source == "Fitted"], tail_df$fitted)
  expect_equal(nrow(.copula_v2_tail_long_data(data.frame())), 0L)
})

test_that("copula contour comparison controls validate public arguments", {
  controls <- .copula_contour_compare_controls(transform = "normal", diff_scale_limit = 0.1)

  expect_identical(controls$transform, "normal")
  expect_equal(controls$diff_scale_limit, 0.1)
  expect_error(.copula_contour_compare_controls("rank", 0.1), "'transform'")
  expect_error(.copula_contour_compare_controls("uniform", 0), "'diff_scale_limit'")
  expect_error(.copula_contour_compare_controls("uniform", c(0.1, 0.2)), "'diff_scale_limit'")
})

test_that("copula contour comparison metric tables are grouped by time pair", {
  grid_df <- data.frame(
    time_pair = rep(c("1-2", "2-3"), each = 4),
    density_emp = c(0.1, 0.2, 0.3, 0.4, 0.2, 0.2, 0.4, 0.4),
    density_fit = c(0.1, 0.25, 0.25, 0.4, 0.1, 0.3, 0.3, 0.5)
  )

  metrics <- .copula_contour_compare_metric_tables(grid_df)

  expect_identical(metrics$summary$time_pair, c("1-2", "2-3"))
  expect_true(all(c("rmse", "mae", "surface_cor") %in% names(metrics$summary)))
  expect_true(all(metrics$overlap$time_pair %in% c("1-2", "2-3")))
})

test_that("copula contour comparison axis labels follow transform scale", {
  uniform <- .copula_contour_compare_axis_labels("uniform")
  normal <- .copula_contour_compare_axis_labels("normal")

  expect_match(deparse(uniform$x), "U\\[t\\]")
  expect_match(deparse(normal$x), "Phi")
})
