test_that("copula diagnostics plot controls are normalized", {
  controls <- .copula_v2_normalize_plot_controls(
    transform = "normal",
    plot1_style = "bins",
    contour_bins = 2.6,
    time_stratified = FALSE,
    plot2_cuts = 5.4,
    tau_ylim = c(-0.2, 0.8),
    tail_thresholds = c(0.2, 0.05, 0.2, 0.8, NA),
    residual_lags = c(3, 1, 1, 0, NA),
    dashboard_ncol = 2.2
  )

  expect_identical(controls$transform, "normal")
  expect_identical(controls$plot1_style, "bins")
  expect_identical(controls$contour_bins, 3L)
  expect_false(controls$time_stratified)
  expect_identical(controls$plot2_cuts, 5L)
  expect_equal(controls$tau_ylim, c(-0.2, 0.8))
  expect_equal(controls$tail_thresholds, c(0.05, 0.2))
  expect_identical(controls$residual_lags, c(1L, 3L))
  expect_identical(controls$dashboard_ncol, 2L)
})

test_that("copula diagnostics plot controls use documented fallbacks", {
  controls <- suppressWarnings(.copula_v2_normalize_plot_controls(
    transform = "uniform",
    plot1_style = "scatter",
    contour_bins = 1,
    time_stratified = TRUE,
    plot2_cuts = 2,
    tau_ylim = NULL,
    tail_thresholds = c(-1, 0, 0.5, Inf, NA),
    residual_lags = c(-1, 0, NA),
    dashboard_ncol = 1
  ))

  expect_equal(controls$tail_thresholds, c(0.05, 0.10, 0.20))
  expect_identical(controls$residual_lags, 1:3)
  expect_null(controls$tau_ylim)
})

test_that("copula diagnostics plot controls reject invalid user controls", {
  base_args <- list(
    transform = "normal",
    plot1_style = "bins",
    contour_bins = 8,
    time_stratified = FALSE,
    plot2_cuts = 10,
    tau_ylim = NULL,
    tail_thresholds = c(0.05, 0.10, 0.20),
    residual_lags = 1:3,
    dashboard_ncol = 3
  )

  expect_error(
    do.call(.copula_v2_normalize_plot_controls, utils::modifyList(base_args, list(transform = "rank"))),
    "'transform' must be either 'uniform' or 'normal'.",
    fixed = TRUE
  )
  expect_error(
    do.call(.copula_v2_normalize_plot_controls, utils::modifyList(base_args, list(plot1_style = "hex"))),
    "'plot1_style' must be either 'bins' or 'scatter'.",
    fixed = TRUE
  )
  expect_error(
    do.call(.copula_v2_normalize_plot_controls, utils::modifyList(base_args, list(contour_bins = 0))),
    "'contour_bins' must be a single finite number >= 1.",
    fixed = TRUE
  )
  expect_error(
    do.call(.copula_v2_normalize_plot_controls, utils::modifyList(base_args, list(time_stratified = NA))),
    "'time_stratified' must be TRUE or FALSE.",
    fixed = TRUE
  )
  expect_error(
    do.call(.copula_v2_normalize_plot_controls, utils::modifyList(base_args, list(plot2_cuts = 1))),
    "'plot2_cuts' must be a single finite number >= 2.",
    fixed = TRUE
  )
  expect_error(
    do.call(.copula_v2_normalize_plot_controls, utils::modifyList(base_args, list(tau_ylim = c(1, 0)))),
    "'tau_ylim' must be NULL or a numeric vector of length 2 with tau_ylim[1] < tau_ylim[2].",
    fixed = TRUE
  )
  expect_error(
    do.call(.copula_v2_normalize_plot_controls, utils::modifyList(base_args, list(dashboard_ncol = 0))),
    "'dashboard_ncol' must be a single positive integer.",
    fixed = TRUE
  )
})
