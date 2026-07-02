test_that("copula diagnostics empirical overlay panel renders finite density grids", {
  pair_data <- data.frame(
    u1 = c(0.2, 0.4, 0.6, 0.8),
    u2 = c(0.3, 0.5, 0.7, 0.9)
  )
  density_grid <- expand.grid(u1 = seq(0.2, 0.8, length.out = 3), u2 = seq(0.2, 0.8, length.out = 3))
  density_grid$density <- 1
  axis_labels <- .copula_v2_axis_labels("uniform")

  p <- .copula_v2_empirical_overlay_panel(
    pair_data_plot = pair_data,
    density_grid = density_grid,
    plot1_style = "scatter",
    contour_bins = 3,
    is_grouped = FALSE,
    copula_family_name = "N",
    x_label = axis_labels$x,
    y_label = axis_labels$y
  )

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Empirical Copula with Fitted Overlay")
  expect_true(is.expression(p$labels$x))
})

test_that("copula diagnostics empirical overlay panel handles grouped and empty density paths", {
  pair_data <- data.frame(
    u1 = c(0.2, 0.4, 0.6, 0.8),
    u2 = c(0.3, 0.5, 0.7, 0.9),
    split_group = c("early", "early", "late", "late")
  )
  density_grid <- expand.grid(u1 = seq(0.2, 0.8, length.out = 3), u2 = seq(0.2, 0.8, length.out = 3))
  density_grid$density <- 1
  density_grid$split_group <- rep(c("early", "late"), length.out = nrow(density_grid))
  axis_labels <- .copula_v2_axis_labels("uniform")

  grouped <- .copula_v2_empirical_overlay_panel(
    pair_data_plot = pair_data,
    density_grid = density_grid,
    plot1_style = "bins",
    contour_bins = 3,
    is_grouped = TRUE,
    copula_family_name = "N",
    x_label = axis_labels$x,
    y_label = axis_labels$y
  )

  empty_density <- density_grid
  empty_density$density <- NA_real_
  empty <- .copula_v2_empirical_overlay_panel(
    pair_data_plot = pair_data,
    density_grid = empty_density,
    plot1_style = "bins",
    contour_bins = 3,
    is_grouped = TRUE,
    copula_family_name = "N",
    x_label = axis_labels$x,
    y_label = axis_labels$y
  )

  expect_s3_class(grouped, "ggplot")
  expect_identical(names(grouped$facet$params$facets), "split_group")
  expect_s3_class(empty, "ggplot")
  expect_equal(empty$labels$title, "Empirical Copula with Fitted Overlay")
})

test_that("copula diagnostics cut summary panel renders finite and empty summaries", {
  quartile_df <- data.frame(
    cut_group = c("C1", "C2", "C3"),
    tau_emp = c(0.1, 0.2, 0.3),
    tau_fit = c(0.15, 0.25, 0.35),
    n_pairs = c(5L, 5L, 5L)
  )

  p <- .copula_v2_cut_summary_panel(
    quartile_df = quartile_df,
    plot2_cuts = 3,
    is_grouped = FALSE,
    tau_ylim = NULL
  )
  empty <- .copula_v2_cut_summary_panel(
    quartile_df = data.frame(),
    plot2_cuts = 3,
    is_grouped = FALSE,
    tau_ylim = NULL
  )

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Observed vs Fitted Correlation by Quantile Bin")
  expect_equal(p$labels$y, "Kendall's tau")
  expect_s3_class(empty, "ggplot")
  expect_equal(empty$labels$title, "Observed vs Fitted Correlation by Quantile Bin")
})

test_that("copula diagnostics cut summary panel supports grouped fixed-scale panels", {
  quartile_df <- data.frame(
    cut_group = c("C1", "C2", "C1", "C2"),
    tau_emp = c(0.1, 0.2, 0.3, 0.4),
    tau_fit = c(0.15, 0.25, 0.35, 0.45),
    n_pairs = c(5L, 5L, 5L, 5L),
    split_group = c("early", "early", "late", "late")
  )

  p <- .copula_v2_cut_summary_panel(
    quartile_df = quartile_df,
    plot2_cuts = 2,
    is_grouped = TRUE,
    tau_ylim = c(-1, 1)
  )

  expect_s3_class(p, "ggplot")
  expect_identical(names(p$facet$params$facets), "split_group")
  expect_equal(p$coordinates$limits$y, c(-1, 1))
})

test_that("copula diagnostics Rosenblatt time panel handles finite and empty data", {
  finite <- data.frame(time = c(1, 1, 2, 2), z = c(-0.5, 0.2, 0.1, 0.6))
  p <- .copula_v2_rosenblatt_time_panel(finite)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Rosenblatt Normal Scores by Time")
  expect_equal(p$labels$y, "Normal score")

  empty <- .copula_v2_rosenblatt_time_panel(data.frame(z = c(NA_real_, Inf)))
  expect_s3_class(empty, "ggplot")
  expect_equal(empty$labels$title, "Rosenblatt Normal Scores by Time")
})

test_that("copula diagnostics Rosenblatt QQ panel handles finite and empty data", {
  p <- .copula_v2_rosenblatt_qq_panel(data.frame(theoretical = c(-1, 1), observed = c(-0.8, 0.9)))

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Rosenblatt Normal QQ")
  expect_equal(p$labels$x, "Theoretical normal quantile")

  empty <- .copula_v2_rosenblatt_qq_panel(data.frame())
  expect_s3_class(empty, "ggplot")
  expect_equal(empty$labels$title, "Rosenblatt Normal QQ")
})

test_that("copula diagnostics Rosenblatt lag panel handles finite and empty data", {
  p <- .copula_v2_rosenblatt_lag_panel(data.frame(z_prev = c(-1, 0, 1), z_curr = c(-0.8, 0.1, 0.7)))

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Rosenblatt Lag Plot")
  expect_true(is.expression(p$labels$x))

  empty <- .copula_v2_rosenblatt_lag_panel(data.frame(z_prev = NA_real_, z_curr = NA_real_))
  expect_s3_class(empty, "ggplot")
  expect_equal(empty$labels$title, "Rosenblatt Lag Plot")
})

test_that("copula diagnostics Kendall panel handles finite and empty data", {
  p <- .copula_v2_kendall_panel(data.frame(fitted = c(0.1, 0.5, 0.9), empirical = c(0.2, 0.4, 0.8)))

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Kendall Function Diagnostic")
  expect_equal(p$labels$x, "Fitted copula probability")

  empty <- .copula_v2_kendall_panel(data.frame())
  expect_s3_class(empty, "ggplot")
  expect_equal(empty$labels$title, "Kendall Function Diagnostic")
})

test_that("copula diagnostics tail panel renders finite tail data", {
  tail_long <- data.frame(
    threshold = c(0.05, 0.10, 0.05, 0.10),
    tail = c("Lower", "Lower", "Upper", "Upper"),
    source = c("Empirical", "Fitted", "Empirical", "Fitted"),
    probability = c(0.01, 0.02, 0.03, 0.04)
  )

  p <- .copula_v2_tail_panel(
    tail_long = tail_long,
    title = "Tail Co-occurrence",
    message_subtitle = "message subtitle",
    empty_message = "empty",
    plot_subtitle = "plot subtitle",
    y_label = "Joint probability"
  )

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Tail Co-occurrence")
  expect_equal(p$labels$y, "Joint probability")
})

test_that("copula diagnostics tail panel renders message plot for empty data", {
  p <- .copula_v2_tail_panel(
    tail_long = data.frame(),
    title = "Conditional Tail Exceedance",
    message_subtitle = "message subtitle",
    empty_message = "No finite conditional tail diagnostics",
    plot_subtitle = "plot subtitle",
    y_label = "Conditional probability",
    ylim = c(0, 1)
  )

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Conditional Tail Exceedance")
})

test_that("copula diagnostics residual lag panel handles finite and empty data", {
  p <- .copula_v2_residual_lag_panel(data.frame(lag = c(1L, 2L), cor_z = c(0.1, -0.2), n_pairs = c(12L, 10L)))

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Residual Dependence by Lag")
  expect_equal(p$labels$y, "Correlation")

  empty <- .copula_v2_residual_lag_panel(data.frame(lag = 1L, cor_z = NA_real_, n_pairs = 0L))
  expect_s3_class(empty, "ggplot")
  expect_equal(empty$labels$title, "Residual Dependence by Lag")
})
