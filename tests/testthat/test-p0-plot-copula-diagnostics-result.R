test_that("copula diagnostics dashboard plot list preserves public names", {
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  plots <- .copula_v2_dashboard_plot_list(
    empirical_overlay = p,
    quartile_correlation = p,
    rosenblatt_by_time = p,
    rosenblatt_qq = p,
    rosenblatt_lag = p,
    kendall_function = p,
    tail_cooccurrence = p,
    conditional_tail_exceedance = p,
    residual_lag_correlation = p
  )

  expect_identical(
    names(plots),
    c(
      "empirical_overlay",
      "quartile_correlation",
      "kendall_function",
      "rosenblatt_by_time",
      "rosenblatt_qq",
      "rosenblatt_lag",
      "tail_cooccurrence",
      "conditional_tail_exceedance",
      "residual_lag_correlation"
    )
  )
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
})

test_that("copula diagnostics dashboard arrangement returns a plot object", {
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  plots <- .copula_v2_dashboard_plot_list(p, p, p, p, p, p, p, p, p)

  dashboard <- .copula_v2_arrange_dashboard(plots, dashboard_ncol = 3)

  expect_true(inherits(dashboard, "ggplot") || inherits(dashboard, "ggarrange"))
})

test_that("copula diagnostics result bundle preserves returned names", {
  result <- .copula_v2_plot_result(
    plots = list(a = "plot"),
    dashboard = "dashboard",
    fit_data = data.frame(a = 1),
    pair_data = data.frame(b = 1),
    pair_data_uniform = data.frame(c = 1),
    rosenblatt = data.frame(d = 1),
    rosenblatt_pairs = data.frame(e = 1),
    quartile_summary = data.frame(f = 1),
    kendall_summary = data.frame(g = 1),
    tail_summary = data.frame(h = 1),
    conditional_tail_summary = data.frame(i = 1),
    residual_lag_summary = data.frame(j = 1)
  )

  expect_identical(
    names(result),
    c(
      "plots",
      "dashboard",
      "fit_data",
      "pair_data",
      "pair_data_uniform",
      "rosenblatt",
      "rosenblatt_pairs",
      "quartile_summary",
      "kendall_summary",
      "tail_summary",
      "conditional_tail_summary",
      "residual_lag_summary"
    )
  )
  expect_identical(result$plots, list(a = "plot"))
  expect_identical(result$dashboard, "dashboard")
})
