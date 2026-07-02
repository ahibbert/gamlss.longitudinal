test_that("plot_dist data helpers build per-time response and pseudo-observation slices", {
  plot_data <- data.frame(
    subject = c("a", "b", "c", "a", "b", "c"),
    time = c(1, 1, 1, 2, 2, 2),
    response = c(2, 1, NA, 4, 6, 5),
    stringsAsFactors = FALSE
  )

  slices <- gamlss.longitudinal:::.plot_dist_margin_slices(plot_data, time_values = 1:2)

  expect_equal(length(slices$margin_data), 2L)
  expect_equal(slices$margin_data[[1]]$response, c(2, 1, NA))
  expect_equal(slices$margin_pseudo[[1]]$u, c(2 / 3, 1 / 3, NA))
  expect_equal(slices$margin_pseudo[[2]]$u, c(1 / 4, 3 / 4, 2 / 4))
})

test_that("plot_dist control helper normalizes overlay choices and validates inputs", {
  fit <- structure(list(), class = "gamlss.longitudinal")

  default_controls <- gamlss.longitudinal:::.plot_dist_resolve_controls(
    fit = NULL,
    offdiag_scale = "pseudo",
    transform = "normal",
    overlay = NULL,
    copula_dist = NULL,
    margin_dist = NULL
  )
  fit_controls <- gamlss.longitudinal:::.plot_dist_resolve_controls(
    fit = fit,
    offdiag_scale = "pseudo",
    transform = "uniform",
    overlay = NULL,
    copula_dist = NULL,
    margin_dist = NULL
  )

  expect_equal(default_controls$overlay, "none")
  expect_equal(fit_controls$overlay, "model")
  expect_error(
    gamlss.longitudinal:::.plot_dist_resolve_controls(
      fit = NULL,
      offdiag_scale = "pseudo",
      transform = "normal",
      overlay = "model",
      copula_dist = NULL,
      margin_dist = NULL
    ),
    "'fit' is required"
  )
  expect_warning(
    gamlss.longitudinal:::.plot_dist_resolve_controls(
      fit = fit,
      offdiag_scale = "response",
      transform = "normal",
      overlay = "model",
      copula_dist = NULL,
      margin_dist = NULL
    ),
    "offdiag_scale = 'pseudo'"
  )
})

test_that("plot_dist data helpers build response-scale off-diagonal data", {
  margin_data <- list(
    data.frame(subject = c("a", "b", "c"), response = c(1, 2, NA)),
    data.frame(subject = c("a", "b", "c"), response = c(3, NA, 5))
  )
  margin_pseudo <- list()

  out <- gamlss.longitudinal:::.plot_dist_offdiag_data(
    i = 1,
    j = 2,
    margin_data = margin_data,
    margin_pseudo = margin_pseudo,
    offdiag_scale = "response",
    transform = "uniform"
  )

  expect_equal(out$input_data, data.frame(X1 = 1, X2 = 3))
})

test_that("plot_dist data helpers build pseudo-observation off-diagonal data", {
  margin_data <- list()
  margin_pseudo <- list(
    data.frame(subject = c("a", "b", "c"), u = c(0.25, 0.5, NA)),
    data.frame(subject = c("a", "b", "c"), u = c(0.75, NA, 0.5))
  )

  uniform_out <- gamlss.longitudinal:::.plot_dist_offdiag_data(
    i = 1,
    j = 2,
    margin_data = margin_data,
    margin_pseudo = margin_pseudo,
    offdiag_scale = "pseudo",
    transform = "uniform"
  )
  normal_out <- gamlss.longitudinal:::.plot_dist_offdiag_data(
    i = 1,
    j = 2,
    margin_data = margin_data,
    margin_pseudo = margin_pseudo,
    offdiag_scale = "pseudo",
    transform = "normal"
  )

  expect_equal(uniform_out$input_data, data.frame(X1 = 0.25, X2 = 0.75))
  expect_equal(
    normal_out$input_data,
    data.frame(X1 = stats::qnorm(0.25), X2 = stats::qnorm(0.75)),
    tolerance = 1e-12
  )
})

test_that("plot_dist diagonal panel helper returns histogram panel", {
  margin_data <- list(data.frame(subject = c("a", "b", "c"), response = c(1, 2, 3)))

  p <- gamlss.longitudinal:::.plot_dist_diagonal_panel(
    i = 1,
    margin_data = margin_data,
    time_values = 1,
    overlay = "none",
    margin_dist = NULL,
    fit = NULL,
    fit_diag_data = NULL,
    grid_n = 20
  )

  layer_classes <- vapply(p$layers, function(layer) class(layer$geom)[1], character(1))
  expect_s3_class(p, "ggplot")
  expect_true("GeomBar" %in% layer_classes)
})

test_that("plot_dist off-diagonal panel helper returns scatter density panel with correlations", {
  margin_data <- list(
    data.frame(subject = c("a", "b", "c"), response = c(1, 2, 3)),
    data.frame(subject = c("a", "b", "c"), response = c(3, 2, 1))
  )

  p <- gamlss.longitudinal:::.plot_dist_offdiag_panel(
    i = 1,
    j = 2,
    margin_data = margin_data,
    margin_pseudo = list(),
    offdiag_scale = "response",
    transform = "uniform",
    show_cor_stats = TRUE,
    overlay = "none",
    copula_spec = NULL,
    fit_pair_data = NULL,
    fit_copula_spec = NULL,
    time_values = 1:2,
    grid_n = 20,
    contour_bins = 4
  )

  layer_classes <- vapply(p$layers, function(layer) class(layer$geom)[1], character(1))
  expect_s3_class(p, "ggplot")
  expect_true("GeomPoint" %in% layer_classes)
  expect_true("GeomDensity2d" %in% layer_classes)
  expect_match(p$labels$subtitle, "Pearson r = -1.000")
})

test_that("plot_dist overlay helper adds selected copula contours", {
  input_data <- data.frame(
    X1 = c(0.2, 0.4, 0.6, 0.8),
    X2 = c(0.3, 0.5, 0.7, 0.9)
  )
  base_plot <- ggplot2::ggplot(input_data, ggplot2::aes(X1, X2)) +
    ggplot2::geom_point()
  copula_spec <- list(family = "N", par = 0, par2 = 0, tau = 0)

  p <- gamlss.longitudinal:::.plot_dist_add_selected_copula_overlay(
    base_plot,
    input_data = input_data,
    copula_spec = copula_spec,
    transform = "uniform",
    grid_n = 4,
    contour_bins = 3
  )

  layer_classes <- vapply(p$layers, function(layer) class(layer$geom)[1], character(1))
  expect_true("GeomContour" %in% layer_classes)
})
