test_that("plot_margin_fit transform helpers preserve response-scale data", {
  obs <- data.frame(response = c(1, 2), split_group = "All")
  density <- data.frame(response = c(1, 2), density = c(0.4, 0.2), split_group = "All")

  expect_equal(
    gamlss.longitudinal:::.plot_margin_transform_observed(obs, "response"),
    obs
  )
  expect_equal(
    gamlss.longitudinal:::.plot_margin_transform_density(density, "response"),
    density
  )
})

test_that("plot_margin_fit transform helpers apply log response scale", {
  obs <- data.frame(response = c(1, exp(1)), split_group = "All")
  density <- data.frame(
    response = c(1, exp(1), -1, NA),
    density = c(0.4, 0.2, 99, 99),
    split_group = "All"
  )

  obs_out <- gamlss.longitudinal:::.plot_margin_transform_observed(obs, "log")
  density_out <- gamlss.longitudinal:::.plot_margin_transform_density(density, "log")

  expect_equal(obs_out$response, c(0, 1), tolerance = 1e-12)
  expect_equal(density_out$response, c(0, 1), tolerance = 1e-12)
  expect_equal(density_out$density, c(0.4, 0.2 * exp(1)), tolerance = 1e-12)
})

test_that("plot_margin_fit observed log transform rejects non-positive responses", {
  obs <- data.frame(response = c(1, 0), split_group = "All")

  expect_error(
    gamlss.longitudinal:::.plot_margin_transform_observed(obs, "log"),
    "requires positive responses",
    fixed = TRUE
  )
})

test_that("plot_margin_fit plot builder returns ggplot with optional facets", {
  obs <- data.frame(
    response = c(1, 2, 3, 4),
    split_group = rep(c("A", "B"), each = 2),
    stringsAsFactors = FALSE
  )
  density <- data.frame(
    response = c(1, 2, 3, 4),
    density = c(0.2, 0.3, 0.3, 0.2),
    split_group = rep(c("A", "B"), each = 2),
    stringsAsFactors = FALSE
  )

  p <- gamlss.longitudinal:::.plot_margin_fit_plot(
    obs = obs,
    density_grid = density,
    margin_dist = gamlss.dist::NO(),
    by_time = TRUE
  )

  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetWrap"))
  expect_equal(p$layers[[1L]]$aes_params$fill, .gl_plot_colours()[["neutral"]])
  expect_true(is.na(p$layers[[1L]]$aes_params$colour))
  expect_equal(p$layers[[2L]]$aes_params$colour, .gl_plot_colours()[["fitted"]])
})

test_that("plot_margin_fit preserves the supplied time-panel order", {
  obs <- data.frame(
    response = 1:4,
    split_group = c("8", "16", "24", "32"),
    stringsAsFactors = FALSE
  )
  density <- data.frame(
    response = 1:4,
    density = rep(0.25, 4),
    split_group = c("8", "16", "24", "32"),
    stringsAsFactors = FALSE
  )

  p <- gamlss.longitudinal:::.plot_margin_fit_plot(
    obs = obs,
    density_grid = density,
    margin_dist = gamlss.dist::NO(),
    by_time = TRUE
  )

  built <- ggplot2::ggplot_build(p)
  expect_equal(as.character(built$layout$layout$split_group), c("8", "16", "24", "32"))

  discrete_p <- gamlss.longitudinal:::.plot_margin_fit_plot(
    obs = obs,
    density_grid = density,
    margin_dist = gamlss.dist::PO(),
    by_time = TRUE
  )
  discrete_built <- ggplot2::ggplot_build(discrete_p)
  expect_equal(
    as.character(discrete_built$layout$layout$split_group),
    c("8", "16", "24", "32")
  )
})

test_that("wide discrete margins share capped bins and conserve displayed probability", {
  obs <- data.frame(
    response = c(0:49, 0:99),
    split_group = c(rep("early", 50), rep("late", 100)),
    stringsAsFactors = FALSE
  )
  density <- data.frame(
    response = rep(0:99, 2),
    density = c(stats::dpois(0:99, 10), stats::dpois(0:99, 20)),
    split_group = rep(c("early", "late"), each = 100),
    stringsAsFactors = FALSE
  )

  display <- gamlss.longitudinal:::.plot_margin_discrete_display_data(
    obs,
    density,
    max_bins = 10
  )

  expect_equal(display$bin_width, 10L)
  expect_lte(length(unique(display$obs$response)), 10L)
  expect_equal(
    as.numeric(tapply(display$obs$probability, display$obs$split_group, sum)),
    c(1, 1)
  )
  expect_equal(
    as.numeric(tapply(display$density$density, display$density$split_group, sum)),
    c(sum(stats::dpois(0:99, 10)), sum(stats::dpois(0:99, 20))),
    tolerance = 1e-12
  )
  expect_identical(
    unique(display$obs$bin_start[display$obs$split_group == "early"]),
    seq(0, 40, by = 10)
  )
})

test_that("small discrete supports retain one bin per integer", {
  obs <- data.frame(response = c(0, 1, 1, 3), split_group = "All")
  density <- data.frame(
    response = 0:3,
    density = c(0.2, 0.4, 0.3, 0.1),
    split_group = "All"
  )

  display <- gamlss.longitudinal:::.plot_margin_discrete_display_data(
    obs,
    density,
    max_bins = 30
  )

  expect_equal(display$bin_width, 1L)
  expect_equal(display$density$response, 0:3)
  expect_equal(display$density$density, density$density)
})

test_that("discrete plots can reserve a separate zero bin", {
  obs <- data.frame(
    response = c(rep(0, 40), 1:99),
    split_group = "All"
  )
  density <- data.frame(
    response = 0:99,
    density = c(0.4, 0.6 * stats::dpois(0:98, 20)),
    split_group = "All"
  )

  display <- gamlss.longitudinal:::.plot_margin_discrete_display_data(
    obs,
    density,
    max_bins = 10,
    separate_zero = TRUE
  )

  zero_obs <- display$obs[display$obs$bin_start == 0, ]
  zero_density <- display$density[display$density$bin_start == 0, ]
  expect_equal(nrow(zero_obs), 1L)
  expect_equal(zero_obs$bin_end, 0)
  expect_equal(zero_obs$bin_width, 1)
  expect_equal(zero_obs$probability, 40 / 139)
  expect_equal(zero_density$density, 0.4)
  expect_equal(nrow(display$density), 10L)
  expect_equal(display$density$bin_start[2L], 1)
  expect_equal(display$density$bin_end[2L], 11)
  expect_equal(sum(display$obs$probability), 1)
  expect_equal(sum(display$density$density), sum(density$density))

  p <- gamlss.longitudinal:::.plot_margin_fit_plot(
    obs = obs,
    density_grid = density,
    margin_dist = gamlss.dist::ZIP(),
    bins = 10,
    separate_zero = TRUE
  )
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  zero_bar <- built$data[[1L]][which.min(built$data[[1L]]$xmin), ]
  expect_equal(zero_bar$xmax - zero_bar$xmin, 1)
})

test_that("plot_margin_fit validates separate_zero", {
  expect_error(
    plot_margin_fit(
      x = 0:3,
      margin_dist = gamlss.dist::PO(),
      separate_zero = NA,
      plot = FALSE
    ),
    "'separate_zero' must be TRUE or FALSE.",
    fixed = TRUE
  )
})

test_that("plot_margin_fit raw helper requires time column for time overlays", {
  expect_error(
    gamlss.longitudinal:::.plot_margin_fit_from_raw(
      x = data.frame(response = 1:4),
      margin_dist = gamlss.dist::NO(),
      response_var = "response",
      by_time = TRUE,
      time_var = "time"
    ),
    "containing 'time_var'",
    fixed = TRUE
  )
})

test_that("plot_margin_fit facets raw data by a supplied factor column", {
  set.seed(2401)
  dat <- data.frame(
    response = stats::rnorm(60),
    time = rep(1:3, 20),
    arm = factor(rep(c("Treatment", "Control"), each = 30),
      levels = c("Control", "Treatment"))
  )

  out <- suppressWarnings(plot_margin_fit(
    data = dat,
    margin_dist = gamlss.dist::NO(),
    response_var = "response",
    by = "arm",
    plot = FALSE
  ))

  expect_identical(out$by, "arm")
  expect_s3_class(out$plot, "ggplot")
  expect_true(inherits(out$plot$facet, "FacetWrap"))
  expect_equal(levels(out$data$split_group), c("Control", "Treatment"))
  expect_setequal(unique(out$density$split_group), c("Control", "Treatment"))

  expect_warning(
    by_wins <- plot_margin_fit(
      data = dat,
      margin_dist = gamlss.dist::NO(),
      response_var = "response",
      time_var = "time",
      by_time = TRUE,
      by = "arm",
      plot = FALSE
    ),
    "using by='arm'",
    fixed = TRUE
  )
  expect_identical(by_wins$by, "arm")
  expect_error(
    plot_margin_fit(
      data = dat,
      margin_dist = gamlss.dist::NO(),
      response_var = "response",
      by = "missing_group",
      plot = FALSE
    ),
    "not found in provided data",
    fixed = TRUE
  )
})
