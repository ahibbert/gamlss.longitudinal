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
