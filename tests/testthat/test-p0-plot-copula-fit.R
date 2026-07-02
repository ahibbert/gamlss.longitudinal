test_that("plot_copula_fit input helper requires copula or fit", {
  expect_error(
    gamlss.longitudinal:::.plot_copula_fit_inputs(
      u1 = c(0.2, 0.4, 0.6),
      u2 = c(0.3, 0.5, 0.7),
      copula_dist = NULL
    ),
    "Supply 'copula_dist' or a fitted 'fit'",
    fixed = TRUE
  )
})

test_that("plot_copula_fit density helper builds ungrouped grids", {
  pair_data <- data.frame(
    u1 = c(0.2, 0.4, 0.6),
    u2 = c(0.3, 0.5, 0.7),
    theta_pair = c(0, 0, 0),
    zeta_pair = c(0, 0, 0)
  )
  spec <- list(family = "N", par = 0, par2 = 0, tau = 0)

  grid <- gamlss.longitudinal:::.plot_copula_fit_density_grid(
    pair_data = pair_data,
    spec = spec,
    grid_n = 4,
    max_pairs_overlay = 10
  )

  expect_equal(nrow(grid), 16L)
  expect_setequal(names(grid), c("u1", "u2", "density"))
  expect_true(all(is.finite(grid$density)))
})

test_that("plot_copula_fit plot builder returns ggplot with optional facets", {
  pair_plot <- data.frame(
    u1 = c(0.2, 0.4, 0.6, 0.8),
    u2 = c(0.3, 0.5, 0.7, 0.9),
    split_group = rep(c("early", "late"), each = 2)
  )
  density_grid <- expand.grid(
    u1 = seq(0.1, 0.9, length.out = 4),
    u2 = seq(0.1, 0.9, length.out = 4)
  )
  density_grid$density <- stats::dnorm(density_grid$u1) * stats::dnorm(density_grid$u2)
  density_grid$split_group <- rep(c("early", "late"), length.out = nrow(density_grid))
  spec <- list(family = "N", par = 0, par2 = 0, tau = 0)

  p <- gamlss.longitudinal:::.plot_copula_fit_plot(
    pair_plot = pair_plot,
    density_grid = density_grid,
    spec = spec,
    transform = "uniform",
    by_time = TRUE
  )

  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetWrap"))
})
