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
  expect_false(p$facet$params$free$x)
  expect_false(p$facet$params$free$y)
  expect_match(p$labels$subtitle, "Empirical density with fitted copula contours", fixed = TRUE)
  expect_identical(p$scales$get_scales("fill")$guide, "none")
})

test_that("plot_copula_fit preserves the supplied facet order", {
  pair_plot <- data.frame(
    u1 = c(0.2, 0.4, 0.6, 0.8),
    u2 = c(0.3, 0.5, 0.7, 0.9),
    split_group = factor(
      rep(c("Treatment", "Control"), each = 2),
      levels = c("Treatment", "Control")
    )
  )
  density_grid <- expand.grid(
    u1 = seq(0.1, 0.9, length.out = 4),
    u2 = seq(0.1, 0.9, length.out = 4),
    split_group = c("Control", "Treatment"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  density_grid$density <- stats::dnorm(density_grid$u1) *
    stats::dnorm(density_grid$u2)
  spec <- list(family = "N", par = 0, par2 = 0, tau = 0)

  p <- gamlss.longitudinal:::.plot_copula_fit_plot(
    pair_plot = pair_plot,
    density_grid = density_grid,
    spec = spec,
    transform = "uniform",
    by_time = TRUE
  )

  built <- ggplot2::ggplot_build(p)
  expect_equal(
    as.character(built$layout$layout$split_group),
    c("Treatment", "Control")
  )
})

test_that("plot_copula_fit facets raw adjacent pairs by a supplied factor column", {
  dat <- simulate_longitudinal_dataset(
    n = 30,
    times = 1:3,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "N",
    copula_params = list(theta = 0.45),
    seed = 2402
  )
  subject_number <- as.integer(as.character(dat$subject))
  dat$arm <- factor(
    ifelse(subject_number <= 15, "Treatment", "Control"),
    levels = c("Control", "Treatment")
  )

  out <- suppressWarnings(plot_copula_fit(
    data = dat,
    copula_dist = "N",
    u_var = "u",
    by = "arm",
    grid_n = 8,
    plot = FALSE
  ))

  expect_identical(out$by, "arm")
  expect_s3_class(out$plot, "ggplot")
  expect_true(inherits(out$plot$facet, "FacetWrap"))
  expect_equal(levels(out$pair_data$split_group), c("Control", "Treatment"))
  expect_setequal(unique(out$density$split_group), c("Control", "Treatment"))
  expect_true(all(c("subject", "time_left", "time_right") %in% names(out$pair_data)))

  expect_warning(
    by_wins <- plot_copula_fit(
      data = dat,
      copula_dist = "N",
      u_var = "u",
      by_time = TRUE,
      by = "arm",
      grid_n = 8,
      plot = FALSE
    ),
    "using by='arm'",
    fixed = TRUE
  )
  expect_identical(by_wins$by, "arm")
})
