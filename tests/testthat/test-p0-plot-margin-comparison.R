test_that("margin comparison display data shares bins and preserves model masses", {
  obs <- data.frame(
    response = c(0:19, 0:39),
    split_group = c(rep("early", 20), rep("late", 40))
  )
  density <- expand.grid(
    response = 0:39,
    split_group = c("early", "late"),
    model = c("PO", "NBI"),
    stringsAsFactors = FALSE
  )
  density$density <- stats::dpois(
    density$response,
    lambda = ifelse(density$model == "PO", 5, 8)
  )

  display <- gamlss.longitudinal:::.plot_margin_comparison_display_data(
    obs,
    density,
    max_bins = 10
  )

  expect_equal(display$bin_width, 4L)
  expect_equal(
    as.numeric(tapply(display$obs$probability, display$obs$split_group, sum)),
    c(1, 1)
  )
  fitted_totals <- aggregate(density ~ split_group + model, display$density, sum)
  expect_equal(
    fitted_totals$density[fitted_totals$model == "PO"],
    rep(sum(stats::dpois(0:39, 5)), 2),
    tolerance = 1e-12
  )
  expect_equal(
    fitted_totals$density[fitted_totals$model == "NBI"],
    rep(sum(stats::dpois(0:39, 8)), 2),
    tolerance = 1e-12
  )
})

test_that("margin comparison KDE is nonnegative and scaled for every time panel", {
  obs <- data.frame(
    response = c(0, 0, 1, 2, 3, 0, 2, 4, 6, 8),
    split_group = rep(c("early", "late"), each = 5)
  )

  kde <- gamlss.longitudinal:::.plot_margin_comparison_kde(
    obs,
    bin_width = 2,
    support = c(lower = 0, upper = 8),
    n = 128
  )

  expect_setequal(unique(kde$split_group), c("early", "late"))
  expect_true(all(is.finite(kde$density)))
  expect_true(all(kde$density >= 0))
  expect_gte(min(kde$response), 0)
})

test_that("margin comparison plot draws bars, fitted steps, KDE, and time facets", {
  obs <- data.frame(
    response = rep(0:3, 2),
    probability = rep(c(0.4, 0.3, 0.2, 0.1), 2),
    split_group = factor(rep(c("8", "16"), each = 4), levels = c("8", "16"))
  )
  density <- expand.grid(
    response = 0:3,
    model = c("PO", "NBI"),
    split_group = factor(c("8", "16"), levels = c("8", "16")),
    stringsAsFactors = FALSE
  )
  density$density <- rep(c(0.4, 0.3, 0.2, 0.1), 4)
  kde <- data.frame(
    response = rep(seq(0, 3, length.out = 20), 2),
    density = rep(seq(0.4, 0.1, length.out = 20), 2),
    split_group = factor(rep(c("8", "16"), each = 20), levels = c("8", "16"))
  )

  p <- gamlss.longitudinal:::.plot_margin_comparison_plot(
    obs,
    density,
    kde,
    bin_width = 1
  )

  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetWrap"))
  expect_equal(length(p$layers), 3L)
  expect_equal(p$labels$colour, "Distribution")

  no_kde <- gamlss.longitudinal:::.plot_margin_comparison_plot(
    obs,
    density,
    kde[0, ],
    bin_width = 1
  )
  expect_equal(length(no_kde$layers), 2L)
  expect_false(grepl("KDE", no_kde$labels$subtitle, fixed = TRUE))
})

test_that("margin comparison validates routes and matched fitted data", {
  expect_error(
    plot_margin_comparison(),
    "exactly one",
    fixed = TRUE
  )
  expect_error(
    plot_margin_comparison(data = data.frame(response = 1:3, time = 1:3), families = "PO"),
    "at least two",
    fixed = TRUE
  )
})
