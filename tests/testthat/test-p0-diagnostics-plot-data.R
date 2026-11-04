test_that("diagnostic QQ helper builds unsplit and split frames", {
  z <- c(0.5, -1, 0)
  no_split <- list(split_by = FALSE, group = NULL)
  split_info <- list(split_by = TRUE, group = factor(c("b", "a", "b")))

  qq <- gamlss.longitudinal:::.gl_qq_plot_frame(z, no_split)
  qq_split <- gamlss.longitudinal:::.gl_qq_plot_frame(z, split_info)

  expect_equal(qq$observed, sort(z))
  expect_equal(qq$theoretical, stats::qnorm(stats::ppoints(length(z))))
  expect_true("split_group" %in% names(qq_split))
  expect_setequal(as.character(unique(qq_split$split_group)), c("a", "b"))
})

test_that("diagnostic worm helpers build bands and split frames", {
  pit <- c(0.1, 0.3, 0.6, 0.9)
  no_split <- list(split_by = FALSE, group = NULL)
  split_info <- list(split_by = TRUE, group = factor(c("a", "a", "b", "b")))

  band <- gamlss.longitudinal:::.gl_worm_band_frame(c(-1, 0, 1), n = 4)
  frames <- gamlss.longitudinal:::.gl_worm_plot_frames(pit, no_split)
  split_frames <- gamlss.longitudinal:::.gl_worm_plot_frames(pit, split_info)

  expect_equal(names(band), c("theoretical", "band_lower", "band_upper"))
  expect_true(all(band$band_lower <= 0))
  expect_true(all(band$band_upper >= 0))
  expect_true(all(c("worm_df", "worm_band") %in% names(frames)))
  expect_equal(nrow(frames$worm_df), length(pit))
  expect_equal(nrow(frames$worm_band), length(pit))
  expect_null(split_frames$worm_band)
  expect_true("split_group" %in% names(split_frames$worm_df))
  expect_setequal(as.character(unique(split_frames$worm_df$split_group)), c("a", "b"))
})

test_that("diagnostic rootogram helper builds observed and expected bin rows", {
  y <- c(-0.5, 0, 0.5)
  params <- list(mu = c(0, 0, 0), sigma = c(1, 1, 1))
  breaks <- c(-1, 0, 1)
  split_info <- list(split_by = TRUE, group = factor(c("left", "left", "right")))

  root <- gamlss.longitudinal:::.gl_rootogram_frame(
    y_i = y,
    params_i = params,
    breaks = breaks,
    family = "NO"
  )
  root_split <- gamlss.longitudinal:::.gl_rootogram_plot_frame(
    y = y,
    params = params,
    breaks = breaks,
    family = "NO",
    split_info = split_info
  )

  expect_equal(root$observed, c(1L, 2L))
  expect_true(all(root$expected >= 0))
  expect_equal(root$midpoint, c(-0.5, 0.5))
  expect_true("split_group" %in% names(root_split))
  expect_setequal(as.character(unique(root_split$split_group)), c("left", "right"))
})
