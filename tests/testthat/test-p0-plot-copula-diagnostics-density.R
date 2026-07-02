test_that("copula diagnostics density grid is prepared for ungrouped overlays", {
  pair_data <- data.frame(
    u1 = c(0.2, 0.4, 0.6),
    u2 = c(0.3, 0.5, 0.7),
    theta_pair = c(0, 0, 0),
    zeta_pair = c(0, 0, 0)
  )

  grid <- .copula_v2_density_grid_for_plot(
    family_num = "N",
    pair_data_uniform = pair_data,
    grid_n = 4,
    max_pairs_overlay = 10,
    transform = "uniform",
    is_grouped = FALSE
  )

  expect_equal(nrow(grid), 16L)
  expect_setequal(names(grid), c("u1", "u2", "density"))
  expect_true(all(is.finite(grid$u1)))
  expect_true(all(is.finite(grid$u2)))
  expect_true(all(is.finite(grid$density)))
})

test_that("copula diagnostics density grid handles grouped overlays", {
  pair_data <- data.frame(
    u1 = c(0.2, 0.4, 0.6, 0.8),
    u2 = c(0.3, 0.5, 0.7, 0.9),
    theta_pair = c(0, 0, 0, 0),
    zeta_pair = c(0, 0, 0, 0),
    split_group = c("early", "early", "late", "late")
  )

  grid <- .copula_v2_density_grid_for_plot(
    family_num = "N",
    pair_data_uniform = pair_data,
    grid_n = 3,
    max_pairs_overlay = 10,
    transform = "uniform",
    is_grouped = TRUE
  )

  expect_equal(nrow(grid), 18L)
  expect_setequal(unique(grid$split_group), c("early", "late"))
  expect_true(all(is.finite(grid$density)))
})

test_that("copula diagnostics density grid can be transformed to normal scale", {
  pair_data <- data.frame(
    u1 = c(0.2, 0.4, 0.6),
    u2 = c(0.3, 0.5, 0.7),
    theta_pair = c(0, 0, 0),
    zeta_pair = c(0, 0, 0)
  )

  uniform_grid <- .copula_v2_density_grid_for_plot(
    family_num = "N",
    pair_data_uniform = pair_data,
    grid_n = 4,
    max_pairs_overlay = 10,
    transform = "uniform",
    is_grouped = FALSE
  )
  normal_grid <- .copula_v2_density_grid_for_plot(
    family_num = "N",
    pair_data_uniform = pair_data,
    grid_n = 4,
    max_pairs_overlay = 10,
    transform = "normal",
    is_grouped = FALSE
  )

  expect_equal(nrow(normal_grid), nrow(uniform_grid))
  expect_true(all(is.finite(normal_grid$u1)))
  expect_true(all(is.finite(normal_grid$u2)))
  expect_true(all(is.finite(normal_grid$density)))
  expect_false(isTRUE(all.equal(normal_grid$u1, uniform_grid$u1)))
})

test_that("copula diagnostics axis labels match the plotting transform", {
  uniform_labels <- .copula_v2_axis_labels("uniform")
  normal_labels <- .copula_v2_axis_labels("normal")

  expect_true(is.expression(uniform_labels$x))
  expect_true(is.expression(uniform_labels$y))
  expect_true(is.expression(normal_labels$x))
  expect_true(is.expression(normal_labels$y))
  expect_match(deparse(uniform_labels$x), "U\\[t\\]")
  expect_match(deparse(normal_labels$x), "Phi")
})

test_that("copula diagnostics cut summaries use fitted-tau ranks", {
  pair_data <- data.frame(
    u1 = c(0.1, 0.2, 0.3, 0.7, 0.8, 0.9),
    u2 = c(0.2, 0.1, 0.4, 0.6, 0.9, 0.8),
    tau_fit = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
  )

  summary <- .copula_v2_cut_summary(pair_data, plot2_cuts = 3)

  expect_equal(nrow(summary), 3L)
  expect_equal(summary$cut_group, c("C1", "C2", "C3"))
  expect_equal(summary$n_pairs, c(2L, 2L, 2L))
  expect_true(all(is.finite(summary$tau_fit)))
})

test_that("copula diagnostics cut summaries handle ties and grouped data", {
  pair_data <- data.frame(
    u1 = c(0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9),
    u2 = c(0.2, 0.3, 0.4, 0.1, 0.7, 0.8, 0.9, 0.6),
    tau_fit = c(0.2, 0.2, NA, Inf, 0.2, 0.2, 0.2, 0.2),
    split_group = c("early", "early", "early", "early", "late", "late", "late", "late")
  )

  summary <- .copula_v2_cut_summary_by_group(
    pair_data_plot = pair_data,
    plot2_cuts = 2,
    is_grouped = TRUE
  )

  expect_setequal(summary$split_group, c("early", "late"))
  expect_true(all(summary$cut_group %in% c("C1", "C2")))
  expect_equal(sum(summary$n_pairs), 6L)
})
