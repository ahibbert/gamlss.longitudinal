test_that("copula time summary print method displays both summary sections", {
  x <- list(
    time_summary = data.frame(time = 1:2, n_obs = c(3L, 3L), theta_fit = c(0.2, 0.3), tau_fit = c(0.1, 0.2)),
    pair_summary = data.frame(time_pair = c("1-2"), n_pairs = 3L, theta_pair = 0.25, tau_pair = 0.15)
  )
  class(x) <- "copula_time_summary"

  printed <- utils::capture.output(out <- print(x))

  expect_identical(out, x)
  expect_true(any(grepl("Copula Dependence Summary", printed)))
  expect_true(any(grepl("Fitted dependence by time", printed)))
  expect_true(any(grepl("Adjacent-pair dependence", printed)))
})

test_that("copula summary aggregation helper selects requested statistic", {
  expect_identical(.gl_copula_summary_agg_fun("median"), stats::median)
  expect_identical(.gl_copula_summary_agg_fun("mean"), mean)
})

test_that("copula time table aggregates and sorts fitted summaries", {
  fit_data <- data.frame(
    time = c(2, 1, 2, 1),
    theta_fit = c(0.5, 0.1, 0.7, 0.3),
    tau_fit = c(0.4, 0.05, 0.6, 0.15),
    zeta_fit = c(4, 2, 6, 8)
  )

  no_zeta <- .gl_copula_time_table(fit_data, mean, has_zeta = FALSE)
  with_zeta <- .gl_copula_time_table(fit_data, stats::median, has_zeta = TRUE)

  expect_identical(no_zeta$time, c(1, 2))
  expect_identical(no_zeta$n_obs, c(2L, 2L))
  expect_equal(no_zeta$theta_fit, c(0.2, 0.6))
  expect_false("zeta_fit" %in% names(no_zeta))
  expect_equal(with_zeta$zeta_fit, c(5, 5))
})

test_that("copula pair table aggregates adjacent-pair summaries", {
  pair_data <- data.frame(
    time_pair = c("1-2", "1-2", "2-3"),
    theta_pair = c(0.2, 0.4, 0.8),
    tau_fit = c(0.1, 0.3, 0.7),
    zeta_pair = c(5, 7, 9)
  )

  pair_summary <- .gl_copula_pair_table(pair_data, mean, has_zeta = TRUE)

  expect_identical(pair_summary$time_pair, c("1-2", "2-3"))
  expect_identical(pair_summary$n_pairs, c(2L, 1L))
  expect_equal(pair_summary$theta_pair, c(0.3, 0.8))
  expect_equal(pair_summary$tau_pair, c(0.2, 0.7))
  expect_equal(pair_summary$zeta_pair, c(6, 9))
})

test_that("copula summary data tidying drops zeta only for one-parameter copulas", {
  fit_data <- data.frame(time = 1, zeta_fit = 3)
  pair_data <- data.frame(time_pair = "1-2", zeta_pair = 4)

  dropped <- .gl_copula_drop_zeta_if_absent(fit_data, pair_data, has_zeta = FALSE)
  kept <- .gl_copula_drop_zeta_if_absent(fit_data, pair_data, has_zeta = TRUE)

  expect_false("zeta_fit" %in% names(dropped$fit_data))
  expect_false("zeta_pair" %in% names(dropped$pair_data))
  expect_true("zeta_fit" %in% names(kept$fit_data))
  expect_true("zeta_pair" %in% names(kept$pair_data))
})

test_that("copula time summary plot method returns component plots", {
  skip_if_not_installed("ggpubr")

  x <- list(
    time_summary = data.frame(
      time = 1:3,
      n_obs = c(4L, 4L, 4L),
      theta_fit = c(0.2, 0.3, 0.35),
      tau_fit = c(0.1, 0.18, 0.2),
      zeta_fit = c(5, 6, 7)
    ),
    pair_summary = data.frame(time_pair = c("1-2", "2-3"), n_pairs = c(4L, 4L), theta_pair = c(0.25, 0.32), tau_pair = c(0.14, 0.19))
  )
  class(x) <- "copula_time_summary"

  plotted <- plot(x, plot = FALSE)

  expect_identical(plotted$summary, x)
  expect_s3_class(plotted$p_theta, "ggplot")
  expect_s3_class(plotted$p_tau, "ggplot")
  expect_s3_class(plotted$p_zeta, "ggplot")
  expect_false(is.null(plotted$dashboard))
})

test_that("copula time summary plot method rejects empty summaries", {
  x <- list(time_summary = data.frame(), pair_summary = data.frame())
  class(x) <- "copula_time_summary"

  expect_error(plot(x, plot = FALSE), "No fitted copula summaries")
})
