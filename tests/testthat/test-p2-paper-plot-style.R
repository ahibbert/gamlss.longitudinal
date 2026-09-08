local_paper_plot_style_env <- function() {
  root <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE)
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "paper", "R", "paper-plot-style.R"), envir = env)
  env
}

test_that("paper figures use a stable semantic colour contract", {
  env <- local_paper_plot_style_env()
  colours <- env$jss_plot_colours()

  expect_identical(unname(colours[c("package", "comparator")]), c("#0072B2", "#6B7280"))
  expect_identical(unname(colours[["truth"]]), "#111827")
  expect_identical(unname(env$jss_method_colours(c("gamlss.longitudinal", "gamlss2"))),
    c("#0072B2", "#6B7280"))
  expect_false(any(tolower(colours) == "#e41a1c"))
})

test_that("paper figure dimensions and theme are centralized", {
  skip_if_not_installed("ggplot2")
  env <- local_paper_plot_style_env()

  expect_identical(env$jss_figure_dimensions("standard"), c(width = 7.1, height = 4.5))
  expect_identical(env$jss_figure_dimensions("shallow"), c(width = 7.1, height = 2.8))
  expect_identical(env$jss_figure_dimensions("dashboard"), c(width = 7.1, height = 5.5))
  theme <- env$theme_jss_paper(show_title = TRUE)
  expect_s3_class(theme, "theme")
  expect_identical(theme$legend.position, "top")
  expect_identical(theme$plot.title$face, "bold")
})

test_that("paper save helper writes a non-empty figure", {
  skip_if_not_installed("ggplot2")
  env <- local_paper_plot_style_env()
  path <- tempfile(fileext = ".png")
  plot <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_line(colour = env$jss_plot_colours()[["package"]]) +
    env$theme_jss_paper()

  saved <- env$jss_save_figure(path, plot, size = "shallow")
  expect_identical(saved, normalizePath(path, winslash = "/", mustWork = TRUE))
  expect_true(file.exists(path))
  expect_gt(file.info(path)$size, 0)
})
