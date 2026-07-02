test_that("smooth-term variable parser handles calls names and malformed text", {
  expect_equal(.plot_smooth_terms_extract_var("s(log(age), bs = \"ps\")"), "log(age)")
  expect_equal(.plot_smooth_terms_extract_var("s(`visit time`)"), "visit time")
  expect_equal(.plot_smooth_terms_extract_var("not a call"), "not a call")
})

test_that("smooth-term x recovery prefers basis attributes and data alignment", {
  B <- matrix(1, nrow = 3, ncol = 1)
  attr(B, "smooth_x") <- c(10, 20, 30)
  attr(B, "smooth_var") <- "basis_age"
  object <- list(model_matrix = list(x = list(mu = matrix(nrow = 3, ncol = 0))))

  basis_info <- .plot_smooth_terms_x_info("mu", "s(age)", B, object)

  expect_equal(basis_info$x, c(10, 20, 30))
  expect_equal(basis_info$x_var, "basis_age")

  B2 <- matrix(1, nrow = 2, ncol = 1, dimnames = list(c("b", "a"), NULL))
  dat <- data.frame(age = c(1, 2, 3), row.names = c("a", "b", "c"))

  row_info <- .plot_smooth_terms_x_info("mu", "s(age)", B2, object, data = dat)

  expect_equal(row_info$x, c(2, 1))
  expect_equal(row_info$x_var, "age")
})

test_that("smooth-term x recovery evaluates transformed covariates and reports fallback", {
  B <- matrix(1, nrow = 3, ncol = 1)
  object <- list(model_matrix = list(x = list(mu = matrix(nrow = 3, ncol = 0))))
  dat <- data.frame(age = c(2, 4, 8))

  expr_info <- .plot_smooth_terms_x_info("mu", "s(log(age))", B, object, data = dat)
  expect_warning(
    fallback_info <- .plot_smooth_terms_x_info("mu", "s(missing)", B, object, fallback_to_index = TRUE),
    "Falling back to index"
  )

  expect_equal(expr_info$x, log(dat$age))
  expect_equal(expr_info$x_var, "log(age)")
  expect_equal(fallback_info$x, 1:3)
  expect_equal(fallback_info$x_var, "index")
  expect_error(
    .plot_smooth_terms_x_info("mu", "s(missing)", B, object, fallback_to_index = FALSE),
    "Could not infer x-axis"
  )
})

test_that("smooth-term index includes only available smooth basis and coefficients", {
  object <- list(
    par_s = list(mu = list(`s(age)` = 1:2, `s(time)` = NULL), sigma = list()),
    model_matrix = list(s = list(mu = list(`s(age)` = matrix(1, 3, 2), `s(time)` = matrix(1, 3, 1))))
  )

  index <- .plot_smooth_terms_index(object)

  expect_equal(length(index), 1L)
  expect_equal(index[[1]]$par_name, "mu")
  expect_equal(index[[1]]$s_name, "s(age)")
})

test_that("smooth-term fit SE uses covariance matrix coefficient SEs and missing fallback", {
  B <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
  smooth_vcov <- diag(c(0.04, 0.09))
  smooth_se <- c(0.2, 0.3)

  from_vcov <- .plot_smooth_terms_fit_se(B, smooth_vcov = smooth_vcov)
  from_se <- .plot_smooth_terms_fit_se(B, smooth_se = smooth_se)
  missing <- .plot_smooth_terms_fit_se(B)

  expect_equal(from_vcov, c(0.2, sqrt(0.13)))
  expect_equal(from_se, c(0.2, sqrt(0.13)))
  expect_true(all(is.na(missing)))
})

test_that("smooth-term plot data builder grids duplicate x values and falls back to sorted observations", {
  x <- c(2, 1, 1, 3)
  fitted <- c(2, 1, 3, 4)
  ci_lower <- fitted - 0.1
  ci_upper <- fitted + 0.1

  gridded <- .plot_smooth_terms_plot_df(x, fitted, ci_lower, ci_upper, even_grid = TRUE, grid_n = 5)
  observed <- .plot_smooth_terms_plot_df(x, fitted, ci_lower, ci_upper, even_grid = FALSE, sort_x = TRUE)

  expect_equal(nrow(gridded), 20L)
  expect_equal(range(gridded$x), c(1, 3))
  expect_equal(observed$x, c(1, 1, 2, 3))
  expect_equal(observed$fitted, c(1, 3, 2, 4))
})
