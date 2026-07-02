test_that("plot term helpers clean expression names and detect factor terms", {
  dat <- data.frame(
    sex = factor(c("A", "B")),
    visit_time = factor(c("1", "2")),
    stringsAsFactors = FALSE
  )

  expect_equal(.plot_terms_clean_expr_name("as.factor(`sex`)"), "sex")
  expect_true(.plot_terms_is_factor_term("as.factor(age)", data_for_terms = NULL))
  expect_true(.plot_terms_is_factor_term("sex", data_for_terms = dat))
  expect_true(.plot_terms_is_factor_term("time_covariate", data_for_terms = dat))
  expect_false(.plot_terms_is_factor_term("age", data_for_terms = dat))
})

test_that("plot term count helper groups factor terms and respects intercept interaction flags", {
  X <- matrix(
    c(
      1, 20, 0, 0,
      1, 30, 1, 30,
      1, 40, 0, 0
    ),
    nrow = 3,
    byrow = TRUE
  )
  colnames(X) <- c("intercept", "age", "sexB", "age:sexB")
  attr(X, "assign") <- c(0L, 1L, 2L, 3L)
  attr(X, "term.labels") <- c("age", "sex", "age:sex")

  object <- list(
    par_s = list(mu = list(`s(age)` = 1:2), sigma = list()),
    model_matrix = list(x = list(mu = X)),
    par = c(mu.intercept = 1, mu.age = 0.1, mu.sexB = 0.2, `mu.age:sexB` = 0.3),
    dataset = data.frame(age = c(20, 30, 40), sex = factor(c("A", "B", "A")))
  )

  default_counts <- .plot_terms_count(object)
  full_counts <- .plot_terms_count(
    object,
    include_intercept = TRUE,
    plot_interactions = TRUE
  )

  expect_equal(default_counts, list(smooth = 1, fixed = 2, total = 3))
  expect_equal(full_counts, list(smooth = 1, fixed = 4, total = 5))
})

test_that("plot term dashboard helpers preserve plot order and layout metadata", {
  smooth_results <- list(plots = list("smooth1", "smooth2"))
  fixed_results <- list(plots = list("fixed1"))

  plots <- .plot_terms_collect_plot_objects(smooth_results, fixed_results)
  default_layout <- .plot_terms_dashboard_layout(plots, ncol = NULL)
  explicit_layout <- .plot_terms_dashboard_layout(plots, ncol = 3)
  paginated_layout <- .plot_terms_dashboard_layout(plots, paginate = TRUE)

  expect_identical(plots, list("smooth1", "smooth2", "fixed1"))
  expect_identical(default_layout$ncol, 2)
  expect_equal(default_layout$nrow, 2)
  expect_identical(explicit_layout$ncol, 3)
  expect_equal(explicit_layout$nrow, 1)
  expect_true(paginated_layout$paginate)
  expect_null(.plot_terms_dashboard_layout(list()))
})
