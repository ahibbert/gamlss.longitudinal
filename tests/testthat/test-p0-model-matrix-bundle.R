test_that("model matrix bundle resolves copula links and fixed matrices", {
  prepared <- gamlss.longitudinal:::.gl_prepare_fit_data(
    dataset = data.frame(
      id = c("a", "a", "b", "b"),
      visit = c(1, 2, 1, 2),
      y = c(1.0, 2.0, 1.5, 2.5),
      x = c(0, 1, 0, 1)
    ),
    time_var = "visit",
    subject_var = "id",
    mu.formula = y ~ x,
    sigma.formula = ~ 1,
    nu.formula = ~ 1,
    tau.formula = ~ 1,
    theta.formula = ~ 1,
    zeta.formula = ~ 1,
    verbose = 0
  )

  bundle <- gamlss.longitudinal:::.gl_build_model_matrix_bundle(
    formulas_int = prepared$formulas_int,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    dataset = prepared$dataset
  )

  expect_named(bundle, c("mm", "copula_link"))
  expect_true(all(c("x", "s") %in% names(bundle$mm)))
  expect_true(all(c("mu", "sigma", "theta") %in% names(bundle$mm$x)))
  expect_equal(nrow(bundle$mm$x$mu), nrow(prepared$dataset))
  expect_true("x" %in% colnames(bundle$mm$x$mu))
  expect_true(is.list(bundle$copula_link))
  expect_true("theta.linkinv" %in% names(bundle$copula_link))
})

test_that("model matrix rank warning helper preserves warning contract", {
  rank_deficient <- list(
    mu = cbind(
      `(Intercept)` = 1,
      x = c(1, 2, 3),
      x_duplicate = c(1, 2, 3)
    )
  )

  expect_warning(
    gamlss.longitudinal:::.gl_warn_rank_deficient_model_matrices(rank_deficient),
    "Fixed-effect model matrix for parameter 'mu' is rank deficient",
    fixed = TRUE
  )
})

test_that("model matrix proxy data preserves likelihood data while filling design inputs", {
  dat <- data.frame(
    response = c(1, NA, 3),
    time = c(1, 2, 3),
    subject = c(1, 1, 1),
    x = c(1, NA, 5),
    group = factor(c("b", NA, "a"), levels = c("a", "b"), ordered = TRUE),
    label = c("z", NA, "z"),
    stringsAsFactors = FALSE
  )

  out <- gamlss.longitudinal:::.gl_build_model_matrix_proxy_dataset(dat)

  expect_equal(out$response, c(1, 2, 3))
  expect_equal(out$x, c(1, 3, 5))
  expect_false(is.ordered(out$group))
  expect_equal(as.character(out$group), c("b", "a", "a"))
  expect_equal(out$label, c("z", "z", "z"))
  expect_identical(dat$response, c(1, NA, 3))
})

test_that("model matrix formula helpers preserve parameter and column naming policy", {
  expect_equal(
    gamlss.longitudinal:::.gl_normalize_time_covariate_colnames(
      c("time_covariate.L", "time_covariate.Q", "time_covariate^4", "x")
    ),
    c("time_covariate.1", "time_covariate.2", "time_covariate.4", "x")
  )

  fml <- gamlss.longitudinal:::.gl_to_response_formula("~ x + z")
  expect_s3_class(fml, "formula")
  expect_equal(deparse(fml), "response ~ x + z")

  expect_equal(
    gamlss.longitudinal:::.gl_model_matrix_included_parameters(gamlss.dist::NO(), "N"),
    c("mu", "sigma", "theta")
  )
  expect_equal(
    gamlss.longitudinal:::.gl_model_matrix_included_parameters(gamlss.dist::NO(), "t"),
    c("mu", "sigma", "theta", "zeta")
  )
})

test_that("gamlss2 matrix sanitization normalizes ordered factors and finite numeric inputs", {
  dat <- data.frame(
    response = c(1, 2, 3),
    x = c(1, Inf, NA),
    group = factor(c("low", "high", NA), levels = c("low", "mid", "high"), ordered = TRUE)
  )

  dropped <- gamlss.longitudinal:::.gl_sanitize_for_gamlss2(
    dat,
    response ~ x + group,
    preserve_factor_levels = FALSE
  )
  kept <- gamlss.longitudinal:::.gl_sanitize_for_gamlss2(
    dat,
    response ~ x + group,
    preserve_factor_levels = TRUE
  )

  expect_true(all(is.finite(dropped$x)))
  expect_false(is.ordered(dropped$group))
  expect_equal(levels(dropped$group), c("low", "high"))
  expect_equal(levels(kept$group), c("low", "mid", "high"))
})

test_that("model matrix helpers select parameter data and fixed columns", {
  dat <- data.frame(
    response = c(1, 2, 3, 4, 5, 6),
    subject = rep(c("a", "b"), each = 3),
    time = rep(1:3, times = 2),
    time_covariate = rep(1:3, times = 2),
    x = c(0, 1, 2, 3, 4, 5)
  )

  expect_equal(nrow(gamlss.longitudinal:::.gl_model_matrix_parameter_dataset(dat, "mu")), 6)
  expect_equal(nrow(gamlss.longitudinal:::.gl_model_matrix_parameter_dataset(dat, "theta")), 4)

  fixed <- gamlss.longitudinal:::.gl_build_fixed_model_matrix(
    response ~ time_covariate + x,
    dat
  )

  expect_equal(colnames(fixed), c("intercept", "time_covariate", "x"))
  expect_equal(attr(fixed, "term.labels"), c("time_covariate", "x"))
  expect_equal(nrow(fixed), nrow(dat))
})

test_that("model matrix smooth helper returns smooth basis metadata", {
  skip_if_not_installed("mgcv")
  dat <- data.frame(
    response = seq_len(12),
    x = seq(0, 1, length.out = 12)
  )
  smooth_eval_env <- new.env(parent = baseenv())
  smooth_eval_env$s <- mgcv::s

  smooth <- gamlss.longitudinal:::.gl_build_smooth_model_matrices(
    response ~ s(x, k = 5),
    dat,
    smooth_eval_env
  )

  expect_true("s(x, k = 5)" %in% names(smooth))
  expect_equal(nrow(smooth[["s(x, k = 5)"]]), nrow(dat))
  expect_equal(attr(smooth[["s(x, k = 5)"]], "smooth_var"), "x")
  expect_equal(attr(smooth[["s(x, k = 5)"]], "smooth_x"), dat$x)
  expect_true(is.matrix(attr(smooth[["s(x, k = 5)"]], "penalty")))
})
