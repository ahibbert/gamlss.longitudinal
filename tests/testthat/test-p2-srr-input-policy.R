test_that("srr input policy rejects empty, list, and matrix-like data", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  dat <- make_fixture_factor_time_interaction(n_subject = 8L)

  expect_error(
    fit_fixture_model(dat[0, , drop = FALSE]),
    "dataset must contain at least one row",
    fixed = TRUE
  )

  dat_list <- dat
  dat_list$notes <- lapply(seq_len(nrow(dat_list)), function(i) list(i))
  expect_error(
    fit_fixture_model(dat_list),
    "unsupported list-column",
    fixed = TRUE
  )

  dat_matrix <- dat
  dat_matrix$matrix_col <- I(matrix(seq_len(2L * nrow(dat_matrix)), ncol = 2L))
  expect_error(
    fit_fixture_model(dat_matrix),
    "unsupported matrix/data-frame column",
    fixed = TRUE
  )
})

test_that("srr input policy rejects non-standard predictor classes", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  dat <- make_fixture_factor_time_interaction(n_subject = 8L)
  class(dat$age) <- c("review_units", class(dat$age))

  expect_error(
    fit_fixture_model(dat),
    "predictor column(s) have unsupported classes: age",
    fixed = TRUE
  )
})

test_that("srr input policy distinguishes response and predictor missingness", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  dat <- make_fixture_factor_time_interaction(n_subject = 8L)

  dat_resp_inf <- dat
  dat_resp_inf$y[1] <- Inf
  expect_error(
    fit_fixture_model(dat_resp_inf),
    "response variable contains NaN or Inf",
    fixed = TRUE
  )

  dat_resp_nan <- dat
  dat_resp_nan$y[1] <- NaN
  expect_error(
    fit_fixture_model(dat_resp_nan),
    "response variable contains NaN or Inf",
    fixed = TRUE
  )

  dat_pred_na <- dat
  dat_pred_na$age[1] <- NA_real_
  expect_error(
    fit_fixture_model(dat_pred_na),
    "predictor column(s) contain missing values",
    fixed = TRUE
  )

  dat_pred_inf <- dat
  dat_pred_inf$age[1] <- Inf
  expect_error(
    fit_fixture_model(dat_pred_inf),
    "predictor column(s) contain NaN or Inf",
    fixed = TRUE
  )

  structural <- make_fixture_with_structural_missing_rows()
  fit <- fit_fixture_model(
    structural,
    include_dlcopdpar = FALSE,
    max_outer_iter = 2,
    theta_formula = "~ 1"
  )
  expect_true(any(is.na(stats::model.frame(fit, type = "expanded")$response)))
  expect_false(any(is.na(stats::model.frame(fit, type = "observed")$response)))
})

test_that("srr input policy documents factor and character-time conversions", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  dat <- make_fixture_factor_time_interaction(n_subject = 8L)
  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = FALSE,
    max_outer_iter = 2,
    theta_formula = "~ 1"
  )
  expect_true(is.ordered(stats::model.frame(fit, type = "submitted")$time_raw))
  expect_false(any(grepl("time_covariate\\.L", colnames(fit$model_matrix$x$mu))))

  dat_chr <- make_fixture_numeric_time(n_subject = 8L)
  dat_chr$time_raw <- as.character(dat_chr$time_raw)
  expect_warning(
    fit_chr <- fit_fixture_model(
      dat_chr,
      include_dlcopdpar = FALSE,
      max_outer_iter = 2,
      theta_formula = "~ 1"
    ),
    "Converted character time variable",
    fixed = TRUE
  )
  expect_type(stats::model.frame(fit_chr, type = "expanded")$time_covariate, "double")

  dat_char_pred <- dat
  dat_char_pred$clinic <- rep(c("A", "B"), length.out = nrow(dat_char_pred))
  expect_warning(
    fit_fixture_model(
      dat_char_pred,
      include_dlcopdpar = FALSE,
      mu_formula = "y ~ clinic + age",
      sigma_formula = "~ 1",
      theta_formula = "~ 1",
      max_outer_iter = 2
    ),
    "Character predictor column",
    fixed = TRUE
  )
})
