test_that("srr error map covers core fitting input errors", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  dat <- make_fixture_factor_time_interaction(n_subject = 8L)

  expect_error(fit_fixture_model(dat[0, , drop = FALSE]), "at least one row", fixed = TRUE)
  expect_error(fit_fixture_model(dat, time_var = "missing_time"), "`time_var = \"missing_time\"` was not found in `data`.", fixed = TRUE)
  expect_error(fit_fixture_model(dat, subject_var = "missing_subject"), "`subject_var = \"missing_subject\"` was not found in `data`.", fixed = TRUE)
  expect_error(fit_fixture_model(dat, mu_formula = "missing_response ~ age"), "Response variable `missing_response` was not found in `data`.", fixed = TRUE)
  expect_error(fit_fixture_model(make_fixture_with_duplicate_subject_time()), "Duplicate subject/time combinations were found.", fixed = TRUE)
  expect_error(fit_fixture_model(make_fixture_missingness_margin_all_missing()), "All response values are missing at time point(s): 3", fixed = TRUE)
  expect_error(fit_fixture_model(make_fixture_missingness_zero_complete_pairs()), "No complete copula pairs are available", fixed = TRUE)
})

test_that("srr error map covers value conversion and class diagnostics", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  dat <- make_fixture_factor_time_interaction(n_subject = 8L)

  dat_bad_time <- dat
  dat_bad_time$time_raw <- as.character(dat_bad_time$time_raw)
  expect_error(
    fit_fixture_model(dat_bad_time),
    "`time_var` must be numeric-like unless supplied as a factor.",
    fixed = TRUE
  )

  dat_list <- dat
  dat_list$bad <- lapply(seq_len(nrow(dat_list)), identity)
  expect_error(fit_fixture_model(dat_list), "unsupported list-column", fixed = TRUE)

  dat_char_time <- make_fixture_numeric_time(n_subject = 8L)
  dat_char_time$time_raw <- as.character(dat_char_time$time_raw)
  expect_warning(
    fit_fixture_model(dat_char_time, include_dlcopdpar = FALSE, max_outer_iter = 2, theta_formula = "~ 1"),
    "Converted character time variable",
    fixed = TRUE
  )
})

test_that("srr error map covers prediction newdata diagnostics", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  dat <- make_fixture_factor_time_interaction(n_subject = 8L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = FALSE, max_outer_iter = 2, theta_formula = "~ 1")

  nd <- dat[seq_len(3), , drop = FALSE]
  nd$gender <- factor("other", levels = c(levels(dat$gender), "other"))
  expect_error(
    suppressWarnings(stats::predict(fit, newdata = nd)),
    "contains level(s) not seen during fitting",
    fixed = TRUE
  )
})
