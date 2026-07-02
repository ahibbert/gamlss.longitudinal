prepare_fit_data <- function(...) {
  gamlss.longitudinal:::.gl_prepare_fit_data(...)
}

test_that("fit panel lookup and grid helpers preserve deterministic ordering", {
  dataset <- data.frame(
    subject = c("b", "a", "a"),
    time = c(2, 1, 2),
    time_covariate = ordered(c("week2", "week1", "week2"), levels = c("week1", "week2")),
    response = c(3, 1, 2)
  )

  time_lookup <- gamlss.longitudinal:::.gl_fit_time_lookup(dataset)
  full_grid <- gamlss.longitudinal:::.gl_fit_full_panel_grid(dataset)

  expect_equal(time_lookup$time, c(1, 2))
  expect_equal(as.character(time_lookup$time_covariate), c("week1", "week2"))
  expect_equal(full_grid$subject, c("a", "b", "a", "b"))
  expect_equal(full_grid$time, c(1, 1, 2, 2))
})

test_that("fit data preparation expands panels and preserves formula time covariate", {
  dataset <- data.frame(
    subject_id = c("b", "a", "a"),
    visit = ordered(c("week2", "week1", "week2"), levels = c("week1", "week2")),
    y = c(3, 1, 2),
    x = c(30, 10, 20)
  )

  prepared <- prepare_fit_data(
    dataset = dataset,
    time_var = "visit",
    subject_var = "subject_id",
    mu.formula = y ~ visit + x,
    sigma.formula = ~ 1,
    nu.formula = ~ 1,
    tau.formula = ~ 1,
    theta.formula = ~ 1,
    zeta.formula = ~ 1,
    verbose = 0
  )

  expect_equal(prepared$response_var, "y")
  expect_equal(prepared$var_map[["visit"]], "time")
  expect_equal(prepared$var_map[["subject_id"]], "subject")
  expect_equal(prepared$var_map[["y"]], "response")

  expanded <- prepared$dataset
  expect_equal(nrow(expanded), 4)
  expect_equal(names(expanded)[1:3], c("subject", "time", "response"))
  expect_equal(expanded$subject, c("a", "a", "b", "b"))
  expect_equal(expanded$time, c(1, 2, 1, 2))
  expect_equal(sum(is.na(expanded$response)), 1)
  expect_true(is.ordered(expanded$time_covariate))

  expect_true("time_covariate" %in% all.vars(prepared$formulas_int$mu))
  expect_false("visit" %in% all.vars(prepared$formulas_int$mu))
  expect_equal(prepared$miss_by_time$n_observed_response, c(1, 2))
  expect_equal(prepared$pair_summary$complete_pairs, 1)
})

test_that("fit data preparation keeps duplicate subject/time error contract", {
  dataset <- data.frame(
    subject_id = c("a", "a"),
    visit = c(1, 1),
    y = c(1, 2)
  )

  expect_error(
    prepare_fit_data(
      dataset = dataset,
      time_var = "visit",
      subject_var = "subject_id",
      mu.formula = y ~ visit,
      sigma.formula = ~ 1,
      nu.formula = ~ 1,
      tau.formula = ~ 1,
      theta.formula = ~ 1,
      zeta.formula = ~ 1,
      verbose = 0
    ),
    "Duplicate subject/time combinations found",
    fixed = TRUE
  )
})

test_that("fit input column normalization records internal names and formula map", {
  dataset <- data.frame(
    id = c("a", "a"),
    visit = ordered(c("week1", "week2"), levels = c("week1", "week2")),
    y = c(1, 2),
    x = c(3, 4)
  )

  normalized <- gamlss.longitudinal:::.gl_normalize_fit_input_columns(
    dataset = dataset,
    time_var = "visit",
    subject_var = "id",
    mu.formula = y ~ visit + x,
    verbose = 0
  )

  expect_equal(normalized$response_var, "y")
  expect_equal(names(normalized$dataset)[1:3], c("subject", "time", "response"))
  expect_equal(normalized$var_map[["visit"]], "time")
  expect_equal(normalized$formula_var_map[["visit"]], "time_covariate")
  expect_true(normalized$time_covariate_is_factor)
  expect_true(normalized$time_covariate_ordered)
  expect_equal(normalized$time_covariate_levels, c("week1", "week2"))
})

test_that("fit input column normalization keeps required-column error contracts", {
  dataset <- data.frame(id = 1, visit = 1, y = 1)

  expect_error(
    gamlss.longitudinal:::.gl_normalize_fit_input_columns(
      dataset,
      time_var = "missing_visit",
      subject_var = "id",
      mu.formula = y ~ visit,
      verbose = 0
    ),
    "time_var='missing_visit' not found",
    fixed = TRUE
  )
  expect_error(
    gamlss.longitudinal:::.gl_normalize_fit_input_columns(
      dataset,
      time_var = "visit",
      subject_var = "missing_id",
      mu.formula = y ~ visit,
      verbose = 0
    ),
    "subject_var='missing_id' not found",
    fixed = TRUE
  )
  expect_error(
    gamlss.longitudinal:::.gl_normalize_fit_input_columns(
      dataset,
      time_var = "visit",
      subject_var = "id",
      mu.formula = missing_y ~ visit,
      verbose = 0
    ),
    "response variable 'missing_y' not found",
    fixed = TRUE
  )
})
