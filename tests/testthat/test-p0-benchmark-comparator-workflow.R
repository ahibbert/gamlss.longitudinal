test_that("benchmark comparator input helper normalizes core arguments", {
  dat <- data.frame(id = c("a", "b"), y = c(1, 2))

  inputs <- .benchmark_standard_models_inputs(
    data = dat,
    formula = y ~ 1,
    subject_var = "id",
    family = "gaussian",
    comparators = c("glm", "glm", "gam"),
    add_subject_re_to_gamm = TRUE,
    add_subject_re_to_gam = NULL,
    fit = NULL,
    fit_name = "primary",
    quantile_prob = "0.8",
    interval_level = "0.9"
  )

  expect_s3_class(inputs$formula, "formula")
  expect_equal(inputs$comparators, c("glm", "gam"))
  expect_equal(inputs$quantile_prob, 0.8)
  expect_equal(inputs$interval_level, 0.9)
  expect_equal(inputs$family$family, "gaussian")
})

test_that("benchmark comparator input helper validates reviewer-facing errors", {
  dat <- data.frame(id = c("a", "b"), y = c(1, 2))
  base_args <- list(
    data = dat,
    formula = y ~ 1,
    subject_var = "id",
    family = "gaussian",
    comparators = "glm",
    add_subject_re_to_gamm = TRUE,
    add_subject_re_to_gam = NULL,
    fit = NULL,
    fit_name = "primary",
    quantile_prob = 0.8,
    interval_level = 0.9
  )

  bad_subject <- base_args
  bad_subject$subject_var <- "missing"
  expect_error(do.call(.benchmark_standard_models_inputs, bad_subject), "subject_var")

  bad_quantile <- base_args
  bad_quantile$quantile_prob <- 1
  expect_error(do.call(.benchmark_standard_models_inputs, bad_quantile), "quantile_prob")

  bad_interval <- base_args
  bad_interval$interval_level <- 0
  expect_error(do.call(.benchmark_standard_models_inputs, bad_interval), "interval_level")

  bad_comparator <- base_args
  bad_comparator$comparators <- "not_a_model"
  expect_error(do.call(.benchmark_standard_models_inputs, bad_comparator), "Unknown comparator")

  duplicate_fit <- base_args
  duplicate_fit$fit <- stats::lm(y ~ 1, data = dat)
  duplicate_fit$fit_name <- "glm"
  expect_error(do.call(.benchmark_standard_models_inputs, duplicate_fit), "must not duplicate")
})

test_that("benchmark comparator input helper preserves deprecated GAM alias behavior", {
  dat <- data.frame(id = c("a", "b"), y = c(1, 2))

  expect_warning(
    inputs <- .benchmark_standard_models_inputs(
      data = dat,
      formula = y ~ 1,
      subject_var = "id",
      family = "gaussian",
      comparators = "gam",
      add_subject_re_to_gamm = TRUE,
      add_subject_re_to_gam = FALSE,
      fit = NULL,
      fit_name = "primary",
      quantile_prob = 0.8,
      interval_level = 0.9
    ),
    "deprecated"
  )

  expect_false(inputs$add_subject_re_to_gamm)
})

test_that("benchmark comparator fit dispatcher runs glm with extra arguments", {
  dat <- data.frame(id = c(1, 1, 2, 2), y = c(0, 1, 0, 1), x = c(0, 1, 0, 1))

  fit <- .benchmark_run_comparator_fit(
    data = dat,
    formula = y ~ x,
    subject_var = "id",
    family = stats::binomial(),
    comparator = "glm",
    correlation = "exchangeable",
    add_subject_re_to_gamm = TRUE,
    extra_args = list(control = stats::glm.control(maxit = 50))
  )

  expect_s3_class(fit, "glm")
  expect_equal(fit$control$maxit, 50)
  expect_equal(stats::formula(fit), y ~ x)
})

test_that("benchmark comparator print method reports core reviewer sections", {
  bench <- structure(
    list(
      formula = y ~ x,
      subject_var = "id",
      family = "gaussian",
      results = data.frame(
        method = "glm",
        success = TRUE,
        mae = 0.1,
        rmse = 0.2,
        stringsAsFactors = FALSE
      ),
      coefficients = data.frame(),
      interpretation = structure(list(), class = "gamlss_longitudinal_benchmark_interpretation")
    ),
    class = "gamlss_longitudinal_benchmark"
  )

  out <- utils::capture.output(returned <- print(bench))

  expect_identical(returned, bench)
  expect_true(any(grepl("Standard Longitudinal Comparator Benchmark", out, fixed = TRUE)))
  expect_true(any(grepl("Benchmark Details", out, fixed = TRUE)))
  expect_true(any(grepl("Estimand note", out, fixed = TRUE)))
})
