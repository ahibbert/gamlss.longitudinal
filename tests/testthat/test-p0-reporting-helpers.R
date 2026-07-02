predict.fake_reporting_fit <- function(object, newdata, type, probs = NULL, q = NULL, direction = NULL, ...) {
  n <- nrow(newdata)
  switch(
    type,
    mean = seq_len(n),
    mu = seq_len(n) + 10,
    median = seq_len(n) + 20,
    quantile = data.frame(
      subject = newdata$subject,
      time = newdata$time,
      response = NA_real_,
      q10 = seq_len(n) + 30,
      q90 = seq_len(n) + 40
    ),
    probability = data.frame(probability = seq_len(n) / 10),
    stop("unexpected prediction type", call. = FALSE)
  )
}
registerS3method("predict", "fake_reporting_fit", predict.fake_reporting_fit, envir = asNamespace("stats"))

test_that("model_spec builds printable fitted-model audit metadata", {
  object <- list(
    response = c(1, NA, Inf),
    response_var = "y",
    subject_var = "id",
    time_var = "visit",
    margin_dist = list(family = c("NO", "Normal"), mu.link = "identity", sigma.link = "log"),
    copula_dist = "N",
    formulas = list(mu = y ~ x, sigma = ~ 1),
    optim_method = "RS",
    convergence = list(converged = TRUE, stop_reason = "tolerance", outer_iterations = 3L, max_outer_iter = 20L),
    vcov_meta = list(precomputed = FALSE),
    vcov = list(hessian_diagnostics = list(rank = 2L))
  )
  class(object) <- "gamlss.longitudinal"

  spec <- model_spec(object)
  printed <- utils::capture.output(out <- print(spec))

  expect_s3_class(spec, "gamlss_longitudinal_model_spec")
  expect_identical(out, spec)
  expect_identical(spec$variables$response, "y")
  expect_identical(spec$distributions$margin, "NO")
  expect_identical(spec$optimisation$method, "RS")
  expect_identical(spec$missingness$n_missing_response, 1L)
  expect_identical(spec$missingness$n_nonfinite_response, 1L)
  expect_equal(spec$margin_links$parameter, c("mu", "sigma"))
  expect_true(any(grepl("GAMLSS Longitudinal Model Specification", printed)))
})

test_that("model_spec rejects non-fitted inputs", {
  expect_error(model_spec(list()), "fitted 'gamlss.longitudinal' object")
})

test_that("reporting_table summarizes predictions overall and by group", {
  object <- structure(list(), class = c("fake_reporting_fit", "gamlss.longitudinal"))
  newdata <- data.frame(
    subject = 1:4,
    time = c(1, 1, 2, 2),
    group = c("a", "a", "b", "b")
  )

  overall <- reporting_table(object, newdata = newdata, probs = c(0.1, 0.9), threshold = 5, direction = "above")
  grouped <- reporting_table(object, newdata = newdata, by = "group", probs = c(0.1, 0.9))

  expect_identical(overall$n, 4L)
  expect_equal(overall$mean, 2.5)
  expect_equal(overall$mu, 12.5)
  expect_equal(overall$median, 22.5)
  expect_equal(overall$q10, 32.5)
  expect_equal(overall$q90, 42.5)
  expect_equal(overall$prob_above_5, 0.25)
  expect_identical(grouped$group, c("a", "b"))
  expect_identical(grouped$n, c(2L, 2L))
  expect_equal(grouped$mean, c(1.5, 3.5))
})

test_that("reporting_table validates inputs before prediction", {
  object <- structure(list(), class = c("fake_reporting_fit", "gamlss.longitudinal"))

  expect_error(reporting_table(list(), data.frame(x = 1)), "fitted 'gamlss.longitudinal' object")
  expect_error(reporting_table(object), "'newdata' is required")
  expect_error(reporting_table(object, data.frame(x = 1), by = "missing"), "'by' column")
})
