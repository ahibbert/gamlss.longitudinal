test_that("model-check statuses describe information and review scope", {
  expect_equal(gamlss.longitudinal:::.gl_check_status_severity("descriptive"), "information")
  expect_equal(gamlss.longitudinal:::.gl_check_status_severity("recorded"), "information")
  expect_equal(gamlss.longitudinal:::.gl_check_status_severity("flagged"), "review")
  expect_equal(gamlss.longitudinal:::.gl_check_status_severity("not_converged"), "concern")
})

test_that("dependence cutoff is optional and user-controlled", {
  expect_null(gamlss.longitudinal:::.gl_validate_dependence_cor_cutoff(NULL))
  expect_equal(gamlss.longitudinal:::.gl_validate_dependence_cor_cutoff("0.25"), 0.25)
  expect_equal(gamlss.longitudinal:::.gl_validate_dependence_cor_cutoff(0.999), 0.999)
  expect_error(gamlss.longitudinal:::.gl_validate_dependence_cor_cutoff(0), "between 0 and 1", fixed = TRUE)
  expect_error(gamlss.longitudinal:::.gl_validate_dependence_cor_cutoff(1), "between 0 and 1", fixed = TRUE)
  expect_error(gamlss.longitudinal:::.gl_validate_dependence_cor_cutoff(c(0.2, 0.3)), "between 0 and 1", fixed = TRUE)
})

test_that("model-check calibration helpers remain descriptive", {
  pit <- c(0.01, 0.05, 0.5, 0.95, 0.99)
  pit_stats <- gamlss.longitudinal:::.gl_pit_calibration_stats(pit)
  tail <- gamlss.longitudinal:::.gl_tail_calibration_stats(pit)

  expect_named(pit_stats, c("n", "mean", "sd", "expected_sd", "ks_p_value"))
  expect_equal(pit_stats$n, length(pit))
  expect_equal(pit_stats$expected_sd, sqrt(1 / 12))
  expect_s3_class(tail$tail_summary, "data.frame")
  expect_named(tail$tail_stats, "tail_ratio_max")
})

test_that("model-check table has no package-defined calibration verdicts", {
  checks <- gamlss.longitudinal:::.gl_check_table(
    summary_obj = list(convergence = list(converged = FALSE)),
    scores = list(),
    pit_stats = list(ks_p_value = 0.001),
    tail_stats = list(tail_ratio_max = 20),
    lag1_cor = 0.4,
    dependence_cor_cutoff = NULL,
    vcov_method = "numderiv"
  )

  expect_equal(nrow(checks), 5L)
  expect_false(any(checks$status %in% c("PASS", "FAIL")))
  expect_identical(checks$status[checks$area == "Convergence"], "not_converged")
  expect_true(all(
    checks$status[checks$area %in% c("Marginal fit", "Tail fit", "Copula fit")] == "descriptive"
  ))
  expect_identical(checks$status[checks$area == "Variance calculation"], "review")
})

test_that("only a supplied dependence threshold creates a flag", {
  make_checks <- function(cutoff) gamlss.longitudinal:::.gl_check_table(
    summary_obj = list(convergence = list(converged = TRUE)),
    scores = list(),
    pit_stats = list(ks_p_value = 0.5),
    tail_stats = list(tail_ratio_max = 1),
    lag1_cor = 0.4,
    dependence_cor_cutoff = cutoff,
    vcov_method = NA_character_
  )

  expect_identical(make_checks(NULL)$status[4], "descriptive")
  expect_identical(make_checks(0.25)$status[4], "flagged")
  expect_identical(make_checks(0.5)$status[4], "not_flagged")
})

test_that("diagnostic summary never claims that model fit passed", {
  expect_identical(
    gamlss.longitudinal:::.gl_basic_checks_result(data.frame(status = c("converged", "descriptive"))),
    "descriptive"
  )
  expect_identical(
    gamlss.longitudinal:::.gl_basic_checks_result(data.frame(status = c("converged", "flagged"))),
    "review"
  )
  expect_identical(
    gamlss.longitudinal:::.gl_basic_checks_result(data.frame(status = c("not_converged", "descriptive"))),
    "not_converged"
  )
})

test_that("randomized PIT helper preserves caller RNG state", {
  set.seed(77)
  before <- .Random.seed
  out1 <- gamlss.longitudinal:::.gl_with_preserved_rng(123, stats::runif(5))
  expect_identical(.Random.seed, before)
  out2 <- gamlss.longitudinal:::.gl_with_preserved_rng(123, stats::runif(5))
  expect_identical(out1, out2)
  expect_identical(.Random.seed, before)
})

test_that("discrete randomized PIT is seeded and reproducible", {
  object <- list(margin_dist = list(family = c("PO", "Poisson")))
  testthat::local_mocked_bindings(
    .gl_fitted_distribution = function(...) list(
      response = c(0, 1, 3),
      params = list(mu = c(1, 1, 1)),
      family = list(family = c("PO", "Poisson"))
    ),
    .gl_call_family_fun = function(kind, family, value, params) {
      stats::ppois(value, lambda = params$mu)
    },
    .package = "gamlss.longitudinal"
  )

  set.seed(99)
  before <- .Random.seed
  one <- gamlss.longitudinal:::.gl_pit(object, randomize = TRUE, seed = 2026)
  two <- gamlss.longitudinal:::.gl_pit(object, randomize = TRUE, seed = 2026)
  expect_identical(.Random.seed, before)
  expect_identical(one$pit, two$pit)
  expect_true(one$randomized)
  expect_identical(one$seed, 2026L)
})

test_that("discrete PIT defaults to seeded randomization and is calibrated under the generating model", {
  set.seed(831)
  y <- stats::rpois(5000, lambda = 2)
  object <- list(margin_dist = list(family = c("PO", "Poisson")))
  testthat::local_mocked_bindings(
    .gl_fitted_distribution = function(...) list(
      response = y,
      params = list(mu = rep(2, length(y))),
      family = list(family = c("PO", "Poisson")),
      subject = seq_along(y), time = rep(1, length(y))
    ),
    .gl_call_family_fun = function(kind, family, value, params) {
      stats::ppois(value, lambda = params$mu)
    },
    .package = "gamlss.longitudinal"
  )

  pit_out <- gamlss.longitudinal:::.gl_pit(object)
  expect_true(pit_out$randomized)
  expect_identical(pit_out$seed, 1L)
  expect_lt(abs(mean(pit_out$pit) - 0.5), 0.02)
  expect_gt(stats::ks.test(pit_out$pit, "punif")$p.value, 0.01)

  nonrandom <- gamlss.longitudinal:::.gl_pit(object, randomize = FALSE)
  expect_false(nonrandom$randomized)
  expect_gt(mean(nonrandom$pit), mean(pit_out$pit))
})
