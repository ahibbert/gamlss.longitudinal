test_that("intermittent gaps error by default with an explicit remedy", {
  dat <- make_fixture_with_structural_missing_rows()

  err <- tryCatch(
    fit_fixture_model(dat, include_dlcopdpar = FALSE),
    error = identity
  )

  expect_s3_class(err, "gamlss.longitudinal_gap_error")
  expect_match(conditionMessage(err), 'missingness = "segment"', fixed = TRUE)
  expect_identical(err$n_subjects_with_gaps, 2L)
  expect_identical(err$n_gap_runs, 2L)
})

test_that("complete data and monotone boundary missingness use the ordinary objective", {
  complete <- data.frame(
    response = 1:6,
    time = rep(1:3, 2),
    subject = rep(c("a", "b"), each = 3)
  )
  dropout <- complete
  dropout$response[dropout$subject == "a" & dropout$time == 3] <- NA_real_
  dropout$response[dropout$subject == "b" & dropout$time == 1] <- NA_real_

  complete_contract <- gamlss.longitudinal:::.gl_missingness_contract(complete)
  dropout_contract <- gamlss.longitudinal:::.gl_missingness_contract(dropout)

  expect_false(complete_contract$has_intermittent_gaps)
  expect_false(dropout_contract$has_intermittent_gaps)
  expect_identical(dropout_contract$objective, "ordinary")
  expect_identical(dropout_contract$n_segments, 2L)
})

test_that("gap detection follows scheduled-time order and records contiguous segments", {
  dat <- data.frame(
    response = c(3, 1, NA, 4, 2, 6),
    time = c(5, 0, 2, 5, 0, 2),
    subject = c("a", "a", "a", "b", "b", "b")
  )

  contract <- gamlss.longitudinal:::.gl_missingness_contract(dat, "segment")

  expect_true(contract$has_intermittent_gaps)
  expect_identical(contract$n_subjects_with_gaps, 1L)
  expect_identical(contract$n_interior_missing, 1L)
  expect_identical(contract$n_gap_runs, 1L)
  expect_identical(contract$n_segments, 3L)
  expect_identical(contract$n_within_segment_pairs, 2L)
  expect_identical(contract$n_between_segment_transitions_omitted, 1L)
  expect_identical(contract$affected_subjects, "a")
  expect_identical(contract$between_segment_assumption, "independent")
})

test_that("copula pair cache never joins observations across a gap", {
  cache <- gamlss.longitudinal:::build_copula_pair_cache(
    response = c(1, NA_real_, 3),
    response_margin = c(1, 2, 3),
    response_subject = c("a", "a", "a")
  )

  expect_identical(cache$segment_id, c(1L, NA_integer_, 2L))
  expect_identical(cache$observed_pair_base, c(FALSE, FALSE))
  expect_false(any(cache$pair_segment1 == cache$pair_segment2, na.rm = TRUE))
})

test_that("segmented fits warn and retain an auditable likelihood contract", {
  dat <- make_fixture_with_structural_missing_rows()
  warnings <- character()

  fit <- withCallingHandlers(
    fit_fixture_model(
      dat,
      missingness = "segment",
      muffle_segment_warning = FALSE,
      include_dlcopdpar = FALSE,
      theta_formula = "~ 1"
    ),
    warning = function(w) {
      if (inherits(w, "gamlss.longitudinal_segment_warning")) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    }
  )

  expect_length(warnings, 1L)
  expect_match(warnings, "treated as independent", fixed = TRUE)
  expect_identical(fit$missingness$objective, "segmented")
  expect_identical(fit$missingness$n_subjects_with_gaps, 2L)
  expect_identical(fit$missingness$n_segments, 18L)
  expect_identical(fit$missingness$n_observed, 46L)
  expect_identical(fit$missingness$n_within_segment_pairs, 28L)
  expect_identical(fit$likelihood_contract$criteria_status, "available_under_segment_independence")

  summary_out <- summary(fit, include_vcov = FALSE)
  expect_true(is.finite(summary_out$fit$AIC))
  expect_true(is.finite(summary_out$fit$BIC))
  expect_identical(summary_out$fit$criteria_status, "provisional_nonconverged")
  expect_identical(
    summary_out$fit$likelihood_contract$criteria_status,
    "available_under_segment_independence"
  )
  expect_match(
    paste(capture.output(print(summary_out)), collapse = "\n"),
    "Between-gap assumption",
    fixed = TRUE
  )

  spec <- model_spec(fit)
  expect_identical(spec$likelihood_contract$objective, "segmented")
  expect_match(
    paste(capture.output(print(spec)), collapse = "\n"),
    "AIC/BIC use this likelihood",
    fixed = TRUE
  )
})
