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

test_that("complete data and monotone dropout use the ordinary objective", {
  complete <- data.frame(
    response = 1:6,
    time = rep(1:3, 2),
    subject = rep(c("a", "b"), each = 3)
  )
  dropout <- complete
  dropout$response[dropout$subject == "a" & dropout$time == 3] <- NA_real_

  complete_contract <- gamlss.longitudinal:::.gl_missingness_contract(complete)
  dropout_contract <- gamlss.longitudinal:::.gl_missingness_contract(dropout)

  expect_false(complete_contract$has_intermittent_gaps)
  expect_false(dropout_contract$has_intermittent_gaps)
  expect_identical(dropout_contract$objective, "ordinary")
  expect_identical(dropout_contract$n_segments, 2L)
})

test_that("leading unobserved visits require an explicit segmented contract", {
  delayed <- data.frame(
    response = c(NA, 2, 3, 4, 5, 6),
    time = rep(1:3, 2),
    subject = rep(c("a", "b"), each = 3)
  )

  err <- tryCatch(
    gamlss.longitudinal:::.gl_missingness_contract(delayed),
    error = identity
  )
  expect_s3_class(err, "gamlss.longitudinal_gap_error")
  expect_match(conditionMessage(err), "leading unobserved scheduled visits", fixed = TRUE)

  contract <- gamlss.longitudinal:::.gl_missingness_contract(delayed, "segment")
  expect_true(contract$has_delayed_entry)
  expect_false(contract$has_intermittent_gaps)
  expect_identical(contract$n_subjects_with_delayed_entry, 1L)
  expect_identical(contract$n_leading_missing, 1L)
  expect_identical(contract$objective, "segmented")
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

test_that("six-row panel fixes contribution and left-row dependence semantics", {
  skip_if_not_installed("gamlss.dist")
  panel <- data.frame(
    response = c(0.2, NA, 0.8, -0.5, 0.1, 0.6),
    time = rep(c(0, 2, 5), 2),
    subject = rep(c("a", "b"), each = 3)
  )
  panel$time_covariate <- panel$time
  cache <- gamlss.longitudinal:::build_copula_pair_cache(
    panel$response, panel$time, panel$subject
  )

  expect_identical(cache$row_id1, c(1L, 4L, 2L, 5L))
  expect_identical(cache$row_id2, c(2L, 5L, 3L, 6L))
  expect_identical(cache$observed_pair_base, c(FALSE, TRUE, FALSE, TRUE))
  expect_identical(cache$theta_index_map, c(1L, 2L, NA_integer_, 3L, 4L, NA_integer_))

  formulas <- list(
    mu = response ~ time_covariate,
    sigma = ~ 1,
    nu = ~ 1,
    tau = ~ 1,
    theta = ~ time_covariate,
    zeta = ~ 1
  )
  copula_link <- gamlss.longitudinal:::get_copula_dist("N")$copula_link
  mm <- do.call(
    gamlss.longitudinal:::create_model_matrices,
    c(formulas, list(
      margin.family = gamlss.dist::NO(), copula.family = "N",
      copula.link = copula_link, dataset = panel,
      quiet_gamlss2 = TRUE, preserve_factor_levels = TRUE
    ))
  )
  expect_equal(mm$x$mu$time_covariate, panel$time)
  expect_equal(mm$x$theta$time_covariate, c(0, 2, 0, 2))

  coefficients <- c(
    mu.intercept = 0.1, mu.time_covariate = 0.2,
    sigma.intercept = 0,
    theta.intercept = 0.1, theta.time_covariate = 0.02
  )
  smooth <- list(mu = NA, sigma = NA, theta = NA)
  eta <- gamlss.longitudinal:::calc_eta(
    coefficients, mm, gamlss.dist::NO(), copula_link, smooth
  )
  lik <- gamlss.longitudinal:::calc_likelihood_minimal(
    eta_inv = eta$eta_inv,
    mm = mm$x,
    margin_dist = gamlss.dist::NO(), copula_dist = "N",
    response = panel$response, response_margin = panel$time,
    response_subject = panel$subject, pair_cache = cache
  )

  expect_identical(unname(lik$contribution_counts[["marginal_included"]]), 5L)
  expect_identical(unname(lik$contribution_counts[["pair_included"]]), 2L)
  expect_equal(lik$copula_par1, eta$eta_inv$theta[c(1L, 3L, 2L, 4L)])
  expect_equal(lik$copula_par1[lik$pair_included], eta$eta_inv$theta[c(3L, 4L)])

  object <- structure(list(
    response = panel$response,
    response_margin = panel$time,
    response_subject = panel$subject,
    margin_dist = gamlss.dist::NO(), copula_dist = "N",
    par = coefficients, par_s = smooth, model_matrix = mm,
    formulas_int = formulas, var_map = list(), dataset = panel,
    convergence = list(converged = TRUE)
  ), class = "gamlss.longitudinal")

  predicted <- predict(object, newdata = panel, type = "mu")
  expect_equal(predicted, unname(eta$eta_inv$mu))

  simulation_data <- gamlss.longitudinal:::.gl_simulation_newdata(object, panel)
  expect_equal(
    simulation_data$fit_data$theta_fit[c(1L, 2L, 4L, 5L)],
    unname(eta$eta_inv$theta)
  )
  expect_true(all(is.na(simulation_data$fit_data$theta_fit[c(3L, 6L)])))
  simulated <- simulate(object, nsim = 2L, seed = 431, newdata = panel)
  expect_identical(dim(simulated), c(6L, 2L))
  expect_true(all(is.finite(as.matrix(simulated))))
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
      theta_formula = "~ 1",
      compute_vcov = TRUE
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
  expect_match(warnings, "model-Hessian inference is unavailable", fixed = TRUE)
  expect_identical(fit$missingness$objective, "segmented")
  expect_identical(fit$missingness$n_subjects_with_gaps, 2L)
  expect_identical(fit$missingness$n_segments, 18L)
  expect_identical(fit$missingness$n_observed, 46L)
  expect_identical(fit$missingness$n_within_segment_pairs, 28L)
  expect_identical(fit$likelihood_contract$criteria_status, "available_under_segment_independence")
  expect_identical(fit$vcov_meta$inference_status, "unavailable_for_segmented_objective")
  expect_identical(fit$vcov_meta$hessian_diagnostics$reason, "segmented_objective")

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
