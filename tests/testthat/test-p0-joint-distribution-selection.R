test_that("select_joint_distribution ranks joint fits and exposes best-fit helpers", {
  skip_if_not_installed("gamlss")

  dat <- simulate_longitudinal_dataset(
    n = 12,
    times = 1:3,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "N",
    copula_params = list(theta = 0.3),
    seed = 1001
  )

  selected <- select_joint_distribution(
    data = dat,
    response_var = "response",
    time_var = "time",
    subject_var = "subject",
    type = "realAll",
    margin_families = "NO",
    copula_families = c("N", "C"),
    progress = FALSE,
    fit_args = list(
      max_outer_iter = 3,
      max_inner_iter = 3,
      outer_stop_crit = 1,
      inner_stop_crit = 1,
      backtracking_max_halves = 0
    )
  )

  expect_s3_class(selected, "joint_distribution_selection")
  expect_equal(sort(unique(selected$margin_family)), "NO")
  expect_equal(sort(unique(selected$copula_family)), c("C", "N"))
  expect_true(any(is.finite(selected$AIC)))
  expect_true(all(is.na(selected$error[is.finite(selected$AIC)])))
  expect_equal(selected$rank[is.finite(selected$AIC)], seq_len(sum(is.finite(selected$AIC))))
  expect_equal(attr(selected, "selected"), paste(selected$margin_family[[1L]], selected$copula_family[[1L]], sep = "+"))

  best <- best_fit(selected)
  best_family <- best_fit_family(selected)
  expect_equal(best$margin_family_name, "NO")
  expect_equal(best$margin_family$family[[1L]], "NO")
  expect_true(best$copula_family %in% c("N", "C"))
  expect_equal(best_family$margin_dist$family[[1L]], "NO")
  expect_equal(best_family$copula_dist, best$copula_family)
  expect_equal(selected$best_fit$copula_family, best$copula_family)
})

test_that("select_joint_distribution can screen with time-intercept margins", {
  skip_if_not_installed("gamlss")

  dat <- simulate_longitudinal_dataset(
    n = 10,
    times = 1:3,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "N",
    copula_params = list(theta = 0.25),
    seed = 1002
  )
  dat$response <- dat$response + rep(c(-0.4, 0.1, 0.5), times = length(unique(dat$subject)))

  selected <- select_joint_distribution(
    data = dat,
    response_var = "response",
    time_var = "time",
    subject_var = "subject",
    type = "realAll",
    margin_families = "NO",
    copula_families = "N",
    time_intercepts = TRUE,
    progress = FALSE,
    fit_args = list(
      max_outer_iter = 2,
      max_inner_iter = 2,
      outer_stop_crit = 1,
      inner_stop_crit = 1,
      backtracking_max_halves = 0
    )
  )

  expect_s3_class(selected, "joint_distribution_selection")
  expect_true(isTRUE(attr(selected, "time_intercepts")))
  expect_equal(attr(selected, "time_var"), "time")
  expect_true(isTRUE(attr(attr(selected, "margin_selection"), "time_intercepts")))
  expect_equal(attr(attr(selected, "margin_selection"), "time_var"), "time")
  expect_equal(selected$margin_family, "NO")
  expect_equal(selected$copula_family, "N")
})

test_that("select_joint_distribution can screen with factor copula time intercepts", {
  skip_if_not_installed("gamlss")

  dat <- simulate_longitudinal_dataset(
    n = 10,
    times = 1:3,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "N",
    copula_params = list(theta = 0.25),
    seed = 1003
  )

  selected <- select_joint_distribution(
    data = dat,
    response_var = "response",
    time_var = "time",
    subject_var = "subject",
    type = "realAll",
    margin_families = "NO",
    copula_families = "N",
    copula_time_intercepts = TRUE,
    progress = FALSE,
    fit_args = list(
      max_outer_iter = 2,
      max_inner_iter = 2,
      outer_stop_crit = 1,
      inner_stop_crit = 1,
      backtracking_max_halves = 0
    )
  )

  expect_s3_class(selected, "joint_distribution_selection")
  expect_true(isTRUE(attr(selected, "copula_time_intercepts")))
  expect_equal(attr(selected, "copula_time_var"), "time")
  expect_equal(selected$margin_family, "NO")
  expect_equal(selected$copula_family, "N")
})

test_that("joint distribution fit formulas preserve time-intercept policies", {
  intercept_only <- gamlss.longitudinal:::.joint_selection_fit_formulas(
    response_var = "response",
    time_var = "visit",
    time_intercepts = FALSE,
    copula_time_intercepts = FALSE
  )
  time_specific <- gamlss.longitudinal:::.joint_selection_fit_formulas(
    response_var = "response",
    time_var = "visit",
    time_intercepts = TRUE,
    copula_time_intercepts = TRUE
  )

  expect_equal(deparse(intercept_only$mu_formula), "response ~ 1")
  expect_equal(deparse(intercept_only$par_formula), "~1")
  expect_equal(deparse(intercept_only$theta_formula), "~1")
  expect_match(deparse(time_specific$mu_formula), "response ~ factor\\(visit\\)")
  expect_match(deparse(time_specific$par_formula), "~factor\\(visit\\)")
  expect_match(deparse(time_specific$theta_formula), "~factor\\(visit\\)")
})

test_that("joint distribution fit row helpers preserve success and failure schemas", {
  fit <- list(convergence = list(converged = TRUE, hit_outer_limit = FALSE))
  fit_metrics <- list(
    logLik = -12,
    AIC = 30,
    BIC = 35,
    model_selection = rbind(
      LogLik = c(marginal = -10, copula = -2, joint = -12),
      AIC = c(marginal = 20, copula = 6, joint = 30),
      BIC = c(marginal = 23, copula = 8, joint = 35),
      EDF = c(marginal = 2, copula = 1, joint = 3)
    )
  )

  success <- gamlss.longitudinal:::.joint_selection_success_row(
    margin_family = "NO",
    copula_family = "N",
    fit = fit,
    fit_metrics = fit_metrics,
    elapsed = 1.5,
    warnings = c("a", "a", "b")
  )
  failure <- gamlss.longitudinal:::.joint_selection_failed_row(
    margin_family = "NO",
    copula_family = "C",
    elapsed = 2,
    warnings = "warn",
    error = "failed"
  )

  expect_named(success, names(failure))
  expect_equal(success$logLik, -12)
  expect_equal(success$EDF, 3)
  expect_true(success$converged)
  expect_equal(success$warnings, "a\nb")
  expect_true(is.na(success$error))
  expect_false(failure$converged)
  expect_equal(failure$error, "failed")
})

test_that("select_joint_distribution rejects zero likelihood resets after nonzero history", {
  fit <- list(
    log_lik_history = cbind(
      marginal = c(-20, -18, 0),
      copula = c(3, 4, 0),
      joint = c(-17, -14, 0)
    )
  )
  fit_metrics <- list(
    logLik = 0,
    AIC = 10,
    BIC = 20,
    model_selection = rbind(
      LogLik = c(marginal = 0, copula = 0, joint = 0),
      AIC = c(marginal = 4, copula = 6, joint = 10),
      BIC = c(marginal = 8, copula = 12, joint = 20),
      EDF = c(marginal = 2, copula = 3, joint = 5)
    )
  )

  expect_match(
    gamlss.longitudinal:::.joint_selection_invalid_fit_reason(fit, fit_metrics),
    "zero final likelihood"
  )

  fit_metrics$model_selection["LogLik", ] <- c(marginal = -20, copula = 3, joint = -17)
  fit_metrics$logLik <- -17
  expect_null(gamlss.longitudinal:::.joint_selection_invalid_fit_reason(fit, fit_metrics))
})

test_that("joint distribution finalizer ranks successes and preserves fit ordering", {
  margin_selection <- data.frame(family = "NO", stringsAsFactors = FALSE)
  attr(margin_selection, "response_type") <- "realAll"

  rows <- list(
    data.frame(
      margin_family = "NO",
      copula_family = "C",
      logLik = -12,
      AIC = 30,
      BIC = 35,
      EDF = 3,
      converged = TRUE,
      hit_outer_limit = FALSE,
      elapsed_sec = 1,
      warnings = "",
      error = NA_character_,
      stringsAsFactors = FALSE
    ),
    data.frame(
      margin_family = "NO",
      copula_family = "N",
      logLik = -10,
      AIC = 25,
      BIC = 31,
      EDF = 3,
      converged = TRUE,
      hit_outer_limit = FALSE,
      elapsed_sec = 2,
      warnings = "",
      error = NA_character_,
      stringsAsFactors = FALSE
    ),
    data.frame(
      margin_family = "NO",
      copula_family = "t",
      logLik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      EDF = NA_real_,
      converged = FALSE,
      hit_outer_limit = NA,
      elapsed_sec = 3,
      warnings = "",
      error = "failed",
      stringsAsFactors = FALSE
    )
  )

  out <- gamlss.longitudinal:::.joint_selection_finalize_result(
    rows = rows,
    fit_store = list("fit-C", "fit-N", "fit-t"),
    n_pairs = 8L,
    criterion = "AIC",
    margin_selection = margin_selection,
    time_intercepts = TRUE,
    time_var = "visit",
    copula_time_intercepts = TRUE,
    keep_fits = TRUE
  )

  expect_s3_class(out, "joint_distribution_selection")
  expect_equal(out$copula_family, c("N", "C", "t"))
  expect_equal(out$rank, c(1L, 2L, NA_integer_))
  expect_equal(out$delta_AIC, c(0, 5, NA))
  expect_equal(out$n_pairs, rep(8L, 3))
  expect_equal(attr(out, "selected"), "NO+N")
  expect_equal(attr(out, "response_type"), "realAll")
  expect_true(isTRUE(attr(out, "time_intercepts")))
  expect_equal(attr(out, "time_var"), "visit")
  expect_true(isTRUE(attr(out, "copula_time_intercepts")))
  expect_equal(attr(out, "copula_time_var"), "visit")
  expect_equal(attr(out, "fits"), list("fit-N", "fit-C", "fit-t"))
})

test_that("select_joint_distribution retains failures and includes t candidates", {
  skip_if_not_installed("gamlss")

  dat <- simulate_longitudinal_dataset(
    n = 8,
    times = 1:3,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "N",
    copula_params = list(theta = 0.2),
    seed = 1004
  )

  messages <- character(0)
  selected <- withCallingHandlers(
    select_joint_distribution(
      data = dat,
      response_var = "response",
      time_var = "time",
      subject_var = "subject",
      type = "realAll",
      margin_families = "NO",
      copula_families = c("N", "t"),
      progress = TRUE,
      fit_args = list(max_elapsed_sec = 1e-9)
    ),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_true(any(grepl("Joint distribution screen", messages, fixed = TRUE)))
  expect_s3_class(selected, "joint_distribution_selection")
  expect_equal(sort(unique(selected$copula_family)), c("N", "t"))
  expect_true(all(!is.na(selected$error)))
  expect_true(all(!is.finite(selected$AIC)))
  expect_true(is.na(attr(selected, "selected")))
})
