test_that("optimizer control constructor is strict and method aware", {
  ctl <- gamlss_longitudinal_control(
    outer_tol = 0.01,
    max_outer_iter = 12,
    rs = list(inner_tol = 0.001, max_inner_iter = 8)
  )
  expect_s3_class(ctl, "gamlss.longitudinal.control")
  expect_equal(ctl$shared$outer_tol, 0.01)
  expect_identical(ctl$rs$max_inner_iter, 8)

  expect_error(gamlss_longitudinal_control(outer_tol = 0), "positive")
  expect_error(gamlss_longitudinal_control(rs = list(1)), "must be named")
  expect_error(gamlss_longitudinal_control(rs = list(foo = 1)), "Unknown")
  expect_error(gamlss_longitudinal_control(rs = structure(list(1, 2), names = c("max_steps", "max_steps"))), "Duplicate")
  expect_error(gamlss_longitudinal_control(rs = list(max_inner_iter = 0)), "positive integer")
  expect_error(gamlss_longitudinal_control(cg = list(max_stall = 0)), "positive integer")

  expect_error(
    gamlss.longitudinal:::.gl_resolve_optimizer_control(
      "RS", gamlss_longitudinal_control(cg = list(max_stall = 2)),
      legacy_values = list(), legacy_supplied = logical()
    ),
    "wrong-method|method ="
  )
})

test_that("legacy inactive-method controls deprecate and are ignored for one release", {
  gamlss.longitudinal:::.gl_reset_optimizer_control_deprecation()
  seen <- NULL
  cg <- withCallingHandlers(
    gamlss.longitudinal:::.gl_resolve_optimizer_control(
      "CG", NULL,
      legacy_values = list(max_inner_iter = 0),
      legacy_supplied = c(max_inner_iter = TRUE)
    ),
    warning = function(w) {
      seen <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(seen, "gamlss.longitudinal_deprecated_optimizer_control")
  expect_identical(cg$rs$max_inner_iter, 100L)

  gamlss.longitudinal:::.gl_reset_optimizer_control_deprecation()
  rs <- suppressWarnings(gamlss.longitudinal:::.gl_resolve_optimizer_control(
    "RS", NULL,
    legacy_values = list(cg_max_stall = 0),
    legacy_supplied = c(cg_max_stall = TRUE)
  ))
  expect_identical(rs$cg$max_stall, 5L)
})

test_that("legacy optimizer controls deprecate once and conflicts error", {
  gamlss.longitudinal:::.gl_reset_optimizer_control_deprecation()
  vals <- list(outer_stop_crit = 0.2)
  supplied <- c(outer_stop_crit = TRUE)
  seen <- NULL
  out <- withCallingHandlers(
    gamlss.longitudinal:::.gl_resolve_optimizer_control("RS", NULL, vals, supplied),
    warning = function(w) {
      seen <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_true(inherits(seen, "gamlss.longitudinal_deprecated_optimizer_control"))
  expect_equal(out$shared$outer_tol, 0.2)
  expect_silent(gamlss.longitudinal:::.gl_resolve_optimizer_control("RS", NULL, vals, supplied))

  expect_error(
    gamlss.longitudinal:::.gl_resolve_optimizer_control(
      "RS", gamlss_longitudinal_control(outer_tol = 0.1), vals, supplied
    ),
    "same setting"
  )
})

test_that("requested and effective controls retain automatic tolerance provenance", {
  ctl <- gamlss.longitudinal:::.gl_resolve_optimizer_control(
    "CG", gamlss_longitudinal_control(cg = list(update_lambda = FALSE)),
    legacy_values = list(), legacy_supplied = logical()
  )
  context <- list(
    outer_stop_crit = 1e-3,
    inner_stop_crit = 1e-4,
    cg_grad_tol_eff = 1e-2,
    cg_step_tol_eff = 1e-5
  )
  effective <- gamlss.longitudinal:::.gl_effective_optimizer_control(ctl, context)
  expect_true(is.na(ctl$shared$outer_tol))
  expect_equal(effective$shared$outer_tol, 1e-3)
  expect_equal(effective$cg$grad_tol, 1e-2)
  expect_equal(effective$cg$step_tol, 1e-5)
})

test_that("warm start uses one optimizer control bundle without legacy arguments", {
  captured <- NULL
  fake_fit <- function(...) {
    captured <<- list(...)
    list(
      par_s = list(mu = list()),
      calc_lik_out_end = list(log_lik = c(marginal = -2, copula = -1, joint = -3)),
      include_dlcopdpar = FALSE
    )
  }
  ctl <- gamlss.longitudinal:::.gl_resolve_optimizer_control(
    "RS", gamlss_longitudinal_control(rs = list(warm_start_joint = TRUE)),
    legacy_values = list(), legacy_supplied = logical()
  )
  out <- gamlss.longitudinal:::.gl_call_joint_warm_start_fit(
    fit_fn = fake_fit,
    dataset_original = data.frame(y = 1, id = 1, time = 1),
    margin_dist = gamlss.dist::NO(), copula_dist = "N",
    time_var = "time", subject_var = "id",
    mu.formula = y ~ 1, sigma.formula = ~1, nu.formula = ~1,
    tau.formula = ~1, theta.formula = ~1, zeta.formula = ~1,
    optimizer_control = ctl, true_val = NA, method = "RS",
    warm_start_joint_iter = 3L, vcov_method = "analytical",
    vcov_numderiv = FALSE, use_Rcpp = FALSE, lambda_start = NA,
    lambda_penalty_K = 2
  )
  expect_null(out$error)
  expect_s3_class(captured$optimizer_control, "gamlss.longitudinal.control")
  expect_identical(captured$optimizer_control$shared$max_outer_iter, 3L)
  expect_false(captured$optimizer_control$rs$warm_start_joint)
  expect_false(any(names(captured) %in% names(gamlss.longitudinal:::.gl_legacy_optimizer_map())))
})

test_that("method-neutral convergence contract uses stable stop reasons", {
  base <- list(
    method = "RS", outer_log_lik_change = 0.01, outer_stop_crit = 0.1,
    outer_only_run_counter = 3L, max_outer_iter = 10L,
    cg_stop_reason = NA_character_, cg_last_grad_inf = NA_real_,
    cg_last_step_l2 = NA_real_, cg_best_raw_loglik = -1,
    cg_best_iteration = 2L, cg_raw_loglik_drop_from_best = 0,
    cg_raw_loglik_drop_tol = 10, cg_gradient_method = "forward",
    cg_zeta_hessian = "analytical", cg_hessian_method = "analytical",
    objective = -10
  )
  rs <- do.call(gamlss.longitudinal:::.gl_build_convergence_info, base)
  expect_true(rs$converged)
  expect_identical(rs$stop_reason, "converged")

  base$outer_log_lik_change <- 1
  base$outer_only_run_counter <- 11L
  limited <- do.call(gamlss.longitudinal:::.gl_build_convergence_info, base)
  expect_false(limited$converged)
  expect_identical(limited$stop_reason, "max_iterations")
  expect_identical(limited$events$code, "max_iterations")

  base$method <- "CG"
  base$cg_stop_reason <- "raw_loglik_deterioration"
  base$cg_last_grad_inf <- 0
  base$cg_last_step_l2 <- 0
  base$cg_grad_tol <- 0.1
  base$cg_step_tol <- 0.1
  deteriorated <- do.call(gamlss.longitudinal:::.gl_build_convergence_info, base)
  expect_identical(deteriorated$stop_reason, "objective_deterioration")

  base$method <- "RS"
  base$cg_stop_reason <- NA_character_
  base$outer_only_run_counter <- 2L
  base$outer_log_lik_change <- 0.01
  base$objective <- NA_real_
  invalid <- do.call(gamlss.longitudinal:::.gl_build_convergence_info, base)
  expect_false(invalid$converged)
  expect_identical(invalid$stop_reason, "invalid_likelihood")

  base$objective <- NULL
  omitted <- do.call(gamlss.longitudinal:::.gl_build_convergence_info, base)
  expect_true(omitted$converged)
})

test_that("nonconvergence conditions gate uncertainty but allow warned points", {
  object <- structure(
    list(convergence = list(converged = FALSE, stop_reason = "max_iterations")),
    class = "gamlss.longitudinal"
  )
  err <- tryCatch(
    gamlss.longitudinal:::.gl_require_converged_fit(object, "inference"),
    error = identity
  )
  expect_true(inherits(err, "gamlss.longitudinal_nonconvergence_error"))
  seen <- NULL
  withCallingHandlers(
    gamlss.longitudinal:::.gl_warn_nonconverged_prediction(object),
    warning = function(w) {
      seen <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_true(inherits(seen, "gamlss.longitudinal_nonconverged_prediction_warning"))
  err <- tryCatch(vcov(object), error = identity)
  expect_true(inherits(err, "gamlss.longitudinal_nonconvergence_error"))
})

test_that("stop_on_convergence controls tolerance stopping without zero tolerances", {
  expect_false(gamlss.longitudinal:::.gl_should_continue_rs_outer_loop(
    FALSE, 0.01, 0.1, 2L, 10L, stop_on_convergence = TRUE
  ))
  expect_true(gamlss.longitudinal:::.gl_should_continue_rs_outer_loop(
    FALSE, 0.01, 0.1, 2L, 10L, stop_on_convergence = FALSE
  ))
  assess <- gamlss.longitudinal:::.gl_assess_cg_stopping(
    0.01, 0.01, 0.01, 0L, 5L, 0, Inf, 0L, FALSE,
    0.1, 0.1, 0.1, stop_on_convergence = FALSE
  )
  expect_true(assess$tolerance_met)
  expect_false(assess$stop_requested)
})

test_that("elapsed budget failures carry the time-limit stop vocabulary", {
  err <- tryCatch(
    gamlss.longitudinal:::.gl_check_elapsed_budget(
      as.POSIXct("2020-01-01", tz = "UTC"), 1,
      now = as.POSIXct("2020-01-01 00:00:02", tz = "UTC")
    ),
    error = identity
  )
  expect_s3_class(err, "gamlss.longitudinal_time_limit_error")
  expect_identical(err$stop_reason, "time_limit")
})

test_that("public RS and CG limits execute exactly one outer iteration", {
  dat <- make_fixture_factor_time_interaction(8L)
  common <- list(
    dataset = dat, margin_dist = gamlss.dist::NO(), copula_dist = "N",
    time_var = "time_raw", subject_var = "id", mu.formula = y ~ time_raw,
    sigma.formula = ~1, theta.formula = ~1, include_dlcopdpar = FALSE,
    compute_vcov = FALSE, verbose = 0
  )
  quiet_fit <- function(method, control) {
    withCallingHandlers(
      do.call(gamlss_longitudinal, c(common, list(
        method = method,
        optimizer_control = control
      ))),
      warning = function(w) invokeRestart("muffleWarning")
    )
  }

  rs <- quiet_fit("RS", gamlss_longitudinal_control(
    outer_tol = 1e-12,
    max_outer_iter = 1L,
    rs = list(
      inner_tol = 1e-12,
      max_inner_iter = 1L,
      warm_start_joint = FALSE
    )
  ))
  expect_equal(rs$convergence$outer_iterations, 1L)
  expect_true(nrow(rs$rs_block_trace) > 0L)
  expect_identical(unique(rs$rs_block_trace$inner_iteration), 1L)
  rs_summary <- summary(rs, include_vcov = FALSE)
  expect_false(rs_summary$convergence$converged)
  expect_identical(rs_summary$fit$criteria_status, "provisional_nonconverged")
  printed <- capture.output(print(rs_summary))
  expect_true(any(grepl("Convergence: NOT CONVERGED", printed, fixed = TRUE)))
  expect_true(any(grepl("PROVISIONAL; nonconverged fit", printed, fixed = TRUE)))

  cg <- quiet_fit("CG", gamlss_longitudinal_control(
    outer_tol = 1e-12,
    max_outer_iter = 1L,
    cg = list(grad_tol = 1e-12, step_tol = 1e-12, update_lambda = FALSE)
  ))
  expect_equal(cg$convergence$outer_iterations, 1L)
})

test_that("public RS and CG time budgets abort with classed errors", {
  dat <- make_fixture_factor_time_interaction(8L)
  common <- list(
    dataset = dat, margin_dist = gamlss.dist::NO(), copula_dist = "N",
    time_var = "time_raw", subject_var = "id", mu.formula = y ~ time_raw,
    sigma.formula = ~1, theta.formula = ~1, include_dlcopdpar = FALSE,
    compute_vcov = FALSE, verbose = 0
  )
  for (method in c("RS", "CG")) {
    control <- if (identical(method, "RS")) {
      gamlss_longitudinal_control(
        max_elapsed_sec = .Machine$double.eps,
        rs = list(warm_start_joint = FALSE)
      )
    } else {
      gamlss_longitudinal_control(max_elapsed_sec = .Machine$double.eps)
    }
    err <- tryCatch(
      do.call(gamlss_longitudinal, c(common, list(
        method = method,
        optimizer_control = control
      ))),
      error = identity
    )
    expect_s3_class(err, "gamlss.longitudinal_time_limit_error")
    expect_identical(err$stop_reason, "time_limit")
  }
})

test_that("legacy and control-object routes are numerically equivalent for RS and CG", {
  dat <- make_fixture_factor_time_interaction(8L)
  common <- list(
    dataset = dat, margin_dist = gamlss.dist::NO(), copula_dist = "N",
    time_var = "time_raw", subject_var = "id", mu.formula = y ~ time_raw,
    sigma.formula = ~1, theta.formula = ~1, include_dlcopdpar = FALSE,
    compute_vcov = FALSE, verbose = 0
  )
  quiet_fit <- function(args) {
    withCallingHandlers(
      do.call(gamlss_longitudinal, args),
      warning = function(w) invokeRestart("muffleWarning")
    )
  }

  gamlss.longitudinal:::.gl_reset_optimizer_control_deprecation()
  rs_old <- quiet_fit(c(common, list(
    method = "RS", outer_stop_crit = 1, inner_stop_crit = 1,
    max_outer_iter = 2L, max_inner_iter = 2L, warm_start_joint = FALSE
  )))
  rs_new <- quiet_fit(c(common, list(
    method = "RS",
    optimizer_control = gamlss_longitudinal_control(
      outer_tol = 1, max_outer_iter = 2L,
      rs = list(inner_tol = 1, max_inner_iter = 2L, warm_start_joint = FALSE)
    )
  )))
  expect_equal(rs_new$par, rs_old$par, tolerance = 1e-12)

  cg_old <- quiet_fit(c(common, list(
    method = "CG", outer_stop_crit = 1, max_outer_iter = 2L,
    cg_grad_tol = 100, cg_step_tol = 100, cg_update_lambda = FALSE
  )))
  cg_new <- quiet_fit(c(common, list(
    method = "CG",
    optimizer_control = gamlss_longitudinal_control(
      outer_tol = 1, max_outer_iter = 2L,
      cg = list(grad_tol = 100, step_tol = 100, update_lambda = FALSE)
    )
  )))
  expect_equal(cg_new$par, cg_old$par, tolerance = 1e-12)
})
