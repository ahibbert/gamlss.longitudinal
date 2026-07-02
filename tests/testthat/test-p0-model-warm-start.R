make_warm_start_args <- function(...) {
  defaults <- list(
    start_from = NA,
    method = "RS",
    include_dlcopdpar = TRUE,
    warm_start_joint = TRUE,
    warm_start_joint_iter = 2L,
    user_supplied_start = FALSE,
    dataset_original = data.frame(response = 1:2, time = 1:2, subject = 1),
    margin_dist = list(family = c("NO", "Normal")),
    copula_dist = "N",
    time_var = "time",
    subject_var = "subject",
    mu.formula = response ~ 1,
    sigma.formula = ~ 1,
    nu.formula = ~ 1,
    tau.formula = ~ 1,
    theta.formula = ~ 1,
    zeta.formula = ~ 1,
    inner_stop_crit = 1e-6,
    outer_stop_crit = 1e-6,
    start_step_size = 0.5,
    step_adjustment = 1,
    max_steps = 5L,
    true_val = NA,
    max_inner_iter = 3L,
    max_negative_outer_streak = 10L,
    max_elapsed_sec = Inf,
    use_backtracking = TRUE,
    backtracking_max_halves = 50L,
    cg_max_stall = 5L,
    cg_max_delta = 0.5,
    cg_armijo_c1 = 1e-4,
    cg_grad_tol = NA,
    cg_step_tol = NA,
    cg_update_lambda = TRUE,
    cg_lambda_update_every = 10L,
    cg_line_search = "best",
    cg_max_line_search_evals = Inf,
    cg_gradient_method = "forward",
    discrete_score_method = "analytical",
    cg_zeta_hessian = "analytical",
    cg_hessian_method = "analytical",
    vcov_method = "analytical",
    vcov_numderiv = FALSE,
    use_Rcpp = FALSE,
    lambda_start = NA,
    lambda_penalty_K = 2,
    verbose = 0,
    fit_fn = function(...) {
      list(
        par = c(`mu.(Intercept)` = 1),
        par_s = list(mu = list(`s(x)` = c(`mu.s(x).1` = 0.1))),
        calc_lik_out_end = list(log_lik = c(marginal = -2, copula = -1, joint = -3))
      )
    }
  )
  utils::modifyList(defaults, list(...))
}

run_warm_start <- function(...) {
  do.call(gamlss.longitudinal:::.gl_run_joint_warm_start, make_warm_start_args(...))
}

run_warm_start_fit_call <- function(...) {
  args <- make_warm_start_args(...)
  formals_needed <- names(formals(gamlss.longitudinal:::.gl_call_joint_warm_start_fit))
  args <- args[intersect(names(args), formals_needed)]
  do.call(gamlss.longitudinal:::.gl_call_joint_warm_start_fit, args)
}

test_that("warm-start policy preserves the previous joint RS conditions", {
  expect_true(gamlss.longitudinal:::.gl_should_run_joint_warm_start(
    method = "RS",
    include_dlcopdpar = TRUE,
    warm_start_joint = TRUE,
    warm_start_joint_iter = 1L,
    user_supplied_start = FALSE
  ))
  expect_false(gamlss.longitudinal:::.gl_should_run_joint_warm_start("CG", TRUE, TRUE, 1L, FALSE))
  expect_false(gamlss.longitudinal:::.gl_should_run_joint_warm_start("RS", FALSE, TRUE, 1L, FALSE))
  expect_false(gamlss.longitudinal:::.gl_should_run_joint_warm_start("RS", TRUE, FALSE, 1L, FALSE))
  expect_false(gamlss.longitudinal:::.gl_should_run_joint_warm_start("RS", TRUE, TRUE, 0L, FALSE))
  expect_false(gamlss.longitudinal:::.gl_should_run_joint_warm_start("RS", TRUE, TRUE, 1L, TRUE))
})

test_that("warm-start helper returns default metadata when skipped", {
  out <- run_warm_start(
    user_supplied_start = TRUE,
    start_from = c(`mu.(Intercept)` = 2),
    fit_fn = function(...) stop("should not run")
  )

  expect_equal(out$start_from, c(`mu.(Intercept)` = 2))
  expect_null(out$warm_start_par_s)
  expect_false(out$warm_start_info$used)
  expect_equal(out$warm_start_info$outer_iter, 0L)
  expect_null(out$warm_start_info$log_lik)
})

test_that("warm-start helper captures fit output warnings and starts", {
  fit_args <- NULL
  warm_fit <- function(...) {
    fit_args <<- list(...)
    cat("warm output line\n")
    warning("warm warning")
    list(
      par = c(`mu.(Intercept)` = 1.2),
      par_s = list(mu = list(`s(x)` = c(`mu.s(x).1` = 0.4))),
      calc_lik_out_end = list(log_lik = c(marginal = -4, copula = -2, joint = -6))
    )
  }

  top_output <- capture.output(
    out <- run_warm_start(verbose = 2, fit_fn = warm_fit)
  )

  expect_equal(out$start_from, c(`mu.(Intercept)` = 1.2))
  expect_equal(out$warm_start_par_s$mu$`s(x)`, c(`mu.s(x).1` = 0.4))
  expect_true(out$warm_start_info$used)
  expect_equal(out$warm_start_info$outer_iter, 2L)
  expect_equal(out$warm_start_info$log_lik, c(marginal = -4, copula = -2, joint = -6))
  expect_true(out$warm_start_info$carries_smooth)
  expect_match(out$warm_start_info$captured_output, "warm output line")
  expect_equal(out$warm_start_info$captured_warnings, "warm warning")
  expect_match(paste(top_output, collapse = "\n"), "Running separate RS warm-start phase", fixed = TRUE)
  expect_match(paste(top_output, collapse = "\n"), "warm output line", fixed = TRUE)

  expect_false(fit_args$include_dlcopdpar)
  expect_false(fit_args$check_dlcopdpar_gradient)
  expect_false(fit_args$warm_start_joint)
  expect_equal(fit_args$warm_start_joint_iter, 0L)
  expect_equal(fit_args$max_outer_iter, 2L)
  expect_false(fit_args$compute_vcov)
})

test_that("warm-start result helper records accepted fit metadata", {
  warm_fit <- list(
    par = c(`mu.(Intercept)` = 1.2),
    par_s = list(mu = list(`s(x)` = c(`mu.s(x).1` = 0.4))),
    calc_lik_out_end = list(log_lik = c(marginal = -4, copula = -2, joint = -6))
  )

  expect_invisible(gamlss.longitudinal:::.gl_validate_joint_warm_start_fit(warm_fit, NULL))
  out <- gamlss.longitudinal:::.gl_joint_warm_start_result(
    warm_fit = warm_fit,
    warm_output = "warm output",
    warm_warnings = "warm warning",
    warm_start_joint_iter = 2L
  )

  expect_equal(out$start_from, c(`mu.(Intercept)` = 1.2))
  expect_equal(out$warm_start_par_s, warm_fit$par_s)
  expect_true(out$warm_start_info$used)
  expect_equal(out$warm_start_info$outer_iter, 2L)
  expect_false(out$warm_start_info$include_dlcopdpar)
  expect_equal(out$warm_start_info$log_lik, c(marginal = -4, copula = -2, joint = -6))
  expect_true(out$warm_start_info$carries_smooth)
  expect_equal(out$warm_start_info$captured_output, "warm output")
  expect_equal(out$warm_start_info$captured_warnings, "warm warning")
})

test_that("warm-start fit-call helper captures output warnings and cold-start arguments", {
  fit_args <- NULL

  out <- run_warm_start_fit_call(
    fit_fn = function(...) {
      fit_args <<- list(...)
      cat("fit call output\n")
      warning("fit call warning")
      list(
        par = c(`mu.(Intercept)` = 2),
        par_s = list(mu = list()),
        calc_lik_out_end = list(log_lik = c(joint = -5))
      )
    }
  )

  expect_equal(out$fit$par, c(`mu.(Intercept)` = 2))
  expect_match(out$output, "fit call output")
  expect_equal(out$warnings, "fit call warning")
  expect_null(out$error)
  expect_false(fit_args$include_dlcopdpar)
  expect_false(fit_args$check_dlcopdpar_gradient)
  expect_false(fit_args$warm_start_joint)
  expect_equal(fit_args$warm_start_joint_iter, 0L)
  expect_equal(fit_args$max_outer_iter, 2L)
  expect_false(fit_args$compute_vcov)
})

test_that("warm-start runner delegates separate fit call", {
  captured <- NULL

  out <- run_warm_start(
    warm_start_fit_call_fn = function(...) {
      captured <<- list(...)
      list(
        fit = list(
          par = c(`mu.(Intercept)` = 3),
          par_s = list(mu = list()),
          calc_lik_out_end = list(log_lik = c(marginal = -1, copula = -2, joint = -3))
        ),
        output = "delegated output",
        warnings = "delegated warning",
        error = NULL
      )
    }
  )

  expect_equal(captured$warm_start_joint_iter, 2L)
  expect_equal(captured$method, "RS")
  expect_equal(out$start_from, c(`mu.(Intercept)` = 3))
  expect_equal(out$warm_start_info$captured_output, "delegated output")
  expect_equal(out$warm_start_info$captured_warnings, "delegated warning")
})

test_that("warm-start helper preserves failure messages", {
  expect_error(
    run_warm_start(fit_fn = function(...) stop("boom")),
    "Separate RS warm-start phase failed: boom",
    fixed = TRUE
  )
  expect_error(
    run_warm_start(fit_fn = function(...) list(foo = list(), calc_lik_out_end = list(log_lik = c(joint = -1)))),
    "Separate RS warm-start phase did not return coefficient starting values",
    fixed = TRUE
  )
})
