test_that("stop criteria validation preserves explicit tolerances", {
  out <- gamlss.longitudinal:::.gl_resolve_stop_criteria(
    inner_stop_crit = 0.01,
    outer_stop_crit = 0.02,
    cg_grad_tol = 0.03,
    cg_step_tol = 0.04,
    method = "RS",
    par_cov = c(`mu.(Intercept)` = 1),
    par_s = list(),
    mm = list(x = list(), s = list()),
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    copula_link = list(),
    dataset = data.frame(response = 1, time = 1, subject = "a"),
    pair_cache = list(),
    margin_eval_cache = list(),
    verbose = 0
  )

  expect_equal(out$inner_stop_crit, 0.01)
  expect_equal(out$outer_stop_crit, 0.02)
  expect_equal(out$cg_grad_tol_eff, 0.03)
  expect_equal(out$cg_step_tol_eff, 0.04)

  expect_error(
    gamlss.longitudinal:::.gl_validate_stop_crit(0, "outer_stop_crit"),
    "outer_stop_crit must be a single positive finite number",
    fixed = TRUE
  )
})

test_that("outer negative streak helper increments and resets", {
  expect_equal(
    gamlss.longitudinal:::.gl_update_outer_negative_streak(
      outer_log_lik_change = -0.1,
      outer_negative_streak = 0L,
      max_negative_outer_streak = 3L
    ),
    1
  )

  expect_equal(
    gamlss.longitudinal:::.gl_update_outer_negative_streak(
      outer_log_lik_change = 0,
      outer_negative_streak = 2L,
      max_negative_outer_streak = 3L
    ),
    0
  )

  expect_equal(
    gamlss.longitudinal:::.gl_update_outer_negative_streak(
      outer_log_lik_change = NA_real_,
      outer_negative_streak = 2L,
      max_negative_outer_streak = 3L
    ),
    0
  )
})

test_that("outer negative streak helper preserves threshold warning and stop", {
  expect_warning(
    expect_error(
      gamlss.longitudinal:::.gl_update_outer_negative_streak(
        outer_log_lik_change = -0.1,
        outer_negative_streak = 2L,
        max_negative_outer_streak = 3L
      ),
      "Optimization stopped after 3 consecutive negative outer log-likelihood changes",
      fixed = TRUE
    ),
    "Optimization stopped after 3 consecutive negative outer log-likelihood changes",
    fixed = TRUE
  )
})

test_that("outer iteration summary helper returns named values and respects verbosity", {
  cat_calls <- character()
  print_calls <- list()

  out_silent <- gamlss.longitudinal:::.gl_print_outer_iteration_summary(
    outer_start_log_lik = -10,
    outer_end_log_lik = -8,
    outer_log_lik_change = 2,
    verbose = 0,
    cat_fn = function(...) cat_calls <<- c(cat_calls, paste0(...)),
    print_fn = function(x, ...) print_calls[[length(print_calls) + 1L]] <<- x
  )

  expect_equal(out_silent, c(`Start LogLik` = -10, `End LogLik` = -8, Change = 2))
  expect_equal(cat_calls, character())
  expect_equal(print_calls, list())

  out_verbose <- gamlss.longitudinal:::.gl_print_outer_iteration_summary(
    outer_start_log_lik = -10,
    outer_end_log_lik = -8,
    outer_log_lik_change = 2,
    verbose = 1,
    cat_fn = function(...) cat_calls <<- c(cat_calls, paste0(...)),
    print_fn = function(x, ...) print_calls[[length(print_calls) + 1L]] <<- x
  )

  expect_equal(out_verbose, c(`Start LogLik` = -10, `End LogLik` = -8, Change = 2))
  expect_equal(tail(cat_calls, 1), "\n")
  expect_equal(print_calls[[1]], c(`Start LogLik` = -10, `End LogLik` = -8, Change = 2))
})

test_that("outer convergence helper reports only when tolerance is satisfied and verbose", {
  cat_calls <- character()
  print_calls <- list()

  not_converged <- gamlss.longitudinal:::.gl_print_outer_convergence(
    outer_log_lik_change = 0.2,
    outer_stop_crit = 0.1,
    verbose = 1,
    cat_fn = function(...) cat_calls <<- c(cat_calls, paste0(...)),
    print_fn = function(x, ...) print_calls[[length(print_calls) + 1L]] <<- x
  )
  expect_false(not_converged)
  expect_equal(cat_calls, character())
  expect_equal(print_calls, list())

  converged_silent <- gamlss.longitudinal:::.gl_print_outer_convergence(
    outer_log_lik_change = 0.05,
    outer_stop_crit = 0.1,
    verbose = 0,
    cat_fn = function(...) cat_calls <<- c(cat_calls, paste0(...)),
    print_fn = function(x, ...) print_calls[[length(print_calls) + 1L]] <<- x
  )
  expect_true(converged_silent)
  expect_equal(cat_calls, character())
  expect_equal(print_calls, list())

  converged_verbose <- gamlss.longitudinal:::.gl_print_outer_convergence(
    outer_log_lik_change = -0.05,
    outer_stop_crit = 0.1,
    verbose = 1,
    cat_fn = function(...) cat_calls <<- c(cat_calls, paste0(...)),
    print_fn = function(x, ...) print_calls[[length(print_calls) + 1L]] <<- x
  )
  expect_true(converged_verbose)
  expect_equal(print_calls[[1]], c(-0.05))
  expect_equal(cat_calls, "\nOUTER CONVERGED")
})

test_that("RS inner-loop continuation helper preserves first-run and tolerance logic", {
  expect_true(gamlss.longitudinal:::.gl_should_continue_rs_inner_loop(
    first_inner_run = TRUE,
    change_log_lik = 0,
    inner_stop_crit = 0.1,
    inner_run_counter = 1L,
    max_inner_iter = 5L
  ))

  expect_true(gamlss.longitudinal:::.gl_should_continue_rs_inner_loop(
    first_inner_run = FALSE,
    change_log_lik = 0.2,
    inner_stop_crit = 0.1,
    inner_run_counter = 2L,
    max_inner_iter = 5L
  ))

  expect_false(gamlss.longitudinal:::.gl_should_continue_rs_inner_loop(
    first_inner_run = FALSE,
    change_log_lik = 0.05,
    inner_stop_crit = 0.1,
    inner_run_counter = 2L,
    max_inner_iter = 5L
  ))

  expect_false(gamlss.longitudinal:::.gl_should_continue_rs_inner_loop(
    first_inner_run = TRUE,
    change_log_lik = 1,
    inner_stop_crit = 0.1,
    inner_run_counter = 5L,
    max_inner_iter = 5L
  ))
})

test_that("RS outer-loop continuation helper preserves first-run and tolerance logic", {
  expect_true(gamlss.longitudinal:::.gl_should_continue_rs_outer_loop(
    first_outer_run = TRUE,
    outer_log_lik_change = 0,
    outer_stop_crit = 0.1,
    outer_only_run_counter = 1L,
    max_outer_iter = 5L
  ))

  expect_true(gamlss.longitudinal:::.gl_should_continue_rs_outer_loop(
    first_outer_run = FALSE,
    outer_log_lik_change = 0.2,
    outer_stop_crit = 0.1,
    outer_only_run_counter = 2L,
    max_outer_iter = 5L
  ))

  expect_false(gamlss.longitudinal:::.gl_should_continue_rs_outer_loop(
    first_outer_run = FALSE,
    outer_log_lik_change = 0.05,
    outer_stop_crit = 0.1,
    outer_only_run_counter = 2L,
    max_outer_iter = 5L
  ))

  expect_false(gamlss.longitudinal:::.gl_should_continue_rs_outer_loop(
    first_outer_run = TRUE,
    outer_log_lik_change = 1,
    outer_stop_crit = 0.1,
    outer_only_run_counter = 5L,
    max_outer_iter = 5L
  ))
})

test_that("RS outer iteration state helper updates schedule and delegates reporting", {
  summary_calls <- list()
  streak_calls <- list()
  convergence_calls <- list()
  calc_lik_out_end <- list(log_lik = c(margin = -3, copula = -2, joint = -4))

  out <- gamlss.longitudinal:::.gl_update_rs_outer_iteration_state(
    calc_lik_out_end = calc_lik_out_end,
    outer_start_log_lik = c(joint = -5),
    outer_only_run_counter = 3L,
    outer_negative_streak = 1L,
    step_adjustment = 0.5,
    max_steps = 2L,
    start_step_size = 0.8,
    max_negative_outer_streak = 4L,
    outer_stop_crit = 0.1,
    verbose = 2,
    summary_fn = function(outer_start_log_lik, outer_end_log_lik, outer_log_lik_change, verbose) {
      summary_calls[[length(summary_calls) + 1L]] <<- list(
        outer_start_log_lik = outer_start_log_lik,
        outer_end_log_lik = outer_end_log_lik,
        outer_log_lik_change = outer_log_lik_change,
        verbose = verbose
      )
      invisible(c(outer_start_log_lik, outer_end_log_lik, outer_log_lik_change))
    },
    streak_fn = function(outer_log_lik_change, outer_negative_streak, max_negative_outer_streak) {
      streak_calls[[length(streak_calls) + 1L]] <<- list(
        outer_log_lik_change = outer_log_lik_change,
        outer_negative_streak = outer_negative_streak,
        max_negative_outer_streak = max_negative_outer_streak
      )
      0L
    },
    convergence_fn = function(outer_log_lik_change, outer_stop_crit, verbose) {
      convergence_calls[[length(convergence_calls) + 1L]] <<- list(
        outer_log_lik_change = outer_log_lik_change,
        outer_stop_crit = outer_stop_crit,
        verbose = verbose
      )
      invisible(FALSE)
    }
  )

  expect_equal(out$step_size, 0.2)
  expect_equal(out$outer_only_run_counter, 4L)
  expect_equal(out$outer_end_log_lik, c(joint = -4))
  expect_equal(out$outer_log_lik_change, c(joint = 1))
  expect_equal(out$outer_negative_streak, 0L)

  expect_equal(length(summary_calls), 1)
  expect_equal(summary_calls[[1]]$outer_start_log_lik, c(joint = -5))
  expect_equal(summary_calls[[1]]$outer_end_log_lik, c(joint = -4))
  expect_equal(summary_calls[[1]]$outer_log_lik_change, c(joint = 1))
  expect_equal(summary_calls[[1]]$verbose, 2)

  expect_equal(length(streak_calls), 1)
  expect_equal(streak_calls[[1]]$outer_log_lik_change, c(joint = 1))
  expect_equal(streak_calls[[1]]$outer_negative_streak, 1L)
  expect_equal(streak_calls[[1]]$max_negative_outer_streak, 4L)

  expect_equal(length(convergence_calls), 1)
  expect_equal(convergence_calls[[1]]$outer_log_lik_change, c(joint = 1))
  expect_equal(convergence_calls[[1]]$outer_stop_crit, 0.1)
  expect_equal(convergence_calls[[1]]$verbose, 2)
})

test_that("automatic stop criteria are derived from initial likelihood scale", {
  prepared <- gamlss.longitudinal:::.gl_prepare_fit_data(
    dataset = data.frame(
      id = c("a", "a", "b", "b"),
      visit = c(1, 2, 1, 2),
      y = c(1.0, 2.0, 1.5, 2.5)
    ),
    time_var = "visit",
    subject_var = "id",
    mu.formula = y ~ 1,
    sigma.formula = ~ 1,
    nu.formula = ~ 1,
    tau.formula = ~ 1,
    theta.formula = ~ 1,
    zeta.formula = ~ 1,
    verbose = 0
  )
  bundle <- gamlss.longitudinal:::.gl_build_model_matrix_bundle(
    formulas_int = prepared$formulas_int,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    dataset = prepared$dataset
  )
  initial_state <- gamlss.longitudinal:::.gl_build_initial_parameter_state(
    start_from = NA,
    mm = bundle$mm,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    dataset = prepared$dataset
  )
  optimizer_state <- gamlss.longitudinal:::.gl_initialize_optimizer_state(
    par_cov = initial_state$par_cov,
    start_step_size = 0.5,
    dataset = prepared$dataset,
    margin_dist = gamlss.dist::NO()
  )

  out_rs <- gamlss.longitudinal:::.gl_resolve_stop_criteria(
    inner_stop_crit = NA,
    outer_stop_crit = NA,
    cg_grad_tol = NA,
    cg_step_tol = NA,
    method = "RS",
    par_cov = initial_state$par_cov,
    par_s = initial_state$par_s,
    mm = bundle$mm,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    copula_link = bundle$copula_link,
    dataset = prepared$dataset,
    pair_cache = optimizer_state$pair_cache,
    margin_eval_cache = optimizer_state$margin_eval_cache,
    verbose = 0
  )
  out_cg <- gamlss.longitudinal:::.gl_resolve_stop_criteria(
    inner_stop_crit = NA,
    outer_stop_crit = NA,
    cg_grad_tol = NA,
    cg_step_tol = NA,
    method = "CG",
    par_cov = initial_state$par_cov,
    par_s = initial_state$par_s,
    mm = bundle$mm,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    copula_link = bundle$copula_link,
    dataset = prepared$dataset,
    pair_cache = optimizer_state$pair_cache,
    margin_eval_cache = optimizer_state$margin_eval_cache,
    verbose = 0
  )

  expect_true(out_rs$outer_stop_crit > 0)
  expect_true(out_rs$inner_stop_crit > 0)
  expect_equal(out_cg$outer_stop_crit, out_rs$outer_stop_crit / 10)
  expect_equal(out_rs$cg_grad_tol_eff, max(1e-3, 10 * out_rs$outer_stop_crit))
  expect_equal(out_rs$cg_step_tol_eff, max(1e-5, 0.1 * out_rs$outer_stop_crit))
})
