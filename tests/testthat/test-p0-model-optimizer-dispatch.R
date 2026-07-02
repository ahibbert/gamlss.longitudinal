make_optimizer_context <- function() {
  list(
    par_cov = c(`mu.(Intercept)` = 0),
    par_s = list(mu = list()),
    df_s = list(mu = 1),
    lambda_s = list(mu = list()),
    rs_design_cache = list(mu = list()),
    first_outer_run = TRUE,
    outer_log_lik_change = 0,
    outer_start_log_lik = -6,
    outer_end_log_lik = -5,
    log_lik_history = matrix(ncol = 3, nrow = 0),
    par_history = matrix(ncol = 1, nrow = 0, dimnames = list(NULL, "mu.(Intercept)")),
    cg_stop_reason = NA_character_,
    cg_last_grad_inf = NA_real_,
    cg_last_step_l2 = NA_real_,
    cg_best_raw_loglik = -Inf,
    cg_best_iteration = NA_integer_,
    cg_raw_loglik_drop_from_best = NA_real_,
    rs_block_trace = list(),
    outer_run_counter = 1L,
    outer_only_run_counter = 1L,
    outer_negative_streak = 0L,
    step_size = 1,
    weights_final = list(),
    pair_cache = list(pair = TRUE),
    margin_eval_cache = list(margin = TRUE),
    rs_calc_eta = function(...) list(),
    inner_stop_crit = 0.01,
    outer_stop_crit = 0.1,
    cg_grad_tol_eff = 0.2,
    cg_step_tol_eff = 0.3
  )
}

make_dispatch_args <- function(method, optimizer_context, cg_runner_fn, rs_runner_fn) {
  mm <- list(x = list(mu = matrix(1, nrow = 2, ncol = 1)), s = list(mu = list()))
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c(1, 1))

  list(
    method = method,
    optimizer_context = optimizer_context,
    mm = mm,
    margin_dist = "NO",
    copula_link = "identity",
    copula_dist = "N",
    dataset = dataset,
    include_dlcopdpar = FALSE,
    check_dlcopdpar_gradient = FALSE,
    discrete_score_method = "analytical",
    cg_gradient_method = "forward",
    cg_hessian_method = "analytical",
    verbose = 0,
    cg_max_delta = 1,
    max_outer_iter = 10L,
    check_elapsed_budget = function(stage) NULL,
    cg_update_lambda = FALSE,
    cg_max_lambda_updates = 0L,
    cg_lambda_update_every = 1L,
    lambda_penalty_K = 2,
    cg_armijo_c1 = 1e-4,
    cg_line_search = "best",
    cg_max_line_search_evals = 5L,
    use_backtracking = FALSE,
    backtracking_max_halves = 0L,
    cg_zeta_hessian = "analytical",
    cg_max_stall = 5L,
    cg_raw_loglik_drop_tol = 10,
    rs_smooth_trust_radius = Inf,
    rs_update_lambda = FALSE,
    plot_results = FALSE,
    true_val = NULL,
    max_inner_iter = 10L,
    max_negative_outer_streak = 5L,
    step_adjustment = 0.5,
    max_steps = 5L,
    start_step_size = 1,
    cg_runner_fn = cg_runner_fn,
    rs_runner_fn = rs_runner_fn
  )
}

test_that("fit optimizer dispatch routes CG state through the CG runner", {
  context <- make_optimizer_context()
  captured <- list()
  cg_state <- list(
    par_cov = c(`mu.(Intercept)` = 1),
    par_s = list(mu = list()),
    lambda_s = list(mu = list()),
    df_s = list(mu = 2),
    calc_lik_out_end = list(log_lik = c(marginal = -4, copula = -1, joint = -5)),
    outer_only_run_counter = 2L,
    outer_start_log_lik = -6,
    outer_end_log_lik = -5,
    outer_log_lik_change = 1,
    log_lik_history = matrix(c(-4, -1, -5), nrow = 1, dimnames = list(NULL, c("marginal", "copula", "joint"))),
    par_history = matrix(1, nrow = 1, dimnames = list(NULL, "mu.(Intercept)")),
    weights_final = list(mu = c(1, 1)),
    cg_stop_reason = "tolerance",
    cg_last_grad_inf = 0.01,
    cg_last_step_l2 = 0.02,
    cg_best_raw_loglik = -5,
    cg_best_iteration = 1L,
    cg_raw_loglik_drop_from_best = 0,
    cg_lambda_trace = data.frame(),
    cg_step_trace = list(data.frame(outer_iteration = 1L))
  )

  args <- make_dispatch_args(
    method = "CG",
    optimizer_context = context,
    cg_runner_fn = function(...) {
      captured$cg <<- list(...)
      cg_state
    },
    rs_runner_fn = function(...) stop("RS runner should not be called", call. = FALSE)
  )

  out <- do.call(gamlss.longitudinal:::.gl_run_fit_optimizer, args)

  expect_equal(captured$cg$par_cov, context$par_cov)
  expect_identical(captured$cg$pair_cache, context$pair_cache)
  expect_equal(captured$cg$outer_stop_crit, context$outer_stop_crit)
  expect_equal(captured$cg$cg_grad_tol_eff, context$cg_grad_tol_eff)
  expect_equal(out$par_cov, cg_state$par_cov)
  expect_equal(out$cg_stop_reason, "tolerance")
  expect_equal(out$rs_block_trace, context$rs_block_trace)
  expect_equal(out$outer_stop_crit, context$outer_stop_crit)
})

test_that("fit optimizer dispatch routes RS state through the RS runner", {
  context <- make_optimizer_context()
  captured <- list()
  rs_state <- list(
    par_cov = c(`mu.(Intercept)` = 2),
    par_s = list(mu = list()),
    lambda_s = list(mu = list()),
    df_s = list(mu = 3),
    calc_lik_out_end = list(log_lik = c(marginal = -3, copula = -1, joint = -4)),
    rs_block_trace = list(data.frame(parameter = "mu")),
    outer_run_counter = 2L,
    outer_only_run_counter = 2L,
    outer_end_log_lik = -4,
    outer_log_lik_change = 2,
    outer_negative_streak = 0L,
    step_size = 0.5,
    weights_final = list(mu = c(0.5, 0.5)),
    log_lik_history = matrix(c(-3, -1, -4), nrow = 1, dimnames = list(NULL, c("marginal", "copula", "joint"))),
    par_history = matrix(2, nrow = 1, dimnames = list(NULL, "mu.(Intercept)")),
    first_outer_run = FALSE,
    outer_start_log_lik = -6
  )

  args <- make_dispatch_args(
    method = "RS",
    optimizer_context = context,
    cg_runner_fn = function(...) stop("CG runner should not be called", call. = FALSE),
    rs_runner_fn = function(...) {
      captured$rs <<- list(...)
      rs_state
    }
  )

  out <- do.call(gamlss.longitudinal:::.gl_run_fit_optimizer, args)

  expect_equal(captured$rs$par_cov, context$par_cov)
  expect_identical(captured$rs$rs_design_cache, context$rs_design_cache)
  expect_identical(captured$rs$margin_eval_cache, context$margin_eval_cache)
  expect_equal(captured$rs$inner_stop_crit, context$inner_stop_crit)
  expect_equal(out$par_cov, rs_state$par_cov)
  expect_equal(out$rs_block_trace, rs_state$rs_block_trace)
  expect_equal(out$cg_stop_reason, context$cg_stop_reason)
  expect_s3_class(out$cg_lambda_trace, "data.frame")
  expect_equal(out$cg_step_trace, list())
  expect_equal(out$outer_stop_crit, context$outer_stop_crit)
})

test_that("prepared fit optimizer argument builder unpacks normalized workflow state", {
  context <- make_optimizer_context()
  mm <- list(x = list(mu = matrix(1, nrow = 2, ncol = 1)), s = list(mu = list()))
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c(1, 1))
  check_elapsed_budget <- function(stage) NULL

  fit_workflow <- list(
    controls = list(
      method = "CG",
      discrete_score_method = "finite",
      cg_gradient_method = "central",
      cg_hessian_method = "auto",
      cg_max_delta = 0.25,
      cg_max_lambda_updates = 3L,
      cg_lambda_update_every = 4L,
      cg_line_search = "first",
      cg_max_line_search_evals = 5L,
      backtracking_max_halves = 6L,
      cg_zeta_hessian = "finite",
      cg_max_stall = 7L,
      cg_raw_loglik_drop_tol = 8,
      rs_smooth_trust_radius = 9
    ),
    fit_data = list(dataset = dataset),
    matrix_bundle = list(mm = mm, copula_link = "identity"),
    step_controls = list(
      start_step_size = 0.5,
      max_steps = 11L,
      step_adjustment = 0.25
    ),
    optimizer_context = context
  )

  args <- gamlss.longitudinal:::.gl_prepared_fit_optimizer_args(
    fit_workflow = fit_workflow,
    margin_dist = "NO",
    copula_dist = "N",
    include_dlcopdpar = TRUE,
    check_dlcopdpar_gradient = TRUE,
    verbose = 2,
    cg_armijo_c1 = 1e-5,
    cg_update_lambda = TRUE,
    lambda_penalty_K = 3,
    rs_update_lambda = FALSE,
    plot_results = TRUE,
    true_val = c(mu = 1),
    max_outer_iter = 12L,
    max_inner_iter = 13L,
    max_negative_outer_streak = 14L,
    use_backtracking = TRUE,
    check_elapsed_budget = check_elapsed_budget
  )

  expect_equal(args$method, "CG")
  expect_identical(args$optimizer_context, context)
  expect_identical(args$mm, mm)
  expect_identical(args$dataset, dataset)
  expect_equal(args$discrete_score_method, "finite")
  expect_equal(args$cg_gradient_method, "central")
  expect_equal(args$cg_hessian_method, "auto")
  expect_equal(args$cg_max_delta, 0.25)
  expect_equal(args$cg_max_lambda_updates, 3L)
  expect_equal(args$cg_lambda_update_every, 4L)
  expect_equal(args$cg_line_search, "first")
  expect_equal(args$cg_max_line_search_evals, 5L)
  expect_equal(args$backtracking_max_halves, 6L)
  expect_equal(args$cg_zeta_hessian, "finite")
  expect_equal(args$cg_max_stall, 7L)
  expect_equal(args$cg_raw_loglik_drop_tol, 8)
  expect_equal(args$rs_smooth_trust_radius, 9)
  expect_equal(args$start_step_size, 0.5)
  expect_equal(args$max_steps, 11L)
  expect_equal(args$step_adjustment, 0.25)
  expect_identical(args$check_elapsed_budget, check_elapsed_budget)
})

test_that("prepared fit optimizer bridge forwards normalized workflow state", {
  context <- make_optimizer_context()
  mm <- list(x = list(mu = matrix(1, nrow = 2, ncol = 1)), s = list(mu = list()))
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c(1, 1))
  check_elapsed_budget <- function(stage) NULL

  fit_workflow <- list(
    controls = list(
      method = "CG",
      discrete_score_method = "finite",
      cg_gradient_method = "central",
      cg_hessian_method = "auto",
      cg_max_delta = 0.25,
      cg_max_lambda_updates = 3L,
      cg_lambda_update_every = 4L,
      cg_line_search = "first",
      cg_max_line_search_evals = 5L,
      backtracking_max_halves = 6L,
      cg_zeta_hessian = "finite",
      cg_max_stall = 7L,
      cg_raw_loglik_drop_tol = 8,
      rs_smooth_trust_radius = 9
    ),
    fit_data = list(dataset = dataset),
    matrix_bundle = list(mm = mm, copula_link = "identity"),
    step_controls = list(
      start_step_size = 0.5,
      max_steps = 11L,
      step_adjustment = 0.25
    ),
    optimizer_context = context
  )

  captured <- list()
  out <- gamlss.longitudinal:::.gl_run_prepared_fit_optimizer(
    fit_workflow = fit_workflow,
    margin_dist = "NO",
    copula_dist = "N",
    include_dlcopdpar = TRUE,
    check_dlcopdpar_gradient = TRUE,
    verbose = 2,
    cg_armijo_c1 = 1e-5,
    cg_update_lambda = TRUE,
    lambda_penalty_K = 3,
    rs_update_lambda = FALSE,
    plot_results = TRUE,
    true_val = c(mu = 1),
    max_outer_iter = 12L,
    max_inner_iter = 13L,
    max_negative_outer_streak = 14L,
    use_backtracking = TRUE,
    check_elapsed_budget = check_elapsed_budget,
    optimizer_fn = function(...) {
      captured$args <<- list(...)
      list(ok = TRUE)
    }
  )

  expect_identical(out, list(ok = TRUE))
  expect_equal(captured$args$method, "CG")
  expect_identical(captured$args$optimizer_context, context)
  expect_identical(captured$args$mm, mm)
  expect_identical(captured$args$dataset, dataset)
  expect_equal(captured$args$discrete_score_method, "finite")
  expect_equal(captured$args$cg_gradient_method, "central")
  expect_equal(captured$args$cg_hessian_method, "auto")
  expect_equal(captured$args$cg_max_delta, 0.25)
  expect_equal(captured$args$cg_max_lambda_updates, 3L)
  expect_equal(captured$args$cg_lambda_update_every, 4L)
  expect_equal(captured$args$cg_line_search, "first")
  expect_equal(captured$args$cg_max_line_search_evals, 5L)
  expect_equal(captured$args$backtracking_max_halves, 6L)
  expect_equal(captured$args$cg_zeta_hessian, "finite")
  expect_equal(captured$args$cg_max_stall, 7L)
  expect_equal(captured$args$cg_raw_loglik_drop_tol, 8)
  expect_equal(captured$args$rs_smooth_trust_radius, 9)
  expect_equal(captured$args$start_step_size, 0.5)
  expect_equal(captured$args$max_steps, 11L)
  expect_equal(captured$args$step_adjustment, 0.25)
  expect_identical(captured$args$check_elapsed_budget, check_elapsed_budget)
})
