test_that("CG runtime-state helper wires runtime helpers into initializer", {
  calls <- character()
  captured <- list()

  record <- function(label) {
    calls <<- c(calls, label)
  }

  mm <- list(x = list(mu = matrix(1, nrow = 2, ncol = 1)), s = list(mu = list()))
  par_cov <- c(`mu.(Intercept)` = 0)
  par_s <- list(mu = list())
  lambda_s <- list(mu = list())
  runtime <- list(
    build_model = function(...) list(),
    build_penalty = function(...) matrix(0, 1, 1),
    evaluate = function(...) list(),
    objective = function(...) 0,
    gradient = function(...) numeric(),
    finite_hessian_block = function(...) matrix(0, 1, 1),
    observed_hessian = function(...) matrix(0, 1, 1),
    smooth_edf_list = function(...) list()
  )
  state <- list(mm_cg = list(), beta_all = par_cov, cg_trust_radius = 1)

  out <- gamlss.longitudinal:::.gl_prepare_cg_runtime_state(
    par_cov = par_cov,
    par_s = par_s,
    mm = mm,
    margin_dist = "NO",
    copula_link = "identity",
    copula_dist = "N",
    dataset = data.frame(response = 1:2),
    pair_cache = list(pair = TRUE),
    margin_eval_cache = list(margin = TRUE),
    cg_gradient_method = "forward",
    cg_hessian_method = "analytical",
    verbose = 0,
    lambda_s = lambda_s,
    cg_max_delta = 0.5,
    runtime_helpers_fn = function(...) {
      record("runtime")
      captured$runtime <<- list(...)
      runtime
    },
    initialize_state_fn = function(...) {
      record("initialize")
      captured$initialize <<- list(...)
      state
    }
  )

  expect_equal(calls, c("runtime", "initialize"))
  expect_identical(captured$runtime$mm, mm)
  expect_equal(captured$runtime$cg_gradient_method, "forward")
  expect_identical(captured$initialize$build_model_fn, runtime$build_model)
  expect_identical(captured$initialize$build_penalty_fn, runtime$build_penalty)
  expect_equal(captured$initialize$cg_max_delta, 0.5)
  expect_identical(out$runtime, runtime)
  expect_identical(out$state, state)
})

test_that("CG optimizer runner wires helper sequence and returns updated state", {
  calls <- character()

  record <- function(label) {
    calls <<- c(calls, label)
  }

  mm <- list(
    x = list(mu = matrix(1, nrow = 2, ncol = 1)),
    s = list(mu = list())
  )
  colnames(mm$x$mu) <- "mu.(Intercept)"

  out <- gamlss.longitudinal:::.gl_run_cg_optimizer(
    par_cov = c(`mu.(Intercept)` = 0),
    par_s = list(mu = list()),
    mm = mm,
    margin_dist = "NO",
    copula_link = "identity",
    copula_dist = "Gaussian",
    dataset = data.frame(response = 1:2, time = 1:2, subject = c(1, 1)),
    pair_cache = list(pair = TRUE),
    margin_eval_cache = list(margin = TRUE),
    cg_gradient_method = "finite",
    cg_hessian_method = "finite",
    verbose = 0,
    lambda_s = list(mu = list()),
    cg_max_delta = 1,
    log_lik_history = matrix(ncol = 3, nrow = 0),
    par_history = matrix(ncol = 1, nrow = 0, dimnames = list(NULL, "mu.(Intercept)")),
    cg_best_raw_loglik = -Inf,
    cg_best_iteration = NA_integer_,
    outer_only_run_counter = 1L,
    max_outer_iter = 10L,
    check_elapsed_budget = function(stage) record(paste("budget", stage)),
    include_dlcopdpar = FALSE,
    cg_update_lambda = FALSE,
    cg_max_lambda_updates = 0L,
    cg_lambda_update_every = 1L,
    cg_step_tol_eff = 1e-4,
    lambda_penalty_K = list(),
    cg_armijo_c1 = 1e-4,
    cg_line_search = "best",
    cg_max_line_search_evals = 5L,
    use_backtracking = FALSE,
    backtracking_max_halves = 0L,
    cg_zeta_hessian = "analytical",
    cg_max_stall = 5L,
    cg_raw_loglik_drop_tol = Inf,
    outer_stop_crit = 1e-4,
    cg_grad_tol_eff = 1e-4,
    df_s = list(mu = 1),
    outer_start_log_lik = 0,
    outer_end_log_lik = 0,
    outer_log_lik_change = 0,
    cg_stop_reason = NA_character_,
    cg_last_grad_inf = NA_real_,
    cg_last_step_l2 = NA_real_,
    cg_raw_loglik_drop_from_best = NA_real_,
    weights_final = list(),
    ensure_hessian_fn = function() record("ensure-hessian"),
    runtime_helpers_fn = function(...) {
      record("runtime-helpers")
      list(
        build_model = function(mm, par_cov, par_s) {
          record("build-model")
          list(mm = mm, beta = par_cov)
        },
        build_penalty = function(beta_names, lambda_current) {
          record("build-penalty")
          matrix(0, nrow = length(beta_names), ncol = length(beta_names),
                 dimnames = list(beta_names, beta_names))
        },
        evaluate = function(beta_vec, mm_cg) list(),
        objective = function(beta_vec, loglik, penalty_current) as.numeric(loglik),
        gradient = function(beta_vec, base_ll, mm_cg) rep(0, length(beta_vec)),
        finite_hessian_block = function(...) matrix(0, 1, 1),
        observed_hessian = function(...) diag(1),
        smooth_edf_list = function(...) list(mu = 1)
      )
    },
    initialize_state_fn = .gl_initialize_cg_optimizer_state,
    iteration_start_fn = function(
      beta_vec,
      mm_cg,
      penalty_current,
      log_lik_history,
      par_history,
      best_raw_loglik,
      best_iteration,
      current_iteration,
      margin_dist,
      copula_dist,
      include_dlcopdpar,
      dataset,
      cg_gradient_method,
      eval_fn,
      objective_fn,
      finite_gradient_fn
    ) {
      record("iteration-start")
      list(
        eval_start = list(calc_lik = list(log_lik = c(marginal = -5, copula = -1, joint = -6))),
        log_lik_history = rbind(log_lik_history, c(marginal = -5, copula = -1, joint = -6)),
        par_history = rbind(par_history, beta_vec),
        outer_start_log_lik = -6,
        best_raw_loglik = -6,
        best_iteration = current_iteration,
        obj_start = -6,
        grad = c(`mu.(Intercept)` = 0.1)
      )
    },
    curvature_line_search_fn = function(
      dataset,
      margin_dist,
      copula_dist,
      mm_cg,
      beta_vec,
      grad_vec,
      lambda_current,
      lambda_trace,
      penalty_current,
      lambda_update_count,
      update_lambda,
      max_lambda_updates,
      lambda_update_every,
      outer_iteration,
      trust_radius,
      max_delta,
      step_tol,
      build_penalty_fn,
      eval_fn,
      edf_fn,
      objective_fn,
      lambda_penalty_K,
      obj_start,
      armijo_c1,
      line_search,
      max_line_search_evals,
      use_backtracking,
      backtracking_max_halves,
      verbose,
      observed_hessian_fn,
      finite_hessian_fn,
      cg_zeta_hessian
    ) {
      record("curvature-line-search")
      list(
        tmp_obj = list(),
        H_obs = diag(1),
        H_zeta_fd = NULL,
        lambda_s = list(mu = list()),
        cg_lambda_trace = data.frame(iteration = integer()),
        penalty_mat = matrix(0, 1, 1, dimnames = list("mu.(Intercept)", "mu.(Intercept)")),
        cg_trust_radius = 1,
        cg_lambda_update_count = 0L,
        lambda_changed = FALSE,
        df_s = list(mu = 1.25),
        g_pen = c(`mu.(Intercept)` = 0.05),
        H_pen = matrix(1, 1, 1),
        candidate_steps = list(0.1),
        line_search_out = list(best = list(), line_eval_count = 1L),
        best = list(beta = c(`mu.(Intercept)` = 0.2), improvement = 1, step_l2 = 0.2),
        line_eval_count = 1L
      )
    },
    line_search_diagnostics_fn = function(
      best,
      eval_start,
      beta_vec,
      par_cov_template,
      par_s_template,
      stall_count,
      trust_radius,
      step_tol,
      max_delta,
      max_stall,
      verbose,
      best_raw_loglik,
      best_iteration,
      outer_iteration,
      outer_start_log_lik,
      obj_start,
      raw_loglik_drop_tol,
      step_trace,
      g_pen,
      trust_radius_start,
      line_eval_count,
      lambda_update_count,
      lambda_changed,
      outer_stop_crit,
      grad_tol,
      stop_on_convergence
    ) {
      record("line-search-diagnostics")
      list(
        step_acceptance = list(has_step = TRUE, accept = TRUE),
        step_state = list(),
        beta = c(`mu.(Intercept)` = 0.2),
        par_cov = c(`mu.(Intercept)` = 0.2),
        par_s = list(mu = list()),
        calc_lik_out_end = list(log_lik = c(marginal = -4, copula = -1, joint = -5)),
        stall_count = 0L,
        trust_radius = 1,
        accepted_improvement = 1,
        best = best,
        prevented_deterioration = FALSE,
        prevented_raw_loglik_drop = NA_real_,
        outer_end_log_lik = -5,
        iteration_diagnostics = list(),
        outer_log_lik_change = 1,
        best_raw_loglik = -5,
        best_iteration = outer_iteration,
        raw_loglik_drop_from_best = 0,
        grad_inf = 0.05,
        step_l2 = 0.2,
        stopping = list(stop_requested = TRUE),
        tolerance_met = TRUE,
        max_stall_hit = FALSE,
        deterioration_hit = FALSE,
        stop_requested = TRUE,
        step_trace = list(data.frame(outer_iteration = outer_iteration))
      )
    },
    stop_request_state_fn = function(
      update_lambda,
      has_smooths,
      lambda_update_count,
      lambda_current,
      lambda_trace,
      penalty_current,
      df_s_current,
      stall_count,
      H_obs_current,
      beta_vec,
      grad_vec,
      mm_cg,
      trust_radius,
      outer_iteration,
      max_delta,
      build_penalty_fn,
      eval_fn,
      edf_fn,
      objective_fn,
      lambda_penalty_K,
      tolerance_met,
      deterioration_hit,
      raw_loglik_drop_from_best,
      verbose
    ) {
      record("stop-request")
      list(
        stop_state = list(converged = TRUE),
        lambda = lambda_current,
        lambda_trace = lambda_trace,
        lambda_update_count = lambda_update_count,
        penalty_mat = penalty_current,
        df_s = df_s_current,
        stall_count = stall_count,
        stop_reason = "tolerance",
        converged = TRUE
      )
    },
    finalize_state_fn = function(
      dataset,
      margin_dist,
      copula_dist,
      mm,
      mm_cg,
      beta_vec,
      lambda_current,
      penalty_current,
      df_s_current,
      observed_hessian_fn,
      build_penalty_fn,
      edf_fn
    ) {
      record("finalize")
      list(
        final_edf = list(final = TRUE),
        penalty_mat = penalty_current,
        df_s = df_s_current,
        weights_final = list(mu = c(1, 1))
      )
    },
    startup_report_fn = function(...) {
      record("startup-report")
    },
    result_fn = function(...) {
      record("result")
      .gl_build_cg_optimizer_result(...)
    }
  )

  expect_equal(
    calls,
    c(
      "ensure-hessian",
      "startup-report",
      "runtime-helpers",
      "build-model",
      "build-penalty",
      "budget CG outer iteration",
      "iteration-start",
      "curvature-line-search",
      "line-search-diagnostics",
      "stop-request",
      "finalize",
      "result"
    )
  )
  expect_equal(out$par_cov, c(`mu.(Intercept)` = 0.2))
  expect_equal(out$calc_lik_out_end$log_lik["joint"], c(joint = -5))
  expect_equal(out$outer_only_run_counter, 2L)
  expect_equal(out$outer_log_lik_change, 1)
  expect_equal(out$cg_stop_reason, "tolerance")
  expect_equal(out$cg_last_grad_inf, 0.05)
  expect_equal(out$cg_last_step_l2, 0.2)
  expect_equal(out$weights_final, list(mu = c(1, 1)))
  expect_equal(nrow(out$log_lik_history), 1)
  expect_equal(nrow(out$par_history), 1)
})

test_that("CG outer-loop state helper preserves iteration state fields", {
  state <- list(
    beta_all = c(beta = 1),
    par_cov = c(mu = 2),
    par_s = list(mu = list()),
    lambda_s = list(mu = list(`s(x)` = 3)),
    df_s = list(mu = 1.5),
    penalty_mat = matrix(0, 1, 1),
    cg_trust_radius = 0.5,
    cg_stall_count = 1L,
    cg_lambda_update_count = 2L,
    cg_lambda_trace = data.frame(iter = 1L),
    cg_step_trace = list(data.frame(step = 1)),
    cg_converged = TRUE,
    calc_lik_out_end = list(log_lik = c(joint = -1)),
    outer_only_run_counter = 3L,
    outer_start_log_lik = -2,
    outer_end_log_lik = -1,
    outer_log_lik_change = 1,
    log_lik_history = matrix(-1, nrow = 1),
    par_history = matrix(2, nrow = 1),
    cg_stop_reason = "tolerance",
    cg_last_grad_inf = 0.01,
    cg_last_step_l2 = 0.02,
    cg_best_raw_loglik = -1,
    cg_best_iteration = 2L,
    cg_raw_loglik_drop_from_best = 0
  )

  expect_identical(.gl_apply_cg_outer_iteration_state(state), state)
})

test_that("CG optimizer result helper preserves returned-state contract", {
  result <- .gl_build_cg_optimizer_result(
    par_cov = c(mu = 1),
    par_s = list(mu = list()),
    lambda_s = list(mu = list()),
    df_s = list(mu = 1),
    calc_lik_out_end = list(log_lik = c(joint = -2)),
    outer_only_run_counter = 3L,
    outer_start_log_lik = -3,
    outer_end_log_lik = -2,
    outer_log_lik_change = 1,
    log_lik_history = matrix(-2, nrow = 1),
    par_history = matrix(1, nrow = 1),
    weights_final = list(mu = c(1, 1)),
    cg_stop_reason = "tolerance",
    cg_last_grad_inf = 0.01,
    cg_last_step_l2 = 0.02,
    cg_best_raw_loglik = -2,
    cg_best_iteration = 2L,
    cg_raw_loglik_drop_from_best = 0,
    cg_lambda_trace = data.frame(),
    cg_step_trace = list(data.frame(step = 1)),
    final_edf = list(mu = 1),
    penalty_mat = matrix(0, nrow = 1, ncol = 1)
  )

  expect_named(
    result,
    c(
      "par_cov",
      "par_s",
      "lambda_s",
      "df_s",
      "calc_lik_out_end",
      "outer_only_run_counter",
      "outer_start_log_lik",
      "outer_end_log_lik",
      "outer_log_lik_change",
      "log_lik_history",
      "par_history",
      "weights_final",
      "cg_stop_reason",
      "cg_last_grad_inf",
      "cg_last_step_l2",
      "cg_best_raw_loglik",
      "cg_best_iteration",
      "cg_raw_loglik_drop_from_best",
      "cg_lambda_trace",
      "cg_step_trace",
      "final_edf",
      "penalty_mat"
    )
  )
  expect_equal(result$cg_stop_reason, "tolerance")
  expect_equal(result$outer_log_lik_change, 1)
  expect_equal(result$weights_final, list(mu = c(1, 1)))
})

test_that("CG optimizer startup report is quiet unless verbose", {
  quiet <- capture.output(
    .gl_report_cg_optimizer_start(
      verbose = 0,
      cg_max_delta = 0.5,
      cg_lambda_update_every = 10L,
      cg_update_lambda = TRUE,
      cg_line_search = "best",
      cg_gradient_method = "forward",
      cg_hessian_method = "analytical"
    )
  )
  expect_length(quiet, 0)

  noisy <- capture.output(
    .gl_report_cg_optimizer_start(
      verbose = 1,
      cg_max_delta = 0.5,
      cg_lambda_update_every = 10L,
      cg_update_lambda = TRUE,
      cg_line_search = "best",
      cg_gradient_method = "forward",
      cg_hessian_method = "analytical"
    )
  )
  expect_match(paste(noisy, collapse = "\n"), "Using optimization method: CG")
  expect_match(paste(noisy, collapse = "\n"), "line_search=best")
  expect_match(paste(noisy, collapse = "\n"), "gradient=forward")
  expect_match(paste(noisy, collapse = "\n"), "hessian=analytical")
})
