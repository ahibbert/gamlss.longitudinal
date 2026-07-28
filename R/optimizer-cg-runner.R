#' Run the CG optimizer loop
#'
#' Owns CG outer loop and stopping. Numerical details stay delegated to the
#'  existing CG helper files, making this the reviewer entry point for 
#' CG optimizer flow.
#' 
#' Key helpers called are:
#' - .gl_ensure_cg_hessian_available (optimizer-cg-hessian-availability.R)
#' - .gl_build_cg_runtime_helpers (optimizer-cg-runtime.R)
#' - .gl_initialize_cg_optimizer_state (optimizer-cg-state-start.R)
#' - .gl_evaluate_cg_iteration_start (optimizer-cg-state-start.R)
#' - .gl_prepare_cg_curvature_line_search_state (optimizer-cg-curvature-search.R)
#' - .gl_apply_cg_line_search_diagnostics_state (optimizer-cg-line-search.R)
#' - .gl_apply_cg_stop_request_state (optimizer-cg-stop-request.R)
#' - .gl_finalize_cg_optimizer_state (optimizer-cg-finalization.R)
#' - .gl_report_cg_optimizer_start (optimizer-cg-reporting.R)
#' - .gl_build_cg_optimizer_result (optimizer-cg-result.R)
#' - .gl_prepare_cg_runtime_state (optimizer-cg-runtime-state.R)
#' - .gl_run_cg_outer_iteration_step (optimizer-cg-outer-iteration.R)
#' - .gl_run_cg_outer_loop (optimizer-cg-outer-loop.R)
#'
#' @noRd
.gl_run_cg_optimizer <- function(
    par_cov,
    par_s,
    mm,
    margin_dist,
    copula_link,
    copula_dist,
    dataset,
    pair_cache,
    margin_eval_cache,
    cg_gradient_method,
    cg_hessian_method,
    verbose,
    lambda_s,
    cg_max_delta,
    log_lik_history,
    par_history,
    cg_best_raw_loglik,
    cg_best_iteration,
    outer_only_run_counter,
    max_outer_iter,
    check_elapsed_budget,
    include_dlcopdpar,
    cg_update_lambda,
    cg_max_lambda_updates,
    cg_lambda_update_every,
    cg_step_tol_eff,
    lambda_penalty_K,
    cg_armijo_c1,
    cg_line_search,
    cg_max_line_search_evals,
    use_backtracking,
    backtracking_max_halves,
    cg_zeta_hessian,
    cg_max_stall,
    cg_raw_loglik_drop_tol,
    outer_stop_crit,
    cg_grad_tol_eff,
    df_s,
    outer_start_log_lik,
    outer_end_log_lik,
    outer_log_lik_change,
    cg_stop_reason,
    cg_last_grad_inf,
    cg_last_step_l2,
    cg_raw_loglik_drop_from_best,
    weights_final,
    ensure_hessian_fn = .gl_ensure_cg_hessian_available,
    runtime_helpers_fn = .gl_build_cg_runtime_helpers,
    initialize_state_fn = .gl_initialize_cg_optimizer_state,
    iteration_start_fn = .gl_evaluate_cg_iteration_start,
    curvature_line_search_fn = .gl_prepare_cg_curvature_line_search_state,
    line_search_diagnostics_fn = .gl_apply_cg_line_search_diagnostics_state,
    stop_request_state_fn = .gl_apply_cg_stop_request_state,
    finalize_state_fn = .gl_finalize_cg_optimizer_state,
    startup_report_fn = .gl_report_cg_optimizer_start,
    result_fn = .gl_build_cg_optimizer_result,
    runtime_state_fn = .gl_prepare_cg_runtime_state,
    outer_iteration_step_fn = .gl_run_cg_outer_iteration_step,
    outer_loop_fn = .gl_run_cg_outer_loop) {
  ensure_hessian_fn()

  startup_report_fn(
    verbose = verbose,
    cg_max_delta = cg_max_delta,
    cg_lambda_update_every = cg_lambda_update_every,
    cg_update_lambda = cg_update_lambda,
    cg_line_search = cg_line_search,
    cg_gradient_method = cg_gradient_method,
    cg_hessian_method = cg_hessian_method
  )

  cg_runtime_state <- runtime_state_fn(
    par_cov = par_cov,
    par_s = par_s,
    mm = mm,
    margin_dist = margin_dist,
    copula_link = copula_link,
    copula_dist = copula_dist,
    dataset = dataset,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache,
    cg_gradient_method = cg_gradient_method,
    cg_hessian_method = cg_hessian_method,
    verbose = verbose,
    lambda_s = lambda_s,
    cg_max_delta = cg_max_delta,
    runtime_helpers_fn = runtime_helpers_fn,
    initialize_state_fn = initialize_state_fn
  )
  cg_runtime <- cg_runtime_state$runtime
  build_cg_model <- cg_runtime$build_model
  build_cg_penalty <- cg_runtime$build_penalty
  cg_eval <- cg_runtime$evaluate
  cg_objective <- cg_runtime$objective
  cg_gradient <- cg_runtime$gradient
  cg_finite_hessian_block <- cg_runtime$finite_hessian_block
  cg_observed_hessian <- cg_runtime$observed_hessian
  cg_smooth_edf_list <- cg_runtime$smooth_edf_list

  cg_state <- cg_runtime_state$state
  mm_cg <- cg_state$mm_cg
  beta_all <- cg_state$beta_all
  penalty_mat <- cg_state$penalty_mat
  cg_trust_radius <- cg_state$cg_trust_radius
  cg_stall_count <- cg_state$cg_stall_count
  cg_converged <- cg_state$cg_converged
  cg_lambda_update_count <- cg_state$cg_lambda_update_count
  cg_has_smooths <- cg_state$cg_has_smooths
  cg_lambda_trace <- cg_state$cg_lambda_trace
  cg_step_trace <- cg_state$cg_step_trace

  cg_outer_loop <- outer_loop_fn(
    dataset = dataset,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    mm_cg = mm_cg,
    beta_all = beta_all,
    penalty_mat = penalty_mat,
    log_lik_history = log_lik_history,
    par_history = par_history,
    cg_best_raw_loglik = cg_best_raw_loglik,
    cg_best_iteration = cg_best_iteration,
    outer_only_run_counter = outer_only_run_counter,
    max_outer_iter = max_outer_iter,
    check_elapsed_budget = check_elapsed_budget,
    include_dlcopdpar = include_dlcopdpar,
    cg_gradient_method = cg_gradient_method,
    cg_eval = cg_eval,
    cg_objective = cg_objective,
    cg_gradient = cg_gradient,
    lambda_s = lambda_s,
    cg_lambda_trace = cg_lambda_trace,
    cg_lambda_update_count = cg_lambda_update_count,
    cg_update_lambda = cg_update_lambda,
    cg_max_lambda_updates = cg_max_lambda_updates,
    cg_lambda_update_every = cg_lambda_update_every,
    cg_trust_radius = cg_trust_radius,
    cg_max_delta = cg_max_delta,
    cg_step_tol_eff = cg_step_tol_eff,
    build_cg_penalty = build_cg_penalty,
    cg_smooth_edf_list = cg_smooth_edf_list,
    lambda_penalty_K = lambda_penalty_K,
    cg_armijo_c1 = cg_armijo_c1,
    cg_line_search = cg_line_search,
    cg_max_line_search_evals = cg_max_line_search_evals,
    use_backtracking = use_backtracking,
    backtracking_max_halves = backtracking_max_halves,
    verbose = verbose,
    cg_observed_hessian = cg_observed_hessian,
    cg_finite_hessian_block = cg_finite_hessian_block,
    cg_zeta_hessian = cg_zeta_hessian,
    par_cov = par_cov,
    par_s = par_s,
    cg_stall_count = cg_stall_count,
    cg_max_stall = cg_max_stall,
    cg_raw_loglik_drop_tol = cg_raw_loglik_drop_tol,
    cg_step_trace = cg_step_trace,
    outer_stop_crit = outer_stop_crit,
    cg_grad_tol_eff = cg_grad_tol_eff,
    df_s = df_s,
    cg_has_smooths = cg_has_smooths,
    cg_converged = cg_converged,
    outer_iteration_step_fn = outer_iteration_step_fn,
    iteration_start_fn = iteration_start_fn,
    curvature_line_search_fn = curvature_line_search_fn,
    line_search_diagnostics_fn = line_search_diagnostics_fn,
    stop_request_state_fn = stop_request_state_fn
  )
  beta_all <- cg_outer_loop$beta_all
  par_cov <- cg_outer_loop$par_cov
  par_s <- cg_outer_loop$par_s
  lambda_s <- cg_outer_loop$lambda_s
  df_s <- cg_outer_loop$df_s
  penalty_mat <- cg_outer_loop$penalty_mat
  calc_lik_out_end <- cg_outer_loop$calc_lik_out_end
  outer_only_run_counter <- cg_outer_loop$outer_only_run_counter
  outer_start_log_lik <- cg_outer_loop$outer_start_log_lik
  outer_end_log_lik <- cg_outer_loop$outer_end_log_lik
  outer_log_lik_change <- cg_outer_loop$outer_log_lik_change
  log_lik_history <- cg_outer_loop$log_lik_history
  par_history <- cg_outer_loop$par_history
  cg_stop_reason <- cg_outer_loop$cg_stop_reason
  cg_last_grad_inf <- cg_outer_loop$cg_last_grad_inf
  cg_last_step_l2 <- cg_outer_loop$cg_last_step_l2
  cg_best_raw_loglik <- cg_outer_loop$cg_best_raw_loglik
  cg_best_iteration <- cg_outer_loop$cg_best_iteration
  cg_raw_loglik_drop_from_best <- cg_outer_loop$cg_raw_loglik_drop_from_best
  cg_lambda_trace <- cg_outer_loop$cg_lambda_trace
  cg_step_trace <- cg_outer_loop$cg_step_trace

  cg_final_state <- finalize_state_fn(
    dataset = dataset,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    mm = mm,
    mm_cg = mm_cg,
    beta_vec = beta_all,
    lambda_current = lambda_s,
    penalty_current = penalty_mat,
    df_s_current = df_s,
    observed_hessian_fn = cg_observed_hessian,
    build_penalty_fn = build_cg_penalty,
    edf_fn = cg_smooth_edf_list
  )
  final_edf <- cg_final_state$final_edf
  penalty_mat <- cg_final_state$penalty_mat
  df_s <- cg_final_state$df_s
  weights_final <- cg_final_state$weights_final

  result_fn(
    par_cov = par_cov,
    par_s = par_s,
    lambda_s = lambda_s,
    df_s = df_s,
    calc_lik_out_end = calc_lik_out_end,
    outer_only_run_counter = outer_only_run_counter,
    outer_start_log_lik = outer_start_log_lik,
    outer_end_log_lik = outer_end_log_lik,
    outer_log_lik_change = outer_log_lik_change,
    log_lik_history = log_lik_history,
    par_history = par_history,
    weights_final = weights_final,
    cg_stop_reason = cg_stop_reason,
    cg_last_grad_inf = cg_last_grad_inf,
    cg_last_step_l2 = cg_last_step_l2,
    cg_best_raw_loglik = cg_best_raw_loglik,
    cg_best_iteration = cg_best_iteration,
    cg_raw_loglik_drop_from_best = cg_raw_loglik_drop_from_best,
    cg_lambda_trace = cg_lambda_trace,
    cg_step_trace = cg_step_trace,
    final_edf = final_edf,
    penalty_mat = penalty_mat
  )
}
