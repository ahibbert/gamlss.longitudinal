#' Run the RS optimizer branch
#'
#' @keywords internal
#' @noRd
.gl_run_fit_rs_optimizer <- function(
    optimizer_context,
    mm,
    margin_dist,
    copula_link,
    copula_dist,
    dataset,
    include_dlcopdpar,
    check_dlcopdpar_gradient,
    discrete_score_method,
    verbose,
    max_outer_iter,
    check_elapsed_budget,
    lambda_penalty_K,
    use_backtracking,
    backtracking_max_halves,
    rs_smooth_trust_radius,
    rs_update_lambda,
    plot_results,
    true_val,
    max_inner_iter,
    max_negative_outer_streak,
    step_adjustment,
    max_steps,
    start_step_size,
    rs_runner_fn = .gl_run_rs_optimizer) {
  optimizer_state <- rs_runner_fn(
    first_outer_run = optimizer_context$first_outer_run,
    outer_log_lik_change = optimizer_context$outer_log_lik_change,
    outer_stop_crit = optimizer_context$outer_stop_crit,
    outer_only_run_counter = optimizer_context$outer_only_run_counter,
    max_outer_iter = max_outer_iter,
    check_elapsed_budget = check_elapsed_budget,
    verbose = verbose,
    mm = mm,
    rs_calc_eta = optimizer_context$rs_calc_eta,
    par_cov = optimizer_context$par_cov,
    par_s = optimizer_context$par_s,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    dataset = dataset,
    pair_cache = optimizer_context$pair_cache,
    margin_eval_cache = optimizer_context$margin_eval_cache,
    outer_start_log_lik = optimizer_context$outer_start_log_lik,
    log_lik_history = optimizer_context$log_lik_history,
    par_history = optimizer_context$par_history,
    copula_link = copula_link,
    discrete_score_method = discrete_score_method,
    include_dlcopdpar = include_dlcopdpar,
    check_dlcopdpar_gradient = check_dlcopdpar_gradient,
    rs_design_cache = optimizer_context$rs_design_cache,
    rs_smooth_trust_radius = rs_smooth_trust_radius,
    lambda_penalty_K = lambda_penalty_K,
    lambda_s = optimizer_context$lambda_s,
    df_s = optimizer_context$df_s,
    step_size = optimizer_context$step_size,
    rs_update_lambda = rs_update_lambda,
    outer_run_counter = optimizer_context$outer_run_counter,
    rs_block_trace = optimizer_context$rs_block_trace,
    use_backtracking = use_backtracking,
    backtracking_max_halves = backtracking_max_halves,
    plot_results = plot_results,
    true_val = true_val,
    inner_stop_crit = optimizer_context$inner_stop_crit,
    max_inner_iter = max_inner_iter,
    outer_negative_streak = optimizer_context$outer_negative_streak,
    step_adjustment = step_adjustment,
    max_steps = max_steps,
    start_step_size = start_step_size,
    max_negative_outer_streak = max_negative_outer_streak,
    weights_final = optimizer_context$weights_final
  )
  optimizer_state$cg_stop_reason <- optimizer_context$cg_stop_reason
  optimizer_state$cg_last_grad_inf <- optimizer_context$cg_last_grad_inf
  optimizer_state$cg_last_step_l2 <- optimizer_context$cg_last_step_l2
  optimizer_state$cg_best_raw_loglik <- optimizer_context$cg_best_raw_loglik
  optimizer_state$cg_best_iteration <- optimizer_context$cg_best_iteration
  optimizer_state$cg_raw_loglik_drop_from_best <- optimizer_context$cg_raw_loglik_drop_from_best
  optimizer_state$cg_lambda_trace <- data.frame()
  optimizer_state$cg_step_trace <- list()
  optimizer_state$outer_stop_crit <- optimizer_context$outer_stop_crit
  optimizer_state$inner_stop_crit <- optimizer_context$inner_stop_crit
  optimizer_state$cg_grad_tol_eff <- optimizer_context$cg_grad_tol_eff
  optimizer_state$cg_step_tol_eff <- optimizer_context$cg_step_tol_eff
  optimizer_state
}
