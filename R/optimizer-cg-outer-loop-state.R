#' Apply one CG outer-iteration result to loop state
#'
#' @noRd
.gl_apply_cg_outer_iteration_state <- function(cg_outer_iteration) {
  list(
    beta_all = cg_outer_iteration$beta_all,
    par_cov = cg_outer_iteration$par_cov,
    par_s = cg_outer_iteration$par_s,
    lambda_s = cg_outer_iteration$lambda_s,
    df_s = cg_outer_iteration$df_s,
    penalty_mat = cg_outer_iteration$penalty_mat,
    cg_trust_radius = cg_outer_iteration$cg_trust_radius,
    cg_stall_count = cg_outer_iteration$cg_stall_count,
    cg_lambda_update_count = cg_outer_iteration$cg_lambda_update_count,
    cg_lambda_trace = cg_outer_iteration$cg_lambda_trace,
    cg_step_trace = cg_outer_iteration$cg_step_trace,
    cg_converged = cg_outer_iteration$cg_converged,
    calc_lik_out_end = cg_outer_iteration$calc_lik_out_end,
    outer_only_run_counter = cg_outer_iteration$outer_only_run_counter,
    outer_start_log_lik = cg_outer_iteration$outer_start_log_lik,
    outer_end_log_lik = cg_outer_iteration$outer_end_log_lik,
    outer_log_lik_change = cg_outer_iteration$outer_log_lik_change,
    log_lik_history = cg_outer_iteration$log_lik_history,
    par_history = cg_outer_iteration$par_history,
    cg_stop_reason = cg_outer_iteration$cg_stop_reason,
    cg_last_grad_inf = cg_outer_iteration$cg_last_grad_inf,
    cg_last_step_l2 = cg_outer_iteration$cg_last_step_l2,
    cg_best_raw_loglik = cg_outer_iteration$cg_best_raw_loglik,
    cg_best_iteration = cg_outer_iteration$cg_best_iteration,
    cg_raw_loglik_drop_from_best = cg_outer_iteration$cg_raw_loglik_drop_from_best
  )
}
