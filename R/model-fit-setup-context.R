#' This function collects all the fit criteria into a list that can be passed to the optimisers.
#'
#' @noRd
.gl_fit_optimizer_context_bundle <- function(
    initial_state,
    rs_design_cache,
    rs_calc_eta,
    optimizer_state,
    stop_criteria) {
  c(
    initial_state,
    list(
      rs_design_cache = rs_design_cache,
      rs_calc_eta = rs_calc_eta,
      first_outer_run = optimizer_state$first_outer_run,
      outer_log_lik_change = optimizer_state$outer_log_lik_change,
      outer_start_log_lik = optimizer_state$outer_start_log_lik,
      outer_end_log_lik = optimizer_state$outer_end_log_lik,
      log_lik_history = optimizer_state$log_lik_history,
      par_history = optimizer_state$par_history,
      cg_stop_reason = optimizer_state$cg_stop_reason,
      cg_last_grad_inf = optimizer_state$cg_last_grad_inf,
      cg_last_step_l2 = optimizer_state$cg_last_step_l2,
      cg_best_raw_loglik = optimizer_state$cg_best_raw_loglik,
      cg_best_iteration = optimizer_state$cg_best_iteration,
      cg_raw_loglik_drop_from_best = optimizer_state$cg_raw_loglik_drop_from_best,
      rs_block_trace = optimizer_state$rs_block_trace,
      outer_run_counter = optimizer_state$outer_run_counter,
      outer_only_run_counter = optimizer_state$outer_only_run_counter,
      outer_negative_streak = optimizer_state$outer_negative_streak,
      step_size = optimizer_state$step_size,
      weights_final = optimizer_state$weights_final,
      pair_cache = optimizer_state$pair_cache,
      margin_eval_cache = optimizer_state$margin_eval_cache,
      inner_stop_crit = stop_criteria$inner_stop_crit,
      outer_stop_crit = stop_criteria$outer_stop_crit,
      cg_grad_tol_eff = stop_criteria$cg_grad_tol_eff,
      cg_step_tol_eff = stop_criteria$cg_step_tol_eff
    )
  )
}
