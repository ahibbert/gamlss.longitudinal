#' Assemble CG line-search diagnostics state
#'
#' @noRd
.gl_build_cg_line_search_diagnostics_state <- function(
    step_acceptance,
    step_state,
    outer_end_log_lik,
    iteration_diagnostics) {
  cg_stopping <- iteration_diagnostics$stopping

  list(
    step_acceptance = step_acceptance,
    step_state = step_state,
    beta = step_state$beta,
    par_cov = step_state$par_cov,
    par_s = step_state$par_s,
    calc_lik_out_end = step_state$calc_lik_out_end,
    stall_count = step_state$stall_count,
    trust_radius = step_state$trust_radius,
    accepted_improvement = step_state$accepted_improvement,
    best = step_state$best,
    prevented_deterioration = step_state$prevented_deterioration,
    prevented_raw_loglik_drop = step_state$prevented_raw_loglik_drop,
    outer_end_log_lik = outer_end_log_lik,
    iteration_diagnostics = iteration_diagnostics,
    outer_log_lik_change = iteration_diagnostics$outer_log_lik_change,
    best_raw_loglik = iteration_diagnostics$best_raw_loglik,
    best_iteration = iteration_diagnostics$best_iteration,
    raw_loglik_drop_from_best = iteration_diagnostics$raw_loglik_drop_from_best,
    grad_inf = iteration_diagnostics$grad_inf,
    step_l2 = iteration_diagnostics$step_l2,
    stopping = cg_stopping,
    tolerance_met = cg_stopping$tolerance_met,
    max_stall_hit = cg_stopping$max_stall_hit,
    deterioration_hit = cg_stopping$deterioration_hit,
    stop_requested = cg_stopping$stop_requested,
    step_trace = iteration_diagnostics$step_trace
  )
}
