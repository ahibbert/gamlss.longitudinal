#' Update CG iteration diagnostics, stopping state, and trace
#'
#' @noRd
.gl_update_cg_iteration_diagnostics <- function(
    step_trace,
    g_pen,
    best,
    best_raw_loglik,
    best_iteration,
    outer_iteration,
    outer_start_log_lik,
    outer_end_log_lik,
    obj_start,
    accepted_improvement,
    trust_radius_start,
    trust_radius_end,
    line_eval_count,
    lambda_update_count,
    lambda_changed,
    stall_count,
    prevented_deterioration,
    prevented_raw_loglik_drop,
    max_stall,
    raw_loglik_drop_tol,
    outer_stop_crit,
    grad_tol,
    step_tol,
    verbose) {
  outer_log_lik_change <- outer_end_log_lik - outer_start_log_lik

  best_state <- .gl_update_cg_best_loglik(
    candidate_loglik = outer_end_log_lik,
    best_raw_loglik = best_raw_loglik,
    best_iteration = best_iteration,
    current_iteration = outer_iteration
  )
  best_raw_loglik <- best_state$best_raw_loglik
  best_iteration <- best_state$best_iteration

  raw_loglik_drop_from_best <- .gl_cg_raw_loglik_drop(
    best_raw_loglik = best_raw_loglik,
    current_loglik = outer_end_log_lik,
    prevented_deterioration = prevented_deterioration,
    prevented_raw_loglik_drop = prevented_raw_loglik_drop
  )

  .gl_print_outer_iteration_summary(
    outer_start_log_lik = outer_start_log_lik,
    outer_end_log_lik = outer_end_log_lik,
    outer_log_lik_change = outer_log_lik_change,
    verbose = verbose
  )

  grad_inf <- max(abs(g_pen), na.rm = TRUE)
  step_l2 <- if (is.null(best)) 0 else best$step_l2

  stopping <- .gl_assess_cg_stopping(
    outer_log_lik_change = outer_log_lik_change,
    grad_inf = grad_inf,
    step_l2 = step_l2,
    stall_count = stall_count,
    max_stall = max_stall,
    raw_loglik_drop_from_best = raw_loglik_drop_from_best,
    raw_loglik_drop_tol = raw_loglik_drop_tol,
    lambda_update_count = lambda_update_count,
    prevented_deterioration = prevented_deterioration,
    outer_stop_crit = outer_stop_crit,
    grad_tol = grad_tol,
    step_tol = step_tol
  )

  trace_row <- .gl_build_cg_step_trace_row(
    outer_iteration = outer_iteration,
    start_loglik = outer_start_log_lik,
    end_loglik = outer_end_log_lik,
    raw_loglik_change = outer_log_lik_change,
    start_penalized_loglik = obj_start,
    accepted_penalized_improvement = accepted_improvement,
    grad_inf = grad_inf,
    step_l2 = step_l2,
    trust_radius_start = trust_radius_start,
    trust_radius_end = trust_radius_end,
    line_search_evals = line_eval_count,
    accepted_step = !is.null(best),
    lambda_update_count = lambda_update_count,
    lambda_changed = lambda_changed,
    stall_count = stall_count,
    tolerance_met = stopping$tolerance_met,
    max_stall_hit = stopping$max_stall_hit,
    raw_deterioration_hit = stopping$deterioration_hit,
    raw_loglik_drop_from_best = raw_loglik_drop_from_best
  )
  step_trace[[length(step_trace) + 1L]] <- trace_row

  list(
    outer_log_lik_change = outer_log_lik_change,
    best_raw_loglik = best_raw_loglik,
    best_iteration = best_iteration,
    raw_loglik_drop_from_best = raw_loglik_drop_from_best,
    grad_inf = grad_inf,
    step_l2 = step_l2,
    stopping = stopping,
    step_trace = step_trace
  )
}
