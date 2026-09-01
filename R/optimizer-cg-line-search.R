#' Evaluate CG line-search candidate steps
#'
#' @noRd
.gl_evaluate_cg_line_search <- function(
    candidate_steps,
    beta_vec,
    mm_cg,
    penalty_current,
    obj_start,
    trust_radius,
    max_delta,
    armijo_c1,
    line_search,
    max_line_search_evals,
    max_backtrack,
    eval_fn,
    objective_fn,
    candidate_eval_fn = .gl_evaluate_cg_line_search_candidate) {
  best <- NULL
  line_eval_count <- 0L
  stop_line_search <- FALSE

  for (delta0 in candidate_steps) {
    if (stop_line_search) break
    for (bt in seq_len(max_backtrack + 1L)) {
      if (line_eval_count >= max_line_search_evals) {
        stop_line_search <- TRUE
        break
      }
      line_eval_count <- line_eval_count + 1L
      candidate <- candidate_eval_fn(
        delta0 = delta0,
        backtrack_index = bt,
        beta_vec = beta_vec,
        mm_cg = mm_cg,
        penalty_current = penalty_current,
        obj_start = obj_start,
        trust_radius = trust_radius,
        max_delta = max_delta,
        eval_fn = eval_fn,
        objective_fn = objective_fn
      )
      if (is.null(candidate)) next
      if (is.finite(candidate$improvement) && candidate$improvement > max(1e-8, armijo_c1 * max(1, abs(obj_start)))) {
        if (is.null(best) || candidate$improvement > best$improvement) {
          best <- candidate
        }
        if (identical(line_search, "first")) {
          stop_line_search <- TRUE
          break
        }
      }
    }
  }

  list(best = best, line_eval_count = line_eval_count)
}

#' Run one CG line-search step selection
#'
#' @noRd
.gl_run_cg_line_search <- function(
    candidate_steps,
    beta_vec,
    mm_cg,
    penalty_current,
    obj_start,
    trust_radius,
    max_delta,
    armijo_c1,
    line_search,
    max_line_search_evals,
    use_backtracking,
    backtracking_max_halves,
    eval_fn,
    objective_fn,
    verbose,
    line_search_fn = .gl_evaluate_cg_line_search) {
  max_backtrack <- if (isTRUE(use_backtracking)) as.integer(backtracking_max_halves) else 0L

  line_search_out <- line_search_fn(
    candidate_steps = candidate_steps,
    beta_vec = beta_vec,
    mm_cg = mm_cg,
    penalty_current = penalty_current,
    obj_start = obj_start,
    trust_radius = trust_radius,
    max_delta = max_delta,
    armijo_c1 = armijo_c1,
    line_search = line_search,
    max_line_search_evals = max_line_search_evals,
    max_backtrack = max_backtrack,
    eval_fn = eval_fn,
    objective_fn = objective_fn
  )

  if (verbose > 1) {
    cat(paste0("\nCG line search likelihood evaluations: ", line_search_out$line_eval_count))
  }

  line_search_out
}

#' Apply CG line-search result and update iteration diagnostics
#'
#' @noRd
.gl_apply_cg_line_search_diagnostics_state <- function(
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
    lambda_update_count,
    step_trace,
    g_pen,
    trust_radius_start,
    line_eval_count,
    lambda_changed,
    outer_stop_crit,
    grad_tol,
    stop_on_convergence = TRUE,
    assess_step_fn = .gl_assess_cg_step_acceptance,
    apply_step_fn = .gl_apply_cg_step_acceptance,
    diagnostics_fn = .gl_update_cg_iteration_diagnostics,
    state_builder_fn = .gl_build_cg_line_search_diagnostics_state) {
  step_acceptance <- assess_step_fn(
    best = best,
    best_raw_loglik = best_raw_loglik,
    outer_start_loglik = outer_start_log_lik,
    raw_loglik_drop_tol = raw_loglik_drop_tol,
    lambda_update_count = lambda_update_count
  )

  step_state <- apply_step_fn(
    step_acceptance = step_acceptance,
    best = best,
    eval_start = eval_start,
    beta_vec = beta_vec,
    par_cov_template = par_cov_template,
    par_s_template = par_s_template,
    stall_count = stall_count,
    trust_radius = trust_radius,
    step_tol = step_tol,
    max_delta = max_delta,
    max_stall = max_stall,
    verbose = verbose
  )

  outer_end_log_lik <- as.numeric(step_state$calc_lik_out_end$log_lik["joint"])

  iteration_diagnostics <- diagnostics_fn(
    step_trace = step_trace,
    g_pen = g_pen,
    best = step_state$best,
    best_raw_loglik = best_raw_loglik,
    best_iteration = best_iteration,
    outer_iteration = outer_iteration,
    outer_start_log_lik = outer_start_log_lik,
    outer_end_log_lik = outer_end_log_lik,
    obj_start = obj_start,
    accepted_improvement = step_state$accepted_improvement,
    trust_radius_start = trust_radius_start,
    trust_radius_end = step_state$trust_radius,
    line_eval_count = line_eval_count,
    lambda_update_count = lambda_update_count,
    lambda_changed = lambda_changed,
    stall_count = step_state$stall_count,
    prevented_deterioration = step_state$prevented_deterioration,
    prevented_raw_loglik_drop = step_state$prevented_raw_loglik_drop,
    max_stall = max_stall,
    raw_loglik_drop_tol = raw_loglik_drop_tol,
    outer_stop_crit = outer_stop_crit,
    grad_tol = grad_tol,
    step_tol = step_tol,
    verbose = verbose,
    stop_on_convergence = stop_on_convergence
  )

  state_builder_fn(
    step_acceptance = step_acceptance,
    step_state = step_state,
    outer_end_log_lik = outer_end_log_lik,
    iteration_diagnostics = iteration_diagnostics
  )
}
