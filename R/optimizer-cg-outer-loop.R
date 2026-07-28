#' Run the CG outer optimizer loop
#'
#' @noRd
.gl_run_cg_outer_loop <- function(
    dataset,
    margin_dist,
    copula_dist,
    mm_cg,
    beta_all,
    penalty_mat,
    log_lik_history,
    par_history,
    cg_best_raw_loglik,
    cg_best_iteration,
    outer_only_run_counter,
    max_outer_iter,
    check_elapsed_budget,
    include_dlcopdpar,
    cg_gradient_method,
    cg_eval,
    cg_objective,
    cg_gradient,
    lambda_s,
    cg_lambda_trace,
    cg_lambda_update_count,
    cg_update_lambda,
    cg_max_lambda_updates,
    cg_lambda_update_every,
    cg_trust_radius,
    cg_max_delta,
    cg_step_tol_eff,
    build_cg_penalty,
    cg_smooth_edf_list,
    lambda_penalty_K,
    cg_armijo_c1,
    cg_line_search,
    cg_max_line_search_evals,
    use_backtracking,
    backtracking_max_halves,
    verbose,
    cg_observed_hessian,
    cg_finite_hessian_block,
    cg_zeta_hessian,
    par_cov,
    par_s,
    cg_stall_count,
    cg_max_stall,
    cg_raw_loglik_drop_tol,
    cg_step_trace,
    outer_stop_crit,
    cg_grad_tol_eff,
    df_s,
    cg_has_smooths,
    cg_converged,
    outer_iteration_step_fn,
    iteration_start_fn,
    curvature_line_search_fn,
    line_search_diagnostics_fn,
    stop_request_state_fn) {
  while (!cg_converged && outer_only_run_counter < max_outer_iter) {
    check_elapsed_budget("CG outer iteration")

    if (verbose > 0) cat(paste("\nOUTER ITERATION:", outer_only_run_counter))

    cg_outer_iteration <- outer_iteration_step_fn(
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
      iteration_start_fn = iteration_start_fn,
      curvature_line_search_fn = curvature_line_search_fn,
      line_search_diagnostics_fn = line_search_diagnostics_fn,
      stop_request_state_fn = stop_request_state_fn
    )
    cg_outer_state <- .gl_apply_cg_outer_iteration_state(cg_outer_iteration)
    beta_all <- cg_outer_state$beta_all
    par_cov <- cg_outer_state$par_cov
    par_s <- cg_outer_state$par_s
    lambda_s <- cg_outer_state$lambda_s
    df_s <- cg_outer_state$df_s
    penalty_mat <- cg_outer_state$penalty_mat
    cg_trust_radius <- cg_outer_state$cg_trust_radius
    cg_stall_count <- cg_outer_state$cg_stall_count
    cg_lambda_update_count <- cg_outer_state$cg_lambda_update_count
    cg_lambda_trace <- cg_outer_state$cg_lambda_trace
    cg_step_trace <- cg_outer_state$cg_step_trace
    cg_converged <- cg_outer_state$cg_converged
    calc_lik_out_end <- cg_outer_state$calc_lik_out_end
    outer_only_run_counter <- cg_outer_state$outer_only_run_counter
    outer_start_log_lik <- cg_outer_state$outer_start_log_lik
    outer_end_log_lik <- cg_outer_state$outer_end_log_lik
    outer_log_lik_change <- cg_outer_state$outer_log_lik_change
    log_lik_history <- cg_outer_state$log_lik_history
    par_history <- cg_outer_state$par_history
    cg_stop_reason <- cg_outer_state$cg_stop_reason
    cg_last_grad_inf <- cg_outer_state$cg_last_grad_inf
    cg_last_step_l2 <- cg_outer_state$cg_last_step_l2
    cg_best_raw_loglik <- cg_outer_state$cg_best_raw_loglik
    cg_best_iteration <- cg_outer_state$cg_best_iteration
    cg_raw_loglik_drop_from_best <- cg_outer_state$cg_raw_loglik_drop_from_best
  }

  list(
    beta_all = beta_all,
    par_cov = par_cov,
    par_s = par_s,
    lambda_s = lambda_s,
    df_s = df_s,
    penalty_mat = penalty_mat,
    cg_trust_radius = cg_trust_radius,
    cg_stall_count = cg_stall_count,
    cg_lambda_update_count = cg_lambda_update_count,
    cg_lambda_trace = cg_lambda_trace,
    cg_step_trace = cg_step_trace,
    cg_converged = cg_converged,
    calc_lik_out_end = calc_lik_out_end,
    outer_only_run_counter = outer_only_run_counter,
    outer_start_log_lik = outer_start_log_lik,
    outer_end_log_lik = outer_end_log_lik,
    outer_log_lik_change = outer_log_lik_change,
    log_lik_history = log_lik_history,
    par_history = par_history,
    cg_stop_reason = cg_stop_reason,
    cg_last_grad_inf = cg_last_grad_inf,
    cg_last_step_l2 = cg_last_step_l2,
    cg_best_raw_loglik = cg_best_raw_loglik,
    cg_best_iteration = cg_best_iteration,
    cg_raw_loglik_drop_from_best = cg_raw_loglik_drop_from_best
  )
}
