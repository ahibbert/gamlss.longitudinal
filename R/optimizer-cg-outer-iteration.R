#' Run one CG outer optimizer iteration
#'
#' @noRd
.gl_run_cg_outer_iteration_step <- function(
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
    stop_on_convergence,
    iteration_start_fn,
    curvature_line_search_fn,
    line_search_diagnostics_fn,
    stop_request_state_fn,
    maybe_stop_request_state_fn = .gl_maybe_apply_cg_stop_request_state) {
  cg_trust_radius_start <- cg_trust_radius

  cg_iteration_start <- iteration_start_fn(
    beta_vec = beta_all,
    mm_cg = mm_cg,
    penalty_current = penalty_mat,
    log_lik_history = log_lik_history,
    par_history = par_history,
    best_raw_loglik = cg_best_raw_loglik,
    best_iteration = cg_best_iteration,
    current_iteration = outer_only_run_counter,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    include_dlcopdpar = include_dlcopdpar,
    dataset = dataset,
    cg_gradient_method = cg_gradient_method,
    eval_fn = cg_eval,
    objective_fn = cg_objective,
    finite_gradient_fn = cg_gradient
  )
  eval_start <- cg_iteration_start$eval_start
  log_lik_history <- cg_iteration_start$log_lik_history
  par_history <- cg_iteration_start$par_history
  outer_start_log_lik <- cg_iteration_start$outer_start_log_lik
  cg_best_raw_loglik <- cg_iteration_start$best_raw_loglik
  cg_best_iteration <- cg_iteration_start$best_iteration
  obj_start <- cg_iteration_start$obj_start
  grad <- cg_iteration_start$grad

  cg_curvature_line_search <- curvature_line_search_fn(
    dataset = dataset,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    mm_cg = mm_cg,
    beta_vec = beta_all,
    grad_vec = grad,
    lambda_current = lambda_s,
    lambda_trace = cg_lambda_trace,
    penalty_current = penalty_mat,
    lambda_update_count = cg_lambda_update_count,
    update_lambda = cg_update_lambda,
    max_lambda_updates = cg_max_lambda_updates,
    lambda_update_every = cg_lambda_update_every,
    outer_iteration = outer_only_run_counter,
    trust_radius = cg_trust_radius,
    max_delta = cg_max_delta,
    step_tol = cg_step_tol_eff,
    build_penalty_fn = build_cg_penalty,
    eval_fn = cg_eval,
    edf_fn = cg_smooth_edf_list,
    objective_fn = cg_objective,
    lambda_penalty_K = lambda_penalty_K,
    obj_start = obj_start,
    armijo_c1 = cg_armijo_c1,
    line_search = cg_line_search,
    max_line_search_evals = cg_max_line_search_evals,
    use_backtracking = use_backtracking,
    backtracking_max_halves = backtracking_max_halves,
    verbose = verbose,
    observed_hessian_fn = cg_observed_hessian,
    finite_hessian_fn = cg_finite_hessian_block,
    cg_zeta_hessian = cg_zeta_hessian
  )
  H_obs <- cg_curvature_line_search$H_obs
  lambda_s <- cg_curvature_line_search$lambda_s
  cg_lambda_trace <- cg_curvature_line_search$cg_lambda_trace
  penalty_mat <- cg_curvature_line_search$penalty_mat
  cg_trust_radius <- cg_curvature_line_search$cg_trust_radius
  cg_lambda_update_count <- cg_curvature_line_search$cg_lambda_update_count
  lambda_changed <- cg_curvature_line_search$lambda_changed
  df_s <- cg_curvature_line_search$df_s
  g_pen <- cg_curvature_line_search$g_pen
  best <- cg_curvature_line_search$best
  line_eval_count <- cg_curvature_line_search$line_eval_count

  cg_line_search_diagnostics <- line_search_diagnostics_fn(
    best = best,
    eval_start = eval_start,
    beta_vec = beta_all,
    par_cov_template = par_cov,
    par_s_template = par_s,
    stall_count = cg_stall_count,
    trust_radius = cg_trust_radius,
    step_tol = cg_step_tol_eff,
    max_delta = cg_max_delta,
    max_stall = cg_max_stall,
    verbose = verbose,
    best_raw_loglik = cg_best_raw_loglik,
    best_iteration = cg_best_iteration,
    outer_iteration = outer_only_run_counter,
    outer_start_log_lik = outer_start_log_lik,
    obj_start = obj_start,
    raw_loglik_drop_tol = cg_raw_loglik_drop_tol,
    step_trace = cg_step_trace,
    g_pen = g_pen,
    trust_radius_start = cg_trust_radius_start,
    line_eval_count = line_eval_count,
    lambda_update_count = cg_lambda_update_count,
    lambda_changed = lambda_changed,
    outer_stop_crit = outer_stop_crit,
    grad_tol = cg_grad_tol_eff,
    stop_on_convergence = stop_on_convergence
  )
  beta_all <- cg_line_search_diagnostics$beta
  par_cov <- cg_line_search_diagnostics$par_cov
  par_s <- cg_line_search_diagnostics$par_s
  calc_lik_out_end <- cg_line_search_diagnostics$calc_lik_out_end
  cg_stall_count <- cg_line_search_diagnostics$stall_count
  cg_trust_radius <- cg_line_search_diagnostics$trust_radius
  outer_end_log_lik <- cg_line_search_diagnostics$outer_end_log_lik
  outer_log_lik_change <- cg_line_search_diagnostics$outer_log_lik_change
  cg_best_raw_loglik <- cg_line_search_diagnostics$best_raw_loglik
  cg_best_iteration <- cg_line_search_diagnostics$best_iteration
  cg_raw_loglik_drop_from_best <- cg_line_search_diagnostics$raw_loglik_drop_from_best
  grad_inf <- cg_line_search_diagnostics$grad_inf
  step_l2 <- cg_line_search_diagnostics$step_l2
  cg_last_grad_inf <- grad_inf
  cg_last_step_l2 <- step_l2
  cg_step_trace <- cg_line_search_diagnostics$step_trace

  cg_tolerance_met <- cg_line_search_diagnostics$tolerance_met
  cg_deterioration_hit <- cg_line_search_diagnostics$deterioration_hit
  cg_stop_requested <- cg_line_search_diagnostics$stop_requested
  cg_stop_state <- maybe_stop_request_state_fn(
    stop_requested = cg_stop_requested,
    cg_update_lambda = cg_update_lambda,
    cg_has_smooths = cg_has_smooths,
    cg_lambda_update_count = cg_lambda_update_count,
    lambda_s = lambda_s,
    cg_lambda_trace = cg_lambda_trace,
    penalty_mat = penalty_mat,
    df_s = df_s,
    cg_stall_count = cg_stall_count,
    H_obs = H_obs,
    beta_all = beta_all,
    grad = grad,
    mm_cg = mm_cg,
    cg_trust_radius = cg_trust_radius,
    outer_only_run_counter = outer_only_run_counter,
    cg_max_delta = cg_max_delta,
    build_cg_penalty = build_cg_penalty,
    cg_eval = cg_eval,
    cg_smooth_edf_list = cg_smooth_edf_list,
    cg_objective = cg_objective,
    lambda_penalty_K = lambda_penalty_K,
    cg_tolerance_met = cg_tolerance_met,
    cg_deterioration_hit = cg_deterioration_hit,
    cg_raw_loglik_drop_from_best = cg_raw_loglik_drop_from_best,
    verbose = verbose,
    stop_request_state_fn = stop_request_state_fn
  )
  lambda_s <- cg_stop_state$lambda_s
  cg_lambda_trace <- cg_stop_state$cg_lambda_trace
  cg_lambda_update_count <- cg_stop_state$cg_lambda_update_count
  penalty_mat <- cg_stop_state$penalty_mat
  df_s <- cg_stop_state$df_s
  cg_stall_count <- cg_stop_state$cg_stall_count
  cg_stop_reason <- cg_stop_state$cg_stop_reason
  cg_converged <- cg_stop_state$cg_converged

  outer_only_run_counter <- outer_only_run_counter + 1L

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
