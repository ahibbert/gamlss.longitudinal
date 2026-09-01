#' Finalize a fit object from a completed model fit
#' 
#' First collects convergence metadata then calls `.gl_finalize_fit_object()` to put together the final fitted object. 
#' 
#' This is the post-optimizer counterpart to `.gl_prepare_fit_workflow()` and keeps `gamlss_longitudinal()` readable as a high-level sequence.
#'
#' @noRd
.gl_finalize_fit_workflow <- function(
    optimizer_state,
    fit_data,
    matrix_bundle,
    margin_dist,
    copula_dist,
    include_dlcopdpar,
    method,
    original_formulas,
    time_var,
    subject_var,
    warm_start_info,
    fit_start_time,
    compute_vcov,
    vcov_numderiv,
    vcov_method,
    verbose,
    max_outer_iter,
    cg_raw_loglik_drop_tol,
    cg_gradient_method,
    cg_zeta_hessian,
    cg_hessian_method,
    optimizer_control_requested = NULL,
    optimizer_control_effective = NULL,
    convergence_fn = .gl_build_convergence_info,
    finalize_fn = .gl_finalize_fit_object) {
  elapsed_sec <- as.numeric(difftime(Sys.time(), fit_start_time, units = "secs"))
  objective <- as.numeric(optimizer_state$calc_lik_out_end$log_lik["joint"])
  convergence_info <- convergence_fn(
    method = method,
    outer_log_lik_change = optimizer_state$outer_log_lik_change,
    outer_stop_crit = optimizer_state$outer_stop_crit,
    outer_only_run_counter = optimizer_state$outer_only_run_counter,
    max_outer_iter = max_outer_iter,
    cg_stop_reason = optimizer_state$cg_stop_reason,
    cg_last_grad_inf = optimizer_state$cg_last_grad_inf,
    cg_last_step_l2 = optimizer_state$cg_last_step_l2,
    cg_best_raw_loglik = optimizer_state$cg_best_raw_loglik,
    cg_best_iteration = optimizer_state$cg_best_iteration,
    cg_raw_loglik_drop_from_best = optimizer_state$cg_raw_loglik_drop_from_best,
    cg_raw_loglik_drop_tol = cg_raw_loglik_drop_tol,
    cg_gradient_method = cg_gradient_method,
    cg_zeta_hessian = cg_zeta_hessian,
    cg_hessian_method = cg_hessian_method,
    objective = objective,
    elapsed_sec = elapsed_sec,
    optimizer_control_requested = optimizer_control_requested,
    optimizer_control_effective = optimizer_control_effective,
    cg_grad_tol = optimizer_state$cg_grad_tol_eff,
    cg_step_tol = optimizer_state$cg_step_tol_eff
  )

  finalize_fn(
    par_cov = optimizer_state$par_cov,
    log_lik_history = optimizer_state$log_lik_history,
    par_history = optimizer_state$par_history,
    calc_lik_out_end = optimizer_state$calc_lik_out_end,
    mm = matrix_bundle$mm,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    include_dlcopdpar = include_dlcopdpar,
    dataset = fit_data$dataset,
    dataset_original = fit_data$dataset_original,
    response_var = fit_data$response_var,
    time_var = time_var,
    subject_var = subject_var,
    missingness_contract = fit_data$missingness_contract,
    formulas = original_formulas,
    formulas_int = fit_data$formulas_int,
    var_map = fit_data$var_map,
    par_s = optimizer_state$par_s,
    lambda_s = optimizer_state$lambda_s,
    df_s = optimizer_state$df_s,
    weights_final = optimizer_state$weights_final,
    method = method,
    warm_start_info = warm_start_info,
    rs_block_trace = optimizer_state$rs_block_trace,
    cg_lambda_trace = optimizer_state$cg_lambda_trace,
    cg_step_trace = optimizer_state$cg_step_trace,
    convergence_info = convergence_info,
    optimizer_control_requested = optimizer_control_requested,
    optimizer_control_effective = optimizer_control_effective,
    fit_start_time = fit_start_time,
    compute_vcov = compute_vcov,
    vcov_numderiv = vcov_numderiv,
    vcov_method = vcov_method,
    verbose = verbose
  )
}
