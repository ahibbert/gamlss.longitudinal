#' Create the final fitted gamlss.longitudinal object from all the components
#'
#' @noRd
.gl_finalize_fit_object <- function(
    par_cov,
    log_lik_history,
    par_history,
    calc_lik_out_end,
    mm,
    margin_dist,
    copula_dist,
    include_dlcopdpar,
    dataset,
    dataset_original,
    response_var,
    time_var,
    subject_var,
    formulas,
    formulas_int,
    var_map,
    par_s,
    lambda_s,
    df_s,
    weights_final,
    method,
    warm_start_info,
    rs_block_trace,
    cg_lambda_trace,
    cg_step_trace,
    convergence_info,
    fit_start_time,
    compute_vcov,
    vcov_numderiv,
    vcov_method,
    verbose,
    optimizer_control_requested = NULL,
    optimizer_control_effective = NULL) {
  if (!isTRUE(convergence_info$converged)) .gl_warn_nonconverged_fit(convergence_info)

  total_fit_time <- as.numeric(difftime(Sys.time(), fit_start_time, units = "secs"))
  aics <- .gl_fit_information_criteria(
    par_cov = par_cov,
    par_s = par_s,
    df_s = df_s,
    calc_lik_out_end = calc_lik_out_end,
    dataset = dataset
  )

  if (verbose > 0) {
    .gl_print_fit_summary(
      par_cov = par_cov,
      aics = aics,
      dataset = dataset,
      margin_dist = margin_dist,
      copula_dist = copula_dist,
      total_fit_time = total_fit_time
    )
  }

  return_list <- .gl_build_fit_object_core(
    par_cov = par_cov,
    log_lik_history = log_lik_history,
    par_history = par_history,
    calc_lik_out_end = calc_lik_out_end,
    mm = mm,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    include_dlcopdpar = include_dlcopdpar,
    dataset = dataset,
    dataset_original = dataset_original,
    response_var = response_var,
    time_var = time_var,
    subject_var = subject_var,
    formulas = formulas,
    formulas_int = formulas_int,
    var_map = var_map,
    par_s = par_s,
    lambda_s = lambda_s,
    df_s = df_s,
    weights_final = weights_final,
    method = method,
    warm_start_info = warm_start_info,
    convergence_info = convergence_info,
    optimizer_control_requested = optimizer_control_requested,
    optimizer_control_effective = optimizer_control_effective
  )

  return_list <- .gl_attach_fit_optimizer_traces(
    return_list = return_list,
    method = method,
    rs_block_trace = rs_block_trace,
    cg_lambda_trace = cg_lambda_trace,
    cg_step_trace = cg_step_trace
  )

  return_list <- .gl_attach_fit_vcov(
    return_list = return_list,
    compute_vcov = isTRUE(compute_vcov) && isTRUE(convergence_info$converged),
    vcov_numderiv = vcov_numderiv,
    vcov_method = vcov_method,
    verbose = verbose
  )

  class(return_list) <- "gamlss.longitudinal"
  return_list
}
