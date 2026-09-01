#' Finalize an already prepared and optimised model fit
#'
#' @noRd
.gl_finalize_prepared_fit <- function(
    optimizer_state,
    fit_workflow,
    margin_dist,
    copula_dist,
    include_dlcopdpar,
    original_formulas,
    time_var,
    subject_var,
    fit_start_time,
    compute_vcov,
    verbose,
    max_outer_iter = NULL,
    finalize_fn = .gl_finalize_fit_workflow) {
  fit_controls <- fit_workflow$controls

  finalize_fn(
    optimizer_state = optimizer_state,
    fit_data = fit_workflow$fit_data,
    matrix_bundle = fit_workflow$matrix_bundle,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    include_dlcopdpar = include_dlcopdpar,
    method = fit_controls$method,
    original_formulas = original_formulas,
    time_var = time_var,
    subject_var = subject_var,
    warm_start_info = fit_workflow$warm_start$warm_start_info,
    fit_start_time = fit_start_time,
    compute_vcov = compute_vcov,
    vcov_numderiv = fit_controls$vcov_numderiv,
    vcov_method = fit_controls$vcov_method,
    verbose = verbose,
    max_outer_iter = fit_workflow$optimizer_control_effective$shared$max_outer_iter,
    cg_raw_loglik_drop_tol = fit_controls$cg_raw_loglik_drop_tol,
    cg_gradient_method = fit_controls$cg_gradient_method,
    cg_zeta_hessian = fit_controls$cg_zeta_hessian,
    cg_hessian_method = fit_controls$cg_hessian_method,
    optimizer_control_requested = fit_workflow$optimizer_control_requested,
    optimizer_control_effective = fit_workflow$optimizer_control_effective
  )
}
