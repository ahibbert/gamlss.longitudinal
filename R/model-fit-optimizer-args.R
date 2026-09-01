#' Build arguments for the selected optimizer from a prepared workflow
#'
#' @noRd
.gl_prepared_fit_optimizer_args <- function(
    fit_workflow,
    margin_dist,
    copula_dist,
    include_dlcopdpar,
    check_dlcopdpar_gradient,
    verbose,
    lambda_penalty_K,
    plot_results,
    true_val,
    check_elapsed_budget) {
  fit_controls <- fit_workflow$controls
  fit_data <- fit_workflow$fit_data
  matrix_bundle <- fit_workflow$matrix_bundle
  step_controls <- fit_workflow$step_controls
  optimizer_control <- fit_workflow$optimizer_control_effective

  list(
    method = fit_controls$method,
    optimizer_control = optimizer_control,
    optimizer_context = fit_workflow$optimizer_context,
    mm = matrix_bundle$mm,
    margin_dist = margin_dist,
    copula_link = matrix_bundle$copula_link,
    copula_dist = copula_dist,
    dataset = fit_data$dataset,
    include_dlcopdpar = include_dlcopdpar,
    check_dlcopdpar_gradient = check_dlcopdpar_gradient,
    discrete_score_method = fit_controls$discrete_score_method,
    cg_gradient_method = fit_controls$cg_gradient_method,
    cg_hessian_method = fit_controls$cg_hessian_method,
    verbose = verbose,
    check_elapsed_budget = check_elapsed_budget,
    lambda_penalty_K = lambda_penalty_K,
    plot_results = plot_results,
    true_val = true_val,
    step_adjustment = step_controls$step_adjustment,
    max_steps = step_controls$max_steps,
    start_step_size = step_controls$start_step_size
  )
}
