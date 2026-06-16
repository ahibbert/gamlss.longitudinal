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
    cg_armijo_c1,
    cg_update_lambda,
    lambda_penalty_K,
    rs_update_lambda,
    plot_results,
    true_val,
    max_outer_iter,
    max_inner_iter,
    max_negative_outer_streak,
    use_backtracking,
    check_elapsed_budget) {
  fit_controls <- fit_workflow$controls
  fit_data <- fit_workflow$fit_data
  matrix_bundle <- fit_workflow$matrix_bundle
  step_controls <- fit_workflow$step_controls

  list(
    method = fit_controls$method,
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
    cg_max_delta = fit_controls$cg_max_delta,
    max_outer_iter = max_outer_iter,
    check_elapsed_budget = check_elapsed_budget,
    cg_update_lambda = cg_update_lambda,
    cg_max_lambda_updates = fit_controls$cg_max_lambda_updates,
    cg_lambda_update_every = fit_controls$cg_lambda_update_every,
    lambda_penalty_K = lambda_penalty_K,
    cg_armijo_c1 = cg_armijo_c1,
    cg_line_search = fit_controls$cg_line_search,
    cg_max_line_search_evals = fit_controls$cg_max_line_search_evals,
    use_backtracking = use_backtracking,
    backtracking_max_halves = fit_controls$backtracking_max_halves,
    cg_zeta_hessian = fit_controls$cg_zeta_hessian,
    cg_max_stall = fit_controls$cg_max_stall,
    cg_raw_loglik_drop_tol = fit_controls$cg_raw_loglik_drop_tol,
    rs_smooth_trust_radius = fit_controls$rs_smooth_trust_radius,
    rs_update_lambda = rs_update_lambda,
    plot_results = plot_results,
    true_val = true_val,
    max_inner_iter = max_inner_iter,
    max_negative_outer_streak = max_negative_outer_streak,
    step_adjustment = step_controls$step_adjustment,
    max_steps = step_controls$max_steps,
    start_step_size = step_controls$start_step_size
  )
}
