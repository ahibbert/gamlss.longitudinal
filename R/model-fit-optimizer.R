#' Run the selected optimizer for an already prepared fitting workflow
#'
#' Bridges `.gl_prepare_fit_workflow()` to `.gl_run_fit_optimizer()` by structuring
#' the inputs like controls, model matrices etc. into a single list of arguments for the optimizer. 
#' This is a bit of an extra step, but it keeps the optimizer function focused on the optimization 
#' logic and not on the setup which is a bit easier for review and testing.
#'
#' @noRd
.gl_run_prepared_fit_optimizer <- function(
    fit_workflow,
    margin_dist,
    copula_dist,
    include_dlcopdpar,
    check_dlcopdpar_gradient,
    verbose,
    lambda_penalty_K,
    plot_results,
    true_val,
    check_elapsed_budget,
    optimizer_fn = .gl_run_fit_optimizer) {
  optimizer_args <- .gl_prepared_fit_optimizer_args(
    fit_workflow = fit_workflow,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    include_dlcopdpar = include_dlcopdpar,
    check_dlcopdpar_gradient = check_dlcopdpar_gradient,
    verbose = verbose,
    lambda_penalty_K = lambda_penalty_K,
    plot_results = plot_results,
    true_val = true_val,
    check_elapsed_budget = check_elapsed_budget
  )

  do.call(optimizer_fn, optimizer_args)
}
