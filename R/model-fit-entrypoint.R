#' Main entry point for model fitting, called from gamlss_longitudinal().
#' 
#' This function orchestrates the entire model fitting process, including:
#' 1. Preprocessing and validation of input data and formulas by calling workflow_fn = .gl_prepare_fit_workflow (model-fit-workflow.R),
#' 2. Execution of the optimization loop by calling optimizer_fn = .gl_run_prepared_fit_optimizer (model-fit-optimizer.R),
#' 3. Finalization of the fit by calling finalize_fn = .gl_finalize_prepared_fit (model-fit-finalize.R)
#' 
#' @noRd
.gl_run_gamlss_longitudinal_entrypoint <- function(
    dataset,
    margin_dist,
    copula_dist,
    time_var,
    subject_var,
    missingness,
    mu.formula,
    sigma.formula,
    nu.formula,
    tau.formula,
    theta.formula,
    zeta.formula,
    include_dlcopdpar,
    check_dlcopdpar_gradient,
    start_from,
    verbose,
    plot_results,
    true_val,
    method,
    optimizer_control,
    compute_vcov,
    vcov_method,
    vcov_numderiv,
    use_Rcpp,
    lambda_start,
    lambda_penalty_K,
    time_fn = Sys.time,
    capability_preflight_fn = .gl_preflight_fit_capabilities,
    margin_normalizer_fn = .normalise_margin_dist_links,
    budget_checker_fn = .gl_build_elapsed_budget_checker,
    workflow_fn = .gl_prepare_fit_workflow,
    optimizer_fn = .gl_run_prepared_fit_optimizer,
    finalize_fn = .gl_finalize_prepared_fit) {
  fit_start_time <- time_fn()

  capability_preflight_fn(
    dataset = dataset,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    mu.formula = mu.formula
  )

  margin_dist <- margin_normalizer_fn(margin_dist)

  check_elapsed_budget <- budget_checker_fn(
    fit_start_time = fit_start_time,
    max_elapsed_sec = optimizer_control$shared$max_elapsed_sec
  )

  fit_workflow <- workflow_fn(
    dataset = dataset,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    time_var = time_var,
    subject_var = subject_var,
    missingness = missingness,
    mu.formula = mu.formula,
    sigma.formula = sigma.formula,
    nu.formula = nu.formula,
    tau.formula = tau.formula,
    theta.formula = theta.formula,
    zeta.formula = zeta.formula,
    include_dlcopdpar = include_dlcopdpar,
    start_from = start_from,
    verbose = verbose,
    true_val = true_val,
    method = method,
    optimizer_control = optimizer_control,
    vcov_method = vcov_method,
    vcov_numderiv = vcov_numderiv,
    use_Rcpp = use_Rcpp,
    lambda_start = lambda_start,
    lambda_penalty_K = lambda_penalty_K
  )

  .gl_warn_segmented_missingness(fit_workflow$fit_data$missingness_contract)

  optimizer_state <- optimizer_fn(
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

  finalize_fn(
    optimizer_state = optimizer_state,
    fit_workflow = fit_workflow,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    include_dlcopdpar = include_dlcopdpar,
    original_formulas = .gl_original_formula_bundle(
      mu.formula = mu.formula,
      sigma.formula = sigma.formula,
      nu.formula = nu.formula,
      tau.formula = tau.formula,
      theta.formula = theta.formula,
      zeta.formula = zeta.formula
    ),
    time_var = time_var,
    subject_var = subject_var,
    fit_start_time = fit_start_time,
    compute_vcov = compute_vcov,
    verbose = verbose
  )
}
