#' Initialises starting parameters and components needed for model fitting
#' 
#' Initializes the coefficients, design matrixes, and stopping criteria for optimiser.
#'
#' @noRd
.gl_initialize_fit_optimizer_context <- function(
    start_from,
    warm_start_par_s,
    mm,
    margin_dist,
    copula_dist,
    dataset,
    lambda_start,
    start_step_size,
    copula_link,
    inner_stop_crit,
    outer_stop_crit,
    cg_grad_tol,
    cg_step_tol,
    method,
    verbose,
    initial_state_fn = .gl_build_initial_parameter_state,
    rs_design_cache_fn = .gl_build_rs_design_cache,
    optimizer_state_fn = .gl_initialize_optimizer_state,
    rs_eta_calculator_fn = .gl_build_rs_eta_calculator,
    stop_criteria_fn = .gl_resolve_stop_criteria) {
  initial_state <- initial_state_fn(
    start_from = start_from,
    warm_start_par_s = warm_start_par_s,
    mm = mm,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    dataset = dataset,
    lambda_start = lambda_start
  )

  par_cov <- initial_state$par_cov
  par_s <- initial_state$par_s
  df_s <- initial_state$df_s
  lambda_s <- initial_state$lambda_s

  rs_design_cache <- rs_design_cache_fn(mm = mm, par_s = par_s)

  optimizer_state <- optimizer_state_fn(
    par_cov = par_cov,
    start_step_size = start_step_size,
    dataset = dataset,
    margin_dist = margin_dist
  )

  rs_calc_eta <- rs_eta_calculator_fn(
    rs_design_cache = rs_design_cache,
    mm = mm,
    margin_dist = margin_dist,
    copula_link = copula_link
  )

  stop_criteria <- stop_criteria_fn(
    inner_stop_crit = inner_stop_crit,
    outer_stop_crit = outer_stop_crit,
    cg_grad_tol = cg_grad_tol,
    cg_step_tol = cg_step_tol,
    method = method,
    par_cov = par_cov,
    par_s = par_s,
    mm = mm,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    copula_link = copula_link,
    dataset = dataset,
    pair_cache = optimizer_state$pair_cache,
    margin_eval_cache = optimizer_state$margin_eval_cache,
    verbose = verbose
  )

  .gl_fit_optimizer_context_bundle(
    initial_state = initial_state,
    rs_design_cache = rs_design_cache,
    rs_calc_eta = rs_calc_eta,
    optimizer_state = optimizer_state,
    stop_criteria = stop_criteria
  )
}
