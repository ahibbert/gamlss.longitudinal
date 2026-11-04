#' Sets off the workflow for the selected optimizer 
#' 
#' @noRd
.gl_run_fit_optimizer <- function(
    method,
    optimizer_control,
    optimizer_context,
    mm,
    margin_dist,
    copula_link,
    copula_dist,
    dataset,
    include_dlcopdpar,
    check_dlcopdpar_gradient,
    discrete_score_method,
    cg_gradient_method,
    cg_hessian_method,
    verbose,
    check_elapsed_budget,
    lambda_penalty_K,
    plot_results,
    true_val,
    step_adjustment,
    max_steps,
    start_step_size,
    cg_runner_fn = .gl_run_cg_optimizer,
    rs_runner_fn = .gl_run_rs_optimizer) {
  shared <- optimizer_control$shared
  if (identical(method, "CG")) {
    cg <- optimizer_control$cg
    return(.gl_run_fit_cg_optimizer(
      optimizer_context = optimizer_context,
      mm = mm,
      margin_dist = margin_dist,
      copula_link = copula_link,
      copula_dist = copula_dist,
      dataset = dataset,
      include_dlcopdpar = include_dlcopdpar,
      cg_gradient_method = cg_gradient_method,
      cg_hessian_method = cg_hessian_method,
      verbose = verbose,
      cg_max_delta = cg$max_delta,
      max_outer_iter = shared$max_outer_iter,
      check_elapsed_budget = check_elapsed_budget,
      cg_update_lambda = cg$update_lambda,
      cg_max_lambda_updates = cg$max_lambda_updates,
      cg_lambda_update_every = cg$lambda_update_every,
      lambda_penalty_K = lambda_penalty_K,
      cg_armijo_c1 = cg$armijo_c1,
      cg_line_search = cg$line_search,
      cg_max_line_search_evals = cg$max_line_search_evals,
      use_backtracking = cg$use_backtracking,
      backtracking_max_halves = cg$backtracking_max_halves,
      cg_zeta_hessian = cg$zeta_hessian,
      cg_max_stall = cg$max_stall,
      cg_raw_loglik_drop_tol = cg$raw_loglik_drop_tol,
      stop_on_convergence = shared$stop_on_convergence,
      cg_runner_fn = cg_runner_fn
    ))
  }

  rs <- optimizer_control$rs
  .gl_run_fit_rs_optimizer(
    optimizer_context = optimizer_context,
    mm = mm,
    margin_dist = margin_dist,
    copula_link = copula_link,
    copula_dist = copula_dist,
    dataset = dataset,
    include_dlcopdpar = include_dlcopdpar,
    check_dlcopdpar_gradient = check_dlcopdpar_gradient,
    discrete_score_method = discrete_score_method,
    verbose = verbose,
    max_outer_iter = shared$max_outer_iter,
    check_elapsed_budget = check_elapsed_budget,
    lambda_penalty_K = lambda_penalty_K,
    use_backtracking = rs$use_backtracking,
    backtracking_max_halves = rs$backtracking_max_halves,
    rs_smooth_trust_radius = rs$smooth_trust_radius,
    rs_update_lambda = rs$update_lambda,
    plot_results = plot_results,
    true_val = true_val,
    max_inner_iter = rs$max_inner_iter,
    max_negative_outer_streak = rs$max_negative_outer_streak,
    stop_on_convergence = shared$stop_on_convergence,
    step_adjustment = step_adjustment,
    max_steps = max_steps,
    start_step_size = start_step_size,
    rs_runner_fn = rs_runner_fn
  )
}
