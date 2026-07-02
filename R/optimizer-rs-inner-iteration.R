#' Run one RS inner likelihood-score-backfitting iteration:
#' calculates likelihood with .gl_evaluate_rs_likelihood_state (see R/optimizer-rs-likelihood.R), 
#' calculates score with .gl_evaluate_rs_parameter_score_state (see R/optimizer-rs-score.R), 
#' and performs backfitting step with .gl_run_rs_backfitting_step (see R/optimizer-rs-backfitting.R).
#'
#' @noRd
.gl_run_rs_inner_iteration_step <- function(
    par_name,
    rs_calc_eta,
    par_cov,
    par_s,
    mm,
    margin_dist,
    copula_dist,
    dataset,
    pair_cache,
    margin_eval_cache,
    current_calc_lik_out = NULL,
    first_outer_run,
    outer_start_log_lik,
    timer,
    timer_start,
    log_lik_history,
    par_history,
    copula_link,
    discrete_score_method,
    include_dlcopdpar,
    check_dlcopdpar_gradient,
    outer_only_run_counter,
    verbose,
    rs_design_cache,
    rs_smooth_trust_radius,
    lambda_penalty_K,
    lambda_s,
    df_s,
    step_size,
    rs_update_lambda,
    inner_run_counter,
    outer_run_counter,
    rs_block_trace,
    run_counter,
    use_backtracking,
    backtracking_max_halves,
    plot_results,
    true_val,
    likelihood_state_fn,
    score_state_fn,
    backfitting_step_fn) {

  likelihood_state_args <- list(
    rs_calc_eta = rs_calc_eta,
    par_cov = par_cov,
    par_s = par_s,
    mm = mm,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    dataset = dataset,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache,
    first_outer_run = first_outer_run,
    outer_start_log_lik = outer_start_log_lik,
    timer = timer,
    timer_start = timer_start,
    log_lik_history = log_lik_history,
    par_history = par_history
  )
  likelihood_state_formals <- names(formals(likelihood_state_fn))
  if ("..." %in% likelihood_state_formals || "par_name" %in% likelihood_state_formals) {
    likelihood_state_args$par_name <- par_name
  }
  if ("..." %in% likelihood_state_formals || "current_calc_lik_out" %in% likelihood_state_formals) {
    likelihood_state_args$current_calc_lik_out <- current_calc_lik_out
  }
  rs_likelihood_state <- do.call(likelihood_state_fn, likelihood_state_args)
  eta <- rs_likelihood_state$eta
  eta_dr <- rs_likelihood_state$eta_dr
  eta_inv <- rs_likelihood_state$eta_inv
  calc_lik_out <- rs_likelihood_state$calc_lik_out
  log_lik <- rs_likelihood_state$log_lik
  margin_deriv <- rs_likelihood_state$margin_deriv
  copula_d <- rs_likelihood_state$copula_d
  Fx_1_2 <- rs_likelihood_state$Fx_1_2
  first_outer_run <- rs_likelihood_state$first_outer_run
  outer_start_log_lik <- rs_likelihood_state$outer_start_log_lik
  timer <- rs_likelihood_state$timer
  log_lik_history <- rs_likelihood_state$log_lik_history
  par_history <- rs_likelihood_state$par_history

  rs_score_state <- score_state_fn(
    par_name = par_name,
    eta = eta,
    eta_inv = eta_inv,
    Fx_1_2 = Fx_1_2,
    copula_dist = copula_dist,
    calc_lik_out = calc_lik_out,
    mm = mm,
    margin_dist = margin_dist,
    dataset = dataset,
    pair_cache = pair_cache,
    discrete_score_method = discrete_score_method,
    include_dlcopdpar = include_dlcopdpar,
    margin_deriv = margin_deriv,
    copula_d = copula_d,
    log_lik = log_lik,
    check_dlcopdpar_gradient = check_dlcopdpar_gradient,
    outer_only_run_counter = outer_only_run_counter,
    verbose = verbose,
    timer = timer,
    timer_start = timer_start
  )
  d1 <- rs_score_state$d1
  timer <- rs_score_state$timer

  rs_backfitting_step <- backfitting_step_fn(
    eta = eta,
    eta_dr = eta_dr,
    d1 = d1,
    par_name = par_name,
    rs_design_cache = rs_design_cache,
    par_cov = par_cov,
    par_s = par_s,
    timer = timer,
    timer_start = timer_start,
    rs_smooth_trust_radius = rs_smooth_trust_radius,
    rs_calc_eta = rs_calc_eta,
    calc_lik_out = calc_lik_out,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache,
    lambda_penalty_K = lambda_penalty_K,
    lambda_s = lambda_s,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    dataset = dataset,
    mm = mm,
    copula_link = copula_link,
    df_s = df_s,
    step_size = step_size,
    rs_update_lambda = rs_update_lambda,
    inner_run_counter = inner_run_counter,
    outer_run_counter = outer_run_counter,
    outer_only_run_counter = outer_only_run_counter,
    verbose = verbose,
    rs_block_trace = rs_block_trace,
    run_counter = run_counter,
    use_backtracking = use_backtracking,
    backtracking_max_halves = backtracking_max_halves,
    plot_results = plot_results,
    log_lik_history = log_lik_history,
    par_history = par_history,
    true_val = true_val
  )

  list(
    score = rs_backfitting_step$score,
    lambda_s = rs_backfitting_step$lambda_s,
    par_cov = rs_backfitting_step$par_cov,
    par_s = rs_backfitting_step$par_s,
    calc_lik_out_end = rs_backfitting_step$calc_lik_out_end,
    df_s = rs_backfitting_step$df_s,
    rs_block_trace = rs_backfitting_step$rs_block_trace,
    change_log_lik = rs_backfitting_step$change_log_lik,
    run_counter = rs_backfitting_step$run_counter,
    outer_run_counter = rs_backfitting_step$outer_run_counter,
    inner_run_counter = rs_backfitting_step$inner_run_counter,
    timer = rs_backfitting_step$timer,
    first_outer_run = first_outer_run,
    outer_start_log_lik = outer_start_log_lik,
    log_lik_history = log_lik_history,
    par_history = par_history
  )
}
