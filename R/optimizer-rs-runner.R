#' Run the RS optimizer loop
#'
#' Owns the RS outer loop (and stopping) and per-parameter inner GLIM loop. The numerical
#' work remains delegated to the likelihood, score, backfitting, and stopping
#' helpers so this file is the reviewer entry point for RS optimizer workflow.
#'
#' Key helpers called are:
#' - .gl_should_continue_rs_outer_loop (optimizer-rs-loop-stopping.R)
#' - .gl_should_continue_rs_inner_loop (optimizer-rs-loop-stopping.R)
#' - .gl_initialize_rs_parameter_state (optimizer-state.R)
#' - .gl_evaluate_rs_iteration_likelihood_state (optimizer-rs-likelihood.R)
#' - .gl_evaluate_rs_parameter_score_state (optimizer-rs-likelihood.R)
#' - .gl_run_rs_inner_backfitting_step (optimizer-rs-backfitting.R)
#' - .gl_update_rs_outer_iteration_state (optimizer-rs-loop-stopping.R)
#' - .gl_run_rs_inner_iteration_step (optimizer-rs-inner-iteration.R)
#' - .gl_run_rs_parameter_iterations (optimizer-rs-parameter-loop.R)
#' - .gl_build_rs_optimizer_result (optimizer-rs-result.R)
#' 
#' @noRd
.gl_run_rs_optimizer <- function(
    first_outer_run,
    outer_log_lik_change,
    outer_stop_crit,
    outer_only_run_counter,
    max_outer_iter,
    check_elapsed_budget,
    verbose,
    mm,
    rs_calc_eta,
    par_cov,
    par_s,
    margin_dist,
    copula_dist,
    dataset,
    pair_cache,
    margin_eval_cache,
    outer_start_log_lik,
    log_lik_history,
    par_history,
    copula_link,
    discrete_score_method,
    include_dlcopdpar,
    check_dlcopdpar_gradient,
    rs_design_cache,
    rs_smooth_trust_radius,
    lambda_penalty_K,
    lambda_s,
    df_s,
    step_size,
    rs_update_lambda,
    outer_run_counter,
    rs_block_trace,
    use_backtracking,
    backtracking_max_halves,
    plot_results,
    true_val,
    inner_stop_crit,
    max_inner_iter,
    outer_negative_streak,
    step_adjustment,
    max_steps,
    start_step_size,
    max_negative_outer_streak,
    weights_final,
    outer_continue_fn = .gl_should_continue_rs_outer_loop,
    inner_continue_fn = .gl_should_continue_rs_inner_loop,
    initialize_parameter_state_fn = .gl_initialize_rs_parameter_state,
    likelihood_state_fn = .gl_evaluate_rs_iteration_likelihood_state,
    score_state_fn = .gl_evaluate_rs_parameter_score_state,
    backfitting_step_fn = .gl_run_rs_inner_backfitting_step,
    outer_iteration_state_fn = .gl_update_rs_outer_iteration_state,
    inner_iteration_step_fn = .gl_run_rs_inner_iteration_step,
    parameter_iterations_fn = .gl_run_rs_parameter_iterations,
    time_fn = Sys.time,
    result_fn = .gl_build_rs_optimizer_result) {

  # While stopping criteria are not met, run outer loop
  while (outer_continue_fn(
    first_outer_run = first_outer_run,
    outer_log_lik_change = outer_log_lik_change,
    outer_stop_crit = outer_stop_crit,
    outer_only_run_counter = outer_only_run_counter,
    max_outer_iter = max_outer_iter
  )) {
    # Check optional elapsed time budget (often used for debugging / stopping long runs proactively) 
    check_elapsed_budget("RS outer iteration")

    if (verbose > 0) cat(paste("\nOUTER ITERATION:", outer_only_run_counter))

    #Per-outer-loop flag for the first iteration in the outer loop, 
    # used to control some logging and caching behavior
    first_outer_run <- TRUE
    
    # Run the RS parameter iterations (inner GLIM loop)
    rs_parameter_iterations <- parameter_iterations_fn(
      mm = mm,
      verbose = verbose,
      rs_calc_eta = rs_calc_eta,
      par_cov = par_cov,
      par_s = par_s,
      margin_dist = margin_dist,
      copula_dist = copula_dist,
      dataset = dataset,
      pair_cache = pair_cache,
      margin_eval_cache = margin_eval_cache,
      first_outer_run = first_outer_run,
      outer_start_log_lik = outer_start_log_lik,
      log_lik_history = log_lik_history,
      par_history = par_history,
      copula_link = copula_link,
      discrete_score_method = discrete_score_method,
      include_dlcopdpar = include_dlcopdpar,
      check_dlcopdpar_gradient = check_dlcopdpar_gradient,
      outer_only_run_counter = outer_only_run_counter,
      rs_design_cache = rs_design_cache,
      rs_smooth_trust_radius = rs_smooth_trust_radius,
      lambda_penalty_K = lambda_penalty_K,
      lambda_s = lambda_s,
      df_s = df_s,
      step_size = step_size,
      rs_update_lambda = rs_update_lambda,
      outer_run_counter = outer_run_counter,
      rs_block_trace = rs_block_trace,
      use_backtracking = use_backtracking,
      backtracking_max_halves = backtracking_max_halves,
      plot_results = plot_results,
      true_val = true_val,
      inner_stop_crit = inner_stop_crit,
      max_inner_iter = max_inner_iter,
      weights_final = weights_final,
      inner_continue_fn = inner_continue_fn,
      initialize_parameter_state_fn = initialize_parameter_state_fn,
      inner_iteration_step_fn = inner_iteration_step_fn,
      likelihood_state_fn = likelihood_state_fn,
      score_state_fn = score_state_fn,
      backfitting_step_fn = backfitting_step_fn,
      time_fn = time_fn
    )
      par_cov <- rs_parameter_iterations$par_cov
      par_s <- rs_parameter_iterations$par_s
      lambda_s <- rs_parameter_iterations$lambda_s
      df_s <- rs_parameter_iterations$df_s
      calc_lik_out_end <- rs_parameter_iterations$calc_lik_out_end
      rs_block_trace <- rs_parameter_iterations$rs_block_trace
      outer_run_counter <- rs_parameter_iterations$outer_run_counter
      weights_final <- rs_parameter_iterations$weights_final
      first_outer_run <- rs_parameter_iterations$first_outer_run
      outer_start_log_lik <- rs_parameter_iterations$outer_start_log_lik
      log_lik_history <- rs_parameter_iterations$log_lik_history
      par_history <- rs_parameter_iterations$par_history

    # Update outer iteration state (step size, counters, stopping criteria)
    outer_state <- outer_iteration_state_fn(
      calc_lik_out_end = calc_lik_out_end,
      outer_start_log_lik = outer_start_log_lik,
      outer_only_run_counter = outer_only_run_counter,
      outer_negative_streak = outer_negative_streak,
      step_adjustment = step_adjustment,
      max_steps = max_steps,
      start_step_size = start_step_size,
      max_negative_outer_streak = max_negative_outer_streak,
      outer_stop_crit = outer_stop_crit,
      verbose = verbose
    )
      step_size <- outer_state$step_size
      outer_only_run_counter <- outer_state$outer_only_run_counter
      outer_end_log_lik <- outer_state$outer_end_log_lik
      outer_log_lik_change <- outer_state$outer_log_lik_change
      outer_negative_streak <- outer_state$outer_negative_streak    
  }

  # Build and return the final RS optimizer output after the outer loop has finished
  result_fn(
    par_cov = par_cov,
    par_s = par_s,
    lambda_s = lambda_s,
    df_s = df_s,
    calc_lik_out_end = calc_lik_out_end,
    rs_block_trace = rs_block_trace,
    outer_run_counter = outer_run_counter,
    outer_only_run_counter = outer_only_run_counter,
    outer_end_log_lik = outer_end_log_lik,
    outer_log_lik_change = outer_log_lik_change,
    outer_negative_streak = outer_negative_streak,
    step_size = step_size,
    weights_final = weights_final,
    log_lik_history = log_lik_history,
    par_history = par_history,
    first_outer_run = first_outer_run,
    outer_start_log_lik = outer_start_log_lik
  )
}
