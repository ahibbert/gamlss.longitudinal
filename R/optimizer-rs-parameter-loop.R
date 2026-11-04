#' Runs the RS inner loop, iterating over each parameter in the model
#' and updating parameter estimates and likelihood etc.
#'
#' @noRd
.gl_run_rs_parameter_iterations <- function(
    mm,
    verbose,
    rs_calc_eta,
    par_cov,
    par_s,
    margin_dist,
    copula_dist,
    dataset,
    pair_cache,
    margin_eval_cache,
    first_outer_run,
    outer_start_log_lik,
    log_lik_history,
    par_history,
    copula_link,
    discrete_score_method,
    include_dlcopdpar,
    check_dlcopdpar_gradient,
    outer_only_run_counter,
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
    weights_final,
    inner_continue_fn,
    initialize_parameter_state_fn,
    inner_iteration_step_fn,
    likelihood_state_fn,
    score_state_fn,
    backfitting_step_fn,
    time_fn,
    check_elapsed_budget = function(stage) invisible(NULL)) {
  block_status <- list()
  lambda_searches <- list()
  lambda_search_failed <- FALSE

  # For each parameter in the model, run the inner iteration (GLIM) loop to update the parameter estimates
  for (par_name in names(mm$x)) {
    if (verbose > 2) {
      cat(paste("\nINNER ITERATION: Parameter:", par_name))
    }

    rs_parameter_state <- initialize_parameter_state_fn()
    first_inner_run <- rs_parameter_state$first_inner_run
    change_log_lik <- rs_parameter_state$change_log_lik
    change_objective <- 0
    terminal_status <- "not_run"
    beta_change_inner <- rs_parameter_state$beta_change_inner
    run_counter <- rs_parameter_state$run_counter
    inner_run_counter <- rs_parameter_state$inner_run_counter

    # INNER ITERATION (GLIM)

    # While inner stopping criteria are not met, run the inner iteration (GLIM) loop for the current parameter
    while (inner_continue_fn(
      first_inner_run = first_inner_run,
      change_log_lik = change_objective,
      inner_stop_crit = inner_stop_crit,
      inner_run_counter = inner_run_counter,
      max_inner_iter = max_inner_iter
    )) {
      timer <- c()
      check_elapsed_budget("RS inner iteration")
      timer_start <- time_fn()

      first_inner_run <- FALSE

      # Run an inner step for the current parameter, updating the parameter estimates and state
      rs_inner_step <- inner_iteration_step_fn(
        par_name = par_name,
        rs_calc_eta = rs_calc_eta,
        par_cov = par_cov,
        par_s = par_s,
        mm = mm,
        margin_dist = margin_dist,
        copula_dist = copula_dist,
        dataset = dataset,
        pair_cache = pair_cache,
        margin_eval_cache = margin_eval_cache,
        current_calc_lik_out = if (exists("calc_lik_out_end", inherits = FALSE)) calc_lik_out_end else NULL,
        first_outer_run = first_outer_run,
        outer_start_log_lik = outer_start_log_lik,
        timer = timer,
        timer_start = timer_start,
        log_lik_history = log_lik_history,
        par_history = par_history,
        copula_link = copula_link,
        discrete_score_method = discrete_score_method,
        include_dlcopdpar = include_dlcopdpar,
        check_dlcopdpar_gradient = check_dlcopdpar_gradient,
        outer_only_run_counter = outer_only_run_counter,
        verbose = verbose,
        rs_design_cache = rs_design_cache,
        rs_smooth_trust_radius = rs_smooth_trust_radius,
        lambda_penalty_K = lambda_penalty_K,
        lambda_s = lambda_s,
        df_s = df_s,
        step_size = step_size,
        rs_update_lambda = rs_update_lambda,
        inner_run_counter = inner_run_counter,
        outer_run_counter = outer_run_counter,
        rs_block_trace = rs_block_trace,
        run_counter = run_counter,
        use_backtracking = use_backtracking,
        backtracking_max_halves = backtracking_max_halves,
        plot_results = plot_results,
        true_val = true_val,
        likelihood_state_fn = likelihood_state_fn,
        score_state_fn = score_state_fn,
        backfitting_step_fn = backfitting_step_fn
      )
      score <- rs_inner_step$score
      lambda_s <- rs_inner_step$lambda_s
      par_cov <- rs_inner_step$par_cov
      par_s <- rs_inner_step$par_s
      calc_lik_out_end <- rs_inner_step$calc_lik_out_end
      df_s <- rs_inner_step$df_s
      rs_block_trace <- rs_inner_step$rs_block_trace
      change_log_lik <- rs_inner_step$change_log_lik
      change_objective <- rs_inner_step$change_objective %||% change_log_lik
      terminal_status <- rs_inner_step$coefficient_step_status %||% "accepted"
      search <- rs_inner_step$lambda_search %||% list(status = "not_searched", searches = list())
      lambda_searches <- c(lambda_searches, search$searches)
      lambda_search_failed <- identical(search$status, "failed")
      run_counter <- rs_inner_step$run_counter
      outer_run_counter <- rs_inner_step$outer_run_counter
      inner_run_counter <- rs_inner_step$inner_run_counter
      timer <- rs_inner_step$timer
      first_outer_run <- rs_inner_step$first_outer_run
      outer_start_log_lik <- rs_inner_step$outer_start_log_lik
      log_lik_history <- rs_inner_step$log_lik_history
      par_history <- rs_inner_step$par_history

      if (lambda_search_failed || identical(terminal_status, "rejected")) break
    }

    weights_final[[par_name]] <- score$w_k
    block_status[[par_name]] <- if (lambda_search_failed) {
      "lambda_search_failed"
    } else if (identical(terminal_status, "rejected")) {
      "rejected"
    } else if (!is.finite(change_objective)) {
      "invalid_objective"
    } else if (inner_run_counter > max_inner_iter && abs(change_objective) > inner_stop_crit) {
      "max_iterations"
    } else {
      "objective_plateau"
    }
    if (lambda_search_failed) break
  }

  list(
    par_cov = par_cov,
    par_s = par_s,
    lambda_s = lambda_s,
    df_s = df_s,
    calc_lik_out_end = calc_lik_out_end,
    rs_block_trace = rs_block_trace,
    outer_run_counter = outer_run_counter,
    weights_final = weights_final,
    first_outer_run = first_outer_run,
    outer_start_log_lik = outer_start_log_lik,
    log_lik_history = log_lik_history,
    par_history = par_history,
    block_status = block_status,
    complete_sweep = length(block_status) == length(mm$x),
    lambda_search_failed = lambda_search_failed,
    lambda_searches = lambda_searches
  )
}
