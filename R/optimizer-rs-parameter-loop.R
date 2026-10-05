#' Runs the RS inner loop, iterating over each parameter in the model
#' and updating parameter estimates and likelihood etc.
#'
#' @noRd
.gl_update_rs_inner_budget <- function(
    current_budget,
    max_inner_iter,
    iterations_used,
    final_loglik_change,
    inner_stop_crit) {
  still_moving_at_cap <- iterations_used >= current_budget &&
    is.finite(final_loglik_change) &&
    abs(final_loglik_change) > inner_stop_crit
  if (!isTRUE(still_moving_at_cap) || current_budget >= max_inner_iter) {
    return(as.integer(current_budget))
  }
  min(
    as.integer(max_inner_iter),
    max(as.integer(current_budget) + 1L, 2L * as.integer(current_budget))
  )
}

#' Adapt one RS block's nominal step after its accepted/rejected trace
#'
#' @noRd
.gl_update_rs_block_step <- function(
    nominal_step,
    block_trace,
    adaptive_step_size,
    previously_rejected = FALSE) {
  if (!isTRUE(adaptive_step_size) || !is.data.frame(block_trace) || !nrow(block_trace)) {
    return(list(step_size = nominal_step, rejected = FALSE))
  }
  accepted <- block_trace$accepted_step_size[
    is.finite(block_trace$accepted_step_size) & block_trace$accepted_step_size > 0
  ]
  rejected <- any(block_trace$rejected %in% TRUE) && !length(accepted)
  if (rejected && !isTRUE(previously_rejected)) {
    return(list(step_size = nominal_step / 2, rejected = TRUE))
  }
  # A block that cannot move even on its one reduced-step, no-backtracking
  # retry has exhausted the adaptive proposal. Do not halve indefinitely or
  # prevent the unchanged outer convergence contract from being evaluated.
  if (rejected) {
    return(list(step_size = nominal_step, rejected = FALSE))
  }
  if (length(accepted) && any(block_trace$backtracking_attempts > 0L)) {
    nominal_step <- min(nominal_step, stats::median(accepted))
  }
  list(step_size = nominal_step, rejected = FALSE)
}

#' Slice newly appended RS trace rows regardless of storage representation
#'
#' @noRd
.gl_rs_new_block_trace <- function(rs_block_trace, start_index) {
  trace_size <- if (is.data.frame(rs_block_trace)) {
    nrow(rs_block_trace)
  } else {
    length(rs_block_trace)
  }
  if (trace_size <= start_index) return(data.frame())
  indices <- seq.int(start_index + 1L, trace_size)
  if (is.data.frame(rs_block_trace)) {
    return(rs_block_trace[indices, , drop = FALSE])
  }
  do.call(rbind, rs_block_trace[indices])
}

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
    rs_adaptive_step_size = FALSE,
    rs_block_step_sizes = stats::setNames(rep(step_size, length(mm$x)), names(mm$x)),
    rs_block_rejected_last = stats::setNames(rep(FALSE, length(mm$x)), names(mm$x)),
    plot_results,
    true_val,
    inner_stop_crit,
    max_inner_iter,
    rs_inner_budgets = stats::setNames(
      rep(as.integer(max_inner_iter), length(mm$x)),
      names(mm$x)
    ),
    weights_final,
    inner_continue_fn,
    initialize_parameter_state_fn,
    inner_iteration_step_fn,
    likelihood_state_fn,
    score_state_fn,
    backfitting_step_fn,
    time_fn) {

  # For each parameter in the model, run the inner iteration (GLIM) loop to update the parameter estimates
  for (par_name in names(mm$x)) {
    block_trace_start <- if (is.data.frame(rs_block_trace)) {
      nrow(rs_block_trace)
    } else {
      length(rs_block_trace)
    }
    if (verbose > 2) {
      cat(paste("\nINNER ITERATION: Parameter:", par_name))
    }

    rs_parameter_state <- initialize_parameter_state_fn()
    first_inner_run <- rs_parameter_state$first_inner_run
    change_log_lik <- rs_parameter_state$change_log_lik
    beta_change_inner <- rs_parameter_state$beta_change_inner
    run_counter <- rs_parameter_state$run_counter
    inner_run_counter <- rs_parameter_state$inner_run_counter
    block_inner_budget <- min(
      as.integer(max_inner_iter),
      as.integer(rs_inner_budgets[[par_name]] %||% max_inner_iter)
    )
    block_step_size <- min(
      as.numeric(step_size),
      as.numeric(rs_block_step_sizes[[par_name]] %||% step_size)
    )
    block_backtracking_halves <- if (
      isTRUE(rs_adaptive_step_size) &&
      isTRUE(rs_block_rejected_last[[par_name]])
    ) 0L else backtracking_max_halves

    # INNER ITERATION (GLIM)

    # While inner stopping criteria are not met, run the inner iteration (GLIM) loop for the current parameter
    while (inner_continue_fn(
      first_inner_run = first_inner_run,
      change_log_lik = change_log_lik,
      inner_stop_crit = inner_stop_crit,
      inner_run_counter = inner_run_counter,
      max_inner_iter = block_inner_budget
    )) {
      timer <- c()
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
        step_size = block_step_size,
        rs_update_lambda = rs_update_lambda,
        inner_run_counter = inner_run_counter,
        outer_run_counter = outer_run_counter,
        rs_block_trace = rs_block_trace,
        run_counter = run_counter,
        use_backtracking = use_backtracking,
        backtracking_max_halves = block_backtracking_halves,
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
      run_counter <- rs_inner_step$run_counter
      outer_run_counter <- rs_inner_step$outer_run_counter
      inner_run_counter <- rs_inner_step$inner_run_counter
      timer <- rs_inner_step$timer
      first_outer_run <- rs_inner_step$first_outer_run
      outer_start_log_lik <- rs_inner_step$outer_start_log_lik
      log_lik_history <- rs_inner_step$log_lik_history
      par_history <- rs_inner_step$par_history

      # print(timer)
    }

    iterations_used <- as.integer(inner_run_counter) - 1L
    rs_inner_budgets[[par_name]] <- .gl_update_rs_inner_budget(
      current_budget = block_inner_budget,
      max_inner_iter = max_inner_iter,
      iterations_used = iterations_used,
      final_loglik_change = change_log_lik,
      inner_stop_crit = inner_stop_crit
    )
    block_trace <- .gl_rs_new_block_trace(rs_block_trace, block_trace_start)
    step_state <- .gl_update_rs_block_step(
      nominal_step = block_step_size,
      block_trace = block_trace,
      adaptive_step_size = rs_adaptive_step_size,
      previously_rejected = isTRUE(rs_block_rejected_last[[par_name]])
    )
    rs_block_step_sizes[[par_name]] <- step_state$step_size
    rs_block_rejected_last[[par_name]] <- step_state$rejected

    weights_final[[par_name]] <- score$w_k
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
    rs_inner_budgets = rs_inner_budgets,
    rs_block_step_sizes = rs_block_step_sizes,
    rs_block_rejected_last = rs_block_rejected_last
  )
}
