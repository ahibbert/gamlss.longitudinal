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
.gl_rs_lambda_update_scheduled <- function(
    outer_iteration,
    update_lambda,
    lambda_update_every) {
  !isTRUE(update_lambda) ||
    lambda_update_every == 1L ||
    (outer_iteration != 1L &&
      ((outer_iteration - 2L) %% lambda_update_every == 0L))
}

#' @noRd
.gl_rs_outer_change_for_stopping <- function(
    objective_change,
    update_lambda,
    lambda_update_every,
    lambda_update_scheduled) {
  if (isTRUE(update_lambda) && lambda_update_every > 1L &&
      !isTRUE(lambda_update_scheduled)) {
    return(Inf)
  }
  objective_change
}

#' @noRd
.gl_rs_outer_change_after_step_adaptation <- function(
    objective_change,
    adaptive_step_size,
    block_rejected_last) {
  if (isTRUE(adaptive_step_size) && any(block_rejected_last %in% TRUE)) {
    return(Inf)
  }
  objective_change
}

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
    rs_lambda_update_every = 1L,
    rs_bfgs_refine = FALSE,
    rs_bfgs_refine_after = 10L,
    rs_bfgs_method = "BFGS",
    rs_bfgs_max_refinements = 1L,
    rs_bfgs_max_iter = 100L,
    rs_bfgs_reltol = 1e-12,
    rs_bfgs_parscale_max = 25,
    rs_bfgs_active_parameters = NULL,
    rs_bfgs_coordinate_rescue = FALSE,
    rs_bfgs_coordinate_threshold = 1.5,
    rs_bfgs_coordinate_radius = 4,
    rs_bfgs_coordinate_max = 10L,
    rs_bfgs_coordinate_tol = 1e-3,
    rs_bfgs_repolish_iter = 50L,
    outer_run_counter,
    rs_block_trace,
    use_backtracking,
    backtracking_max_halves,
    rs_adaptive_step_size = FALSE,
    plot_results,
    true_val,
    inner_stop_crit,
    max_inner_iter,
    initial_inner_iter = max_inner_iter,
    outer_negative_streak,
    step_adjustment,
    max_steps,
    start_step_size,
    max_negative_outer_streak,
    weights_final,
    stop_on_convergence = TRUE,
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

  rs_inner_budgets <- stats::setNames(
    rep(as.integer(initial_inner_iter), length(mm$x)),
    names(mm$x)
  )
  rs_block_step_sizes <- stats::setNames(rep(step_size, length(mm$x)), names(mm$x))
  rs_block_rejected_last <- stats::setNames(rep(FALSE, length(mm$x)), names(mm$x))
  rs_bfgs_refinements <- 0L
  rs_bfgs_diagnostics <- list(enabled = isTRUE(rs_bfgs_refine), used = FALSE)

  # While stopping criteria are not met, run outer loop
  while (outer_continue_fn(
    first_outer_run = first_outer_run,
    outer_log_lik_change = outer_log_lik_change,
    outer_stop_crit = outer_stop_crit,
    outer_only_run_counter = outer_only_run_counter,
    max_outer_iter = max_outer_iter,
    stop_on_convergence = stop_on_convergence
  )) {
    # Check optional elapsed time budget (often used for debugging / stopping long runs proactively) 
    check_elapsed_budget("RS outer iteration")

    .gl_progress(
      verbose,
      "Outer iteration ", outer_only_run_counter, ".",
      prefix = "gamlss_longitudinal"
    )

    #Per-outer-loop flag for the first iteration in the outer loop, 
    # used to control some logging and caching behavior
    first_outer_run <- TRUE
    lambda_update_scheduled <- .gl_rs_lambda_update_scheduled(
      outer_iteration = outer_only_run_counter,
      update_lambda = rs_update_lambda,
      lambda_update_every = rs_lambda_update_every
    )

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
      rs_update_lambda = isTRUE(rs_update_lambda) && isTRUE(lambda_update_scheduled),
      outer_run_counter = outer_run_counter,
      rs_block_trace = rs_block_trace,
      use_backtracking = use_backtracking,
      backtracking_max_halves = backtracking_max_halves,
      rs_adaptive_step_size = rs_adaptive_step_size,
      rs_block_step_sizes = rs_block_step_sizes,
      rs_block_rejected_last = rs_block_rejected_last,
      plot_results = plot_results,
      true_val = true_val,
      inner_stop_crit = inner_stop_crit,
      max_inner_iter = max_inner_iter,
      rs_inner_budgets = rs_inner_budgets,
      weights_final = weights_final,
      check_elapsed_budget = check_elapsed_budget,
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
      rs_inner_budgets <- rs_parameter_iterations$rs_inner_budgets %||%
        rs_inner_budgets
      rs_block_step_sizes <- rs_parameter_iterations$rs_block_step_sizes %||%
        rs_block_step_sizes
      rs_block_rejected_last <- rs_parameter_iterations$rs_block_rejected_last %||%
        rs_block_rejected_last

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
      outer_objective_change <- outer_log_lik_change
      outer_negative_streak <- outer_state$outer_negative_streak    
      # Do not declare convergence between scheduled smoothing-parameter
      # updates. The next scheduled cycle must verify that re-optimising
      # lambda no longer causes a material objective change.
      outer_log_lik_change <- .gl_rs_outer_change_for_stopping(
        objective_change = outer_log_lik_change,
        update_lambda = rs_update_lambda,
        lambda_update_every = rs_lambda_update_every,
        lambda_update_scheduled = lambda_update_scheduled
      )
      # A rejected adaptive block has just reduced its nominal step for the
      # next outer cycle. Do not declare convergence before that smaller step
      # has actually been tried; otherwise a zero accepted change can stop the
      # fit at the very point where adaptation was meant to recover it.
      outer_log_lik_change <- .gl_rs_outer_change_after_step_adaptation(
        objective_change = outer_log_lik_change,
        adaptive_step_size = rs_adaptive_step_size,
        block_rejected_last = rs_block_rejected_last
      )

      completed_outer <- outer_only_run_counter - 1L
      if (.gl_should_run_rs_bfgs_refinement(
        completed_outer = completed_outer,
        enabled = rs_bfgs_refine,
        refine_after = rs_bfgs_refine_after,
        already_refined = rs_bfgs_refinements >= rs_bfgs_max_refinements,
        refinements_completed = rs_bfgs_refinements,
        max_refinements = rs_bfgs_max_refinements,
        lambda_update_scheduled = lambda_update_scheduled,
        objective_change = outer_objective_change,
        outer_tol = outer_stop_crit,
        prior_refinement_incomplete = isTRUE(
          (rs_bfgs_diagnostics$convergence %||% 0L) != 0L
        )
      )) {
        refinement <- .gl_run_rs_bfgs_refinement(
          mm = mm,
          par_cov = par_cov,
          par_s = par_s,
          lambda_s = lambda_s,
          margin_dist = margin_dist,
          copula_link = copula_link,
          copula_dist = copula_dist,
          dataset = dataset,
          pair_cache = pair_cache,
          margin_eval_cache = margin_eval_cache,
          include_dlcopdpar = include_dlcopdpar,
          method = rs_bfgs_method,
          max_iter = rs_bfgs_max_iter,
          reltol = rs_bfgs_reltol,
          parscale_max = rs_bfgs_parscale_max,
          active_parameters = rs_bfgs_active_parameters,
          coordinate_rescue = rs_bfgs_coordinate_rescue,
          coordinate_threshold = rs_bfgs_coordinate_threshold,
          coordinate_radius = rs_bfgs_coordinate_radius,
          coordinate_max = rs_bfgs_coordinate_max,
          coordinate_tol = rs_bfgs_coordinate_tol,
          repolish_iter = rs_bfgs_repolish_iter,
          check_elapsed_budget = check_elapsed_budget,
          verbose = verbose
        )
        par_cov <- refinement$par_cov
        par_s <- refinement$par_s
        rs_bfgs_refinements <- rs_bfgs_refinements + 1L
        prior_runs <- rs_bfgs_diagnostics$runs %||% list()
        rs_bfgs_diagnostics <- refinement$diagnostics
        rs_bfgs_diagnostics$runs <- c(prior_runs, list(refinement$diagnostics))
        rs_bfgs_diagnostics$refinements <- rs_bfgs_refinements
        if (isTRUE(refinement$diagnostics$accepted)) {
          # Require subsequent ordinary RS cycles to establish the unchanged
          # raw-likelihood and scheduled-lambda convergence contract.
          outer_log_lik_change <- Inf
        }
      }
  }

  # Build and return the final RS optimizer output after the outer loop has finished
  result <- result_fn(
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
  result$rs_inner_budgets <- rs_inner_budgets
  result$rs_block_step_sizes <- rs_block_step_sizes
  result$rs_block_rejected_last <- rs_block_rejected_last
  result$rs_bfgs_diagnostics <- rs_bfgs_diagnostics
  result
}
