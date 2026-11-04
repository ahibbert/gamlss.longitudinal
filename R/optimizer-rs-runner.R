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
    result_fn = .gl_build_rs_optimizer_result,
    final_weights_fn = .gl_refresh_rs_final_weights) {


  has_smooths <- any(vapply(rs_design_cache, function(x) length(x$smooth_penalty_meta) > 0L, logical(1)))
  n_smooths <- sum(vapply(rs_design_cache, function(x) length(x$smooth_penalty_meta), integer(1)))
  phase <- if (has_smooths && isTRUE(rs_update_lambda)) "selection" else "fixed_lambda"
  selection_stop_reason <- if (identical(phase, "selection")) NA_character_ else "not_requested"
  lambda_epoch <- 0L
  raw_negative_streak <- 0L
  objective_change <- NA_real_
  objective <- NA_real_
  settled <- FALSE
  stop_reason <- "max_iterations"
  detail_reason <- "outer_iteration_limit"
  sweep_trace <- list()
  lambda_searches <- list()
  block_status <- list()
  comparable <- FALSE
  objective_start <- NA_real_
  outer_end_log_lik <- NA_real_
  calc_lik_out_end <- list(
    log_lik = c(marginal = NA_real_, copula = NA_real_, joint = NA_real_), valid = FALSE)
  final_edf_refreshed <- FALSE

  # Both phases use the same counters and elapsed budget. The outer helper is
  # retained as the iteration-limit seam; protocol stopping is assessed below.
  while (outer_continue_fn(
    first_outer_run = TRUE, outer_log_lik_change = outer_log_lik_change,
    outer_stop_crit = outer_stop_crit, outer_only_run_counter = outer_only_run_counter,
    max_outer_iter = max_outer_iter, stop_on_convergence = FALSE
  )) {
    check_elapsed_budget("RS outer iteration")
    .gl_progress(verbose, "Outer iteration ", outer_only_run_counter, " (", phase, ").",
                 prefix = "gamlss_longitudinal")
    first_outer_run <- TRUE
    lambda_start_sweep <- lambda_s
    par_s_start_sweep <- par_s
    trace_start <- length(rs_block_trace)
    phase_this_sweep <- phase
    rs_parameter_iterations <- parameter_iterations_fn(
      mm = mm, verbose = verbose, rs_calc_eta = rs_calc_eta, par_cov = par_cov,
      par_s = par_s, margin_dist = margin_dist, copula_dist = copula_dist,
      dataset = dataset, pair_cache = pair_cache, margin_eval_cache = margin_eval_cache,
      first_outer_run = first_outer_run, outer_start_log_lik = outer_start_log_lik,
      log_lik_history = log_lik_history, par_history = par_history, copula_link = copula_link,
      discrete_score_method = discrete_score_method, include_dlcopdpar = include_dlcopdpar,
      check_dlcopdpar_gradient = check_dlcopdpar_gradient,
      outer_only_run_counter = outer_only_run_counter, rs_design_cache = rs_design_cache,
      rs_smooth_trust_radius = rs_smooth_trust_radius, lambda_penalty_K = lambda_penalty_K,
      lambda_s = lambda_s, df_s = df_s, step_size = step_size,
      rs_update_lambda = identical(phase, "selection"), outer_run_counter = outer_run_counter,
      rs_block_trace = rs_block_trace, use_backtracking = use_backtracking,
      backtracking_max_halves = backtracking_max_halves, plot_results = plot_results,
      true_val = true_val, inner_stop_crit = inner_stop_crit, max_inner_iter = max_inner_iter,
      weights_final = weights_final, inner_continue_fn = inner_continue_fn,
      initialize_parameter_state_fn = initialize_parameter_state_fn,
      inner_iteration_step_fn = inner_iteration_step_fn, likelihood_state_fn = likelihood_state_fn,
      score_state_fn = score_state_fn, backfitting_step_fn = backfitting_step_fn,
      time_fn = time_fn, check_elapsed_budget = check_elapsed_budget
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
    block_status <- rs_parameter_iterations$block_status %||% list()
    searches <- rs_parameter_iterations$lambda_searches %||% list()
    lambda_searches <- c(lambda_searches, searches)
    if (length(rs_block_trace) > trace_start) {
      for (i in seq.int(trace_start + 1L, length(rs_block_trace))) {
        rs_block_trace[[i]]$phase <- phase_this_sweep
        rs_block_trace[[i]]$lambda_epoch_before <- lambda_epoch
        if (isTRUE(rs_block_trace[[i]]$lambda_changed)) lambda_epoch <- lambda_epoch + 1L
        rs_block_trace[[i]]$lambda_epoch <- lambda_epoch
      }
    }
    lambda_changed <- !.gl_rs_same_lambda(lambda_start_sweep, lambda_s)
    comparable <- identical(phase_this_sweep, "fixed_lambda") && (!has_smooths || !lambda_changed)
    objective <- .gl_rs_objective(calc_lik_out_end, par_s, lambda_s, rs_design_cache)
    objective_start <- if (comparable) as.numeric(outer_start_log_lik) -
      .gl_rs_penalty_value(par_s_start_sweep, lambda_start_sweep, rs_design_cache) else NA_real_
    objective_change <- if (comparable) objective - objective_start else NA_real_
    outer_state <- outer_iteration_state_fn(
      calc_lik_out_end = calc_lik_out_end, outer_start_log_lik = outer_start_log_lik,
      outer_only_run_counter = outer_only_run_counter, outer_negative_streak = outer_negative_streak,
      step_adjustment = step_adjustment, max_steps = max_steps, start_step_size = start_step_size,
      max_negative_outer_streak = max_negative_outer_streak,
      outer_stop_crit = outer_stop_crit, verbose = verbose,
      objective_change = objective_change, objective_kind = "conditional_penalized",
      comparable = comparable, report_convergence = FALSE
    )
    step_size <- outer_state$step_size
    outer_only_run_counter <- outer_state$outer_only_run_counter
    outer_end_log_lik <- outer_state$outer_end_log_lik
    outer_log_lik_change <- outer_state$outer_log_lik_change
    outer_negative_streak <- outer_state$outer_negative_streak
    raw_negative_streak <- if (is.finite(outer_log_lik_change) && outer_log_lik_change < 0)
      raw_negative_streak + 1L else 0L
    finite_state <- is.finite(objective) && is.finite(as.numeric(outer_end_log_lik)) &&
      !identical(calc_lik_out_end$valid, FALSE)
    complete <- isTRUE(rs_parameter_iterations$complete_sweep)
    unresolved <- unlist(block_status, use.names = FALSE)
    rejected <- any(unresolved == "rejected")
    inner_limit <- any(unresolved == "max_iterations")
    failed_search <- isTRUE(rs_parameter_iterations$lambda_search_failed)
    plateau <- comparable && finite_state && is.finite(objective_change) &&
      abs(objective_change) <= outer_stop_crit
    settled <- plateau && objective_change >= 0 && complete && !rejected && !inner_limit &&
      !any(unresolved %in% c("invalid_objective", "lambda_search_failed"))
    sweep_trace[[length(sweep_trace) + 1L]] <- list(
      iteration = outer_only_run_counter - 1L, phase = phase_this_sweep,
      lambda_epoch = lambda_epoch, raw_change = as.numeric(outer_log_lik_change),
      objective = objective, objective_change = objective_change, comparable = comparable,
      complete = complete, block_status = block_status, objective_plateau = settled
    )
    check_elapsed_budget("RS completed sweep")
    if (failed_search) {
      stop_reason <- "numerical_failure"; detail_reason <- "lambda_search_failed"; settled <- FALSE
      break
    }
    if (!finite_state || (identical(phase_this_sweep, "fixed_lambda") && !comparable)) {
      stop_reason <- if (identical(calc_lik_out_end$valid, FALSE) ||
        !is.finite(as.numeric(outer_end_log_lik))) "invalid_likelihood" else "numerical_failure"
      detail_reason <- if (!comparable) "fixed_lambda_changed" else "invalid_final_state"
      settled <- FALSE
      break
    }
    if (identical(phase_this_sweep, "selection")) {
      eligible_search_sweep <- complete && length(searches) >= n_smooths &&
        all(vapply(searches, function(x) identical(x$status, "selected"), logical(1)))
      if (eligible_search_sweep && is.finite(outer_log_lik_change) &&
          abs(outer_log_lik_change) <= outer_stop_crit) {
        phase <- "fixed_lambda"
        selection_stop_reason <- "raw_likelihood_plateau"
        outer_negative_streak <- 0L
      }
      next
    }
    if (plateau && rejected) {
      stop_reason <- "max_stall"; detail_reason <- "terminal_block_rejected"
      break
    }
    if (plateau && inner_limit) {
      stop_reason <- "max_iterations"; detail_reason <- "terminal_inner_iteration_limit"
      break
    }
    if (settled && isTRUE(stop_on_convergence)) break
  }
  if (settled) {
    stop_reason <- "converged"; detail_reason <- "conditional_objective_plateau"
  } else if (identical(stop_reason, "max_iterations") &&
             any(unlist(block_status, use.names = FALSE) == "max_iterations")) {
    detail_reason <- "terminal_inner_iteration_limit"
  }
  # Recompute working weights at the common final q, not the successive
  # pre-update block states. No coefficient or smoothing update is made here.
  if (has_smooths && length(sweep_trace) > 0L && is.finite(objective) &&
      is.finite(as.numeric(outer_end_log_lik)) && !identical(calc_lik_out_end$valid, FALSE)) {
    check_elapsed_budget("RS final working weights")
    final_state <- final_weights_fn(
      par_cov = par_cov, par_s = par_s, mm = mm, rs_calc_eta = rs_calc_eta,
      margin_dist = margin_dist, copula_dist = copula_dist, dataset = dataset,
      pair_cache = pair_cache, margin_eval_cache = margin_eval_cache,
      copula_link = copula_link, discrete_score_method = discrete_score_method,
      include_dlcopdpar = include_dlcopdpar, check_elapsed_budget = check_elapsed_budget
    )
    weights_final <- final_state$weights
    calc_lik_out_end <- final_state$likelihood
    outer_end_log_lik <- calc_lik_out_end$log_lik["joint"]
    objective <- .gl_rs_objective(calc_lik_out_end, par_s, lambda_s, rs_design_cache)
    outer_log_lik_change <- outer_end_log_lik - outer_start_log_lik
    if (comparable) {
      objective_change <- objective - objective_start
      if (settled && (!is.finite(objective_change) || objective_change < 0 ||
                      abs(objective_change) > outer_stop_crit)) {
        settled <- FALSE
        stop_reason <- "numerical_failure"
        detail_reason <- "final_refresh_changed_objective_contract"
      }
    }
    df_s <- .gl_refresh_rs_smooth_edf(df_s, mm$s, lambda_s, weights_final)
    final_edf_refreshed <- TRUE
    last <- length(sweep_trace)
    sweep_trace[[last]]$objective <- objective
    sweep_trace[[last]]$objective_change <- objective_change
    sweep_trace[[last]]$raw_change <- as.numeric(outer_log_lik_change)
    sweep_trace[[last]]$objective_plateau <- settled
    check_elapsed_budget("RS final EDF")
  }
  out <- result_fn(
    par_cov = par_cov, par_s = par_s, lambda_s = lambda_s, df_s = df_s,
    calc_lik_out_end = calc_lik_out_end, rs_block_trace = rs_block_trace,
    outer_run_counter = outer_run_counter, outer_only_run_counter = outer_only_run_counter,
    outer_end_log_lik = outer_end_log_lik, outer_log_lik_change = outer_log_lik_change,
    outer_negative_streak = outer_negative_streak, step_size = step_size,
    weights_final = weights_final, log_lik_history = log_lik_history, par_history = par_history,
    first_outer_run = first_outer_run, outer_start_log_lik = outer_start_log_lik
  )
  out$rs_protocol <- list(
    version = 1L, target = "joint_loglik_minus_half_quadratic_penalty",
    phase = phase, objective = objective, objective_change = objective_change,
    objective_tolerance = outer_stop_crit, objective_plateau = settled,
    stop_reason = stop_reason, detail_reason = detail_reason,
    selection_stop_reason = selection_stop_reason, lambda_converged = NA,
    selection_status = if (!has_smooths || !isTRUE(rs_update_lambda)) "not_requested" else
      if (identical(detail_reason, "lambda_search_failed")) "failed" else
      if (identical(selection_stop_reason, "raw_likelihood_plateau")) "ended_by_raw_plateau" else "incomplete",
    settling_status = if (settled) "objective_plateau" else
      if (any(vapply(sweep_trace, function(x) identical(x$phase, "fixed_lambda"), logical(1))))
        "incomplete" else "not_started",
    lambda_epoch = lambda_epoch, selected_lambda = lambda_s,
    lambda_searches = lambda_searches, block_status = block_status,
    raw_negative_streak = raw_negative_streak, sweeps = sweep_trace,
    stationarity_status = "not_checked", final_edf_refreshed = final_edf_refreshed
  )
  out
}
