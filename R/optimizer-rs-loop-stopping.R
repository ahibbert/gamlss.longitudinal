.gl_update_outer_negative_streak <- function(
    outer_log_lik_change,
    outer_negative_streak,
    max_negative_outer_streak) {
  if (is.finite(outer_log_lik_change) && outer_log_lik_change < 0) {
    outer_negative_streak <- outer_negative_streak + 1
  } else {
    outer_negative_streak <- 0
  }

  if (outer_negative_streak >= max_negative_outer_streak) {
    msg <- paste0(
      "Optimization stopped after ", max_negative_outer_streak, " consecutive negative outer log-likelihood changes. ",
      "We believe the model may be misspecified and the likelihood may be malformed. ",
      "Try different starting parameters or covariate combinations. Other options include switching between joint and separate optimisation. In general, joint optimisation provides more stable convergence.",
      "Alternatively, you can increase the max_negative_outer_streak parameter to allow more negative changes before stopping, but we recommend investigating the cause of the consecutive negative changes in likelihood."
    )

    warning(msg, call. = FALSE)
    stop(msg, call. = FALSE)
  }

  outer_negative_streak
}

#' Determine whether the RS inner loop should continue
#' (stopping criteria not met and max iterations not exceeded)
#'
#' @noRd
.gl_should_continue_rs_inner_loop <- function(
    first_inner_run,
    change_log_lik,
    inner_stop_crit,
    inner_run_counter,
    max_inner_iter) {
  (isTRUE(first_inner_run) || abs(change_log_lik) > inner_stop_crit) &&
    inner_run_counter < max_inner_iter
}

#' Determine whether the RS outer loop should continue
#' (stopping criteria not met and max iterations not exceeded)
#'
#' @noRd
.gl_should_continue_rs_outer_loop <- function(
    first_outer_run,
    outer_log_lik_change,
    outer_stop_crit,
    outer_only_run_counter,
    max_outer_iter) {
  (isTRUE(first_outer_run) || abs(outer_log_lik_change) > outer_stop_crit) &&
    outer_only_run_counter < max_outer_iter
}

.gl_update_rs_outer_iteration_state <- function(
    calc_lik_out_end,
    outer_start_log_lik,
    outer_only_run_counter,
    outer_negative_streak,
    step_adjustment,
    max_steps,
    start_step_size,
    max_negative_outer_streak,
    outer_stop_crit,
    verbose,
    summary_fn = .gl_print_outer_iteration_summary,
    streak_fn = .gl_update_outer_negative_streak,
    convergence_fn = .gl_print_outer_convergence) {

  step_size <- (step_adjustment^min(outer_only_run_counter, max_steps)) * start_step_size
  outer_only_run_counter <- outer_only_run_counter + 1
  outer_end_log_lik <- calc_lik_out_end$log_lik["joint"]
  outer_log_lik_change <- outer_end_log_lik - outer_start_log_lik

  summary_fn(
    outer_start_log_lik = outer_start_log_lik,
    outer_end_log_lik = outer_end_log_lik,
    outer_log_lik_change = outer_log_lik_change,
    verbose = verbose
  )

  outer_negative_streak <- streak_fn(
    outer_log_lik_change = outer_log_lik_change,
    outer_negative_streak = outer_negative_streak,
    max_negative_outer_streak = max_negative_outer_streak
  )

  convergence_fn(
    outer_log_lik_change = outer_log_lik_change,
    outer_stop_crit = outer_stop_crit,
    verbose = verbose
  )

  list(
    step_size = step_size,
    outer_only_run_counter = outer_only_run_counter,
    outer_end_log_lik = outer_end_log_lik,
    outer_log_lik_change = outer_log_lik_change,
    outer_negative_streak = outer_negative_streak
  )
}