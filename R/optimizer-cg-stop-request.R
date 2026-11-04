.gl_handle_cg_stop_request <- function(
    update_lambda,
    has_smooths,
    lambda_update_count,
    lambda_current,
    lambda_trace,
    H_obs_current,
    beta_vec,
    grad_vec,
    mm_cg,
    trust_radius,
    outer_iteration,
    max_delta,
    build_penalty_fn,
    eval_fn,
    edf_fn,
    objective_fn,
    lambda_penalty_K,
    tolerance_met,
    deterioration_hit,
    raw_loglik_drop_from_best,
    verbose,
    update_lambdas_fn = .gl_update_cg_lambdas,
    delay_fn = .gl_should_delay_cg_convergence_for_lambda_update,
    stop_reason_fn = .gl_cg_stop_reason) {
  if (delay_fn(
    update_lambda = update_lambda,
    has_smooths = has_smooths,
    lambda_update_count = lambda_update_count
  )) {
    lambda_update <- update_lambdas_fn(
      lambda_current = lambda_current,
      lambda_trace = lambda_trace,
      H_obs_current = H_obs_current,
      beta_vec = beta_vec,
      grad_vec = grad_vec,
      mm_cg = mm_cg,
      trust_radius = trust_radius,
      outer_iteration = outer_iteration,
      update_lambda = update_lambda,
      max_delta = max_delta,
      build_penalty_fn = build_penalty_fn,
      eval_fn = eval_fn,
      edf_fn = edf_fn,
      objective_fn = objective_fn,
      lambda_penalty_K = lambda_penalty_K,
      verbose = verbose
    )
    lambda_current <- lambda_update$lambda
    lambda_trace <- lambda_update$lambda_trace
    penalty_mat <- build_penalty_fn(names(beta_vec), lambda_current)
    df_s <- edf_fn(H_obs_current, penalty_mat, names(beta_vec))

    if (verbose > 0) {
      cat("\nCG convergence delayed for first smoother lambda update")
    }

    return(list(
      lambda = lambda_current,
      lambda_trace = lambda_trace,
      penalty_mat = penalty_mat,
      df_s = df_s,
      lambda_update_count = lambda_update_count + 1L,
      stall_count_reset = TRUE,
      stop_reason = NA_character_,
      converged = FALSE
    ))
  }

  stop_reason <- stop_reason_fn(
    tolerance_met = tolerance_met,
    deterioration_hit = deterioration_hit
  )

  if (identical(stop_reason, "tolerance") && verbose > 0) {
    cat("\nOUTER CONVERGED")
  }

  if (identical(stop_reason, "raw_loglik_deterioration") && verbose > 0) {
    cat(paste0(
      "\nCG stopped after raw log-likelihood dropped ",
      signif(raw_loglik_drop_from_best, 5),
      " below best seen value."
    ))
  }

  list(
    lambda = lambda_current,
    lambda_trace = lambda_trace,
    penalty_mat = NULL,
    df_s = NULL,
    lambda_update_count = lambda_update_count,
    stall_count_reset = FALSE,
    stop_reason = stop_reason,
    converged = TRUE
  )
}

#' Apply a requested CG stop to optimizer state
#'
#' @noRd
.gl_apply_cg_stop_request_state <- function(
    update_lambda,
    has_smooths,
    lambda_update_count,
    lambda_current,
    lambda_trace,
    penalty_current,
    df_s_current,
    stall_count,
    H_obs_current,
    beta_vec,
    grad_vec,
    mm_cg,
    trust_radius,
    outer_iteration,
    max_delta,
    build_penalty_fn,
    eval_fn,
    edf_fn,
    objective_fn,
    lambda_penalty_K,
    tolerance_met,
    deterioration_hit,
    raw_loglik_drop_from_best,
    verbose,
    stop_request_fn = .gl_handle_cg_stop_request) {
  stop_state <- stop_request_fn(
    update_lambda = update_lambda,
    has_smooths = has_smooths,
    lambda_update_count = lambda_update_count,
    lambda_current = lambda_current,
    lambda_trace = lambda_trace,
    H_obs_current = H_obs_current,
    beta_vec = beta_vec,
    grad_vec = grad_vec,
    mm_cg = mm_cg,
    trust_radius = trust_radius,
    outer_iteration = outer_iteration,
    max_delta = max_delta,
    build_penalty_fn = build_penalty_fn,
    eval_fn = eval_fn,
    edf_fn = edf_fn,
    objective_fn = objective_fn,
    lambda_penalty_K = lambda_penalty_K,
    tolerance_met = tolerance_met,
    deterioration_hit = deterioration_hit,
    raw_loglik_drop_from_best = raw_loglik_drop_from_best,
    verbose = verbose
  )

  list(
    stop_state = stop_state,
    lambda = stop_state$lambda,
    lambda_trace = stop_state$lambda_trace,
    lambda_update_count = stop_state$lambda_update_count,
    penalty_mat = if (!is.null(stop_state$penalty_mat)) stop_state$penalty_mat else penalty_current,
    df_s = if (!is.null(stop_state$df_s)) stop_state$df_s else df_s_current,
    stall_count = if (isTRUE(stop_state$stall_count_reset)) 0L else stall_count,
    stop_reason = stop_state$stop_reason,
    converged = stop_state$converged
  )
}

#' Apply CG stop state only when the iteration requested a stop
#'
#' @noRd
.gl_maybe_apply_cg_stop_request_state <- function(
    stop_requested,
    cg_update_lambda,
    cg_has_smooths,
    cg_lambda_update_count,
    lambda_s,
    cg_lambda_trace,
    penalty_mat,
    df_s,
    cg_stall_count,
    H_obs,
    beta_all,
    grad,
    mm_cg,
    cg_trust_radius,
    outer_only_run_counter,
    cg_max_delta,
    build_cg_penalty,
    cg_eval,
    cg_smooth_edf_list,
    cg_objective,
    lambda_penalty_K,
    cg_tolerance_met,
    cg_deterioration_hit,
    cg_raw_loglik_drop_from_best,
    verbose,
    stop_request_state_fn = .gl_apply_cg_stop_request_state) {
  if (!isTRUE(stop_requested)) {
    return(list(
      lambda_s = lambda_s,
      cg_lambda_trace = cg_lambda_trace,
      cg_lambda_update_count = cg_lambda_update_count,
      penalty_mat = penalty_mat,
      df_s = df_s,
      cg_stall_count = cg_stall_count,
      cg_stop_reason = NA_character_,
      cg_converged = FALSE
    ))
  }

  cg_stop_state <- stop_request_state_fn(
    update_lambda = cg_update_lambda,
    has_smooths = cg_has_smooths,
    lambda_update_count = cg_lambda_update_count,
    lambda_current = lambda_s,
    lambda_trace = cg_lambda_trace,
    penalty_current = penalty_mat,
    df_s_current = df_s,
    stall_count = cg_stall_count,
    H_obs_current = H_obs,
    beta_vec = beta_all,
    grad_vec = grad,
    mm_cg = mm_cg,
    trust_radius = cg_trust_radius,
    outer_iteration = outer_only_run_counter,
    max_delta = cg_max_delta,
    build_penalty_fn = build_cg_penalty,
    eval_fn = cg_eval,
    edf_fn = cg_smooth_edf_list,
    objective_fn = cg_objective,
    lambda_penalty_K = lambda_penalty_K,
    tolerance_met = cg_tolerance_met,
    deterioration_hit = cg_deterioration_hit,
    raw_loglik_drop_from_best = cg_raw_loglik_drop_from_best,
    verbose = verbose
  )

  list(
    lambda_s = cg_stop_state$lambda,
    cg_lambda_trace = cg_stop_state$lambda_trace,
    cg_lambda_update_count = cg_stop_state$lambda_update_count,
    penalty_mat = cg_stop_state$penalty_mat,
    df_s = cg_stop_state$df_s,
    cg_stall_count = cg_stop_state$stall_count,
    cg_stop_reason = cg_stop_state$stop_reason,
    cg_converged = cg_stop_state$converged
  )
}
