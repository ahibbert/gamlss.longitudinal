#' Put together convergence information for a fitted model
#'
#' @noRd
.gl_build_convergence_info <- function(
    method,
    outer_log_lik_change,
    outer_stop_crit,
    outer_only_run_counter,
    max_outer_iter,
    cg_stop_reason,
    cg_last_grad_inf,
    cg_last_step_l2,
    cg_best_raw_loglik,
    cg_best_iteration,
    cg_raw_loglik_drop_from_best,
    cg_raw_loglik_drop_tol,
    cg_gradient_method,
    cg_zeta_hessian,
    cg_hessian_method,
    objective = NULL,
    elapsed_sec = NA_real_,
    optimizer_control_requested = NULL,
    optimizer_control_effective = NULL,
    cg_grad_tol = NA_real_,
    cg_step_tol = NA_real_) {
  objective_supplied <- !missing(objective) && !is.null(objective)
  objective_ok <- !objective_supplied ||
    (length(objective) == 1L && is.finite(objective))
  objective_contract <- is.finite(outer_log_lik_change) &&
    is.finite(outer_stop_crit) && abs(outer_log_lik_change) <= outer_stop_crit
  converged <- objective_ok && objective_contract

  if (identical(method, "CG")) {
    gradient_contract <- is.finite(cg_last_grad_inf) && is.finite(cg_grad_tol) && cg_last_grad_inf <= cg_grad_tol
    step_contract <- is.finite(cg_last_step_l2) && is.finite(cg_step_tol) && cg_last_step_l2 <= cg_step_tol
    converged <- objective_ok && objective_contract && gradient_contract && step_contract &&
      !cg_stop_reason %in% c("max_stall", "raw_loglik_deterioration")
  } else {
    gradient_contract <- NA
    step_contract <- NA
  }

  hit_outer_limit <- outer_only_run_counter > max_outer_iter && !isTRUE(converged)
  safeguard_reason <- if (identical(cg_stop_reason, "max_stall")) {
    "max_stall"
  } else if (identical(cg_stop_reason, "raw_loglik_deterioration")) {
    "objective_deterioration"
  } else {
    NA_character_
  }
  termination_reason <- if (!objective_ok) {
    "invalid_likelihood"
  } else if (!is.na(safeguard_reason)) {
    safeguard_reason
  } else if (isTRUE(converged)) {
    "converged"
  } else if (outer_only_run_counter > max_outer_iter) {
    "max_iterations"
  } else {
    "numerical_failure"
  }
  events <- data.frame(
    level = character(), code = character(), message = character(),
    iteration = integer(), stringsAsFactors = FALSE
  )
  if (!isTRUE(converged)) {
    events <- rbind(events, data.frame(
      level = "warning",
      code = termination_reason,
      message = paste0("Optimizer returned without satisfying the ", method, " convergence contract."),
      iteration = max(0L, outer_only_run_counter - 1L),
      stringsAsFactors = FALSE
    ))
  }

  list(
    converged = isTRUE(converged),
    hit_outer_limit = isTRUE(hit_outer_limit),
    hit_max_stall = isTRUE(identical(cg_stop_reason, "max_stall")),
    hit_raw_loglik_deterioration = isTRUE(identical(cg_stop_reason, "raw_loglik_deterioration")),
    stop_reason = termination_reason,
    termination_reason = termination_reason,
    objective = if (objective_supplied) as.numeric(objective)[1L] else NA_real_,
    logLik = if (objective_supplied) as.numeric(objective)[1L] else NA_real_,
    objective_change = as.numeric(outer_log_lik_change),
    objective_tolerance = as.numeric(outer_stop_crit),
    objective_criterion_met = isTRUE(objective_contract),
    gradient_norm = if (identical(method, "CG")) as.numeric(cg_last_grad_inf) else NA_real_,
    gradient_tolerance = if (identical(method, "CG")) as.numeric(cg_grad_tol) else NA_real_,
    gradient_criterion_met = if (identical(method, "CG")) isTRUE(gradient_contract) else NA,
    step_norm = if (identical(method, "CG")) as.numeric(cg_last_step_l2) else NA_real_,
    step_tolerance = if (identical(method, "CG")) as.numeric(cg_step_tol) else NA_real_,
    step_criterion_met = if (identical(method, "CG")) isTRUE(step_contract) else NA,
    grad_inf = as.numeric(cg_last_grad_inf),
    step_l2 = as.numeric(cg_last_step_l2),
    best_raw_loglik = as.numeric(cg_best_raw_loglik),
    best_raw_loglik_iteration = as.integer(cg_best_iteration),
    raw_loglik_drop_from_best = as.numeric(cg_raw_loglik_drop_from_best),
    raw_loglik_drop_tol = as.numeric(cg_raw_loglik_drop_tol),
    outer_iterations = max(0L, outer_only_run_counter - 1L),
    max_outer_iter = max_outer_iter,
    outer_log_lik_change = as.numeric(outer_log_lik_change),
    outer_stop_crit = outer_stop_crit,
    elapsed_sec = as.numeric(elapsed_sec),
    requested_controls = optimizer_control_requested,
    effective_controls = optimizer_control_effective,
    events = events,
    warnings = events[events$level == "warning", , drop = FALSE],
    method = method,
    cg_gradient_method = if (identical(method, "CG")) cg_gradient_method else NA_character_,
    cg_zeta_hessian = if (identical(method, "CG")) cg_zeta_hessian else NA_character_,
    cg_hessian_method = if (identical(method, "CG")) cg_hessian_method else NA_character_
  )
}

#' Extract information criteria from fitted model to pass to the final fitted object
#'
#' @noRd
.gl_fit_information_criteria <- function(par_cov, par_s, df_s, calc_lik_out_end, dataset) {
  p_cop <- par_cov[grepl("theta", names(par_cov)) | grepl("zeta", names(par_cov))]
  p_mar <- par_cov[!(grepl("theta", names(par_cov)) | grepl("zeta", names(par_cov)))]

  df_s_total <- df_s_cop_total <- df_s_margin_total <- 0

  for (par_name in names(par_s)) {
    if (par_name %in% c("theta", "zeta")) {
      df_s_cop_total <- df_s_cop_total + sum(unlist(df_s[[par_name]]))
    } else {
      df_s_margin_total <- df_s_margin_total + sum(unlist(df_s[[par_name]]))
    }

    df_s_total <- df_s_total + sum(unlist(df_s[[par_name]]))
  }

  n_observed <- sum(is.finite(dataset$response))
  if (n_observed < 1L) {
    stop("Information criteria require at least one finite observed response.", call. = FALSE)
  }

  aics <- rbind(
    t(calc_lik_out_end$log_lik),
    t(-calc_lik_out_end$log_lik * 2) + 2 * c(
      length(p_mar) + df_s_margin_total,
      length(p_cop) + df_s_cop_total,
      length(par_cov) + df_s_total
    ),
    t(-calc_lik_out_end$log_lik * 2) + c(
      length(p_mar) + df_s_margin_total,
      length(p_cop) + df_s_cop_total,
      length(par_cov) + df_s_total
    ) * log(n_observed),
    t(c(
      length(p_mar) + df_s_margin_total,
      length(p_cop) + df_s_cop_total,
      length(par_cov) + df_s_total
    ))
  )

  rownames(aics) <- c("LogLik", "AIC", "BIC", "EDF")
  aics
}
