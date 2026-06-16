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
    cg_hessian_method) {
  converged <- is.finite(outer_log_lik_change) && abs(outer_log_lik_change) <= outer_stop_crit

  if (identical(method, "CG")) {
    converged <- identical(cg_stop_reason, "tolerance")
  }

  hit_outer_limit <- outer_only_run_counter >= max_outer_iter && !isTRUE(converged)

  list(
    converged = isTRUE(converged),
    hit_outer_limit = isTRUE(hit_outer_limit),
    hit_max_stall = isTRUE(identical(cg_stop_reason, "max_stall")),
    hit_raw_loglik_deterioration = isTRUE(identical(cg_stop_reason, "raw_loglik_deterioration")),
    stop_reason = if (identical(method, "CG")) cg_stop_reason else if (isTRUE(converged)) "tolerance" else NA_character_,
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
    ) * log(nrow(dataset)),
    t(c(
      length(p_mar) + df_s_margin_total,
      length(p_cop) + df_s_cop_total,
      length(par_cov) + df_s_total
    ))
  )

  rownames(aics) <- c("LogLik", "AIC", "BIC", "EDF")
  aics
}
