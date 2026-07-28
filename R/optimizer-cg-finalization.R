#' Refresh final CG smooth effective degrees of freedom
#'
#' @noRd
.gl_refresh_final_cg_smooth_edf <- function(
    dataset,
    margin_dist,
    copula_dist,
    mm_cg,
    beta_vec,
    lambda_current,
    penalty_current,
    df_s_current,
    observed_hessian_fn,
    build_penalty_fn,
    edf_fn) {
  final_obj <- .gl_build_cg_hessian_object(
    dataset = dataset,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    mm_cg = mm_cg,
    beta_vec = beta_vec
  )

  final_H <- tryCatch(
    observed_hessian_fn(final_obj, beta_vec, mm_cg, context = "final smooth EDF update"),
    error = function(e) NULL
  )

  if (!is.null(final_H)) {
    penalty_current <- build_penalty_fn(names(beta_vec), lambda_current)
    df_s_current <- edf_fn(final_H, penalty_current, names(beta_vec))
  }

  list(
    final_H = final_H,
    penalty_mat = penalty_current,
    df_s = df_s_current
  )
}

#' Finalize CG optimizer state before fitted-object assembly
#'
#' @noRd
.gl_finalize_cg_optimizer_state <- function(
    dataset,
    margin_dist,
    copula_dist,
    mm,
    mm_cg,
    beta_vec,
    lambda_current,
    penalty_current,
    df_s_current,
    observed_hessian_fn,
    build_penalty_fn,
    edf_fn,
    final_edf_fn = .gl_refresh_final_cg_smooth_edf,
    weights_fn = .gl_initialize_unit_weights) {
  final_edf <- final_edf_fn(
    dataset = dataset,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    mm_cg = mm_cg,
    beta_vec = beta_vec,
    lambda_current = lambda_current,
    penalty_current = penalty_current,
    df_s_current = df_s_current,
    observed_hessian_fn = observed_hessian_fn,
    build_penalty_fn = build_penalty_fn,
    edf_fn = edf_fn
  )

  list(
    final_edf = final_edf,
    penalty_mat = final_edf$penalty_mat,
    df_s = final_edf$df_s,
    weights_final = weights_fn(mm)
  )
}

#' Optionally replace CG zeta Hessian block with finite differences
#'
#' @noRd
.gl_apply_cg_zeta_hessian_override <- function(
    H_obs,
    beta_vec,
    mm_cg,
    zeta_hessian,
    finite_hessian_fn,
    verbose) {
  H_zeta_fd <- NULL

  if (identical(zeta_hessian, "finite")) {
    zeta_names <- grep("^zeta\\.", names(beta_vec), value = TRUE)

    if (length(zeta_names) > 0L) {
      H_zeta_fd <- finite_hessian_fn(beta_vec, zeta_names, mm_cg)

      if (all(is.finite(H_zeta_fd))) {
        H_obs[zeta_names, zeta_names] <- 0.5 * (H_zeta_fd + t(H_zeta_fd))
      } else if (verbose > 0) {
        cat("\nCG finite zeta Hessian skipped because the block was not finite.")
      }
    }
  }

  list(H_obs = H_obs, H_zeta_fd = H_zeta_fd)
}

#' Build one CG step trace row
#'
#' @noRd
.gl_build_cg_step_trace_row <- function(
    outer_iteration,
    start_loglik,
    end_loglik,
    raw_loglik_change,
    start_penalized_loglik,
    accepted_penalized_improvement,
    grad_inf,
    step_l2,
    trust_radius_start,
    trust_radius_end,
    line_search_evals,
    accepted_step,
    lambda_update_count,
    lambda_changed,
    stall_count,
    tolerance_met,
    max_stall_hit,
    raw_deterioration_hit,
    raw_loglik_drop_from_best) {
  data.frame(
    outer_iteration = as.integer(outer_iteration),
    start_logLik = as.numeric(start_loglik),
    end_logLik = as.numeric(end_loglik),
    raw_logLik_change = as.numeric(raw_loglik_change),
    start_penalized_logLik = as.numeric(start_penalized_loglik),
    accepted_penalized_improvement = as.numeric(accepted_penalized_improvement),
    grad_inf = as.numeric(grad_inf),
    step_l2 = as.numeric(step_l2),
    trust_radius_start = as.numeric(trust_radius_start),
    trust_radius_end = as.numeric(trust_radius_end),
    line_search_evals = as.integer(line_search_evals),
    accepted_step = isTRUE(accepted_step),
    lambda_update_count = as.integer(lambda_update_count),
    lambda_changed = isTRUE(lambda_changed),
    stall_count = as.integer(stall_count),
    tolerance_met = isTRUE(tolerance_met),
    max_stall_hit = isTRUE(max_stall_hit),
    raw_deterioration_hit = isTRUE(raw_deterioration_hit),
    raw_loglik_drop_from_best = as.numeric(raw_loglik_drop_from_best),
    row.names = NULL
  )
}
