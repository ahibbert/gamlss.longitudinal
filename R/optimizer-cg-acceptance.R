#' Assess whether to accept a CG line-search step
#'
#' @noRd
.gl_assess_cg_step_acceptance <- function(
    best,
    best_raw_loglik,
    outer_start_loglik,
    raw_loglik_drop_tol,
    lambda_update_count) {
  if (is.null(best)) {
    return(list(
      has_step = FALSE,
      accept = FALSE,
      prevented_deterioration = FALSE,
      prevented_raw_loglik_drop = NA_real_
    ))
  }

  prospective_best_raw_loglik <- max(best_raw_loglik, outer_start_loglik, na.rm = TRUE)
  prospective_raw_loglik_drop <- prospective_best_raw_loglik - best$eval$loglik
  prevented_deterioration <- is.finite(raw_loglik_drop_tol) &&
    lambda_update_count > 0L &&
    is.finite(prospective_raw_loglik_drop) &&
    prospective_raw_loglik_drop >= raw_loglik_drop_tol
  prevented_deterioration <- isTRUE(prevented_deterioration)

  list(
    has_step = TRUE,
    accept = !prevented_deterioration,
    prevented_deterioration = prevented_deterioration,
    prevented_raw_loglik_drop = if (prevented_deterioration) prospective_raw_loglik_drop else NA_real_
  )
}

#' Apply a CG line-search acceptance decision
#'
#' @noRd
.gl_apply_cg_step_acceptance <- function(
    step_acceptance,
    best,
    eval_start,
    beta_vec,
    par_cov_template,
    par_s_template,
    stall_count,
    trust_radius,
    step_tol,
    max_delta,
    max_stall,
    verbose) {
  accepted_improvement <- NA_real_
  accepted_best <- best
  calc_lik_out_end <- eval_start$calc_lik
  par_cov <- par_cov_template
  par_s <- par_s_template

  if (!isTRUE(step_acceptance$has_step)) {
    stall_count <- stall_count + 1L
    trust_radius <- .gl_shrink_cg_trust_radius(
      trust_radius = trust_radius,
      step_tol = step_tol
    )
    accepted_best <- NULL

    if (verbose > 0) {
      cat(paste0("\nCG step rejected (stall ", stall_count, "/", max_stall, ")\n"))
    }
  } else if (!isTRUE(step_acceptance$accept)) {
    accepted_best <- NULL
  } else {
    beta_vec <- best$beta
    unpacked <- .gl_unpack_cg_beta(
      beta_vec,
      par_cov_template = par_cov_template,
      par_s_template = par_s_template
    )
    par_cov <- unpacked$par_cov
    par_s <- unpacked$par_s
    calc_lik_out_end <- best$eval$calc_lik
    stall_count <- 0L
    accepted_improvement <- best$improvement
    trust_radius <- .gl_expand_cg_trust_radius(
      trust_radius = trust_radius,
      step_l2 = best$step_l2,
      step_tol = step_tol,
      max_delta = max_delta
    )
  }

  list(
    beta = beta_vec,
    par_cov = par_cov,
    par_s = par_s,
    calc_lik_out_end = calc_lik_out_end,
    stall_count = stall_count,
    trust_radius = trust_radius,
    accepted_improvement = accepted_improvement,
    best = accepted_best,
    prevented_deterioration = step_acceptance$prevented_deterioration,
    prevented_raw_loglik_drop = step_acceptance$prevented_raw_loglik_drop
  )
}
