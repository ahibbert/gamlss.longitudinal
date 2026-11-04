#' Update CG smoothing parameters and lambda trace
#'
#' @noRd
.gl_update_cg_lambdas <- function(
    lambda_current,
    lambda_trace,
    H_obs_current,
    beta_vec,
    grad_vec,
    mm_cg,
    trust_radius,
    outer_iteration,
    update_lambda,
    max_delta,
    build_penalty_fn,
    eval_fn,
    edf_fn,
    objective_fn,
    lambda_penalty_K,
    verbose) {
  lambda_new <- lambda_current

  if (!isTRUE(update_lambda)) {
    return(list(lambda = lambda_new, lambda_trace = lambda_trace))
  }

  for (pn in names(lambda_new)) {
    if (length(lambda_new[[pn]]) == 0) next

    for (sn in names(lambda_new[[pn]])) {
      lambda0 <- as.numeric(lambda_new[[pn]][[sn]])
      lambda0 <- .gl_cg_lambda_base(lambda0)
      candidates <- .gl_cg_lambda_candidates(lambda0)

      lambda_scores <- .gl_score_cg_lambda_candidates(
        candidates = candidates,
        parameter = pn,
        smooth = sn,
        lambda_current = lambda_new,
        H_obs_current = H_obs_current,
        beta_vec = beta_vec,
        grad_vec = grad_vec,
        mm_cg = mm_cg,
        trust_radius = trust_radius,
        max_delta = max_delta,
        build_penalty_fn = build_penalty_fn,
        eval_fn = eval_fn,
        edf_fn = edf_fn,
        objective_fn = objective_fn,
        lambda_penalty_K = lambda_penalty_K
      )

      best <- lambda_scores$best

      if (length(best) == 1 && is.finite(lambda_scores$penalized_loglik[best])) {
        trace_rows <- .gl_build_cg_lambda_trace_rows(
          outer_iteration = outer_iteration,
          parameter = pn,
          smooth = sn,
          lambda_before = lambda0,
          candidates = candidates,
          lambda_scores = lambda_scores,
          best = best
        )

        lambda_trace <- rbind(lambda_trace, trace_rows)
        lambda_new[[pn]][[sn]] <- candidates[best]

        if (verbose > 1) {
          cat(paste0(
            "\nCG lambda update for ", pn, " - ", sn, ": ",
            signif(lambda0, 4), " -> ", signif(candidates[best], 4)
          ))
        }
      }
    }
  }

  list(lambda = lambda_new, lambda_trace = lambda_trace)
}

#' Update CG smoothing parameters on the scheduled iterations
#'
#' @noRd
.gl_maybe_update_cg_lambdas_on_schedule <- function(
    lambda_current,
    lambda_trace,
    penalty_current,
    lambda_update_count,
    update_lambda,
    max_lambda_updates,
    lambda_update_every,
    outer_iteration,
    H_obs_current,
    beta_vec,
    grad_vec,
    mm_cg,
    trust_radius,
    max_delta,
    step_tol,
    build_penalty_fn,
    eval_fn,
    edf_fn,
    objective_fn,
    lambda_penalty_K,
    verbose,
    update_lambdas_fn = .gl_update_cg_lambdas,
    shrink_trust_radius_fn = .gl_shrink_cg_trust_radius,
    lambdas_changed_fn = .gl_cg_lambdas_changed) {
  lambda_changed <- FALSE

  should_update <- isTRUE(update_lambda) &&
    outer_iteration > 1L &&
    lambda_update_count < max_lambda_updates &&
    (outer_iteration %% lambda_update_every == 0L)

  if (!isTRUE(should_update)) {
    return(list(
      lambda = lambda_current,
      lambda_trace = lambda_trace,
      penalty_mat = penalty_current,
      trust_radius = trust_radius,
      lambda_update_count = lambda_update_count,
      lambda_changed = lambda_changed
    ))
  }

  lambda_before <- lambda_current

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
  penalty_current <- build_penalty_fn(names(beta_vec), lambda_current)

  lambda_changed <- lambdas_changed_fn(
    lambda_before = lambda_before,
    lambda_after = lambda_current
  )

  if (isTRUE(lambda_changed)) {
    trust_radius <- shrink_trust_radius_fn(
      trust_radius = trust_radius,
      step_tol = step_tol
    )

    if (verbose > 0) {
      cat(paste0("\nCG trust radius shrunk after lambda update to ", signif(trust_radius, 4)))
    }
  }

  list(
    lambda = lambda_current,
    lambda_trace = lambda_trace,
    penalty_mat = penalty_current,
    trust_radius = trust_radius,
    lambda_update_count = lambda_update_count + 1L,
    lambda_changed = lambda_changed
  )
}
