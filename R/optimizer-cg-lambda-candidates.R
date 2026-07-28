#' Candidate smoothing parameters for CG lambda update
#'
#' @noRd
.gl_cg_lambda_base <- function(lambda0) {
  lambda0 <- as.numeric(lambda0)
  if (!is.finite(lambda0) || lambda0 <= 0) lambda0 <- 1
  lambda0
}

#' Candidate smoothing parameters for CG lambda update
#'
#' @noRd
.gl_cg_lambda_candidates <- function(lambda0) {
  lambda0 <- .gl_cg_lambda_base(lambda0)
  unique(pmax(0.01, pmin(1e6, lambda0 * c(0.1, 0.25, 0.5, 1, 2, 4, 10))))
}

#' Score candidate smoothing parameters for one CG smooth
#'
#' @noRd
.gl_score_cg_lambda_candidates <- function(
    candidates,
    parameter,
    smooth,
    lambda_current,
    H_obs_current,
    beta_vec,
    grad_vec,
    mm_cg,
    trust_radius,
    max_delta,
    build_penalty_fn,
    eval_fn,
    edf_fn,
    objective_fn,
    lambda_penalty_K) {
  gaic_score <- rep(Inf, length(candidates))
  penalty_value <- rep(NA_real_, length(candidates))
  penalized_loglik <- rep(NA_real_, length(candidates))
  raw_loglik <- rep(NA_real_, length(candidates))
  edf_values <- rep(NA_real_, length(candidates))

  for (jj in seq_along(candidates)) {
    lambda_try <- lambda_current
    lambda_try[[parameter]][[smooth]] <- candidates[jj]
    P_try <- build_penalty_fn(names(beta_vec), lambda_try)
    g_try <- grad_vec - as.numeric(P_try %*% beta_vec)
    H_try <- H_obs_current - P_try
    delta <- tryCatch(-as.numeric(.solve_linear_system(H_try, g_try)), error = function(e) NULL)
    if (is.null(delta) || !all(is.finite(delta))) next
    delta <- .gl_limit_cg_step(delta, trust_radius = trust_radius, max_delta = max_delta)
    beta_try <- beta_vec + delta
    eval_try <- eval_fn(beta_try, mm_cg)
    if (is.null(eval_try) || !is.finite(eval_try$loglik)) next
    edf_try <- sum(unlist(edf_fn(H_obs_current, P_try, names(beta_vec))), na.rm = TRUE)
    penalty_try <- sum(as.numeric(beta_try) * as.numeric(P_try %*% beta_try))
    raw_loglik[jj] <- eval_try$loglik
    edf_values[jj] <- edf_try
    penalty_value[jj] <- penalty_try
    penalized_loglik[jj] <- objective_fn(beta_try, eval_try$loglik, P_try)
    gaic_score[jj] <- -2 * eval_try$loglik + lambda_penalty_K * edf_try
  }

  list(
    best = which.max(penalized_loglik),
    raw_loglik = raw_loglik,
    edf_values = edf_values,
    penalty_value = penalty_value,
    penalized_loglik = penalized_loglik,
    gaic_score = gaic_score
  )
}
