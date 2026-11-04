#' Initialize CG optimizer state after runtime helpers are built
#'
#' @noRd
.gl_initialize_cg_optimizer_state <- function(
    mm,
    par_cov,
    par_s,
    lambda_s,
    cg_max_delta,
    build_model_fn,
    build_penalty_fn) {
  cg_aug <- build_model_fn(mm, par_cov, par_s)
  beta_all <- cg_aug$beta

  list(
    mm_cg = cg_aug$mm,
    beta_all = beta_all,
    penalty_mat = build_penalty_fn(names(beta_all), lambda_s),
    cg_trust_radius = as.numeric(cg_max_delta),
    cg_stall_count = 0L,
    cg_converged = FALSE,
    cg_lambda_update_count = 0L,
    cg_has_smooths = length(unlist(lambda_s, use.names = FALSE)) > 0L,
    cg_lambda_trace = data.frame(),
    cg_step_trace = list()
  )
}

#' Evaluate and initialize one CG outer iteration
#'
#' @noRd
.gl_evaluate_cg_iteration_start <- function(
    beta_vec,
    mm_cg,
    penalty_current,
    log_lik_history,
    par_history,
    best_raw_loglik,
    best_iteration,
    current_iteration,
    margin_dist,
    copula_dist,
    include_dlcopdpar,
    dataset,
    cg_gradient_method,
    eval_fn,
    objective_fn,
    finite_gradient_fn,
    history_fn = .gl_append_optimizer_history,
    best_fn = .gl_update_cg_best_loglik,
    gradient_fn = .gl_compute_cg_gradient) {
  eval_start <- eval_fn(beta_vec, mm_cg)

  if (is.null(eval_start) || !is.finite(eval_start$loglik)) {
    stop("CG failed: current likelihood is not finite.")
  }

  history_state <- history_fn(
    log_lik_history = log_lik_history,
    par_history = par_history,
    calc_lik_out = eval_start$calc_lik,
    par_cov = eval_start$par_cov
  )

  outer_start_log_lik <- eval_start$loglik

  best_state <- best_fn(
    candidate_loglik = outer_start_log_lik,
    best_raw_loglik = best_raw_loglik,
    best_iteration = best_iteration,
    current_iteration = current_iteration
  )

  obj_start <- objective_fn(beta_vec, outer_start_log_lik, penalty_current)

  grad <- gradient_fn(
    beta_vec = beta_vec,
    mm_cg = mm_cg,
    eval_start = eval_start,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    include_dlcopdpar = include_dlcopdpar,
    dataset = dataset,
    cg_gradient_method = cg_gradient_method,
    finite_gradient_fn = finite_gradient_fn
  )

  list(
    eval_start = eval_start,
    log_lik_history = history_state$log_lik_history,
    par_history = history_state$par_history,
    outer_start_log_lik = outer_start_log_lik,
    best_raw_loglik = best_state$best_raw_loglik,
    best_iteration = best_state$best_iteration,
    obj_start = obj_start,
    grad = grad
  )
}
