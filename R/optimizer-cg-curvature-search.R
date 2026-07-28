#' Prepare CG curvature, lambda updates, and line-search state
#'
#' @noRd
.gl_prepare_cg_curvature_line_search_state <- function(
    dataset,
    margin_dist,
    copula_dist,
    mm_cg,
    beta_vec,
    grad_vec,
    lambda_current,
    lambda_trace,
    penalty_current,
    lambda_update_count,
    update_lambda,
    max_lambda_updates,
    lambda_update_every,
    outer_iteration,
    trust_radius,
    max_delta,
    step_tol,
    build_penalty_fn,
    eval_fn,
    edf_fn,
    objective_fn,
    lambda_penalty_K,
    cg_zeta_hessian,
    finite_hessian_fn,
    obj_start,
    armijo_c1,
    line_search,
    max_line_search_evals,
    use_backtracking,
    backtracking_max_halves,
    verbose,
    observed_hessian_fn,
    hessian_object_fn = .gl_build_cg_hessian_object,
    zeta_hessian_fn = .gl_apply_cg_zeta_hessian_override,
    lambda_schedule_fn = .gl_maybe_update_cg_lambdas_on_schedule,
    candidate_steps_fn = .gl_build_cg_candidate_steps,
    line_search_fn = .gl_run_cg_line_search) {
  tmp_obj <- hessian_object_fn(
    dataset = dataset,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    mm_cg = mm_cg,
    beta_vec = beta_vec
  )

  H_obs <- observed_hessian_fn(
    tmp_obj,
    beta_vec,
    mm_cg,
    context = paste0("outer iteration ", outer_iteration)
  )

  lambda_changed <- FALSE

  zeta_hessian_state <- zeta_hessian_fn(
    H_obs = H_obs,
    beta_vec = beta_vec,
    mm_cg = mm_cg,
    zeta_hessian = cg_zeta_hessian,
    finite_hessian_fn = finite_hessian_fn,
    verbose = verbose
  )
  H_obs <- zeta_hessian_state$H_obs
  H_zeta_fd <- zeta_hessian_state$H_zeta_fd

  lambda_schedule <- lambda_schedule_fn(
    lambda_current = lambda_current,
    lambda_trace = lambda_trace,
    penalty_current = penalty_current,
    lambda_update_count = lambda_update_count,
    update_lambda = update_lambda,
    max_lambda_updates = max_lambda_updates,
    lambda_update_every = lambda_update_every,
    outer_iteration = outer_iteration,
    H_obs_current = H_obs,
    beta_vec = beta_vec,
    grad_vec = grad_vec,
    mm_cg = mm_cg,
    trust_radius = trust_radius,
    max_delta = max_delta,
    step_tol = step_tol,
    build_penalty_fn = build_penalty_fn,
    eval_fn = eval_fn,
    edf_fn = edf_fn,
    objective_fn = objective_fn,
    lambda_penalty_K = lambda_penalty_K,
    verbose = verbose
  )
  lambda_current <- lambda_schedule$lambda
  lambda_trace <- lambda_schedule$lambda_trace
  penalty_current <- lambda_schedule$penalty_mat
  trust_radius <- lambda_schedule$trust_radius
  lambda_update_count <- lambda_schedule$lambda_update_count
  lambda_changed <- lambda_schedule$lambda_changed

  df_s <- edf_fn(H_obs, penalty_current, names(beta_vec))

  g_pen <- grad_vec - as.numeric(penalty_current %*% beta_vec)
  H_pen <- H_obs - penalty_current

  candidate_steps <- candidate_steps_fn(
    g_pen = g_pen,
    H_pen = H_pen,
    trust_radius = trust_radius
  )

  line_search_out <- line_search_fn(
    candidate_steps = candidate_steps,
    beta_vec = beta_vec,
    mm_cg = mm_cg,
    penalty_current = penalty_current,
    obj_start = obj_start,
    trust_radius = trust_radius,
    max_delta = max_delta,
    armijo_c1 = armijo_c1,
    line_search = line_search,
    max_line_search_evals = max_line_search_evals,
    use_backtracking = use_backtracking,
    backtracking_max_halves = backtracking_max_halves,
    eval_fn = eval_fn,
    objective_fn = objective_fn,
    verbose = verbose
  )

  list(
    tmp_obj = tmp_obj,
    H_obs = H_obs,
    H_zeta_fd = H_zeta_fd,
    lambda_s = lambda_current,
    cg_lambda_trace = lambda_trace,
    penalty_mat = penalty_current,
    cg_trust_radius = trust_radius,
    cg_lambda_update_count = lambda_update_count,
    lambda_changed = lambda_changed,
    df_s = df_s,
    g_pen = g_pen,
    H_pen = H_pen,
    candidate_steps = candidate_steps,
    line_search_out = line_search_out,
    best = line_search_out$best,
    line_eval_count = line_search_out$line_eval_count
  )
}
