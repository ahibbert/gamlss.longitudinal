#' Evaluate one CG line-search candidate
#'
#' @noRd
.gl_evaluate_cg_line_search_candidate <- function(
    delta0,
    backtrack_index,
    beta_vec,
    mm_cg,
    penalty_current,
    obj_start,
    trust_radius,
    max_delta,
    eval_fn,
    objective_fn) {
  delta <- delta0 / (2^(backtrack_index - 1L))
  delta <- .gl_limit_cg_step(delta, trust_radius = trust_radius, max_delta = max_delta)
  beta_try <- beta_vec + delta
  eval_try <- eval_fn(beta_try, mm_cg)

  if (is.null(eval_try) || !is.finite(eval_try$loglik)) {
    return(NULL)
  }

  obj_try <- objective_fn(beta_try, eval_try$loglik, penalty_current)
  improvement <- obj_try - obj_start

  list(
    beta = beta_try,
    eval = eval_try,
    improvement = improvement,
    step_l2 = sqrt(sum(delta^2))
  )
}
