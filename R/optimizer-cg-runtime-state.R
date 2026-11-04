#' Prepare CG runtime helpers and initial optimizer state
#'
#' @param par_cov,par_s Current fixed and smooth coefficients.
#' @param mm Model-matrix bundle.
#' @param margin_dist,copula_link,copula_dist Model family settings.
#' @param dataset Prepared fitting data.
#' @param pair_cache,margin_eval_cache Likelihood caches.
#' @param cg_gradient_method,cg_hessian_method CG derivative methods.
#' @param verbose Verbosity level.
#' @param lambda_s Current smooth penalties.
#' @param cg_max_delta Initial CG trust-region radius.
#' @param runtime_helpers_fn Dependency-injected runtime helper factory.
#' @param initialize_state_fn Dependency-injected state initializer.
#' @return List with `runtime` helper functions and initial `state`.
#' @noRd
.gl_prepare_cg_runtime_state <- function(
    par_cov,
    par_s,
    mm,
    margin_dist,
    copula_link,
    copula_dist,
    dataset,
    pair_cache,
    margin_eval_cache,
    cg_gradient_method,
    cg_hessian_method,
    verbose,
    lambda_s,
    cg_max_delta,
    runtime_helpers_fn = .gl_build_cg_runtime_helpers,
    initialize_state_fn = .gl_initialize_cg_optimizer_state) {
  cg_runtime <- runtime_helpers_fn(
    par_cov = par_cov,
    par_s = par_s,
    mm = mm,
    margin_dist = margin_dist,
    copula_link = copula_link,
    copula_dist = copula_dist,
    dataset = dataset,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache,
    cg_gradient_method = cg_gradient_method,
    cg_hessian_method = cg_hessian_method,
    verbose = verbose
  )

  cg_state <- initialize_state_fn(
    mm = mm,
    par_cov = par_cov,
    par_s = par_s,
    lambda_s = lambda_s,
    cg_max_delta = cg_max_delta,
    build_model_fn = cg_runtime$build_model,
    build_penalty_fn = cg_runtime$build_penalty
  )

  list(runtime = cg_runtime, state = cg_state)
}
