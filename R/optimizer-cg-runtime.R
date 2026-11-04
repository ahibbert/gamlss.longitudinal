#' Build CG runtime helper closures for one fit
#'
#' @noRd
.gl_build_cg_runtime_helpers <- function(
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
    calc_eta_fn = calc_eta,
    likelihood_fn = calc_likelihood_minimal,
    analytical_hessian_fn = calc_analytical_hessian) {
  cg_eval <- function(beta_vec, mm_cg) {
    .gl_evaluate_cg_beta(
      beta_vec = beta_vec,
      mm_cg = mm_cg,
      par_cov_template = par_cov,
      par_s_template = par_s,
      margin_dist = margin_dist,
      copula_link = copula_link,
      copula_dist = copula_dist,
      dataset = dataset,
      pair_cache = pair_cache,
      margin_eval_cache = margin_eval_cache,
      calc_eta_fn = calc_eta_fn,
      likelihood_fn = likelihood_fn
    )
  }

  cg_finite_hessian_block <- function(beta_vec, block_names, mm_cg, h = 1e-4) {
    .gl_cg_finite_hessian_block(
      beta_vec = beta_vec,
      block_names = block_names,
      mm_cg = mm_cg,
      eval_fn = cg_eval,
      h = h
    )
  }

  list(
    build_model = .gl_build_cg_model,
    build_penalty = function(beta_names, lambda_current) {
      .gl_build_cg_penalty(beta_names, lambda_current, par_s = par_s, mm = mm)
    },
    evaluate = cg_eval,
    objective = .gl_cg_objective,
    gradient = function(beta_vec, base_ll, mm_cg) {
      .gl_cg_finite_gradient(
        beta_vec = beta_vec,
        base_ll = base_ll,
        mm_cg = mm_cg,
        eval_fn = cg_eval,
        gradient_method = cg_gradient_method
      )
    },
    finite_hessian_block = cg_finite_hessian_block,
    observed_hessian = function(tmp_obj, beta_vec, mm_cg, context = "CG iteration") {
      .gl_cg_observed_hessian(
        tmp_obj = tmp_obj,
        beta_vec = beta_vec,
        mm_cg = mm_cg,
        analytical_fn = function(obj) analytical_hessian_fn(obj, progress = FALSE),
        finite_fn = cg_finite_hessian_block,
        hessian_method = cg_hessian_method,
        verbose = verbose,
        context = context
      )
    },
    smooth_edf_list = function(H_obs_current, penalty_current, beta_names) {
      .gl_cg_smooth_edf_list(
        H_obs_current = H_obs_current,
        penalty_current = penalty_current,
        beta_names = beta_names,
        par_s = par_s
      )
    }
  )
}
