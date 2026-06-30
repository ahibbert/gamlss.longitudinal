#' Build an RS backfitting runner bound to current derivative state
#'
#' @noRd
.gl_build_rs_backfitting_runner <- function(
    design_info,
    w_k_vec,
    z_k,
    rs_smooth_trust_radius,
    rs_calc_eta,
    calc_lik_out,
    pair_cache,
    margin_eval_cache,
    backfitting_iteration_fn = .gl_rs_backfitting_iteration) {
  function(par_s,
           par_cov,
           beta_start,
           lambda_s,
           first_inner_run,
           K,
           margin_dist,
           copula_dist,
           dataset,
           mm,
           copula_link,
           df_s,
           step_size,
           par_name) {
    backfitting_iteration_fn(
      par_s = par_s,
      par_cov = par_cov,
      beta_start = beta_start,
      lambda_s = lambda_s,
      K = K,
      margin_dist = margin_dist,
      copula_dist = copula_dist,
      dataset = dataset,
      mm = mm,
      df_s = df_s,
      step_size = step_size,
      par_name = par_name,
      design_info = design_info,
      w_k_vec = w_k_vec,
      z_k = z_k,
      rs_smooth_trust_radius = rs_smooth_trust_radius,
      rs_calc_eta = rs_calc_eta,
      calc_lik_out = calc_lik_out,
      pair_cache = pair_cache,
      margin_eval_cache = margin_eval_cache
    )
  }
}
