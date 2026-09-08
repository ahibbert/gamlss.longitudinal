#' @keywords internal
#' @noRd
.calc_likelihood_update_copula <- function(eta_inv, base_lik, copula_dist, pair_cache) {
  copula_eval <- .gl_likelihood_evaluate_copula_pairs(
    eta_inv = eta_inv,
    pair_cache = pair_cache,
    copula_dist = copula_dist,
    margin_p = base_lik$margin_p,
    margin_d = base_lik$margin_d,
    likelihood_type = base_lik$likelihood_type,
    margin_log_d = base_lik$margin_log_d %||% log(base_lik$margin_d),
    margin_p_lower = base_lik$margin_p_lower,
    margin_log_survival = base_lik$margin_log_survival,
    margin_log_survival_lower = base_lik$margin_log_survival_lower,
    Fx_1_2 = base_lik$Fx_1_2,
    Fx_1_2_lower = base_lik$Fx_1_2_lower
  )
  pair_complete <- copula_eval$pair_complete
  par1 <- copula_eval$par1
  par2 <- copula_eval$par2
  copula_d <- copula_eval$copula_d
  copula_log_d <- copula_eval$copula_log_d
  copula_rect_prob <- copula_eval$copula_rect_prob

  margin_log_d <- base_lik$margin_log_d %||% log(base_lik$margin_d)
  margin_included <- base_lik$margin_included %||% !is.na(base_lik$margin_d)
  combined <- .gl_likelihood_combine_components(
    margin_log_d = margin_log_d,
    margin_included = margin_included,
    copula_log_d = copula_log_d,
    pair_included = pair_complete,
    pair_input_valid = copula_eval$pair_input_valid
  )

  base_lik$log_lik <- combined$log_lik
  base_lik$copula_d <- copula_d
  base_lik$copula_log_d <- copula_log_d
  base_lik$copula_rect_prob <- copula_rect_prob
  base_lik$pair_complete <- pair_complete
  base_lik$pair_included <- combined$pair_included
  base_lik$pair_input_valid <- combined$pair_input_valid
  base_lik$pair_contribution_valid <- combined$pair_contribution_valid
  base_lik$valid <- combined$valid
  base_lik$failure <- combined$failure
  base_lik$contribution_counts <- combined$contribution_counts
  base_lik$copula_par1 <- par1
  base_lik$copula_par2 <- par2
  base_lik
}
