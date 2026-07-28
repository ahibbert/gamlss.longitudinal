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
    margin_p_lower = base_lik$margin_p_lower,
    Fx_1_2 = base_lik$Fx_1_2,
    Fx_1_2_lower = base_lik$Fx_1_2_lower
  )
  pair_complete <- copula_eval$pair_complete
  par1 <- copula_eval$par1
  par2 <- copula_eval$par2
  copula_d <- copula_eval$copula_d
  copula_rect_prob <- copula_eval$copula_rect_prob

  copula_loglik_terms <- log(copula_d[pair_complete])
  copula_loglik_terms <- copula_loglik_terms[is.finite(copula_loglik_terms)]
  copula_loglik <- sum(copula_loglik_terms)
  marginal_loglik <- as.numeric(base_lik$log_lik["marginal"])
  log_lik <- c(marginal_loglik, copula_loglik, marginal_loglik + copula_loglik)
  names(log_lik) <- c("marginal", "copula", "joint")

  base_lik$log_lik <- log_lik
  base_lik$copula_d <- copula_d
  base_lik$copula_rect_prob <- copula_rect_prob
  base_lik$pair_complete <- pair_complete
  base_lik$copula_par1 <- par1
  base_lik$copula_par2 <- par2
  base_lik
}
