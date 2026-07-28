#' Calculate minimal joint likelihood components
#'
#' @param pair_cache Optional cache built by build_copula_pair_cache to reuse pair indexing across repeated likelihood calls.
#'
#' @noRd
calc_likelihood_minimal <- function(eta_inv, mm, margin_dist, copula_dist, calc_d2 = FALSE, response, response_margin, response_subject, penalize_smooth = FALSE, par_s = NA, pair_cache = NULL, margin_eval_cache = NULL, calc_margin_deriv = TRUE, margin_deriv_names = NULL) {
  # Setup input matrix of response and parameters
  # response=dataset$response; response_subject=dataset$subject; response_margin=dataset$time; dataset=NA
  if (is.null(pair_cache)) {
    pair_cache <- build_copula_pair_cache(response, response_margin, response_subject)
  }
  if (is.null(margin_eval_cache) || !identical(margin_eval_cache$calc_d2, calc_d2)) {
    margin_eval_cache <- .build_margin_eval_cache(margin_dist, calc_d2 = calc_d2)
  }

  margin_names <- pair_cache$margin_names
  num_margins <- pair_cache$num_margins
  n_obs <- pair_cache$n_obs
  discrete_margin <- .is_discrete_margin(margin_dist)

  order_margin <- cbind(response_margin, response_subject)
  colnames(order_margin) <- c("time", "subject")

  margin_eval <- .gl_likelihood_evaluate_margins(
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = margin_dist,
    margin_eval_cache = margin_eval_cache,
    response = response,
    discrete_margin = discrete_margin,
    calc_margin_deriv = calc_margin_deriv,
    margin_deriv_names = margin_deriv_names
  )
  margin_deriv <- margin_eval$margin_deriv
  margin_p <- margin_eval$margin_p
  margin_p_lower <- margin_eval$margin_p_lower
  margin_d <- margin_eval$margin_d
  likelihood_type <- margin_eval$likelihood_type

  ################ COPULA DERIVATIVES
  # First calculate margin F(x1), F(x2) as inputs to copula
  copula_eval <- .gl_likelihood_evaluate_copula_pairs(
    eta_inv = eta_inv,
    pair_cache = pair_cache,
    copula_dist = copula_dist,
    margin_p = margin_p,
    margin_d = margin_d,
    likelihood_type = likelihood_type,
    margin_p_lower = margin_p_lower
  )
  row_id1 <- copula_eval$row_id1
  row_id2 <- copula_eval$row_id2
  order_copula <- copula_eval$order_copula
  Fx_1_2 <- copula_eval$Fx_1_2
  Fx_1_2_lower <- copula_eval$Fx_1_2_lower
  pair_complete <- copula_eval$pair_complete
  par1 <- copula_eval$par1
  par2 <- copula_eval$par2
  copula_d <- copula_eval$copula_d
  copula_rect_prob <- copula_eval$copula_rect_prob

  ######## COMBINE MARGINS AND COPULA DERVIATIVES

  margin_loglik_terms <- log(margin_d[!is.na(margin_d)])
  margin_loglik_terms <- margin_loglik_terms[is.finite(margin_loglik_terms)]
  copula_loglik_terms <- log(copula_d[pair_complete])
  copula_loglik_terms <- copula_loglik_terms[is.finite(copula_loglik_terms)]

  log_lik <- c(sum(margin_loglik_terms), sum(copula_loglik_terms), sum(margin_loglik_terms) + sum(copula_loglik_terms))
  names(log_lik) <- c("marginal", "copula", "joint")

  copula_p <- rep(NA_real_, length(copula_d))

  return_list <- list(log_lik, margin_d, copula_d, margin_p, copula_p, Fx_1_2, order_copula, margin_deriv, pair_complete, par1, par2, row_id1, row_id2, pair_cache$theta_index_map, likelihood_type, margin_p_lower, Fx_1_2_lower, copula_rect_prob)
  names(return_list) <- c("log_lik", "margin_d", "copula_d", "margin_p", "copula_p", "Fx_1_2", "order_copula", "margin_deriv", "pair_complete", "copula_par1", "copula_par2", "copula_row_id1", "copula_row_id2", "copula_theta_index_map", "likelihood_type", "margin_p_lower", "Fx_1_2_lower", "copula_rect_prob")
  return(return_list)
}
