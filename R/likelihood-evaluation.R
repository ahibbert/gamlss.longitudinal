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
  response_inclusion <- !is.na(response)
  if (!is.null(pair_cache$observed_margin_base) &&
      !identical(as.logical(pair_cache$observed_margin_base), response_inclusion)) {
    stop(
      "'pair_cache' was built for a different response missingness pattern.",
      call. = FALSE
    )
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
  margin_log_survival <- margin_eval$margin_log_survival
  margin_log_survival_lower <- margin_eval$margin_log_survival_lower
  margin_d <- margin_eval$margin_d
  margin_log_d <- margin_eval$margin_log_d
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
    margin_log_d = margin_log_d,
    margin_p_lower = margin_p_lower,
    margin_log_survival = margin_log_survival,
    margin_log_survival_lower = margin_log_survival_lower
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
  copula_log_d <- copula_eval$copula_log_d
  copula_rect_prob <- copula_eval$copula_rect_prob

  ######## COMBINE MARGINS AND COPULA DERVIATIVES

  margin_included <- pair_cache$observed_margin_base
  if (is.null(margin_included)) margin_included <- response_inclusion
  combined <- .gl_likelihood_combine_components(
    margin_log_d = margin_log_d,
    margin_included = margin_included,
    copula_log_d = copula_log_d,
    pair_included = pair_complete,
    pair_input_valid = copula_eval$pair_input_valid
  )

  copula_p <- rep(NA_real_, length(copula_d))

  list(
    log_lik = combined$log_lik,
    margin_d = margin_d,
    copula_d = copula_d,
    margin_p = margin_p,
    copula_p = copula_p,
    Fx_1_2 = Fx_1_2,
    order_copula = order_copula,
    margin_deriv = margin_deriv,
    pair_complete = pair_complete,
    copula_par1 = par1,
    copula_par2 = par2,
    copula_row_id1 = row_id1,
    copula_row_id2 = row_id2,
    copula_theta_index_map = pair_cache$theta_index_map,
    likelihood_type = likelihood_type,
    margin_p_lower = margin_p_lower,
    Fx_1_2_lower = Fx_1_2_lower,
    copula_rect_prob = copula_rect_prob,
    # Keep all legacy fields above in their historical positional order.
    margin_log_d = margin_log_d,
    copula_log_d = copula_log_d,
    margin_included = combined$margin_included,
    pair_included = combined$pair_included,
    margin_contribution_valid = combined$margin_contribution_valid,
    pair_input_valid = combined$pair_input_valid,
    pair_contribution_valid = combined$pair_contribution_valid,
    valid = combined$valid,
    failure = combined$failure,
    contribution_counts = combined$contribution_counts,
    margin_log_survival = margin_log_survival,
    margin_log_survival_lower = margin_log_survival_lower
  )
}
