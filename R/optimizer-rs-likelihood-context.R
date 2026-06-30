#' Evaluate and unpack the RS likelihood context
#'
#' @noRd
.gl_rs_likelihood_context <- function(
    eta_inv,
    mm,
    margin_dist,
    copula_dist,
    dataset,
    pair_cache,
    margin_eval_cache,
    likelihood_fn = calc_likelihood_minimal) {
  calc_lik_out <- likelihood_fn(
    eta_inv,
    mm = mm$x,
    margin_dist,
    copula_dist,
    calc_d2 = FALSE,
    response = dataset$response,
    response_margin = dataset$time,
    response_subject = dataset$subject,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache
  )

  list(
    calc_lik_out = calc_lik_out,
    log_lik = calc_lik_out$log_lik,
    margin_d = calc_lik_out$margin_d,
    margin_p = calc_lik_out$margin_p,
    margin_deriv = calc_lik_out$margin_deriv,
    copula_d = calc_lik_out$copula_d,
    copula_p = calc_lik_out$copula_p,
    Fx_1_2 = calc_lik_out$Fx_1_2,
    order_copula = calc_lik_out$order_copula
  )
}

#' Precompute RS scores for discrete rectangle likelihoods
#'
#' @noRd
.gl_rs_discrete_scores <- function(
    calc_lik_out,
    eta_inv,
    mm,
    margin_dist,
    copula_dist,
    dataset,
    pair_cache,
    discrete_score_method,
    score_fn = .calc_discrete_rectangle_scores) {
  if (!identical(calc_lik_out$likelihood_type, "discrete_rectangle")) {
    return(NULL)
  }

  score_fn(
    eta_inv,
    mm$x,
    margin_dist,
    copula_dist,
    dataset$response,
    dataset$time,
    dataset$subject,
    pair_cache = pair_cache,
    calc_lik = calc_lik_out,
    method = discrete_score_method
  )
}

#' Build RS copula-derivative context from current likelihood values
#'
#' @noRd
.gl_rs_copula_derivative_context <- function(
    eta_inv,
    Fx_1_2,
    copula_dist,
    calc_lik_out,
    derivative_fn = calc_copula_derivatives) {
  Fx_1_2[Fx_1_2 > 1] <- 1
  Fx_1_2[Fx_1_2 < 0] <- 0

  copula_derivatives <- derivative_fn(
    eta_inv,
    Fx_1_2,
    copula_dist,
    par1 = calc_lik_out$copula_par1,
    par2 = calc_lik_out$copula_par2,
    pair_complete = calc_lik_out$pair_complete
  )

  list(
    Fx_1_2 = Fx_1_2,
    copula_derivatives = copula_derivatives,
    dldth = copula_derivatives$dldth,
    dcdth = copula_derivatives$dcdth,
    dcdu1 = copula_derivatives$dcdu1,
    dcdu2 = copula_derivatives$dcdu2,
    dldz = if ("zeta" %in% names(eta_inv)) copula_derivatives$dldz else NULL,
    dcdz = if ("zeta" %in% names(eta_inv)) copula_derivatives$dcdz else NULL
  )
}
