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
    current_calc_lik_out = NULL,
    par_name = NULL,
    likelihood_fn = calc_likelihood_minimal,
    copula_likelihood_update_fn = .calc_likelihood_update_copula) {

  likelihood_args <- list(
    eta_inv,
    mm = mm$x,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    calc_d2 = FALSE,
    response = dataset$response,
    response_margin = dataset$time,
    response_subject = dataset$subject,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache
  )
  margin_deriv_names <- .gl_rs_margin_derivatives_for_parameter(par_name)
  if (par_name %in% c("mu", "sigma", "nu", "tau") &&
      !is.null(current_calc_lik_out) &&
      isTRUE(getOption("gamlss.longitudinal.fast_margin_score_lik", TRUE))) {
    margin_deriv_input <- .gl_likelihood_margin_derivative_input(
      eta_inv = eta_inv,
      mm = mm$x,
      margin_dist = margin_dist,
      response = dataset$response
    )
    margin_deriv <- .gl_likelihood_evaluate_margin_derivatives(
      margin_deriv_input = margin_deriv_input,
      margin_eval_cache = margin_eval_cache,
      response = dataset$response,
      calc_margin_deriv = TRUE,
      margin_deriv_names = margin_deriv_names
    )
    calc_lik_out <- current_calc_lik_out
    calc_lik_out$margin_deriv <- margin_deriv
  } else if (par_name %in% c("theta", "zeta") &&
      !is.null(current_calc_lik_out) &&
      isTRUE(getOption("gamlss.longitudinal.fast_copula_lik", TRUE))) {
    calc_lik_out <- copula_likelihood_update_fn(
      eta_inv = eta_inv,
      base_lik = current_calc_lik_out,
      copula_dist = copula_dist,
      pair_cache = pair_cache
    )
  } else {
    likelihood_formals <- names(formals(likelihood_fn))
    if (!is.null(margin_deriv_names) &&
        ("margin_deriv_names" %in% likelihood_formals || "..." %in% likelihood_formals)) {
      likelihood_args$margin_deriv_names <- margin_deriv_names
    }
    calc_lik_out <- do.call(likelihood_fn, likelihood_args)
  }

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

#' Select the marginal distribution derivatives needed for one RS parameter block
#'
#' @noRd
.gl_rs_margin_derivatives_for_parameter <- function(
    par_name,
    margin_params = c("mu", "sigma", "nu", "tau")) {
  if (is.null(par_name)) {
    return(NULL)
  }
  margin_deriv_subnames <- c(mu = "m", sigma = "d", nu = "v", tau = "t")
  if (par_name %in% margin_params && par_name %in% names(margin_deriv_subnames)) {
    return(paste0("dld", margin_deriv_subnames[[par_name]]))
  }
  character(0)
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
    par_name = NULL,
    include_dlcopdpar = TRUE,
    derivative_fn = calc_copula_derivatives) {
  Fx_1_2[Fx_1_2 > 1] <- 1
  Fx_1_2[Fx_1_2 < 0] <- 0

  derivative_args <- list(
    eta_inv,
    Fx_1_2,
    copula_dist,
    par1 = calc_lik_out$copula_par1,
    par2 = calc_lik_out$copula_par2,
    pair_complete = calc_lik_out$pair_complete
  )
  requested_derivatives <- .gl_rs_copula_derivatives_for_parameter(
    par_name = par_name,
    include_dlcopdpar = include_dlcopdpar,
    has_zeta = "zeta" %in% names(eta_inv)
  )
  derivative_formals <- names(formals(derivative_fn))
  if (!is.null(requested_derivatives) &&
      ("derivatives" %in% derivative_formals || "..." %in% derivative_formals)) {
    derivative_args$derivatives <- requested_derivatives
  }
  copula_derivatives <- do.call(derivative_fn, derivative_args)

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

#' Select the copula derivatives needed for one RS parameter block
#'
#' @noRd
.gl_rs_copula_derivatives_for_parameter <- function(
    par_name,
    include_dlcopdpar,
    has_zeta,
    margin_params = c("mu", "sigma", "nu", "tau")) {
  if (is.null(par_name)) {
    return(NULL)
  }
  if (par_name %in% margin_params) {
    if (isTRUE(include_dlcopdpar)) {
      return(c("dcdu1", "dcdu2"))
    }
    return(character(0))
  }
  if (identical(par_name, "theta")) {
    return("dldth")
  }
  if (identical(par_name, "zeta") && isTRUE(has_zeta)) {
    return("dldz")
  }
  NULL
}
