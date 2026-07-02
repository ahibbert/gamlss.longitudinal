#' Collate the score for a single parameter for the copula and margin score components.
#' Copula score is assembled in .gl_rs_copula_parameter_score (optimizer-rs-score-copula.R) 
#' and margin score is assembled in .gl_rs_margin_parameter_score (optimizer-rs-score-margin.R).
#'
#' @noRd
.gl_rs_parameter_score <- function(
    par_name,
    discrete_scores,
    include_dlcopdpar,
    eta,
    eta_inv,
    response,
    margin_deriv,
    margin_dist,
    copula_dist,
    dataset,
    mm,
    calc_lik_out,
    copula_derivatives,
    dcdu1,
    dcdu2,
    copula_d,
    log_lik,
    pair_cache,
    check_dlcopdpar_gradient,
    outer_only_run_counter,
    verbose,
    margin_params = c("mu", "sigma", "nu", "tau"),
    copula_score_fn = .gl_rs_copula_parameter_score,
    margin_score_fn = .gl_rs_margin_parameter_score) {
  if (
    !is.null(discrete_scores) &&
      (!par_name %in% margin_params || isTRUE(include_dlcopdpar))
  ) {
    d1 <- as.matrix(discrete_scores[[par_name]])
    colnames(d1) <- paste0("dld", par_name)
    return(list(
      d1 = d1,
      d1_m = NULL,
      d1_cop = NULL,
      path = "discrete"
    ))
  }

  if (!par_name %in% margin_params) {
    return(list(
      d1 = copula_score_fn(
        par_name = par_name,
        eta = eta,
        response = response,
        calc_lik_out = calc_lik_out,
        copula_derivatives = copula_derivatives
      ),
      d1_m = NULL,
      d1_cop = NULL,
      path = "copula"
    ))
  }

  margin_score <- margin_score_fn(
    par_name = par_name,
    margin_deriv = margin_deriv,
    include_dlcopdpar = include_dlcopdpar,
    eta = eta,
    eta_inv = eta_inv,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    dataset = dataset,
    mm = mm,
    calc_lik_out = calc_lik_out,
    dcdu1 = dcdu1,
    dcdu2 = dcdu2,
    copula_d = copula_d,
    log_lik = log_lik,
    pair_cache = pair_cache,
    check_dlcopdpar_gradient = check_dlcopdpar_gradient,
    outer_only_run_counter = outer_only_run_counter,
    verbose = verbose
  )
  margin_score$path <- "margin"
  margin_score
}
