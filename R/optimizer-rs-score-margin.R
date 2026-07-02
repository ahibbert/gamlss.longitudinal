#' Return score for all margin parameters 
#' (first derivatives of log-likelihood w.r.t. each margin parameter)
#' The actual derivatives are calculated in the .gl_rs_margin_parameter_score function.
#'
#' @noRd
.gl_rs_margin_parameter_score <- function(
    par_name,
    margin_deriv,
    include_dlcopdpar,
    eta,
    eta_inv,
    margin_dist,
    copula_dist,
    dataset,
    mm,
    calc_lik_out,
    dcdu1,
    dcdu2,
    copula_d,
    log_lik,
    pair_cache,
    check_dlcopdpar_gradient,
    outer_only_run_counter,
    verbose,
    fx_deriv_fn = calc_Fx_derivatives,
    dlcopdpar_fn = .calc_dlcopdpar_indexed,
    gradient_check_fn = check_dlcopdpar_gradient_margin_score) {
  margin_deriv_subnames <- c("m", "d", "v", "t")
  names(margin_deriv_subnames) <- c("mu", "sigma", "nu", "tau")

  d1 <- as.matrix(margin_deriv[grepl(paste("dld", margin_deriv_subnames[par_name], sep = ""), names(margin_deriv))][[1]])
  colnames(d1) <- paste("dld", par_name, sep = "")

  d1_m <- d1
  d1_cop <- d1 * 0

  if (include_dlcopdpar == TRUE) {
    nd_impact_F <- fx_deriv_fn(
      eta_inv,
      mm$x,
      margin_dist,
      response = dataset$response,
      par_names = par_name
    )

    d1_cop <- dlcopdpar_fn(
      row_id1 = calc_lik_out$copula_row_id1,
      row_id2 = calc_lik_out$copula_row_id2,
      dcdu1 = dcdu1,
      dcdu2 = dcdu2,
      copula_d = copula_d,
      F_nd = nd_impact_F[[par_name]],
      n_obs = length(dataset$response),
      pair_complete = calc_lik_out$pair_complete
    )

    d1_m <- d1
    d1 <- d1_m + d1_cop

    if (check_dlcopdpar_gradient && outer_only_run_counter == 1) {
      gradient_check <- gradient_check_fn(
        eta = eta,
        eta_inv = eta_inv,
        par_name = par_name,
        margin_dist = margin_dist,
        copula_dist = copula_dist,
        dataset = dataset,
        mm = mm$x,
        pair_cache = pair_cache,
        d1 = d1,
        base_loglik = log_lik["joint"],
        verbose = verbose
      )

      if (isTRUE(gradient_check$warned)) {
        warning(gradient_check$message, call. = FALSE)
      }
    }
  }

  d1 <- d1[, grepl(par_name, colnames(d1))]

  list(
    d1 = d1,
    d1_m = d1_m,
    d1_cop = d1_cop
  )
}
