#' CG penalized objective
#'
#' @noRd
.gl_cg_objective <- function(beta_vec, loglik, penalty_current) {
  as.numeric(loglik) - 0.5 * sum(as.numeric(beta_vec) * as.numeric(penalty_current %*% beta_vec))
}

#' Evaluate one CG coefficient vector
#'
#' @noRd
.gl_evaluate_cg_beta <- function(
    beta_vec,
    mm_cg,
    par_cov_template,
    par_s_template,
    margin_dist,
    copula_link,
    copula_dist,
    dataset,
    pair_cache,
    margin_eval_cache,
    calc_eta_fn = calc_eta,
    likelihood_fn = calc_likelihood_minimal) {
  unpacked <- .gl_unpack_cg_beta(
    beta_vec,
    par_cov_template = par_cov_template,
    par_s_template = par_s_template
  )

  eta_out <- calc_eta_fn(
    beta_vec,
    mm_cg,
    margin_dist,
    copula_link,
    par_s = setNames(lapply(names(mm_cg$x), function(x) list()), names(mm_cg$x))
  )
  eta_inv <- eta_out$eta_inv

  if (any(!is.finite(unlist(eta_inv, use.names = FALSE)))) {
    return(NULL)
  }

  if (copula_dist %in% c("N", "t") &&
    "theta" %in% names(eta_inv) &&
    any(abs(eta_inv$theta) >= 0.999, na.rm = TRUE)) {
    return(NULL)
  }

  positive_names <- intersect(names(eta_inv), c("sigma", "tau", "zeta"))
  for (pn in positive_names) {
    if (any(eta_inv[[pn]] <= 1e-8, na.rm = TRUE)) {
      return(NULL)
    }
  }

  lik <- tryCatch(likelihood_fn(
    eta_inv,
    mm = mm_cg$x,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    calc_d2 = FALSE,
    response = dataset$response,
    response_margin = dataset$time,
    response_subject = dataset$subject,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache
  ), error = function(e) NULL)

  if (is.null(lik) || identical(lik$valid, FALSE) ||
      !is.finite(as.numeric(lik$log_lik["joint"]))) {
    return(NULL)
  }

  list(
    loglik = as.numeric(lik$log_lik["joint"]),
    calc_lik = lik,
    eta_out = eta_out,
    par_cov = unpacked$par_cov,
    par_s = unpacked$par_s
  )
}
