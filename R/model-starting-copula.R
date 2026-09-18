.select_t_copula_zeta_start <- function(dataset, margin_dist, copula_dist, margin_par, theta_start) {
  # Select the finite grid candidate with the largest joint likelihood. Merely
  # taking the first finite value can pin approximately Gaussian data against
  # the lower degrees-of-freedom boundary, where link-scale steps are tiny.

  # Note theta is selected earlier using kendalls tau to copula parameter transform

  zeta_grid <- c(2.05, 2.2, 2.5, 3, 4, 5, 8, 12, 20, 35)
  fallback_zeta <- 3
  param_names <- c(names(margin_dist$parameters), get_copula_dist(copula_dist)$parameters)
  mm_stub <- as.list(setNames(rep(1, length(param_names)), param_names))
  pair_cache <- build_copula_pair_cache(dataset$response, dataset$time, dataset$subject)
  base_eta_inv <- c(as.list(margin_par), list(theta = as.numeric(theta_start)[1]))
  best_zeta <- fallback_zeta
  best_log_lik <- -Inf

  for (candidate_zeta in zeta_grid) {
    eta_inv <- base_eta_inv
    eta_inv$zeta <- candidate_zeta
    candidate_fit <- tryCatch(
      calc_likelihood_minimal(
        eta_inv = eta_inv,
        mm = mm_stub,
        margin_dist = margin_dist,
        copula_dist = copula_dist,
        calc_d2 = FALSE,
        response = dataset$response,
        response_margin = dataset$time,
        response_subject = dataset$subject,
        pair_cache = pair_cache,
        calc_margin_deriv = FALSE
      ),
      error = function(e) NULL
    )
    candidate_log_lik <- if (is.null(candidate_fit)) {
      NA_real_
    } else {
      as.numeric(candidate_fit$log_lik["joint"])
    }
    if (is.finite(candidate_log_lik) && candidate_log_lik > best_log_lik) {
      best_log_lik <- candidate_log_lik
      best_zeta <- candidate_zeta
    }
  }
  best_zeta
}
