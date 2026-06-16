.select_t_copula_zeta_start <- function(dataset, margin_dist, copula_dist, margin_par, theta_start) {
  # Grid is ordered low-to-high: return the first candidate that yields a finite
  # joint log-likelihood. Starting as low as possible avoids the optimizer being
  # trapped at high df where the link-scale step sizes collapse (known issue).
  # This will still be optimised over during fit, so the staring fit just needs to be at a point where the likelihood is reasonably shaped.

  # Note theta is selected earlier using kendalls tau to copula parameter transform

  zeta_grid <- c(2.05, 2.2, 2.5, 3, 4, 5, 8, 12, 20, 35)
  fallback_zeta <- 3
  param_names <- c(names(margin_dist$parameters), get_copula_dist(copula_dist)$parameters)
  mm_stub <- as.list(setNames(rep(1, length(param_names)), param_names))
  pair_cache <- build_copula_pair_cache(dataset$response, dataset$time, dataset$subject)
  base_eta_inv <- c(as.list(margin_par), list(theta = as.numeric(theta_start)[1]))

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
        pair_cache = pair_cache
      ),
      error = function(e) NULL
    )
    if (!is.null(candidate_fit) && is.finite(candidate_fit$log_lik["joint"])) {
      return(candidate_zeta)
    }
  }
  fallback_zeta
}
