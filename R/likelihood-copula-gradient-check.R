check_dlcopdpar_gradient_margin_score <- function(
    eta,
    eta_inv,
    par_name,
    margin_dist,
    copula_dist,
    dataset,
    mm,
    pair_cache,
    d1,
    base_loglik,
    verbose = FALSE) {
  if (!par_name %in% c("mu", "sigma", "nu", "tau")) {
    return(list(warned = FALSE, message = NULL))
  }

  eta_vec <- as.numeric(eta[[par_name]])
  score_vec <- as.numeric(d1)
  finite_idx <- which(is.finite(eta_vec) & is.finite(score_vec))
  if (length(finite_idx) == 0) {
    return(list(warned = FALSE, message = NULL))
  }

  probe_candidates <- unique(round(seq(1, length(finite_idx), length.out = min(3, length(finite_idx)))))
  probe_idx <- finite_idx[probe_candidates]
  eps <- 1e-6
  tolerance <- 1e-3
  diffs <- numeric(length(probe_idx))

  linkinv_fun <- margin_dist[[paste(par_name, ".linkinv", sep = "")]]
  if (is.null(linkinv_fun)) {
    return(list(warned = FALSE, message = NULL))
  }

  for (k in seq_along(probe_idx)) {
    idx <- probe_idx[k]
    eta_pert <- eta
    eta_pert[[par_name]][idx] <- eta_pert[[par_name]][idx] + eps

    eta_inv_pert <- eta_inv
    eta_inv_pert[[par_name]] <- linkinv_fun(eta_pert[[par_name]])

    lik_pert <- calc_likelihood_minimal(
      eta_inv_pert,
      mm = mm,
      margin_dist = margin_dist,
      copula_dist = copula_dist,
      calc_d2 = FALSE,
      response = dataset$response,
      response_margin = dataset$time,
      response_subject = dataset$subject,
      pair_cache = pair_cache,
      calc_margin_deriv = FALSE
    )$log_lik["joint"]

    finite_diff <- (lik_pert - base_loglik) / eps
    diffs[k] <- finite_diff - score_vec[idx]
  }

  max_abs_diff <- max(abs(diffs), na.rm = TRUE)
  if (!is.finite(max_abs_diff)) {
    return(list(warned = FALSE, message = NULL))
  }

  message_text <- paste0(
    "DLCOPDGRAD check for ", par_name,
    ": max abs(score - finite_diff) = ", signif(max_abs_diff, 4),
    " over ", length(probe_idx), " probe row(s)."
  )

  if (max_abs_diff > tolerance) {
    return(list(warned = TRUE, message = paste0(message_text, " Potential score mismatch detected.")))
  }

  if (!is.null(verbose) && verbose > 1) {
    message(message_text)
  }

  return(list(warned = FALSE, message = message_text))
}
