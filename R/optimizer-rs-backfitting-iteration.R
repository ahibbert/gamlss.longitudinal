#' Compute one RS backfitting proposal
#'
#' Builds the penalized weighted least-squares update for one parameter block,
#' refreshes fitted values, and evaluates the resulting likelihood and GAIC.
#'
#' @noRd
.gl_rs_backfitting_iteration <- function(
    par_s,
    par_cov,
    beta_start,
    lambda_s,
    K,
    margin_dist,
    copula_dist,
    dataset,
    mm,
    df_s,
    step_size,
    par_name,
    design_info,
    w_k_vec,
    z_k,
    rs_smooth_trust_radius,
    rs_calc_eta,
    calc_lik_out,
    pair_cache = NULL,
    margin_eval_cache = NULL,
    likelihood_fn = calc_likelihood_minimal,
    copula_likelihood_update_fn = .calc_likelihood_update_copula) {
  X <- design_info$X
  fixed_names <- design_info$fixed_names
  smooth_penalty_meta <- design_info$smooth_penalty_meta

  # Setup penalty matrix and compute smooth EDFs for each smooth term
  pen_mat <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  if (length(par_s[[par_name]]) > 0) {
    for (s_name in names(smooth_penalty_meta)) {
      meta <- smooth_penalty_meta[[s_name]]
      B <- meta$B
      idx <- meta$idx
      S <- meta$S_base

      pen_mat[idx, idx] <- S * lambda_s[[par_name]][[s_name]]

      BtWB_s <- t(B) %*% (B * as.vector(w_k_vec))
      df_s[[par_name]][[s_name]] <- sum(.solve_linear_system(BtWB_s + pen_mat[idx, idx]) * BtWB_s)
    }
  }

  # Compute the penalized weighted least-squares update for the parameter block
  XtWX <- t(X) %*% (X * w_k_vec)
  XtWz <- t(X) %*% (z_k * w_k_vec)
  beta_update <- as.vector(.solve_linear_system(XtWX + pen_mat, XtWz))
  beta_new <- beta_start * (1 - step_size) + step_size * beta_update

  # Apply trust region constraint to smooth parameter updates if specified
  if (is.finite(rs_smooth_trust_radius) && length(par_s[[par_name]]) > 0) {
    for (s_name in names(smooth_penalty_meta)) {
      idx <- smooth_penalty_meta[[s_name]]$idx
      delta_s <- beta_new[idx] - beta_start[idx]
      delta_norm <- sqrt(sum(delta_s^2))
      if (is.finite(delta_norm) && delta_norm > rs_smooth_trust_radius) {
        beta_new[idx] <- beta_start[idx] + delta_s * (rs_smooth_trust_radius / delta_norm)
      }
    }
  }

  # Update parameter vectors and evaluate the resulting likelihood and GAIC
  temp_par_cov_new <- beta_new[fixed_names]
  par_cov_new <- c(temp_par_cov_new, par_cov[!names(par_cov) %in% names(temp_par_cov_new)])

  temp_par_s_new <- beta_new[!names(beta_new) %in% names(par_cov_new)]
  par_s_new <- par_s
  for (s_name in names(par_s[[par_name]])) {
    smooth_col_names <- colnames(X)[smooth_penalty_meta[[s_name]]$idx]
    par_s_new[[par_name]][[s_name]] <- temp_par_s_new[smooth_col_names]
  }

  eta_out <- rs_calc_eta(par_cov_current = par_cov_new, par_s_current = par_s_new)
  eta_inv <- eta_out$eta_inv

  if (par_name %in% c("theta", "zeta") && isTRUE(getOption("gamlss.longitudinal.fast_copula_lik", TRUE))) {
    calc_lik_out_end <- copula_likelihood_update_fn(
      eta_inv = eta_inv,
      base_lik = calc_lik_out,
      copula_dist = copula_dist,
      pair_cache = pair_cache
    )
  } else {
    calc_lik_out_end <- likelihood_fn(
      eta_inv,
      mm = mm$x,
      margin_dist,
      copula_dist,
      calc_d2 = FALSE,
      response = dataset$response,
      response_margin = dataset$time,
      response_subject = dataset$subject,
      pair_cache = pair_cache,
      margin_eval_cache = margin_eval_cache,
      calc_margin_deriv = FALSE
    )
  }

  GAIC_lambda_k <- -2 * calc_lik_out_end$log_lik["joint"] + K * sum(unlist(df_s[[par_name]]))

  return_list <- list(par_cov_new, par_s_new, calc_lik_out_end, GAIC_lambda_k, df_s)
  names(return_list) <- c("par_cov", "par_s", "calc_lik_out_end", "GAIC_lambda_k", "df_s")
  return_list
}
