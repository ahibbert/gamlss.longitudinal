#' Evaluate pair-level copula likelihood inputs
#'
#' @noRd
.gl_likelihood_evaluate_copula_pairs <- function(
    eta_inv,
    pair_cache,
    copula_dist,
    margin_p,
    margin_d,
    likelihood_type,
    margin_p_lower = NULL,
    Fx_1_2 = NULL,
    Fx_1_2_lower = NULL) {
  row_id1 <- pair_cache$row_id1
  row_id2 <- pair_cache$row_id2
  n_obs <- pair_cache$n_obs

  if (is.null(Fx_1_2)) {
    Fx_1_2 <- pair_cache$Fx_1_2_template
    if (is.null(Fx_1_2)) {
      Fx_1_2 <- matrix(NA_real_, nrow = length(row_id1), ncol = 2)
      colnames(Fx_1_2) <- c("u1", "u2")
    }
    if (length(row_id1) > 0) {
      Fx_1_2[, 1] <- margin_p[row_id1]
      Fx_1_2[, 2] <- margin_p[row_id2]
    }
  }

  if (identical(likelihood_type, "discrete_rectangle") && is.null(Fx_1_2_lower)) {
    Fx_1_2_lower <- Fx_1_2
    if (length(row_id1) > 0) {
      Fx_1_2_lower[, 1] <- margin_p_lower[row_id1]
      Fx_1_2_lower[, 2] <- margin_p_lower[row_id2]
    }
  }

  pair_complete <- pair_cache$observed_pair_base &
    is.finite(Fx_1_2[, 1]) & is.finite(Fx_1_2[, 2])
  if (identical(likelihood_type, "discrete_rectangle")) {
    pair_complete <- pair_complete &
      is.finite(Fx_1_2_lower[, 1]) & is.finite(Fx_1_2_lower[, 2]) &
      is.finite(margin_d[row_id1]) & is.finite(margin_d[row_id2]) &
      margin_d[row_id1] > 0 & margin_d[row_id2] > 0
  }

  par1 <- rep(NA_real_, length(row_id1))
  par2 <- rep(NA_real_, length(row_id1))
  if (length(row_id1) > 0) {
    theta_len <- length(eta_inv[["theta"]])
    if (theta_len == n_obs) {
      theta_idx <- row_id1
    } else {
      theta_idx <- pair_cache$theta_index_map[row_id1]
    }

    par1 <- eta_inv[["theta"]][theta_idx]
    if ("zeta" %in% names(eta_inv)) {
      par2 <- eta_inv[["zeta"]][theta_idx]
    } else {
      par2 <- rep(0, length(par1))
    }
  }

  pair_complete <- pair_complete & is.finite(par1) & is.finite(par2)

  Fx_eval <- Fx_1_2
  if (nrow(Fx_eval) > 0) {
    Fx_eval[!is.finite(Fx_eval)] <- 0.5
    Fx_eval[Fx_eval > 1] <- 1
    Fx_eval[Fx_eval < 0] <- 0
  }

  par1_eval <- par1
  par2_eval <- par2
  par1_eval[!is.finite(par1_eval)] <- 0
  par2_eval[!is.finite(par2_eval)] <- 0

  if (copula_dist == "C") {
    par1_eval[par1_eval >= 28] <- 27.9
  }

  copula_rect_prob <- NULL
  if (length(par1_eval) == 0) {
    copula_d <- numeric(0)
    copula_rect_prob <- numeric(0)
  } else if (identical(likelihood_type, "discrete_rectangle")) {
    Fx_lower_eval <- Fx_1_2_lower
    Fx_lower_eval[!is.finite(Fx_lower_eval)] <- 0.5
    Fx_lower_eval[Fx_lower_eval > 1] <- 1
    Fx_lower_eval[Fx_lower_eval < 0] <- 0
    copula_rect_prob <- .copula_rectangle_prob(
      Fx_eval[, 1], Fx_eval[, 2], Fx_lower_eval[, 1], Fx_lower_eval[, 2],
      family = copula_dist, par = par1_eval, par2 = par2_eval
    )
    denom <- margin_d[row_id1] * margin_d[row_id2]
    copula_d <- copula_rect_prob / denom
  } else {
    copula_d <- .copula_pdf(
      Fx_eval[, 1], Fx_eval[, 2],
      family = copula_dist, par = par1_eval, par2 = par2_eval
    )
  }

  if (length(copula_d) > 0) {
    copula_d[!is.finite(copula_d) | copula_d <= 0] <- 1
    copula_d[!pair_complete] <- 1
  }

  list(
    Fx_1_2 = Fx_1_2,
    Fx_1_2_lower = Fx_1_2_lower,
    order_copula = pair_cache$order_copula,
    pair_complete = pair_complete,
    par1 = par1,
    par2 = par2,
    row_id1 = row_id1,
    row_id2 = row_id2,
    copula_d = copula_d,
    copula_rect_prob = copula_rect_prob
  )
}
