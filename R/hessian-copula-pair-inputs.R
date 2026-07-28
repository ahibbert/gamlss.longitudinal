#' @noRd
.copula_hessian_pair_inputs <- function(eta_inv, pair_cache, copula_dist) {
  n_obs <- length(eta_inv[["mu"]])
  row_id1 <- pair_cache$row_id1
  row_id2 <- pair_cache$row_id2
  n_pairs <- length(row_id1)
  has_zeta <- "zeta" %in% names(eta_inv)

  theta_idx <- if (length(eta_inv[["theta"]]) == n_obs) {
    row_id1
  } else {
    pair_cache$theta_index_map[row_id1]
  }

  par1 <- eta_inv[["theta"]][theta_idx]
  par2 <- if (has_zeta) eta_inv[["zeta"]][theta_idx] else rep(0, n_pairs)

  margin_p <- eta_inv[["margin_p_cache"]]
  u1 <- pmax(pmin(margin_p[row_id1], 1 - 1e-10), 1e-10)
  u2 <- pmax(pmin(margin_p[row_id2], 1 - 1e-10), 1e-10)

  par1_e <- par1
  par2_e <- par2
  par1_e[!is.finite(par1_e)] <- 0
  par2_e[!is.finite(par2_e)] <- 0
  if (copula_dist == "C") par1_e[par1_e >= 28] <- 27.9

  pair_ok <- pair_cache$observed_pair_base &
    is.finite(u1) & is.finite(u2) &
    is.finite(par1_e) & is.finite(par2_e)

  list(
    n_obs = n_obs,
    row_id1 = row_id1,
    row_id2 = row_id2,
    n_pairs = n_pairs,
    has_zeta = has_zeta,
    u1 = u1,
    u2 = u2,
    par1_e = par1_e,
    par2_e = par2_e,
    pair_ok = pair_ok
  )
}
