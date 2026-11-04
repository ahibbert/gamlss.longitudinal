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
  raw_u1 <- margin_p[row_id1]
  raw_u2 <- margin_p[row_id2]
  u1 <- .copula_clamp01(raw_u1)
  u2 <- .copula_clamp01(raw_u2)

  included <- pair_cache$observed_pair_base
  invalid <- !is.finite(raw_u1) | !is.finite(raw_u2) |
    raw_u1 < 0 | raw_u1 > 1 | raw_u2 < 0 | raw_u2 > 1 |
    !is.finite(par1) | !is.finite(par2)
  if (any(included & invalid)) {
    stop("Invalid probability or parameter for an included copula pair.", call. = FALSE)
  }
  pair_ok <- included

  if (identical(copula_dist, "G") && !is.null(eta_inv$theta_link_cache) &&
      any(pair_ok & eta_inv$theta_link_cache[theta_idx] == log(16))) {
    stop("Dependence is at a nondifferentiable Gumbel link boundary.", call. = FALSE)
  }

  lower <- 1e-12
  upper <- 1 - 1e-12
  if (any(pair_ok & (raw_u1 == lower | raw_u1 == upper |
      raw_u2 == lower | raw_u2 == upper))) {
    stop("A marginal CDF is at a nondifferentiable copula probability clamp boundary.", call. = FALSE)
  }
  u1_active <- is.finite(raw_u1) & raw_u1 > lower & raw_u1 < upper
  u2_active <- is.finite(raw_u2) & raw_u2 > lower & raw_u2 < upper
  theta_active <- zeta_active <- rep(TRUE, n_pairs)
  if (identical(.copula_backend(), "native")) {
    if (copula_dist %in% c("N", "t")) {
      if (any(pair_ok & abs(par1) == 0.999999)) {
        stop("Correlation is at a nondifferentiable native copula clamp boundary.", call. = FALSE)
      }
      theta_active <- is.finite(par1) & abs(par1) < 0.999999
    }
    if (has_zeta) {
      if (any(pair_ok & par2 == 2.000001)) {
        stop("Student-t df is at a nondifferentiable native lower clamp boundary.", call. = FALSE)
      }
      zeta_active <- is.finite(par2) & par2 > 2.000001 & par2 < 1e6
    }
  }

  par1_e <- par1
  par2_e <- par2
  par1_e[!pair_ok] <- 0
  par2_e[!pair_ok] <- 0

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
    pair_ok = pair_ok,
    u1_active = u1_active,
    u2_active = u2_active,
    theta_active = theta_active,
    zeta_active = zeta_active
  )
}
