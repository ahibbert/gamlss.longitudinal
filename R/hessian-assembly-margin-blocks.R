# Compare distribution-parameter values, independent of coefficient names.
.hessian_assembly_margin_margin_block <- function(
    pa, pb, xa, xb, dra, drb, margin_d1l, margin_d2l, copula_hess,
    eta_dr, eta_d2, pair_ok, id1_ok, id2_ok, margin_mu_eta_d2 = NULL) {
  if (!is.null(margin_mu_eta_d2) && identical(unname(pa), "mu") &&
      identical(unname(pb), "mu")) {
    n <- length(dra)
    check_row <- function(x) {
      if (!is.numeric(x) || !is.null(dim(x)) || length(x) != n || any(!is.finite(x)))
        stop("BCPE/log mu assembly requires finite, equally sized row vectors.", call. = FALSE)
      x
    }
    # Replace the entire marginal chain rule, never subtract a family score.
    # The copula still uses its original natural derivatives and inverse link.
    check_row(xa); check_row(xb); check_row(dra); check_row(drb)
    d2l_obs <- check_row(margin_mu_eta_d2)
    cop_d2 <- copula_hess$cop_d2l_margin$mu$mu
    cop_d1 <- copula_hess$cop_d1l_margin$mu
    if (!is.null(cop_d2)) d2l_obs <- d2l_obs + check_row(cop_d2) * dra * drb
    if (!is.null(cop_d1)) d2l_obs <- d2l_obs + check_row(cop_d1) * check_row(eta_d2$mu)
    check_row(d2l_obs)
    row_contrib <- xa * xb * d2l_obs
    if (any(!is.finite(row_contrib)))
      stop("BCPE/log mu row products are not finitely representable.", call. = FALSE)
    H_ab <- sum(row_contrib)
    if (!is.finite(H_ab))
      stop("BCPE/log mu row sum is not finitely representable.", call. = FALSE)
  } else {
    d2l_obs <- .hessian_assembly_margin_d2_obs(pa, pb, margin_d2l, copula_hess, eta_dr) * dra * drb
    if (identical(unname(pa), unname(pb)) && pa %in% names(eta_d2)) {
      d2l_obs <- d2l_obs +
        .hessian_assembly_margin_d1_obs(pa, margin_d1l, copula_hess, eta_dr) * as.numeric(eta_d2[[pa]])
    }
    H_ab <- sum(xa * xb * d2l_obs, na.rm = TRUE)
  }
  if (pa %in% names(copula_hess$cross_pair_contribs) &&
    pb %in% names(copula_hess$cross_pair_contribs[[pa]])) {
    cross_ab <- copula_hess$cross_pair_contribs[[pa]][[pb]][pair_ok]
    cross_ba <- copula_hess$cross_pair_contribs[[pb]][[pa]][pair_ok]
    H_ab <- H_ab + sum(
      xa[id1_ok] * dra[id1_ok] * cross_ab * xb[id2_ok] * drb[id2_ok] +
        xa[id2_ok] * dra[id2_ok] * cross_ba * xb[id1_ok] * drb[id1_ok],
      na.rm = TRUE
    )
  }
  H_ab
}
