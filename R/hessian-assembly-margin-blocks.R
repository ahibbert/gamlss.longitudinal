#' Assemble margin-by-margin Hessian blocks
#'
#' Combines observation-level marginal curvature with same-pair copula
#' contributions and cross-pair copula terms.
#'
#' @keywords internal
#' @noRd
.hessian_assembly_margin_margin_block <- function(
    pa,
    pb,
    xa,
    xb,
    dra,
    drb,
    margin_d1l,
    margin_d2l,
    copula_hess,
    eta_dr,
    eta_d2,
    pair_ok,
    id1_ok,
    id2_ok) {
  d2l_obs <- .hessian_assembly_margin_d2_obs(pa, pb, margin_d2l, copula_hess, eta_dr) * dra * drb

  if (identical(pa, pb) && pa %in% names(eta_d2)) {
    d2l_obs <- d2l_obs +
      .hessian_assembly_margin_d1_obs(pa, margin_d1l, copula_hess, eta_dr) * as.numeric(eta_d2[[pa]])
  }

  H_ab <- sum(xa * xb * d2l_obs, na.rm = TRUE)

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
