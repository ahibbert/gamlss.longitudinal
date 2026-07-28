#' Margin-side copula Hessian contribution helpers
#'
#' These helpers accumulate pair-level copula derivative quantities back to
#' margin-parameter structures used by covariate Hessian assembly.
#'
#' @noRd
.copula_hessian_margin_terms <- function(
    margin_pars, n_obs, row_id1, row_id2, pair_ok,
    dF_list, d2F_list, d2F_cross, derivs, has_zeta) {
  n_pairs <- length(row_id1)
  accumulators <- .copula_hessian_margin_accumulators(margin_pars, n_obs, n_pairs, has_zeta)
  cop_d1l_margin <- accumulators$cop_d1l_margin
  cop_d2l_margin <- accumulators$cop_d2l_margin
  cop_d2l_margin_theta_u1 <- accumulators$cop_d2l_margin_theta_u1
  cop_d2l_margin_theta_u2 <- accumulators$cop_d2l_margin_theta_u2
  cop_d2l_margin_zeta_u1 <- accumulators$cop_d2l_margin_zeta_u1
  cop_d2l_margin_zeta_u2 <- accumulators$cop_d2l_margin_zeta_u2

  for (k in seq_along(row_id1)) {
    if (!pair_ok[k]) next
    i1 <- row_id1[k]
    i2 <- row_id2[k]
    cv <- derivs$c_val[k]

    for (pi_idx in seq_along(margin_pars)) {
      pni <- margin_pars[pi_idx]
      dFi1 <- dF_list[[pni]][i1]
      dFi2 <- dF_list[[pni]][i2]
      d2Fi1 <- d2F_list[[pni]][i1]
      d2Fi2 <- d2F_list[[pni]][i2]

      d1_pair <- .copula_hessian_margin_d1_pair(derivs, k, dFi1, dFi2, cv)
      cop_d1l_margin[[pni]][i1] <- cop_d1l_margin[[pni]][i1] + d1_pair$i1
      cop_d1l_margin[[pni]][i2] <- cop_d1l_margin[[pni]][i2] + d1_pair$i2

      diag_pair <- .copula_hessian_margin_diag_pair(derivs, k, dFi1, dFi2, d2Fi1, d2Fi2, cv)
      cop_d2l_margin[[pni]][[pni]][i1] <- cop_d2l_margin[[pni]][[pni]][i1] + diag_pair$i1
      cop_d2l_margin[[pni]][[pni]][i2] <- cop_d2l_margin[[pni]][[pni]][i2] + diag_pair$i2

      for (pj_idx in seq_along(margin_pars)) {
        pnj <- margin_pars[pj_idx]
        if (pj_idx <= pi_idx) next
        dFj1 <- dF_list[[pnj]][i1]
        dFj2 <- dF_list[[pnj]][i2]
        d2Fij_cross1 <- if (!is.null(d2F_cross[[pni]][[pnj]])) {
          d2F_cross[[pni]][[pnj]][i1]
        } else {
          0
        }
        d2Fij_cross2 <- if (!is.null(d2F_cross[[pni]][[pnj]])) {
          d2F_cross[[pni]][[pnj]][i2]
        } else {
          0
        }

        cross_pair <- .copula_hessian_margin_cross_pair(
          derivs, k, dFi1, dFi2, dFj1, dFj2, d2Fij_cross1, d2Fij_cross2, cv
        )

        cop_d2l_margin[[pni]][[pnj]][i1] <- cop_d2l_margin[[pni]][[pnj]][i1] + cross_pair$i1
        cop_d2l_margin[[pnj]][[pni]][i1] <- cop_d2l_margin[[pnj]][[pni]][i1] + cross_pair$i1
        cop_d2l_margin[[pni]][[pnj]][i2] <- cop_d2l_margin[[pni]][[pnj]][i2] + cross_pair$i2
        cop_d2l_margin[[pnj]][[pni]][i2] <- cop_d2l_margin[[pnj]][[pni]][i2] + cross_pair$i2
      }

      theta_pair <- .copula_hessian_margin_theta_pair(derivs, k, dFi1, dFi2, cv)
      cop_d2l_margin_theta_u1[[pni]][k] <- theta_pair$i1
      cop_d2l_margin_theta_u2[[pni]][k] <- theta_pair$i2

      if (has_zeta) {
        zeta_pair <- .copula_hessian_margin_zeta_pair(derivs, k, dFi1, dFi2, cv)
        cop_d2l_margin_zeta_u1[[pni]][k] <- zeta_pair$i1
        cop_d2l_margin_zeta_u2[[pni]][k] <- zeta_pair$i2
      }
    }
  }

  list(
    cop_d1l_margin = cop_d1l_margin,
    cop_d2l_margin = cop_d2l_margin,
    cop_d2l_margin_theta_u1 = cop_d2l_margin_theta_u1,
    cop_d2l_margin_theta_u2 = cop_d2l_margin_theta_u2,
    cop_d2l_margin_zeta_u1 = cop_d2l_margin_zeta_u1,
    cop_d2l_margin_zeta_u2 = cop_d2l_margin_zeta_u2
  )
}

#' @noRd
.copula_hessian_cross_pair_terms <- function(
    margin_pars, row_id1, row_id2, pair_ok, dF_list, derivs) {
  cross_pair_contributions <- vector("list", length(margin_pars))
  names(cross_pair_contributions) <- margin_pars
  for (pn in margin_pars) {
    cross_pair_contributions[[pn]] <- vector("list", length(margin_pars))
    names(cross_pair_contributions[[pn]]) <- margin_pars
    for (pn2 in margin_pars) {
      dF_pn_at_id1 <- dF_list[[pn]][row_id1]
      dF_pn2_at_id2 <- dF_list[[pn2]][row_id2]
      cross_vals <- derivs$d2cdu1u2 * dF_pn_at_id1 * dF_pn2_at_id2 / derivs$c_val -
        (derivs$dcdu1 * dF_pn_at_id1 / derivs$c_val) *
          (derivs$dcdu2 * dF_pn2_at_id2 / derivs$c_val)
      cross_vals[!pair_ok] <- 0
      cross_vals[!is.finite(cross_vals)] <- 0
      cross_pair_contributions[[pn]][[pn2]] <- cross_vals
    }
  }

  cross_pair_contributions
}
