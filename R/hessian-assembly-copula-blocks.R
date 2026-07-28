#' @keywords internal
#' @noRd
.hessian_assembly_copula_parameter_block <- function(
    pa,
    pb,
    xa,
    xb,
    dra,
    drb,
    copula_hess,
    pair_cache,
    row_id1,
    pair_ok) {
  th_idx <- pair_cache$theta_index_map[row_id1[pair_ok]]
  valid <- !is.na(th_idx)

  if (pa == "theta" && pb == "theta") {
    d2_pairs <- copula_hess$cop_d2l_theta[pair_ok]
    return(sum(
      xa[th_idx[valid]] * xb[th_idx[valid]] * d2_pairs[valid],
      na.rm = TRUE
    ))
  }

  if (!is.null(copula_hess$cop_d2l_zeta) && pa == "zeta" && pb == "zeta") {
    d2_pairs <- copula_hess$cop_d2l_zeta[pair_ok]
    return(sum(
      xa[th_idx[valid]] * xb[th_idx[valid]] * d2_pairs[valid],
      na.rm = TRUE
    ))
  }

  if (
    !is.null(copula_hess$cop_d2l_thetazeta) &&
      ((pa == "theta" && pb == "zeta") || (pa == "zeta" && pb == "theta"))
  ) {
    d2_pairs <- copula_hess$cop_d2l_thetazeta[pair_ok]
    return(sum(
      xa[th_idx[valid]] * xb[th_idx[valid]] * d2_pairs[valid] *
        dra[th_idx[valid]] * drb[th_idx[valid]],
      na.rm = TRUE
    ))
  }

  NULL
}
