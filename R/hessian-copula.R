#' Assemble per-observation Hessian contributions from the copula
#'
#' For each pair (i,j) in a subject's longitudinal sequence, the copula
#' log-likelihood is log c(F(y_i), F(y_j); theta).  Its second derivative
#' w.r.t. natural parameters of the margin and copula is assembled here.
#'
#' Returns a list with named elements:
#'
#'   cop_d2l_margin  - matrix (n_obs x n_margin_par x n_margin_par) flattened
#'                     as `list[[par1]][[par2]]`, each of length n_obs
#'   cop_d2l_theta   - length n_obs vector: d2 log c / d theta^2
#'   cop_d2l_zeta    - length n_obs vector or NULL
#'   cop_d2l_thetazeta - length n_obs vector or NULL
#'   cop_d2l_margin_theta - parameter-named list of length n_obs vectors:
#'                     d2logc/(dpar dtheta)
#'
#' @keywords internal
#' @noRd
.calc_copula_hessian_contributions <- function(
    eta_inv, pair_cache, copula_dist,
    dF_list,
    d2F_list,
    d2F_cross) {
  n_obs <- length(eta_inv[["mu"]])
  margin_pars <- names(dF_list)

  pair_inputs <- .copula_hessian_pair_inputs(
    eta_inv = eta_inv,
    pair_cache = pair_cache,
    copula_dist = copula_dist
  )
  row_id1 <- pair_inputs$row_id1
  row_id2 <- pair_inputs$row_id2
  has_zeta <- pair_inputs$has_zeta

  if (pair_inputs$n_pairs == 0) {
    return(.copula_hessian_empty_result(n_obs, margin_pars, has_zeta))
  }

  derivs <- .copula_hessian_pair_derivatives(
    u1 = pair_inputs$u1,
    u2 = pair_inputs$u2,
    fam_num = copula_dist,
    par1_e = pair_inputs$par1_e,
    par2_e = pair_inputs$par2_e,
    pair_ok = pair_inputs$pair_ok,
    has_zeta = has_zeta
  )

  parameter_curvature <- .copula_hessian_parameter_curvature(
    derivs = derivs,
    par1_e = pair_inputs$par1_e,
    par2_e = pair_inputs$par2_e,
    pair_ok = pair_inputs$pair_ok,
    has_zeta = has_zeta
  )

  margin_terms <- .copula_hessian_margin_terms(
    margin_pars = margin_pars,
    n_obs = n_obs,
    row_id1 = row_id1,
    row_id2 = row_id2,
    pair_ok = pair_inputs$pair_ok,
    dF_list = dF_list,
    d2F_list = d2F_list,
    d2F_cross = d2F_cross,
    derivs = derivs,
    has_zeta = has_zeta
  )

  cross_pair_contributions <- .copula_hessian_cross_pair_terms(
    margin_pars = margin_pars,
    row_id1 = row_id1,
    row_id2 = row_id2,
    pair_ok = pair_inputs$pair_ok,
    dF_list = dF_list,
    derivs = derivs
  )

  list(
    cop_d1l_margin = margin_terms$cop_d1l_margin,
    cop_d2l_margin = margin_terms$cop_d2l_margin,
    cop_d2l_theta = parameter_curvature$cop_d2l_theta_pair,
    cop_d2l_zeta = parameter_curvature$cop_d2l_zeta_pair,
    cop_d2l_thetazeta = parameter_curvature$cop_d2l_thetazeta_pair,
    cop_d2l_margin_theta_u1 = margin_terms$cop_d2l_margin_theta_u1,
    cop_d2l_margin_theta_u2 = margin_terms$cop_d2l_margin_theta_u2,
    cop_d2l_margin_zeta_u1 = margin_terms$cop_d2l_margin_zeta_u1,
    cop_d2l_margin_zeta_u2 = margin_terms$cop_d2l_margin_zeta_u2,
    cross_pair_contribs = cross_pair_contributions,
    row_id1 = row_id1,
    row_id2 = row_id2,
    pair_ok = pair_inputs$pair_ok
  )
}
