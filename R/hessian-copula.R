#' Assemble per-observation Hessian contributions from the copula

#'

#' For each pair (i,j) in a subject's longitudinal sequence, the copula

#' log-likelihood is log c(F(y_i), F(y_j); theta).  Its second derivative

#' w.r.t. natural parameters of the margin and copula is assembled here.

#'

#' Returns a list with named elements:

#'   cop_d2l_margin  - matrix (n_obs x n_margin_par x n_margin_par) flattened

#'                     as `list[[par1]][[par2]]`, each of length n_obs

#'   cop_d2l_theta   - length n_obs vector: d2 log c / d theta^2

#'   cop_d2l_zeta    - length n_obs vector or NULL

#'   cop_d2l_thetazeta - length n_obs vector or NULL

#'   cop_d2l_margin_theta - parameter-named list of length n_obs vectors:
#'                     d2logc/(dpar dtheta)

#'

#' @keywords internal

.calc_copula_hessian_contributions <- function(

    eta_inv, pair_cache, copula_dist,

    dF_list,    # output of .calc_dFdpar

    d2F_list,   # output of .calc_d2Fdpar2

    d2F_cross   # output of .calc_d2Fdpar_cross

) {

  n_obs <- length(eta_inv[["mu"]])

  margin_pars <- names(dF_list)

  has_zeta <- "zeta" %in% names(eta_inv)


  row_id1 <- pair_cache$row_id1

  row_id2 <- pair_cache$row_id2

  pair_complete <- pair_cache$observed_pair_base

  n_pairs <- length(row_id1)


  # Retrieve copula parameters for each pair

  theta_idx <- if (length(eta_inv[["theta"]]) == n_obs) row_id1 else pair_cache$theta_index_map[row_id1]

  par1 <- eta_inv[["theta"]][theta_idx]

  par2 <- if (has_zeta) eta_inv[["zeta"]][theta_idx] else rep(0, n_pairs)


  margin_p <- eta_inv[["margin_p_cache"]]  # pre-computed by caller


  u1 <- pmax(pmin(margin_p[row_id1], 1 - 1e-10), 1e-10)

  u2 <- pmax(pmin(margin_p[row_id2], 1 - 1e-10), 1e-10)


  # Clamp copula parameters

  par1_e <- par1; par2_e <- par2

  par1_e[!is.finite(par1_e)] <- 0

  par2_e[!is.finite(par2_e)] <- 0

  if (copula_dist == "C") par1_e[par1_e >= 28] <- 27.9


  fam_num <- copula_dist


  if (n_pairs == 0) {

    empty <- rep(0, n_obs)

    cop_d2l_margin <- vector("list", length(margin_pars))

    names(cop_d2l_margin) <- margin_pars

    for (pn1 in margin_pars) {

      cop_d2l_margin[[pn1]] <- vector("list", length(margin_pars))

      names(cop_d2l_margin[[pn1]]) <- margin_pars

      for (pn2 in margin_pars) cop_d2l_margin[[pn1]][[pn2]] <- empty

    }

    return(list(

      cop_d1l_margin = setNames(lapply(margin_pars, function(p) empty), margin_pars),

      cop_d2l_margin = cop_d2l_margin,

      cop_d2l_theta = empty,

      cop_d2l_zeta = if (has_zeta) empty else NULL,

      cop_d2l_thetazeta = if (has_zeta) empty else NULL,

      cop_d2l_margin_theta = setNames(lapply(margin_pars, function(p) empty), margin_pars),

      cop_d2l_margin_zeta = if (has_zeta) setNames(lapply(margin_pars, function(p) empty), margin_pars) else NULL

    ))

  }


  pair_ok <- pair_complete & is.finite(u1) & is.finite(u2) & is.finite(par1_e) & is.finite(par2_e)


  # --- Copula density and derivatives (all vectorised, one pass) ---

  c_val <- .copula_pdf(u1, u2, family = fam_num, par = par1_e, par2 = par2_e)

  c_val[!is.finite(c_val) | c_val <= 0] <- 1


  dcdu1  <- .copula_deriv(u1, u2, fam_num, par1_e, par2_e, deriv = "u1",   log = FALSE)

  dcdu2  <- .copula_deriv(u1, u2, fam_num, par1_e, par2_e, deriv = "u2",   log = FALSE)

  dcdth  <- .copula_deriv(u1, u2, fam_num, par1_e, par2_e, deriv = "par",  log = FALSE)

  dldth  <- .copula_deriv(u1, u2, fam_num, par1_e, par2_e, deriv = "par",  log = TRUE)


  d2cdu1_2   <- .copula_deriv2(u1, u2, fam_num, par1_e, par2_e, deriv = "u1")

  d2cdu2_2   <- .copula_deriv2(u1, u2, fam_num, par1_e, par2_e, deriv = "u2")

  d2cdth2    <- .copula_deriv2(u1, u2, fam_num, par1_e, par2_e, deriv = "par")

  # d2c/(du1 du2) not implemented in BiCopDeriv2; use 4-point FD on copula density only

  hu <- 1e-5

  u1c <- pmax(pmin(u1, 1 - hu), hu); u2c <- pmax(pmin(u2, 1 - hu), hu)

  c_pp <- .copula_pdf(u1c + hu, u2c + hu, fam_num, par1_e, par2_e)

  c_pm <- .copula_pdf(u1c + hu, u2c - hu, fam_num, par1_e, par2_e)

  c_mp <- .copula_pdf(u1c - hu, u2c + hu, fam_num, par1_e, par2_e)

  c_mm <- .copula_pdf(u1c - hu, u2c - hu, fam_num, par1_e, par2_e)

  d2cdu1u2 <- (c_pp - c_pm - c_mp + c_mm) / (4 * hu^2)

  d2cdu1u2[!is.finite(d2cdu1u2)] <- 0

  # d2c/(du1 dtheta) and d2c/(du2 dtheta): FD on BiCopDeriv w.r.t. par

  hth <- 1e-5

  dcdu1_p <- .copula_deriv(u1, u2, fam_num, par1_e + hth, par2_e, deriv = "u1", log = FALSE)

  dcdu1_m <- .copula_deriv(u1, u2, fam_num, par1_e - hth, par2_e, deriv = "u1", log = FALSE)

  d2cdthu1 <- (dcdu1_p - dcdu1_m) / (2 * hth)

  d2cdthu1[!is.finite(d2cdthu1)] <- 0

  dcdu2_p <- .copula_deriv(u1, u2, fam_num, par1_e + hth, par2_e, deriv = "u2", log = FALSE)

  dcdu2_m <- .copula_deriv(u1, u2, fam_num, par1_e - hth, par2_e, deriv = "u2", log = FALSE)

  d2cdthu2 <- (dcdu2_p - dcdu2_m) / (2 * hth)

  d2cdthu2[!is.finite(d2cdthu2)] <- 0


  # Zero out incomplete pairs

  .z <- function(x) { x[!pair_ok] <- 0; x[!is.finite(x)] <- 0; x }

  dcdu1 <- .z(dcdu1); dcdu2 <- .z(dcdu2); dcdth <- .z(dcdth); dldth <- .z(dldth)

  d2cdu1_2 <- .z(d2cdu1_2); d2cdu2_2 <- .z(d2cdu2_2); d2cdth2 <- .z(d2cdth2)

  d2cdu1u2 <- .z(d2cdu1u2); d2cdthu1 <- .z(d2cdthu1); d2cdthu2 <- .z(d2cdthu2)

  c_val[!pair_ok] <- 1


  # d2 log c / d theta^2 = (c * d2c/dth^2 - (dc/dth)^2) / c^2

  d2ldth2_pair <- (c_val * d2cdth2 - dcdth^2) / c_val^2

  d2ldth2_pair <- .z(d2ldth2_pair)


  # Zeta (t-copula second parameter)

  if (has_zeta) {

    dcdz   <- .copula_deriv(u1, u2, fam_num, par1_e, par2_e, deriv = "par2", log = FALSE)

    d2cdz2 <- .copula_deriv2(u1, u2, fam_num, par1_e, par2_e, deriv = "par2")

    d2cdthdz <- .copula_deriv2(u1, u2, fam_num, par1_e, par2_e, deriv = "par1par2")

    hz <- 1e-5 * pmax(1, abs(par2_e))

    par2_p <- par2_e + hz

    par2_m <- pmax(par2_e - hz, 2 + 1e-8)

    hz_eff <- par2_p - par2_m

    dcdu1_zp <- .copula_deriv(u1, u2, fam_num, par1_e, par2_p, deriv = "u1", log = FALSE)

    dcdu1_zm <- .copula_deriv(u1, u2, fam_num, par1_e, par2_m, deriv = "u1", log = FALSE)

    d2cdzu1 <- (dcdu1_zp - dcdu1_zm) / hz_eff

    dcdu2_zp <- .copula_deriv(u1, u2, fam_num, par1_e, par2_p, deriv = "u2", log = FALSE)

    dcdu2_zm <- .copula_deriv(u1, u2, fam_num, par1_e, par2_m, deriv = "u2", log = FALSE)

    d2cdzu2 <- (dcdu2_zp - dcdu2_zm) / hz_eff

    dcdz <- .z(dcdz); d2cdz2 <- .z(d2cdz2); d2cdthdz <- .z(d2cdthdz)

    d2cdzu1 <- .z(d2cdzu1); d2cdzu2 <- .z(d2cdzu2)

    dldz_pair <- .z(dcdz / c_val)

    d2ldz2_pair <- (c_val * d2cdz2 - dcdz^2) / c_val^2

    d2ldthdz_pair <- (d2cdthdz * c_val - dcdth * dcdz) / c_val^2

    d2ldz2_pair <- .z(d2ldz2_pair)

    d2ldthdz_pair <- .z(d2ldthdz_pair)

  }


  # --- Scatter pair-level contributions back to observation level ---

  # Each observation may appear as row_id1 (u1 role) or row_id2 (u2 role)

  # in potentially multiple pairs.


  # Helper: scatter pair vector to obs level

  .scatter <- function(pair_vec_id1, pair_vec_id2) {

    out <- numeric(n_obs)

    for (k in seq_along(row_id1)) {

      if (!pair_ok[k]) next

      out[row_id1[k]] <- out[row_id1[k]] + pair_vec_id1[k]

      out[row_id2[k]] <- out[row_id2[k]] + pair_vec_id2[k]

    }

    out

  }


  # Theta Hessian: theta is a per-pair parameter. Each pair contributes once.

  # We accumulate to row_id1 only (the "from" observation), which carries the

  # theta design row. The assembly step uses X_theta[row_id1] to get the

  # correct design column, so we must keep contributions pair-indexed.

  #

  # Full chain rule (Fisher-z link: theta = tanh(eta)):

  #   d2 logc / d eta^2 = d2 logc / d theta^2 * (d theta/d eta)^2

  #                     + d logc / d theta    * (d2 theta/d eta^2)

  # where d theta/d eta = 1 - theta^2  (dr)

  #       d2 theta/d eta^2 = -2*theta*(1-theta^2)

  # Store in eta space so assembly does NOT multiply by dr^2 again.

  dr_pair        <- (1 - par1_e^2)                     # d theta/d eta per pair

  d2_linkinv_par <- -2 * par1_e * (1 - par1_e^2)       # d2 theta/d eta^2 per pair

  cop_d2l_theta_pair <- d2ldth2_pair * dr_pair^2 + dldth * d2_linkinv_par

  cop_d2l_theta_pair <- .z(cop_d2l_theta_pair)         # zero invalid pairs

  cop_d2l_zeta_pair   <- NULL

  cop_d2l_thetazeta_pair <- NULL

  if (has_zeta) {

    # t-copula zeta uses log_2plus_inv(eta) = exp(eta) + 2.

    # Store zeta curvature in eta space so assembly does NOT multiply by dr^2 again.

    zeta_dr_pair <- pmax(par2_e - 2, 0)

    zeta_d2_linkinv_pair <- zeta_dr_pair

    cop_d2l_zeta_pair <- d2ldz2_pair * zeta_dr_pair^2 + dldz_pair * zeta_d2_linkinv_pair

    cop_d2l_zeta_pair <- .z(cop_d2l_zeta_pair)

    cop_d2l_thetazeta_pair <- d2ldthdz_pair

  }


  # --- Margin parameter contributions ---

  # For obs i appearing as u1 in a pair:

  #   d log c / d par_i = (dc/du1) * (dF/dpar)_i / c

  #   d2 log c / d par_i^2 =

  #       (d2c/du1^2) * (dF/dpar)_i^2 / c

  #     + (dc/du1) * (d2F/dpar^2)_i / c

  #     - (dc/du1)^2 * (dF/dpar)_i^2 / c^2

  #   d2 log c / (d par_i d par_j) [i in u1, j in u1]:

  #       (d2c/du1^2) * dF1/dpi * dF1/dpj / c

  #     + (dc/du1) * d2F/(dpi dpj) / c

  #     - (dc/du1 * dF1/dpi/c) * (dc/du1 * dF1/dpj/c)

  #   cross u1 x u2 (obs i in u1, obs j in u2, different subjects, skip for now):

  #       (d2c/du1du2) * dF1/dpi * dF2/dpj / c  -- but i != j so zero for same obs

  #

  # For obs i appearing as u2: symmetric with u1 <-> u2.

  #

  # Cross pair i (u1) vs j (u2) where i == j (same observation in both pair slots):

  # This would require a subject to have the same time twice - impossible by data checks.

  # So diagonal obs-level hessian from copula is sum over all pairs the obs appears in.


  cop_d2l_margin <- vector("list", length(margin_pars))

  names(cop_d2l_margin) <- margin_pars

  for (pn in margin_pars) {

    cop_d2l_margin[[pn]] <- vector("list", length(margin_pars))

    names(cop_d2l_margin[[pn]]) <- margin_pars

    for (pn2 in margin_pars) cop_d2l_margin[[pn]][[pn2]] <- numeric(n_obs)

  }

  cop_d1l_margin <- setNames(lapply(margin_pars, function(p) numeric(n_obs)), margin_pars)


  # Pair-indexed marginÃƒÆ’Ã¢â‚¬â€theta cross-terms: d2logc/(d par_pn_i1 d theta_k) for u1 and u2 roles.

  # Stored separately to allow correct theta-row-indexed assembly.

  n_pairs <- length(row_id1)

  cop_d2l_margin_theta_u1 <- setNames(lapply(margin_pars, function(p) numeric(n_pairs)), margin_pars)

  cop_d2l_margin_theta_u2 <- setNames(lapply(margin_pars, function(p) numeric(n_pairs)), margin_pars)

  cop_d2l_margin_zeta_u1 <- if (has_zeta) setNames(lapply(margin_pars, function(p) numeric(n_pairs)), margin_pars) else NULL

  cop_d2l_margin_zeta_u2 <- if (has_zeta) setNames(lapply(margin_pars, function(p) numeric(n_pairs)), margin_pars) else NULL


  for (k in seq_along(row_id1)) {

    if (!pair_ok[k]) next

    i1 <- row_id1[k]; i2 <- row_id2[k]

    cv <- c_val[k]


    for (pi_idx in seq_along(margin_pars)) {

      pni <- margin_pars[pi_idx]

      dFi1 <- dF_list[[pni]][i1]   # obs i1 role u1

      dFi2 <- dF_list[[pni]][i2]   # obs i2 role u2

      d2Fi1 <- d2F_list[[pni]][i1]

      d2Fi2 <- d2F_list[[pni]][i2]


      # Contribution of obs i1 (u1 slot) to its own diagonal Hessian entry

      d1_i1 <- dcdu1[k] * dFi1 / cv

      d1_i2 <- dcdu2[k] * dFi2 / cv

      d1_i1[!is.finite(d1_i1)] <- 0

      d1_i2[!is.finite(d1_i2)] <- 0

      cop_d1l_margin[[pni]][i1] <- cop_d1l_margin[[pni]][i1] + d1_i1

      cop_d1l_margin[[pni]][i2] <- cop_d1l_margin[[pni]][i2] + d1_i2


      diag_i1 <- (d2cdu1_2[k] * dFi1^2 + dcdu1[k] * d2Fi1) / cv -

                 (dcdu1[k] * dFi1 / cv)^2

      # Contribution of obs i2 (u2 slot) to its own diagonal Hessian entry

      diag_i2 <- (d2cdu2_2[k] * dFi2^2 + dcdu2[k] * d2Fi2) / cv -

                 (dcdu2[k] * dFi2 / cv)^2


      diag_i1[!is.finite(diag_i1)] <- 0

      diag_i2[!is.finite(diag_i2)] <- 0


      cop_d2l_margin[[pni]][[pni]][i1] <- cop_d2l_margin[[pni]][[pni]][i1] + diag_i1

      cop_d2l_margin[[pni]][[pni]][i2] <- cop_d2l_margin[[pni]][[pni]][i2] + diag_i2


      # Cross u1 x u2: d2 log c / (d par_i1 d par_i2) where i1 != i2

      # This is a cross-observation term. It contributes to H[beta_{pni,i1}, beta_{pni,i2}]

      # but since our design is X_a beta_a (same beta vector), the Hessian element

      # H[a,b] = sum_k X_a[k] X_b[k] d2l/d eta_a d eta_b ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â here k is observation.

      # Cross-obs copula second derivatives contribute off-diagonal blocks involving

      # different observations and are handled separately in .assemble_covariate_hessian.

      # We store the cross-pair contributions in a separate structure.


      for (pj_idx in seq_along(margin_pars)) {

        pnj <- margin_pars[pj_idx]

        if (pj_idx <= pi_idx) next  # upper triangle only, fill both below

        dFj1 <- dF_list[[pnj]][i1]

        dFj2 <- dF_list[[pnj]][i2]

        d2Fij_cross1 <- if (!is.null(d2F_cross[[pni]][[pnj]])) d2F_cross[[pni]][[pnj]][i1] else 0

        d2Fij_cross2 <- if (!is.null(d2F_cross[[pni]][[pnj]])) d2F_cross[[pni]][[pnj]][i2] else 0


        cross_i1 <- (d2cdu1_2[k] * dFi1 * dFj1 + dcdu1[k] * d2Fij_cross1) / cv -

                    (dcdu1[k] * dFi1 / cv) * (dcdu1[k] * dFj1 / cv)

        cross_i2 <- (d2cdu2_2[k] * dFi2 * dFj2 + dcdu2[k] * d2Fij_cross2) / cv -

                    (dcdu2[k] * dFi2 / cv) * (dcdu2[k] * dFj2 / cv)


        cross_i1[!is.finite(cross_i1)] <- 0

        cross_i2[!is.finite(cross_i2)] <- 0


        cop_d2l_margin[[pni]][[pnj]][i1] <- cop_d2l_margin[[pni]][[pnj]][i1] + cross_i1

        cop_d2l_margin[[pnj]][[pni]][i1] <- cop_d2l_margin[[pnj]][[pni]][i1] + cross_i1

        cop_d2l_margin[[pni]][[pnj]][i2] <- cop_d2l_margin[[pni]][[pnj]][i2] + cross_i2

        cop_d2l_margin[[pnj]][[pni]][i2] <- cop_d2l_margin[[pnj]][[pni]][i2] + cross_i2

      }


      # d2 log c / (d par_margin d theta): chain rule

      # = (d2c / (du1 d theta) * dF/dpar) / c - (dc/du1 * dF/dpar / c) * (dc/dtheta / c)

      # for obs in u1 role; symmetric for u2

      d2logc_dpar_dtheta_i1 <- d2cdthu1[k] * dFi1 / cv - (dcdu1[k] * dFi1 / cv) * (dcdth[k] / cv)

      d2logc_dpar_dtheta_i2 <- d2cdthu2[k] * dFi2 / cv - (dcdu2[k] * dFi2 / cv) * (dcdth[k] / cv)

      d2logc_dpar_dtheta_i1[!is.finite(d2logc_dpar_dtheta_i1)] <- 0

      d2logc_dpar_dtheta_i2[!is.finite(d2logc_dpar_dtheta_i2)] <- 0


      cop_d2l_margin_theta_u1[[pni]][k] <- d2logc_dpar_dtheta_i1

      cop_d2l_margin_theta_u2[[pni]][k] <- d2logc_dpar_dtheta_i2


      if (has_zeta) {

        d2logc_dpar_dzeta_i1 <- d2cdzu1[k] * dFi1 / cv - (dcdu1[k] * dFi1 / cv) * (dcdz[k] / cv)

        d2logc_dpar_dzeta_i2 <- d2cdzu2[k] * dFi2 / cv - (dcdu2[k] * dFi2 / cv) * (dcdz[k] / cv)

        d2logc_dpar_dzeta_i1[!is.finite(d2logc_dpar_dzeta_i1)] <- 0

        d2logc_dpar_dzeta_i2[!is.finite(d2logc_dpar_dzeta_i2)] <- 0


        cop_d2l_margin_zeta_u1[[pni]][k] <- d2logc_dpar_dzeta_i1

        cop_d2l_margin_zeta_u2[[pni]][k] <- d2logc_dpar_dzeta_i2

      }

    }

  }


  # Store cross-pair (i1, i2) contributions for off-diagonal assembly

  # These are d2 log c(F_i1, F_i2) / (d par_i1_pn1 d par_i2_pn2)

  # = (d2c/du1du2) * dF_i1/dpn1 * dF_i2/dpn2 / c  (for pn1 == pn2 case the main one)

  # Returned as a list of pair-indexed vectors for use in .assemble_covariate_hessian.

  cross_pair_contributions <- vector("list", length(margin_pars))

  names(cross_pair_contributions) <- margin_pars

  for (pn in margin_pars) {

    cross_pair_contributions[[pn]] <- vector("list", length(margin_pars))

    names(cross_pair_contributions[[pn]]) <- margin_pars

    for (pn2 in margin_pars) {

      # d2 log c / (d pn at i1) (d pn2 at i2)

      # = d2c/(du1 du2) * dF_i1/dpn * dF_i2/dpn2 / c - (dc/du1*dF_i1/dpn/c)*(dc/du2*dF_i2/dpn2/c)

      dF_pn_at_id1  <- dF_list[[pn]][row_id1]

      dF_pn2_at_id2 <- dF_list[[pn2]][row_id2]

      cross_vals <- d2cdu1u2 * dF_pn_at_id1 * dF_pn2_at_id2 / c_val -

                    (dcdu1 * dF_pn_at_id1 / c_val) * (dcdu2 * dF_pn2_at_id2 / c_val)

      cross_vals[!pair_ok] <- 0

      cross_vals[!is.finite(cross_vals)] <- 0

      cross_pair_contributions[[pn]][[pn2]] <- cross_vals

    }

  }


  list(

    cop_d1l_margin          = cop_d1l_margin,

    cop_d2l_margin         = cop_d2l_margin,

    cop_d2l_theta          = cop_d2l_theta_pair,       # pair-indexed

    cop_d2l_zeta           = cop_d2l_zeta_pair,        # pair-indexed or NULL

    cop_d2l_thetazeta      = cop_d2l_thetazeta_pair,   # pair-indexed or NULL

    cop_d2l_margin_theta_u1 = cop_d2l_margin_theta_u1,  # pair-indexed, u1 role

    cop_d2l_margin_theta_u2 = cop_d2l_margin_theta_u2,  # pair-indexed, u2 role

    cop_d2l_margin_zeta_u1 = cop_d2l_margin_zeta_u1,    # pair-indexed, u1 role or NULL

    cop_d2l_margin_zeta_u2 = cop_d2l_margin_zeta_u2,    # pair-indexed, u2 role or NULL

    cross_pair_contribs    = cross_pair_contributions,

    row_id1 = row_id1, row_id2 = row_id2, pair_ok = pair_ok

  )

}



# ------------------------------------------------------------

# 3.  Assemble covariate-level Hessian

# ------------------------------------------------------------


