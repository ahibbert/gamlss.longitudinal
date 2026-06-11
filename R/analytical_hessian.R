# analytical_hessian.R
# ============================================================
# Analytical / semi-analytical Hessian for gamlss.longitudinal
# ============================================================
# This file is sourced from common_functions.R when
# vcov.gamlss.longitudinal(..., method = "analytical") is called.
#
# The numerical path (calc_true_SE_numderiv_only_covariates) remains the
# source of truth and is used by default.
#
# Design summary
# --------------
# The covariate-level Hessian H[a,b] = sum_i d2l_i / (d beta_a d beta_b)
# is assembled in two steps:
#
#   Step 1 – per-pair / per-observation second derivatives of the joint
#             log-likelihood w.r.t. natural parameters (mu, sigma, ...,
#             theta, zeta). These use analytical formulas from gamlss
#             family objects and .copula_deriv2, plus a cheap
#             per-observation finite difference on F(y) only (not the full
#             likelihood).
#
#   Step 2 – chain rule through the link functions and design matrices:
#             H[a,b] = X_a' diag(d2l / d eta_a d eta_b) X_b
#             where d eta / d beta = X.
#
# Supported copulas: anything supported by .copula_deriv2.
# Supported margins: any gamlss family that exposes d2ldm2 etc.

#' Finite-difference second derivative of an inverse link
#'
#' Used in the chain-rule step of the analytical Hessian assembly when family
#' objects do not expose closed-form inverse-link curvature.
#'
#' @noRd
.calc_linkinv_second_derivative <- function(eta, linkinv_fun, h = 1e-5) {
  eta <- as.numeric(eta)
  step <- h * pmax(abs(eta), 1)
  plus <- linkinv_fun(eta + step)
  base <- linkinv_fun(eta)
  minus <- linkinv_fun(eta - step)
  out <- (plus - 2 * base + minus) / (step^2)
  out[!is.finite(out)] <- 0
  as.numeric(out)
}

#' Collect inverse-link second derivatives for all fitted parameters
#'
#' Returns per-parameter vectors for margin and copula inverse-link curvature.
#'
#' @noRd
.calc_eta_d2_linkinv <- function(eta, margin_dist, copula_link, h = 1e-5) {
  out <- vector("list", length(eta))
  names(out) <- names(eta)
  for (par_name in names(eta)) {
    linkinv_fun <- NULL
    if (par_name %in% names(margin_dist$parameters)) {
      linkinv_fun <- margin_dist[[paste(par_name, ".linkinv", sep = "")]]
    } else if (par_name %in% c("theta", "zeta")) {
      linkinv_fun <- copula_link[[paste(par_name, ".linkinv", sep = "")]]
    }
    out[[par_name]] <- if (is.function(linkinv_fun)) {
      .calc_linkinv_second_derivative(eta[[par_name]], linkinv_fun, h = h)
    } else {
      rep(0, length(eta[[par_name]]))
    }
  }
  out
}

#' Choose a finite-difference step on the natural parameter scale
#'
#' Identity-linked parameters use the base step; other links scale by parameter
#' magnitude to reduce avoidable cancellation.
#'
#' @noRd
.natural_fd_step <- function(par, par_name, margin_dist, h) {
  link_name <- margin_dist[[paste(par_name, "link", sep = ".")]]
  if (is.character(link_name) && identical(link_name[1], "identity")) {
    return(rep(h, length(par)))
  }
  h * pmax(abs(par), 1)
}

#' Warn when GG Hessian curvature is likely unstable near nu equals zero
#'
#' @noRd
.warn_gg_near_zero_nu <- function(margin_dist, eta_inv, threshold = 0.05) {
  family <- margin_dist$family[1]
  if (!identical(family, "GG") || !("nu" %in% names(eta_inv))) {
    return(invisible(FALSE))
  }

  nu <- as.numeric(eta_inv[["nu"]])
  min_abs_nu <- suppressWarnings(min(abs(nu[is.finite(nu)]), na.rm = TRUE))
  if (!is.finite(min_abs_nu) || min_abs_nu >= threshold) {
    return(invisible(FALSE))
  }

  warning(
    sprintf(
      paste(
        "Analytical Hessian for GG may be numerically unstable because fitted",
        "nu is close to 0 (min |nu| = %.4g). Consider vcov(..., method =",
        "\"numderiv\") if standard errors are important."
      ),
      min_abs_nu
    ),
    call. = FALSE
  )
  invisible(TRUE)
}

#' Warn when zero-heavy discrete margins may produce delicate Hessians
#'
#' @noRd
.warn_zero_heavy_discrete_hessian <- function(margin_dist, response, zero_threshold = 0.35) {
  family <- margin_dist$family[1]
  zero_inflated_families <- c(
    "ZIP", "ZIP2", "ZAP",
    "ZINBI", "ZINBII", "ZINBF",
    "ZAGA", "ZAIG", "ZALG"
  )
  count_families <- c(
    zero_inflated_families,
    "PO", "PIG", "NBI", "NBII", "DEL", "SICHEL", "SI", "DPO", "DNO"
  )
  if (!(family %in% count_families)) {
    return(invisible(FALSE))
  }

  observed <- response[is.finite(response)]
  if (!length(observed)) {
    return(invisible(FALSE))
  }

  zero_fraction <- mean(observed == 0)
  is_zero_inflated_family <- family %in% zero_inflated_families
  is_zero_heavy_data <- is.finite(zero_fraction) && zero_fraction >= zero_threshold
  if (!is_zero_inflated_family && !is_zero_heavy_data) {
    return(invisible(FALSE))
  }

  warning(
    sprintf(
      paste(
        "Analytical Hessian for zero-heavy discrete margins may be numerically",
        "delicate, especially zero-inflation/shape curvature (family = %s,",
        "zero fraction = %.3f). Consider vcov(..., method = \"numderiv\") if",
        "standard errors are important."
      ),
      family,
      zero_fraction
    ),
    call. = FALSE
  )
  invisible(TRUE)
}


# ------------------------------------------------------------
# 1.  Per-observation CDF first & second derivatives
# ------------------------------------------------------------

#' First derivative of F(y) w.r.t. natural parameters (per observation)
#'
#' Returns a list named by margin parameter (mu, sigma, nu, tau).
#' Each element is a numeric vector of length n_obs.
#' Uses a 2-point central finite difference on the CDF only (O(2 * q) CDF
#' evaluations, not full-likelihood evaluations).
#'
#' @noRd
.calc_dFdpar <- function(eta_inv, mm, margin_dist, response, h = 1e-4) {
  margin_pars <- names(eta_inv)[names(eta_inv) %in% c("mu", "sigma", "nu", "tau")]

  # Build argument list for the p-function
  pfun_name <- paste0("p", margin_dist$family[1])
  pfun <- tryCatch(eval(parse(text = pfun_name)), error = function(e) NULL)
  if (is.null(pfun)) stop("Cannot find CDF function: ", pfun_name)

  args_base <- list()
  args_base[["q"]] <- response
  for (pn in c("mu", "sigma", "nu", "tau")) {
    if (pn %in% names(eta_inv)) args_base[[pn]] <- eta_inv[[pn]]
  }
  valid_args_base <- args_base[names(args_base) %in% names(formals(pfun))]

  result <- vector("list", length(margin_pars))
  names(result) <- margin_pars

  for (pn in margin_pars) {
    args_p <- valid_args_base
    args_m <- valid_args_base
    hp <- .natural_fd_step(as.numeric(args_p[[pn]]), pn, margin_dist, h)
    args_p[[pn]] <- args_p[[pn]] + hp
    args_m[[pn]] <- args_m[[pn]] - hp
    Fp <- do.call(pfun, args_p)
    Fm <- do.call(pfun, args_m)
    dF <- (Fp - Fm) / (2 * hp)
    dF[!is.finite(dF)] <- 0
    result[[pn]] <- as.numeric(dF)
  }
  result
}


#' Second diagonal derivative d2F/dpar^2 (per observation)
#'
#' @keywords internal
.calc_d2Fdpar2 <- function(eta_inv, mm, margin_dist, response, h = 1e-4) {
  margin_pars <- names(eta_inv)[names(eta_inv) %in% c("mu", "sigma", "nu", "tau")]

  pfun_name <- paste0("p", margin_dist$family[1])
  pfun <- tryCatch(eval(parse(text = pfun_name)), error = function(e) NULL)
  if (is.null(pfun)) stop("Cannot find CDF function: ", pfun_name)

  args_base <- list()
  args_base[["q"]] <- response
  for (pn in c("mu", "sigma", "nu", "tau")) {
    if (pn %in% names(eta_inv)) args_base[[pn]] <- eta_inv[[pn]]
  }
  valid_args_base <- args_base[names(args_base) %in% names(formals(pfun))]
  F0 <- do.call(pfun, valid_args_base)
  F0[!is.finite(F0)] <- NA

  result <- vector("list", length(margin_pars))
  names(result) <- margin_pars

  for (pn in margin_pars) {
    args_p <- valid_args_base
    args_m <- valid_args_base
    hp <- .natural_fd_step(as.numeric(args_p[[pn]]), pn, margin_dist, h)
    args_p[[pn]] <- args_p[[pn]] + hp
    args_m[[pn]] <- args_m[[pn]] - hp
    Fp <- do.call(pfun, args_p)
    Fm <- do.call(pfun, args_m)
    d2F <- (Fp + Fm - 2 * F0) / (hp^2)
    d2F[!is.finite(d2F)] <- 0
    result[[pn]] <- as.numeric(d2F)
  }
  result
}


#' Cross second derivative d2F/(dpar1 dpar2) (per observation)
#'
#' @keywords internal
.calc_d2Fdpar_cross <- function(eta_inv, mm, margin_dist, response, h = 1e-4) {
  margin_pars <- names(eta_inv)[names(eta_inv) %in% c("mu", "sigma", "nu", "tau")]
  n_par <- length(margin_pars)

  pfun_name <- paste0("p", margin_dist$family[1])
  pfun <- tryCatch(eval(parse(text = pfun_name)), error = function(e) NULL)
  if (is.null(pfun)) stop("Cannot find CDF function: ", pfun_name)

  args_base <- list()
  args_base[["q"]] <- response
  for (pn in c("mu", "sigma", "nu", "tau")) {
    if (pn %in% names(eta_inv)) args_base[[pn]] <- eta_inv[[pn]]
  }
  valid_args_base <- args_base[names(args_base) %in% names(formals(pfun))]

  # result[[par1]][[par2]] = d2F/(d par1 d par2)
  result <- vector("list", n_par)
  names(result) <- margin_pars
  for (pn in margin_pars) {
    result[[pn]] <- vector("list", n_par)
    names(result[[pn]]) <- margin_pars
    result[[pn]][[pn]] <- NULL  # diagonal handled by .calc_d2Fdpar2
  }

  if (n_par < 2) return(result)

  for (i in seq_len(n_par - 1)) {
    pn1 <- margin_pars[i]
    for (j in (i + 1):n_par) {
      pn2 <- margin_pars[j]
      h1 <- .natural_fd_step(as.numeric(valid_args_base[[pn1]]), pn1, margin_dist, h)
      h2 <- .natural_fd_step(as.numeric(valid_args_base[[pn2]]), pn2, margin_dist, h)
      # 4-point mixed-derivative formula
      args_pp <- valid_args_base; args_pp[[pn1]] <- args_pp[[pn1]] + h1; args_pp[[pn2]] <- args_pp[[pn2]] + h2
      args_pm <- valid_args_base; args_pm[[pn1]] <- args_pm[[pn1]] + h1; args_pm[[pn2]] <- args_pm[[pn2]] - h2
      args_mp <- valid_args_base; args_mp[[pn1]] <- args_mp[[pn1]] - h1; args_mp[[pn2]] <- args_mp[[pn2]] + h2
      args_mm <- valid_args_base; args_mm[[pn1]] <- args_mm[[pn1]] - h1; args_mm[[pn2]] <- args_mm[[pn2]] - h2

      Fpp <- do.call(pfun, args_pp)
      Fpm <- do.call(pfun, args_pm)
      Fmp <- do.call(pfun, args_mp)
      Fmm <- do.call(pfun, args_mm)
      d2F <- (Fpp - Fpm - Fmp + Fmm) / (4 * h1 * h2)
      d2F[!is.finite(d2F)] <- 0
      result[[pn1]][[pn2]] <- as.numeric(d2F)
      result[[pn2]][[pn1]] <- as.numeric(d2F)
    }
  }
  result
}


# ------------------------------------------------------------
# 2.  Per-pair copula Hessian contributions
# ------------------------------------------------------------

#' Assemble per-observation Hessian contributions from the copula
#'
#' For each pair (i,j) in a subject's longitudinal sequence, the copula
#' log-likelihood is log c(F(y_i), F(y_j); theta).  Its second derivative
#' w.r.t. natural parameters of the margin and copula is assembled here.
#'
#' Returns a list with named elements:
#'   cop_d2l_margin  – matrix (n_obs x n_margin_par x n_margin_par) flattened
#'                     as `list[[par1]][[par2]]`, each of length n_obs
#'   cop_d2l_theta   – length n_obs vector: d2 log c / d theta^2
#'   cop_d2l_zeta    – length n_obs vector or NULL
#'   cop_d2l_thetazeta – length n_obs vector or NULL
#'   cop_d2l_margin_theta – `list[[par_name]]` length n_obs: d2logc/(dpar dtheta)
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

  # Pair-indexed margin×theta cross-terms: d2logc/(d par_pn_i1 d theta_k) for u1 and u2 roles.
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
      # H[a,b] = sum_k X_a[k] X_b[k] d2l/d eta_a d eta_b — here k is observation.
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

#' Assemble the full covariate Hessian from per-observation second derivatives
#'
#' `H[a,b] = X_a' diag(w_ab) X_b`, where `w_ab[i] = d2l/(deta_a deta_b)[i]`.
#'
#' @keywords internal
.assemble_covariate_hessian <- function(
    object,
    margin_d1l,           # list[[pn]]: per-obs d log margin / d par
    margin_d2l,           # list[[pn1]][[pn2]]: per-obs d2 log margin / d par1 d par2
    copula_hess,          # output of .calc_copula_hessian_contributions
    eta_dr,               # link function derivatives d par/d eta, named list
    eta_d2,               # link inverse second derivatives d2 par / d eta2
    mm,                   # model matrix list (mm$x, mm$s)
    pair_cache            # pair_cache for theta_index_map
) {
  par_names  <- names(object$par)
  n_par      <- length(par_names)
  H          <- matrix(0, nrow = n_par, ncol = n_par,
                       dimnames = list(par_names, par_names))

  mm_x <- mm$x
  # Map each covariate coefficient to its distribution parameter block
  # e.g. "mu.(Intercept)", "mu.time", "sigma.(Intercept)", ...
  # par_names has format "param.term" or just "param" for intercept.
  param_of <- function(cn) {
    # parameter block is everything before the first dot
    m <- regexpr("^[^.]+", cn)
    regmatches(cn, m)
  }

  param_block <- sapply(par_names, param_of)

  # Helper: get design matrix column for covariate coefficient a
  get_x_col <- function(a_name) {
    blk <- param_block[a_name]
    if (!blk %in% names(mm_x)) return(NULL)
    X <- mm_x[[blk]]
    # column name within the block is the part after "blk."
    col_name <- sub(paste0("^", blk, "\\."), "", a_name)
    if (col_name == a_name) col_name <- "(Intercept)"  # no dot = intercept
    if (!col_name %in% colnames(X)) {
      # Smooth columns can be appended to mm$x as numbered basis columns while
      # the coefficient names retain the full smooth label, e.g.
      # theta.s(s1, k = 10).3 -> design column "3".
      basis_col <- sub("^.*\\.([0-9]+)$", "\\1", col_name)
      if (!identical(basis_col, col_name) && basis_col %in% colnames(X)) {
        col_name <- basis_col
      }
    }
    if (!col_name %in% colnames(X)) return(NULL)
    X[, col_name, drop = TRUE]
  }

  # Helper: per-obs margin second derivative (margin x margin only)
  d2l_par_margin_obs <- function(pn1, pn2) {
    m_pars <- names(margin_d2l)
    n_obs  <- length(eta_dr[[if (pn1 %in% names(eta_dr)) pn1 else pn2]])
    val    <- numeric(n_obs)
    if (pn1 %in% m_pars && pn2 %in% m_pars) {
      if (!is.null(margin_d2l[[pn1]][[pn2]]))
        val <- val + margin_d2l[[pn1]][[pn2]]
      if (!is.null(copula_hess$cop_d2l_margin[[pn1]][[pn2]]))
        val <- val + copula_hess$cop_d2l_margin[[pn1]][[pn2]]
    }
    val
  }

  d1l_par_margin_obs <- function(pn) {
    n_obs <- length(eta_dr[[pn]])
    val <- numeric(n_obs)
    if (pn %in% names(margin_d1l)) val <- val + margin_d1l[[pn]]
    if (!is.null(copula_hess$cop_d1l_margin) && pn %in% names(copula_hess$cop_d1l_margin)) {
      val <- val + copula_hess$cop_d1l_margin[[pn]]
    }
    val
  }

  row_id1 <- copula_hess$row_id1
  row_id2 <- copula_hess$row_id2
  pair_ok <- copula_hess$pair_ok
  id1_ok  <- row_id1[pair_ok]
  id2_ok  <- row_id2[pair_ok]

  for (a_idx in seq_len(n_par)) {
    a_name <- par_names[a_idx]
    pa     <- param_block[a_name]
    xa     <- get_x_col(a_name)
    if (is.null(xa)) next
    dra    <- as.numeric(eta_dr[[pa]])

    for (b_idx in seq_len(a_idx)) {  # lower triangle + diagonal
      b_name <- par_names[b_idx]
      pb     <- param_block[b_name]
      xb     <- get_x_col(b_name)
      if (is.null(xb)) next
      drb    <- as.numeric(eta_dr[[pb]])

      H_ab <- 0

      # Note: theta/zeta design matrices are indexed by theta_index_map[row_id1]
      # (one row per time-1..T-1 observation = same as mm$x$theta row order).
      # Margin design matrices are obs-indexed (n_obs rows).
      if (pa == "theta" && pb == "theta") {
        # cop_d2l_theta is already in eta space (curvature correction applied);
        # do NOT multiply by dra*drb here.
        th_idx <- pair_cache$theta_index_map[row_id1[pair_ok]]
        valid  <- !is.na(th_idx)
        d2_pairs <- copula_hess$cop_d2l_theta[pair_ok]
        H_ab <- sum(xa[th_idx[valid]] * xb[th_idx[valid]] * d2_pairs[valid],
                    na.rm = TRUE)

      } else if (!is.null(copula_hess$cop_d2l_zeta) &&
                 pa == "zeta" && pb == "zeta") {
        # cop_d2l_zeta is already in eta space (curvature correction applied);
        # do NOT multiply by dra*drb here.
        th_idx <- pair_cache$theta_index_map[row_id1[pair_ok]]
        valid  <- !is.na(th_idx)
        d2_pairs <- copula_hess$cop_d2l_zeta[pair_ok]
        H_ab <- sum(xa[th_idx[valid]] * xb[th_idx[valid]] * d2_pairs[valid],
                    na.rm = TRUE)

      } else if (!is.null(copula_hess$cop_d2l_thetazeta) &&
                 ((pa == "theta" && pb == "zeta") ||
                  (pa == "zeta"  && pb == "theta"))) {
        th_idx <- pair_cache$theta_index_map[row_id1[pair_ok]]
        valid  <- !is.na(th_idx)
        d2_pairs <- copula_hess$cop_d2l_thetazeta[pair_ok]
        H_ab <- sum(xa[th_idx[valid]] * xb[th_idx[valid]] * d2_pairs[valid] *
                    dra[th_idx[valid]] * drb[th_idx[valid]], na.rm = TRUE)

      } else if ((pa %in% names(margin_d2l) && pb == "theta") ||
                 (pb %in% names(margin_d2l) && pa == "theta")) {
        # Margin x theta cross block: pair-indexed to avoid design-matrix length mismatch.
        # H[sigma.a, theta.b] = sum_k {
        #   X_sigma[i1_k] * dr_sigma[i1_k] * d2logc/(d par_i1 d theta_k) * dr_theta[th_k] * X_theta[th_k]
        # + X_sigma[i2_k] * dr_sigma[i2_k] * d2logc/(d par_i2 d theta_k) * dr_theta[th_k] * X_theta[th_k]
        # }
        mp  <- if (pa %in% names(margin_d2l)) pa else pb  # margin parameter name
        xm  <- if (pa %in% names(margin_d2l)) xa else xb  # margin design column (n_obs)
        drm <- if (pa %in% names(margin_d2l)) dra else drb # margin eta_dr (n_obs)
        # theta design column and eta_dr (n_theta_rows); index via th_idx
        xth <- if (pb == "theta") xb else xa
        drth <- if (pb == "theta") drb else dra

        th_idx <- pair_cache$theta_index_map[row_id1[pair_ok]]
        valid  <- !is.na(th_idx)
        th_v   <- th_idx[valid]
        id1_v  <- id1_ok[valid]
        id2_v  <- id2_ok[valid]

        u1_vals <- copula_hess$cop_d2l_margin_theta_u1[[mp]][pair_ok][valid]
        u2_vals <- copula_hess$cop_d2l_margin_theta_u2[[mp]][pair_ok][valid]

        H_ab <- sum(
          xm[id1_v] * drm[id1_v] * u1_vals * drth[th_v] * xth[th_v] +
          xm[id2_v] * drm[id2_v] * u2_vals * drth[th_v] * xth[th_v],
          na.rm = TRUE
        )

      } else if (!is.null(copula_hess$cop_d2l_margin_zeta_u1) &&
                 ((pa %in% names(margin_d2l) && pb == "zeta") ||
                  (pb %in% names(margin_d2l) && pa == "zeta"))) {
        # Margin x zeta cross block. This is analogous to margin x theta,
        # but zeta uses the log_2plus inverse link, so include dzeta/deta.
        mp  <- if (pa %in% names(margin_d2l)) pa else pb
        xm  <- if (pa %in% names(margin_d2l)) xa else xb
        drm <- if (pa %in% names(margin_d2l)) dra else drb
        xze <- if (pb == "zeta") xb else xa
        drze <- if (pb == "zeta") drb else dra

        th_idx <- pair_cache$theta_index_map[row_id1[pair_ok]]
        valid  <- !is.na(th_idx)
        th_v   <- th_idx[valid]
        id1_v  <- id1_ok[valid]
        id2_v  <- id2_ok[valid]

        u1_vals <- copula_hess$cop_d2l_margin_zeta_u1[[mp]][pair_ok][valid]
        u2_vals <- copula_hess$cop_d2l_margin_zeta_u2[[mp]][pair_ok][valid]

        H_ab <- sum(
          xm[id1_v] * drm[id1_v] * u1_vals * drze[th_v] * xze[th_v] +
          xm[id2_v] * drm[id2_v] * u2_vals * drze[th_v] * xze[th_v],
          na.rm = TRUE
        )

      } else {
        # Margin x margin: obs-level
        d2l_obs <- d2l_par_margin_obs(pa, pb) * dra * drb
        if (identical(pa, pb) && pa %in% names(eta_d2)) {
          d2l_obs <- d2l_obs + d1l_par_margin_obs(pa) * as.numeric(eta_d2[[pa]])
        }
        H_ab <- sum(xa * xb * d2l_obs, na.rm = TRUE)

        # Cross-pair copula contributions (obs i1 in u1 role, obs i2 in u2 role)
        if (pa %in% names(copula_hess$cross_pair_contribs) &&
            pb %in% names(copula_hess$cross_pair_contribs[[pa]])) {
          # Term 1: d2logc / (d pa at i1) (d pb at i2)
          cross_ab <- copula_hess$cross_pair_contribs[[pa]][[pb]][pair_ok]
          # Term 2: d2logc / (d pa at i2) (d pb at i1) = cross_pair_contribs[[pb]][[pa]]
          cross_ba <- copula_hess$cross_pair_contribs[[pb]][[pa]][pair_ok]
          H_ab <- H_ab + sum(
            xa[id1_ok] * dra[id1_ok] * cross_ab * xb[id2_ok] * drb[id2_ok] +
            xa[id2_ok] * dra[id2_ok] * cross_ba * xb[id1_ok] * drb[id1_ok],
            na.rm = TRUE
          )
        }
      }

      H[a_idx, b_idx] <- H_ab
      H[b_idx, a_idx] <- H_ab
    }
  }
  H
}


# ------------------------------------------------------------
# 4.  Margin per-observation second derivative (observed, via FD on log-density)
# ------------------------------------------------------------

#' Compute observed per-observation margin log-density second derivatives
#' via finite differences in the natural (pre-link) parameter space.
#'
#' This replaces the gamlss family d2ld*2 functions which return the
#' expected Fisher information (y-independent) rather than the observed
#' second derivative.
#'
#' Returns a `list[[pn1]][[pn2]]` of length n_obs.
#' @keywords internal
.calc_margin_d2l_fd <- function(eta_inv, mm, margin_dist, response, h = 1e-4) {
  margin_pars <- names(mm$x)[names(mm$x) %in% c("mu", "sigma", "nu", "tau")]
  n_obs <- length(response)
  obs_valid <- is.finite(response)

  dfun_name <- paste0("d", margin_dist$family[1])
  dfun <- eval(parse(text = dfun_name))

  # Build base argument list
  base_args <- list(x = response, log = TRUE)
  for (pn in margin_pars) {
    base_args[[pn]] <- as.numeric(eta_inv[[pn]])
  }

  # log-density at base parameters (all observations)
  ll0 <- do.call(dfun, base_args)
  ll0[!obs_valid] <- 0

  result <- vector("list", length(margin_pars))
  names(result) <- margin_pars
  for (pni in margin_pars) {
    result[[pni]] <- vector("list", length(margin_pars))
    names(result[[pni]]) <- margin_pars
  }

  for (i in seq_along(margin_pars)) {
    pni <- margin_pars[i]
    par_i <- as.numeric(eta_inv[[pni]])
    hi <- .natural_fd_step(par_i, pni, margin_dist, h)

    # Log-density at +hi and -hi for pni
    args_ip <- base_args; args_ip[[pni]] <- par_i + hi
    args_im <- base_args; args_im[[pni]] <- par_i - hi
    ll_ip <- do.call(dfun, args_ip); ll_ip[!obs_valid] <- 0
    ll_im <- do.call(dfun, args_im); ll_im[!obs_valid] <- 0

    # Diagonal: d2l/d pni^2 via 3-point stencil
    d2_diag <- (ll_ip - 2 * ll0 + ll_im) / hi^2
    d2_diag[!is.finite(d2_diag)] <- 0
    result[[pni]][[pni]] <- d2_diag

    for (j in seq_along(margin_pars)) {
      if (j <= i) next  # fill upper triangle, then mirror
      pnj <- margin_pars[j]
      par_j <- as.numeric(eta_inv[[pnj]])
      hj <- .natural_fd_step(par_j, pnj, margin_dist, h)

      args_jp <- base_args; args_jp[[pnj]] <- par_j + hj
      args_jm <- base_args; args_jm[[pnj]] <- par_j - hj

      args_ipjp <- args_jp; args_ipjp[[pni]] <- par_i + hi
      args_ipjm <- args_jm; args_ipjm[[pni]] <- par_i + hi
      args_imjp <- args_jp; args_imjp[[pni]] <- par_i - hi
      args_imjm <- args_jm; args_imjm[[pni]] <- par_i - hi

      ll_pp <- do.call(dfun, args_ipjp); ll_pp[!obs_valid] <- 0
      ll_pm <- do.call(dfun, args_ipjm); ll_pm[!obs_valid] <- 0
      ll_mp <- do.call(dfun, args_imjp); ll_mp[!obs_valid] <- 0
      ll_mm <- do.call(dfun, args_imjm); ll_mm[!obs_valid] <- 0

      d2_cross <- (ll_pp - ll_pm - ll_mp + ll_mm) / (4 * hi * hj)
      d2_cross[!is.finite(d2_cross)] <- 0
      result[[pni]][[pnj]] <- d2_cross
      result[[pnj]][[pni]] <- d2_cross  # symmetric
    }
  }
  result
}

#' Extract per-observation margin log-likelihood second derivatives
#'
#' Returns a `list[[pn1]][[pn2]]`, each of length n_obs, sourced from the
#' gamlss family d2ldm2, d2ldmdd, ... functions already called in
#' calc_likelihood_minimal.
#'
#' @keywords internal
.get_margin_d1l <- function(margin_deriv, mm) {
  margin_pars <- names(mm$x)[names(mm$x) %in% c("mu", "sigma", "nu", "tau")]
  subs <- c(mu = "m", sigma = "d", nu = "v", tau = "t")
  result <- vector("list", length(margin_pars))
  names(result) <- margin_pars
  for (pn in margin_pars) {
    key <- paste0("dld", subs[[pn]])
    if (key %in% names(margin_deriv)) {
      result[[pn]] <- as.numeric(margin_deriv[[key]])
    } else {
      n <- if (length(margin_deriv)) length(margin_deriv[[1]]) else nrow(mm$x[[pn]])
      result[[pn]] <- numeric(n)
    }
  }
  result
}

.get_margin_d2l <- function(margin_deriv, mm) {
  margin_pars <- names(mm$x)[names(mm$x) %in% c("mu", "sigma", "nu", "tau")]
  n_par <- length(margin_pars)

  # Map gamlss derivative name endings to (par1, par2) pairs
  subs <- c(mu = "m", sigma = "d", nu = "v", tau = "t")

  result <- vector("list", n_par)
  names(result) <- margin_pars
  for (pn in margin_pars) {
    result[[pn]] <- vector("list", n_par)
    names(result[[pn]]) <- margin_pars
    for (pn2 in margin_pars) result[[pn]][[pn2]] <- NULL
  }

  for (i in seq_along(margin_pars)) {
    pni <- margin_pars[i]
    si  <- subs[pni]
    for (j in seq_along(margin_pars)) {
      pnj <- margin_pars[j]
      sj  <- subs[pnj]
      if (i == j) {
        key <- paste0("d2ld", si, "2")
      } else {
        # gamlss uses alphabetical order of the two sub-letters
        key1 <- paste0("d2ld", si, "d", sj)
        key2 <- paste0("d2ld", sj, "d", si)
        key  <- if (key1 %in% names(margin_deriv)) key1 else
                if (key2 %in% names(margin_deriv)) key2 else NA_character_
      }
      if (!is.na(key) && key %in% names(margin_deriv)) {
        result[[pni]][[pnj]] <- as.numeric(margin_deriv[[key]])
      } else {
        # Not provided by this family: assume zero cross term
        n <- length(margin_deriv[[1]])
        result[[pni]][[pnj]] <- numeric(n)
      }
    }
  }
  result
}


# ------------------------------------------------------------
# 5.  Top-level analytical Hessian function
# ------------------------------------------------------------

#' Compute analytical (semi-analytical) Hessian for gamlss.longitudinal
#'
#' This is the main entry point called by vcov.gamlss.longitudinal when
#' method = "analytical".
#'
#' @param object A fitted gamlss.longitudinal object.
#' @param progress Logical; show a progress indicator.
#' @param h Step size for CDF finite differences (default 1e-4).
#'
#' @return A square matrix of dimension length(object$par) x length(object$par)
#'   representing the Hessian of the joint log-likelihood evaluated at the MLE.
#'
#' @keywords internal
#' @noRd
calc_analytical_hessian <- function(object, progress = interactive(), h = 1e-4) {
  progress <- isTRUE(progress)

  response         <- object$response
  response_margin  <- object$response_margin
  response_subject <- object$response_subject
  margin_dist      <- object$margin_dist
  copula_dist      <- object$copula_dist
  copula_link      <- get_copula_dist(copula_dist)$copula_link
  mm               <- object$model_matrix
  par_cov          <- object$par
  par_s            <- object$par_s

  if (progress) cat("Analytical Hessian: computing eta and likelihood components...\n")

  eta_out  <- calc_eta(par_cov, mm, margin_dist, copula_link, par_s = par_s)
  eta_inv  <- eta_out$eta_inv
  eta_dr   <- eta_out$eta_dr
  eta_d2   <- .calc_eta_d2_linkinv(eta_out$eta, margin_dist, copula_link)
  .warn_gg_near_zero_nu(margin_dist, eta_inv)
  .warn_zero_heavy_discrete_hessian(margin_dist, response)

  pair_cache <- build_copula_pair_cache(response, response_margin, response_subject)

  # Full likelihood call to get margin_deriv (analytical d1 and d2 from gamlss family)
  calc_lik <- calc_likelihood_minimal(
    eta_inv, mm = mm$x, margin_dist, copula_dist,
    calc_d2 = TRUE,
    response = response, response_margin = response_margin,
    response_subject = response_subject,
    pair_cache = pair_cache
  )

  # Attach margin_p to eta_inv for use inside copula hessian function
  eta_inv[["margin_p_cache"]] <- calc_lik$margin_p

  if (progress) cat("Analytical Hessian: computing CDF derivatives...\n")

  # Evaluate only on observed (non-NA) response rows; set to 0 elsewhere
  dF      <- .calc_dFdpar(eta_inv, mm, margin_dist, response, h = h)
  d2F     <- .calc_d2Fdpar2(eta_inv, mm, margin_dist, response, h = h)
  d2F_x   <- .calc_d2Fdpar_cross(eta_inv, mm, margin_dist, response, h = h)

  if (progress) cat("Analytical Hessian: assembling copula Hessian contributions...\n")

  cop_hess <- .calc_copula_hessian_contributions(
    eta_inv, pair_cache, copula_dist, dF, d2F, d2F_x
  )

  if (progress) cat("Analytical Hessian: extracting margin second derivatives...\n")

  margin_d1l <- .get_margin_d1l(calc_lik$margin_deriv, mm)
  margin_d2l <- .calc_margin_d2l_fd(eta_inv, mm, margin_dist, response, h = h)

  if (progress) cat("Analytical Hessian: assembling covariate Hessian matrix...\n")

  H <- .assemble_covariate_hessian(object, margin_d1l, margin_d2l, cop_hess, eta_dr, eta_d2, mm, pair_cache)

  if (progress) cat("Analytical Hessian: done.\n")

  H
}
