#' First derivative of F(y) w.r.t. natural parameters (per observation)

#'

#' Returns a list named by margin parameter (mu, sigma, nu, tau).

#' Each element is a numeric vector of length n_obs.

#' Uses a 2-point central finite difference on the CDF only (O(2 * q) CDF

#' evaluations, not full-likelihood evaluations).

#'

#' @noRd

.calc_dFdpar <- function(eta_inv, mm, margin_dist, response, h = 1e-4) {
  cdf_setup <- .hessian_margin_cdf_setup(eta_inv, margin_dist, response)
  margin_pars <- cdf_setup$margin_pars
  pfun <- cdf_setup$pfun
  valid_args_base <- cdf_setup$args_base

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
  cdf_setup <- .hessian_margin_cdf_setup(eta_inv, margin_dist, response)
  margin_pars <- cdf_setup$margin_pars
  pfun <- cdf_setup$pfun
  valid_args_base <- cdf_setup$args_base

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
  cdf_setup <- .hessian_margin_cdf_setup(eta_inv, margin_dist, response)
  margin_pars <- cdf_setup$margin_pars
  pfun <- cdf_setup$pfun
  valid_args_base <- cdf_setup$args_base

  n_par <- length(margin_pars)

  # result[[par1]][[par2]] = d2F/(d par1 d par2)

  result <- vector("list", n_par)

  names(result) <- margin_pars

  for (pn in margin_pars) {
    result[[pn]] <- vector("list", n_par)

    names(result[[pn]]) <- margin_pars

    result[[pn]][[pn]] <- NULL # diagonal handled by .calc_d2Fdpar2
  }

  if (n_par < 2) {
    return(result)
  }

  for (i in seq_len(n_par - 1)) {
    pn1 <- margin_pars[i]

    for (j in (i + 1):n_par) {
      pn2 <- margin_pars[j]

      h1 <- .natural_fd_step(as.numeric(valid_args_base[[pn1]]), pn1, margin_dist, h)

      h2 <- .natural_fd_step(as.numeric(valid_args_base[[pn2]]), pn2, margin_dist, h)

      # 4-point mixed-derivative formula

      args_pp <- valid_args_base
      args_pp[[pn1]] <- args_pp[[pn1]] + h1
      args_pp[[pn2]] <- args_pp[[pn2]] + h2

      args_pm <- valid_args_base
      args_pm[[pn1]] <- args_pm[[pn1]] + h1
      args_pm[[pn2]] <- args_pm[[pn2]] - h2

      args_mp <- valid_args_base
      args_mp[[pn1]] <- args_mp[[pn1]] - h1
      args_mp[[pn2]] <- args_mp[[pn2]] + h2

      args_mm <- valid_args_base
      args_mm[[pn1]] <- args_mm[[pn1]] - h1
      args_mm[[pn2]] <- args_mm[[pn2]] - h2

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
