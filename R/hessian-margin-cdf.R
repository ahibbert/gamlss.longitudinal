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

  result <- vector("list", length(margin_pars))

  names(result) <- margin_pars

  for (pn in margin_pars) {
    eta_p <- eta_m <- eta_inv
    hp <- .natural_fd_step(as.numeric(eta_inv[[pn]]), pn, margin_dist, h)
    eta_p[[pn]] <- eta_p[[pn]] + hp
    eta_m[[pn]] <- eta_m[[pn]] - hp
    Fp <- calc_F_x(eta_p, mm, margin_dist, response)
    Fm <- calc_F_x(eta_m, mm, margin_dist, response)
    result[[pn]] <- .gl_validate_margin_cdf_derivative(
      (Fp - Fm) / (2 * hp), response, margin_dist,
      context = "Analytical-Hessian first CDF derivative", parameter = pn
    )
  }

  result
}

#' Second diagonal derivative d2F/dpar^2 (per observation)

#'

#' @keywords internal

.calc_d2Fdpar2 <- function(eta_inv, mm, margin_dist, response, h = 1e-4) {
  cdf_setup <- .hessian_margin_cdf_setup(eta_inv, margin_dist, response)
  margin_pars <- cdf_setup$margin_pars
  F0 <- calc_F_x(eta_inv, mm, margin_dist, response)

  result <- vector("list", length(margin_pars))

  names(result) <- margin_pars

  for (pn in margin_pars) {
    eta_p <- eta_m <- eta_inv
    hp <- .natural_fd_step(as.numeric(eta_inv[[pn]]), pn, margin_dist, h)
    eta_p[[pn]] <- eta_p[[pn]] + hp
    eta_m[[pn]] <- eta_m[[pn]] - hp
    Fp <- calc_F_x(eta_p, mm, margin_dist, response)
    Fm <- calc_F_x(eta_m, mm, margin_dist, response)
    result[[pn]] <- .gl_validate_margin_cdf_derivative(
      (Fp + Fm - 2 * F0) / (hp^2), response, margin_dist,
      context = "Analytical-Hessian second CDF derivative", parameter = pn
    )
  }

  result
}

#' Cross second derivative d2F/(dpar1 dpar2) (per observation)

#'

#' @keywords internal

.calc_d2Fdpar_cross <- function(eta_inv, mm, margin_dist, response, h = 1e-4) {
  cdf_setup <- .hessian_margin_cdf_setup(eta_inv, margin_dist, response)
  margin_pars <- cdf_setup$margin_pars

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

      h1 <- .natural_fd_step(as.numeric(eta_inv[[pn1]]), pn1, margin_dist, h)

      h2 <- .natural_fd_step(as.numeric(eta_inv[[pn2]]), pn2, margin_dist, h)

      # 4-point mixed-derivative formula

      eta_pp <- eta_pm <- eta_mp <- eta_mm <- eta_inv
      eta_pp[[pn1]] <- eta_pp[[pn1]] + h1
      eta_pp[[pn2]] <- eta_pp[[pn2]] + h2
      eta_pm[[pn1]] <- eta_pm[[pn1]] + h1
      eta_pm[[pn2]] <- eta_pm[[pn2]] - h2
      eta_mp[[pn1]] <- eta_mp[[pn1]] - h1
      eta_mp[[pn2]] <- eta_mp[[pn2]] + h2
      eta_mm[[pn1]] <- eta_mm[[pn1]] - h1
      eta_mm[[pn2]] <- eta_mm[[pn2]] - h2

      Fpp <- calc_F_x(eta_pp, mm, margin_dist, response)
      Fpm <- calc_F_x(eta_pm, mm, margin_dist, response)
      Fmp <- calc_F_x(eta_mp, mm, margin_dist, response)
      Fmm <- calc_F_x(eta_mm, mm, margin_dist, response)

      d2F <- .gl_validate_margin_cdf_derivative(
        (Fpp - Fpm - Fmp + Fmm) / (4 * h1 * h2),
        response, margin_dist,
        context = "Analytical-Hessian cross CDF derivative",
        parameter = paste(pn1, pn2, sep = "/")
      )

      result[[pn1]][[pn2]] <- as.numeric(d2F)

      result[[pn2]][[pn1]] <- as.numeric(d2F)
    }
  }

  result
}

# ------------------------------------------------------------

# 2.  Per-pair copula Hessian contributions

# ------------------------------------------------------------
