# hessian-*.R

# ============================================================

# Analytical / semi-analytical Hessian for gamlss.longitudinal

# ============================================================

# This file is used by model-vcov.R when

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

#   Step 1 â€“ per-pair / per-observation second derivatives of the joint

#             log-likelihood w.r.t. natural parameters (mu, sigma, ...,

#             theta, zeta). These use analytical formulas from gamlss

#             family objects and .copula_deriv2, plus a cheap

#             per-observation finite difference on F(y) only (not the full

#             likelihood).

#

#   Step 2 â€“ chain rule through the link functions and design matrices:

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


