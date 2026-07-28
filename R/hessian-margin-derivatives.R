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
#' @noRd

.calc_margin_d2l_fd <- function(eta_inv, mm, margin_dist, response, h = 1e-4) {
  density_setup <- .hessian_margin_density_setup(eta_inv, mm, margin_dist, response)

  margin_pars <- density_setup$margin_pars

  n_obs <- length(response)

  obs_valid <- is.finite(response)

  dfun <- density_setup$dfun

  base_args <- density_setup$base_args

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

    args_ip <- base_args
    args_ip[[pni]] <- par_i + hi

    args_im <- base_args
    args_im[[pni]] <- par_i - hi

    ll_ip <- do.call(dfun, args_ip)
    ll_ip[!obs_valid] <- 0

    ll_im <- do.call(dfun, args_im)
    ll_im[!obs_valid] <- 0

    # Diagonal: d2l/d pni^2 via 3-point stencil

    d2_diag <- (ll_ip - 2 * ll0 + ll_im) / hi^2

    d2_diag[!is.finite(d2_diag)] <- 0

    result[[pni]][[pni]] <- d2_diag

    for (j in seq_along(margin_pars)) {
      if (j <= i) next # fill upper triangle, then mirror

      pnj <- margin_pars[j]

      par_j <- as.numeric(eta_inv[[pnj]])

      hj <- .natural_fd_step(par_j, pnj, margin_dist, h)

      args_jp <- base_args
      args_jp[[pnj]] <- par_j + hj

      args_jm <- base_args
      args_jm[[pnj]] <- par_j - hj

      args_ipjp <- args_jp
      args_ipjp[[pni]] <- par_i + hi

      args_ipjm <- args_jm
      args_ipjm[[pni]] <- par_i + hi

      args_imjp <- args_jp
      args_imjp[[pni]] <- par_i - hi

      args_imjm <- args_jm
      args_imjm[[pni]] <- par_i - hi

      ll_pp <- do.call(dfun, args_ipjp)
      ll_pp[!obs_valid] <- 0

      ll_pm <- do.call(dfun, args_ipjm)
      ll_pm[!obs_valid] <- 0

      ll_mp <- do.call(dfun, args_imjp)
      ll_mp[!obs_valid] <- 0

      ll_mm <- do.call(dfun, args_imjm)
      ll_mm[!obs_valid] <- 0

      d2_cross <- (ll_pp - ll_pm - ll_mp + ll_mm) / (4 * hi * hj)

      d2_cross[!is.finite(d2_cross)] <- 0

      result[[pni]][[pnj]] <- d2_cross

      result[[pnj]][[pni]] <- d2_cross # symmetric
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
#' @noRd

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

    si <- subs[pni]

    for (j in seq_along(margin_pars)) {
      pnj <- margin_pars[j]

      sj <- subs[pnj]

      if (i == j) {
        key <- paste0("d2ld", si, "2")
      } else {
        # gamlss uses alphabetical order of the two sub-letters

        key1 <- paste0("d2ld", si, "d", sj)

        key2 <- paste0("d2ld", sj, "d", si)

        key <- if (key1 %in% names(margin_deriv)) key1 else if (key2 %in% names(margin_deriv)) key2 else NA_character_
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
