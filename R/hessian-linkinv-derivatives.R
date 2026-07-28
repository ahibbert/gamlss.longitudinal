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
