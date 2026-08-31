#' Apply distribution-specific vcov method preflight checks
#'
#' @noRd
.gl_vcov_apply_margin_preflight <- function(method, method_used, margin_dist, eta_inv) {
  if (identical(method, "analytical") &&
    identical(as.character(margin_dist$family[1]), "GG") &&
    "nu" %in% names(eta_inv)) {
    nu_abs_min <- suppressWarnings(min(abs(as.numeric(eta_inv$nu[is.finite(eta_inv$nu)])), na.rm = TRUE))
    if (is.finite(nu_abs_min) && nu_abs_min < 0.06) {
      warning(
        sprintf(
          paste(
            "Analytical Hessian for GG may be numerically unstable because fitted",
            "nu is close to 0 (min |nu| = %.4g); falling back to numerical Hessian."
          ),
          nu_abs_min
        ),
        call. = FALSE
      )
      method <- "numderiv"
      method_used <- "numderiv"
    }
  }

  list(method = method, method_used = method_used)
}

#' Apply likelihood-specific vcov method preflight checks
#'
#' @noRd
.gl_vcov_apply_likelihood_preflight <- function(method, method_used, calc_lik_out, response) {
  if (method %in% c("analytical", "analytical_only") &&
    identical(calc_lik_out$likelihood_type, "discrete_rectangle")) {
    zero_fraction <- mean(response == 0, na.rm = TRUE)
    if (is.finite(zero_fraction) && zero_fraction >= 0.35) {
      warning(
        paste0(
          "Analytical Hessian for zero-heavy discrete margins may be numerically delicate; ",
          sprintf("zero fraction = %.3f.", zero_fraction)
        ),
        call. = FALSE
      )
    }
  }

  list(method = method, method_used = method_used)
}
