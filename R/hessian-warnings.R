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
