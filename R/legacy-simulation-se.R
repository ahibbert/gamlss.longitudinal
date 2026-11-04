#' Legacy simulation standard-error helper
#'
#' Retained for compatibility and internal historical support. This helper is
#' not part of the core reviewer path for fitting, prediction, or diagnostics.
#'
#' @noRd
NULL

#' @keywords internal
#' @noRd
bvt_norm_true_SE_B0_Bt <- function(sigma_x, sigma_y, rho, n, d) {
  # sigma_x=2
  # sigma_y=2
  # rho=.75

  # ((1)/(1-rho^2))
  # (1/(sigma_x^2))
  # (1/(sigma_y^2))
  # (1/(sigma_x))
  # (1/(sigma_y))
  # (rho/(sigma_x*sigma_y))

  hessian <- ((-1) / (1 - (rho^2))) * matrix(c(
    (1 / (sigma_x^2)) + (1 / (sigma_y^2)) - 2 * (rho / (sigma_x * sigma_y)),
    ((rho / (sigma_x * sigma_y)) - (1 / (sigma_y^2))),
    ((rho / (sigma_x * sigma_y)) - (1 / (sigma_y^2))),
    (1 / (sigma_y^2))
  ), nrow = 2)

  vcov_matrix <- -solve(hessian)

  true_SE <- (diag(vcov_matrix))
  names(true_SE) <- c("B0", "Bt")

  return(true_SE)
}
