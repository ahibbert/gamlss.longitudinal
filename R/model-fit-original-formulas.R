#' Return original user formulas for fitted-object metadata
#'
#' @param mu.formula,sigma.formula,nu.formula,tau.formula Margin formulas.
#' @param theta.formula,zeta.formula Copula formulas.
#' @return Named list of original formulas, preserving user-facing names.
#' @noRd
.gl_original_formula_bundle <- function(
    mu.formula,
    sigma.formula,
    nu.formula,
    tau.formula,
    theta.formula,
    zeta.formula) {
  list(
    mu = mu.formula,
    sigma = sigma.formula,
    nu = nu.formula,
    tau = tau.formula,
    theta = theta.formula,
    zeta = zeta.formula
  )
}
