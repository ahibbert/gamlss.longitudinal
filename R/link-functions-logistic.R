#### LINK FUNCTIONS ####
#' Logistic-style link helpers
#'
#' These helpers provide simple link, inverse-link, and derivative functions
#' used by `gamlss.longitudinal` parameter transformations.
#'
#' @name link_functions
#' @noRd
NULL

logit <- function(x) {
  return(log(x / (1 - x)))
}
logit_inv <- function(x) {
  return(
    if (all(is.nan(exp(x) / (1 + exp(x))))) {
      return(1)
    } else {
      return(exp(x) / (1 + exp(x)))
    }
  )
}
dlogit <- function(x) {
  return(1 / (x - (x^2)))
}
logit28 <- function(x) {
  return(log(x / (28 - x)))
}
logit28_inv <- function(x) {
  return(
    if (is.nan(exp(x) / (1 + exp(x)))) {
      return(1)
    } else {
      return(28 * exp(x) / (1 + exp(x)))
    }
  )
}
dlogit28 <- function(x) {
  return(1 / (28 * x - (x^2)))
}

dlogit_inv <- function(x) {
  return(exp(x) / ((1 + exp(x))^2))
}
