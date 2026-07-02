#' Copula link helpers
#'
#' @rdname link_functions
#' @noRd
NULL

fisher_z <- function(x) {
  x <- pmin(pmax(x, -0.999999), 0.999999)
  return(atanh(x))
}

fisher_z_inv <- function(x) {
  return(tanh(x))
}

dfisher_z_inv <- function(x) {
  y <- tanh(x)
  return(1 - y^2)
}

gumbel_linkfun <- function(x) {
  x <- pmin(x, 17)
  return(log(x - 1))
}

gumbel_linkinv <- function(x) {
  y <- exp(x) + 1
  y[y > 17] <- 17
  return(y)
}

dgumbel_linkinv <- function(x) {
  y <- exp(x)
  y[y >= 16] <- 0
  return(y)
}
