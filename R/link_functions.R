
#### LINK FUNCTIONS ####
#' Link helpers for bounded copula and distribution parameters
#'
#' These helpers provide simple link, inverse-link, and derivative functions
#' used by `gamlss.longitudinal` parameter transformations.
#'
#' @name link_functions
#' @noRd
NULL

logit <- function(x) {
  return(log(x/(1-x)))
}
logit_inv <- function(x) {
  return(
    if(all(is.nan(exp(x)/(1+exp(x))))) {
      return(1)
    } else {
      return(exp(x)/(1+exp(x)))
    }
  )
}
dlogit <- function(x) {
  return(1/(x-(x^2)))
}
logit28 <- function(x) {
  return(log(x/(28-x)))
}
logit28_inv <- function(x) {
  return(
    if(is.nan(exp(x)/(1+exp(x)))) {
      return(1)
    } else {
      return(28*exp(x)/(1+exp(x)))
    }
  )
}
dlogit28 <- function(x) {
  return(1/(28*x-(x^2)))
}

log_2plus <- function(x) {
  return(
    log(x-2)
  )
}
log_2plus_inv <- function(x) {
  y=exp(x)+2
  #Adjust for error when close to two.
  close_to_boundary <- is.finite(y) & abs(y - 2) <= .Machine$double.eps
  y[close_to_boundary]=y[close_to_boundary]+0.00001
  return(y)
}
dlog_2plus <- function(x) {
  return(1/(x-2))
}
dlog <-function(x) {
  return(1/x)
}
dlog_inv <-function(x) {
  return(exp(x))
}
dlogit_inv <- function(x) {
  return(exp(x)/((1+exp(x))^2))
}

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

dlog_2plus_inv <- function(x) {
  return(exp(x))
}

log_1plus <- function(x) {
  return(
    log(x-1)
  )
}

log_1plus_inv <- function(x) {
  y = exp(x) + 1
  close_to_boundary <- is.finite(y) & abs(y - 1) <= .Machine$double.eps
  y[close_to_boundary] = y[close_to_boundary] + 0.00001
  return(y)
}

dlog_1plus <- function(x) {
  return(1/(x-1))
}

dlog_1plus_inv <- function(x) {
  return(exp(x))
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

