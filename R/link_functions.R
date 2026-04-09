
#### LINK FUNCTIONS ####
#' @export
logit <- function(x) {
  return(log(x/(1-x)))
}
#' @export
logit_inv <- function(x) {
  return(
    if(all(is.nan(exp(x)/(1+exp(x))))) {
      return(1)
    } else {
      return(exp(x)/(1+exp(x)))
    }
  )
}
#' @export
dlogit <- function(x) {
  return(1/(x-(x^2)))
}
#' @export
logit28 <- function(x) {
  return(log(x/(28-x)))
}
#' @export
logit28_inv <- function(x) {
  return(
    if(is.nan(exp(x)/(1+exp(x)))) {
      return(1)
    } else {
      return(28*exp(x)/(1+exp(x)))
    }
  )
}
#' @export
dlogit28 <- function(x) {
  return(1/(28*x-(x^2)))
}

#' @export
log_2plus <- function(x) {
  return(
    log(x-2)
  )
}
#' @export
log_2plus_inv <- function(x) {
  y=exp(x)+2
  #Adjust for error when close to two.
  y[y==2]=y[y==2]+0.00001
  return(y)
}
#' @export
dlog_2plus <- function(x) {
  return(-1/(x-2))
}
#' @export
dlog <-function(x) {
  return(1/x)
}
#' @export
dlog_inv <-function(x) {
  return(exp(x))
}
#' @export
dlogit_inv <- function(x) {
  return(exp(x)/((1+exp(x))^2))
}

#' @export
fisher_z <- function(x) {
  x <- pmin(pmax(x, -0.999999), 0.999999)
  return(atanh(x))
}

#' @export
fisher_z_inv <- function(x) {
  return(tanh(x))
}

#' @export
dfisher_z_inv <- function(x) {
  y <- tanh(x)
  return(1 - y^2)
}

#' @export
dlog_2plus_inv <- function(x) {
  return(exp(x))
}

