#' Positive-bound link helpers
#'
#' @rdname link_functions
#' @noRd
NULL

log_2plus <- function(x) {
  return(
    log(x - 2)
  )
}
log_2plus_inv <- function(x) {
  y <- exp(x) + 2
  # Adjust for error when close to two.
  close_to_boundary <- is.finite(y) & abs(y - 2) <= .Machine$double.eps
  y[close_to_boundary] <- y[close_to_boundary] + 0.00001
  return(y)
}
dlog_2plus <- function(x) {
  return(1 / (x - 2))
}
dlog <- function(x) {
  return(1 / x)
}
dlog_inv <- function(x) {
  return(exp(x))
}

dlog_2plus_inv <- function(x) {
  return(exp(x))
}

log_1plus <- function(x) {
  return(
    log(x - 1)
  )
}

log_1plus_inv <- function(x) {
  y <- exp(x) + 1
  close_to_boundary <- is.finite(y) & abs(y - 1) <= .Machine$double.eps
  y[close_to_boundary] <- y[close_to_boundary] + 0.00001
  return(y)
}

dlog_1plus <- function(x) {
  return(1 / (x - 1))
}

dlog_1plus_inv <- function(x) {
  return(exp(x))
}
