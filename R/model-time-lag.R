normalize_lag_time <- function(time) {
  if (is.factor(time)) {
    time <- as.character(time)
  }

  if (is.character(time)) {
    time_numeric <- suppressWarnings(as.numeric(time))

    if (anyNA(time_numeric)) {
      stop("ERROR: time must be numeric or numeric-like when use_dlcopdpar=TRUE.")
    }

    time <- time_numeric
  }

  time - 1
}
