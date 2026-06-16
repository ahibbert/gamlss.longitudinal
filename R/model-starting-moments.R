.starting_moment_skewness <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) < 3) {
    return(0)
  }

  s <- stats::sd(x)

  if (!is.finite(s) || s <= 0) {
    return(0)
  }

  mean(((x - mean(x)) / s)^3)
}

.starting_moment_kurtosis <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) < 4) {
    return(3)
  }

  s <- stats::sd(x)

  if (!is.finite(s) || s <= 0) {
    return(3)
  }

  mean(((x - mean(x)) / s)^4)
}
