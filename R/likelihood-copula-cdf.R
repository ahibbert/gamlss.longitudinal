#' First derivative of a copula CDF with respect to the first margin
#'
#' Wrapper around the h-function used by likelihood and derivative assembly.
#'
#' @noRd
.copula_cdf_du1 <- function(u1, u2, family, par, par2 = 0) {
  .copula_hfunc1(u1, u2, family = family, par = par, par2 = par2)
}

#' First derivative of a copula CDF with respect to the second margin
#'
#' Uses symmetry of the h-function implementation.
#'
#' @noRd
.copula_cdf_du2 <- function(u1, u2, family, par, par2 = 0) {
  .copula_hfunc1(u2, u1, family = family, par = par, par2 = par2)
}

#' Cached Gaussian copula CDF
#'
#' Evaluates the Gaussian copula CDF with boundary shortcuts and a vectorised
#' bivariate-normal backend.
#'
#' @noRd
.copula_gaussian_cdf_cached <- function(u1, u2, par) {
  vals <- .copula_recycle(as.numeric(u1), as.numeric(u2), .copula_gaussian_rho(par))
  u1 <- vals[[1]]
  u2 <- vals[[2]]
  rho <- vals[[3]]
  out <- rep(NA_real_, length(u1))

  finite <- is.finite(u1) & is.finite(u2) & is.finite(rho)
  if (!any(finite)) {
    return(out)
  }

  u1 <- pmin(pmax(u1, 0), 1)
  u2 <- pmin(pmax(u2, 0), 1)

  zero <- finite & (u1 <= 0 | u2 <= 0)
  out[zero] <- 0

  one1 <- finite & !zero & u1 >= 1
  out[one1] <- u2[one1]

  one2 <- finite & !zero & !one1 & u2 >= 1
  out[one2] <- u1[one2]

  independent <- finite & !zero & !one1 & !one2 & abs(rho) <= 1e-12
  out[independent] <- u1[independent] * u2[independent]

  remaining <- finite & !zero & !one1 & !one2 & !independent
  if (any(remaining)) {
    out[remaining] <- pbivnorm::pbivnorm(
      stats::qnorm(u1[remaining]),
      stats::qnorm(u2[remaining]),
      rho[remaining]
    )
  }

  pmin(pmax(as.numeric(out), 0), 1)
}

#' Cached Gaussian rectangle probability
#'
#' Computes discrete/count likelihood rectangle masses from four cached
#' Gaussian CDF evaluations.
#'
#' @noRd
.copula_gaussian_rectangle_prob_cached <- function(u1, u2, l1, l2, par, floor_probability = TRUE) {
  vals <- .copula_recycle(as.numeric(u1), as.numeric(u2), as.numeric(l1), as.numeric(l2), .copula_gaussian_rho(par))
  n <- length(vals[[1]])
  if (n == 0L) {
    return(numeric(0))
  }
  cdf <- .copula_gaussian_cdf_cached(
    c(vals[[1]], vals[[3]], vals[[1]], vals[[3]]),
    c(vals[[2]], vals[[2]], vals[[4]], vals[[4]]),
    rep(vals[[5]], 4L)
  )
  rect <- cdf[seq_len(n)] -
    cdf[n + seq_len(n)] -
    cdf[2L * n + seq_len(n)] +
    cdf[3L * n + seq_len(n)]
  rect <- as.numeric(rect)
  if (isTRUE(floor_probability)) pmax(rect, 1e-300) else rect
}

#' Copula rectangle probability
#'
#' Returns positive lower-bounded rectangle probabilities for discrete margins,
#' using the cached Gaussian path where possible and generic inclusion-exclusion
#' for other copulas.
#'
#' @noRd
.copula_rectangle_prob <- function(u1, u2, l1, l2, family, par, par2 = 0,
                                   floor_probability = TRUE) {
  family <- .copula_family_code(family)
  if (identical(family, "N")) {
    return(.copula_gaussian_rectangle_prob_cached(
      u1, u2, l1, l2, par, floor_probability = floor_probability
    ))
  }
  if (identical(family, "t")) {
    return(.copula_t_rectangle_prob_conditional(
      u1, u2, l1, l2, par, par2,
      floor_probability = floor_probability
    ))
  }
  rect <- .copula_cdf(u1, u2, family = family, par = par, par2 = par2) -
    .copula_cdf(l1, u2, family = family, par = par, par2 = par2) -
    .copula_cdf(u1, l2, family = family, par = par, par2 = par2) +
    .copula_cdf(l1, l2, family = family, par = par, par2 = par2)
  rect <- as.numeric(rect)
  if (isTRUE(floor_probability)) pmax(rect, 1e-300) else rect
}
