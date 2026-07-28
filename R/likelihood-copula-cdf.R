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
#' Evaluates the Gaussian copula CDF with boundary shortcuts and duplicated-row
#' caching. Uses `mvtnorm` when available, otherwise falls back to `VineCopula`.
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
    key <- paste(signif(u1[remaining], 15), signif(u2[remaining], 15), signif(rho[remaining], 15), sep = "\r")
    unique_key <- !duplicated(key)
    u1_unique <- u1[remaining][unique_key]
    u2_unique <- u2[remaining][unique_key]
    rho_unique <- rho[remaining][unique_key]
    if (requireNamespace("mvtnorm", quietly = TRUE)) {
      z1 <- stats::qnorm(u1_unique)
      z2 <- stats::qnorm(u2_unique)
      unique_vals <- vapply(seq_along(rho_unique), function(ii) {
        as.numeric(mvtnorm::pmvnorm(
          upper = c(z1[[ii]], z2[[ii]]),
          corr = matrix(c(1, rho_unique[[ii]], rho_unique[[ii]], 1), 2, 2)
        ))
      }, numeric(1))
    } else {
      .copula_require_vinecopula("the Gaussian copula CDF backend")
      unique_vals <- VineCopula::BiCopCDF(
        u1_unique,
        u2_unique,
        family = 1,
        par = rho_unique
      )
    }
    out[remaining] <- unique_vals[match(key, key[unique_key])]
  }

  pmin(pmax(as.numeric(out), 0), 1)
}

#' Cached Gaussian rectangle probability
#'
#' Computes discrete/count likelihood rectangle masses from four cached
#' Gaussian CDF evaluations.
#'
#' @noRd
.copula_gaussian_rectangle_prob_cached <- function(u1, u2, l1, l2, par) {
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
  pmax(as.numeric(rect), 1e-300)
}

#' Copula rectangle probability
#'
#' Returns positive lower-bounded rectangle probabilities for discrete margins,
#' using the cached Gaussian path where possible and generic inclusion-exclusion
#' for other copulas.
#'
#' @noRd
.copula_rectangle_prob <- function(u1, u2, l1, l2, family, par, par2 = 0) {
  family <- .copula_family_code(family)
  if (identical(family, "N")) {
    return(.copula_gaussian_rectangle_prob_cached(u1, u2, l1, l2, par))
  }
  rect <- .copula_cdf(u1, u2, family = family, par = par, par2 = par2) -
    .copula_cdf(l1, u2, family = family, par = par, par2 = par2) -
    .copula_cdf(u1, l2, family = family, par = par, par2 = par2) +
    .copula_cdf(l1, l2, family = family, par = par, par2 = par2)
  pmax(as.numeric(rect), 1e-300)
}
