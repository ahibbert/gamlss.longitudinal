.copula_cdf <- function(u1, u2, family, par, par2 = 0) {
  family <- .copula_family_code(family)

  if (.copula_backend() == "native" && family == "N" && all(abs(.copula_gaussian_rho(par)) <= 1e-12)) {
    vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), par)

    return(vals[[1]] * vals[[2]])
  }

  if (.copula_backend() == "native" && family == "C") {
    return(.copula_clayton_cdf(u1, u2, par))
  }

  if (.copula_backend() == "native" && family == "F") {
    return(.copula_frank_cdf(u1, u2, par))
  }

  if (.copula_backend() == "native" && family == "G") {
    return(.copula_gumbel_cdf(u1, u2, par))
  }

  if (.copula_backend() == "native" && family == "J") {
    return(.copula_joe_cdf(u1, u2, par))
  }

  if (.copula_backend() == "native" && family == "t") {
    return(.copula_t_cdf(u1, u2, par, par2))
  }

  .copula_require_vinecopula("the delegated VineCopula CDF backend")

  VineCopula::BiCopCDF(

    u1,
    u2,
    family = .copula_family_number(family),
    par = par,
    par2 = par2
  )
}
