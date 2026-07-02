.copula_hfunc1 <- function(u1, u2, family, par, par2 = 0) {
  family <- .copula_family_code(family)

  if (.copula_backend() == "native" && family == "N") {
    return(.copula_gaussian_hfunc1(u1, u2, par))
  }

  if (.copula_backend() == "native" && family == "C") {
    return(.copula_clayton_hfunc1(u1, u2, par))
  }

  if (.copula_backend() == "native" && family == "F") {
    return(.copula_frank_hfunc1(u1, u2, par))
  }

  if (.copula_backend() == "native" && family == "G") {
    return(.copula_gumbel_hfunc1(u1, u2, par))
  }

  if (.copula_backend() == "native" && family == "J") {
    return(.copula_joe_hfunc1(u1, u2, par))
  }

  if (.copula_backend() == "native" && family == "t") {
    return(.copula_t_hfunc1(u1, u2, par, par2))
  }

  .copula_require_vinecopula("the delegated VineCopula h-function backend")

  VineCopula::BiCopHfunc1(

    u1,
    u2,
    family = .copula_family_number(family),
    par = par,
    par2 = par2
  )
}
