.copula_pdf <- function(u1, u2, family, par, par2 = 0) {
  family <- .copula_family_code(family)

  if (.copula_backend() == "native" && family == "N") {
    return(.copula_gaussian_pdf(u1, u2, par))
  }

  if (.copula_backend() == "native" && family == "C") {
    return(.copula_clayton_pdf(u1, u2, par))
  }

  if (.copula_backend() == "native" && family == "F") {
    return(.copula_frank_pdf(u1, u2, par))
  }

  if (.copula_backend() == "native" && family == "G") {
    return(.copula_gumbel_pdf(u1, u2, par))
  }

  if (.copula_backend() == "native" && family == "J") {
    return(.copula_joe_pdf(u1, u2, par))
  }

  if (.copula_backend() == "native" && family == "t") {
    return(.copula_t_pdf(u1, u2, par, par2))
  }

  .copula_require_vinecopula("the delegated VineCopula PDF backend")

  VineCopula::BiCopPDF(

    u1,
    u2,
    family = .copula_family_number(family),
    par = par,
    par2 = par2
  )
}

#' Evaluate a copula log density without avoidable density-scale underflow
#'
#' @noRd
.copula_logpdf <- function(u1, u2, family, par, par2 = 0) {
  family <- .copula_family_code(family)

  if (.copula_backend() == "native" && family == "N") {
    return(.copula_gaussian_logpdf(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "C") {
    return(.copula_clayton_logpdf(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "t") {
    return(.copula_t_logpdf(u1, u2, par, par2))
  }

  log(.copula_pdf(u1, u2, family = family, par = par, par2 = par2))
}
