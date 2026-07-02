.copula_tau_to_par <- function(family, tau) {
  family <- .copula_family_code(family)

  if (.copula_backend() == "native" && family == "N") {
    return(sin(pi * tau / 2))
  }

  if (.copula_backend() == "native" && family == "C") {
    tau <- pmin(pmax(as.numeric(tau), 0), 1 - 1e-12)

    return(2 * tau / (1 - tau))
  }

  if (.copula_backend() == "native" && family == "F") {
    return(.copula_frank_tau_to_par(tau))
  }

  if (.copula_backend() == "native" && family == "G") {
    tau <- pmin(pmax(as.numeric(tau), 0), 1 - 1e-12)

    return(1 / (1 - tau))
  }

  if (.copula_backend() == "native" && family == "J") {
    return(.copula_joe_tau_to_par(tau))
  }

  if (.copula_backend() == "native" && family == "t") {
    return(sin(pi * tau / 2))
  }

  .copula_require_vinecopula("the delegated VineCopula tau-to-parameter backend")

  VineCopula::BiCopTau2Par(
    family = .copula_family_number(family),
    tau = tau
  )
}


.copula_par_to_tau <- function(family, par, par2 = 0) {
  family <- .copula_family_code(family)

  if (.copula_backend() == "native" && family == "N") {
    return(2 / pi * asin(.copula_gaussian_rho(par)))
  }

  if (.copula_backend() == "native" && family == "C") {
    theta <- .copula_clayton_theta(par)

    return(theta / (theta + 2))
  }

  if (.copula_backend() == "native" && family == "F") {
    return(.copula_frank_par_to_tau(par))
  }

  if (.copula_backend() == "native" && family == "G") {
    theta <- .copula_gumbel_theta(par)

    return(1 - 1 / theta)
  }

  if (.copula_backend() == "native" && family == "J") {
    return(.copula_joe_par_to_tau(par))
  }

  if (.copula_backend() == "native" && family == "t") {
    return(2 / pi * asin(.copula_gaussian_rho(par)))
  }

  .copula_require_vinecopula("the delegated VineCopula parameter-to-tau backend")

  VineCopula::BiCopPar2Tau(
    family = .copula_family_number(family),
    par = par,
    par2 = par2
  )
}
