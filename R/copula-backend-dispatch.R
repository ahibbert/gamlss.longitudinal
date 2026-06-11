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


.copula_deriv <- function(u1, u2, family, par, par2 = 0, deriv, log = FALSE) {

  family <- .copula_family_code(family)

  if (.copula_backend() == "native" && family == "N" && deriv %in% c("u1", "u2", "par")) {

    return(.copula_gaussian_deriv(u1, u2, par, deriv = deriv, log = log))

  }

  if (.copula_backend() == "native" && family == "C" && deriv %in% c("u1", "u2", "par")) {

    return(.copula_clayton_deriv(u1, u2, par, deriv = deriv, log = log))

  }

  if (.copula_backend() == "native" && family == "G" && deriv %in% c("u1", "u2", "par")) {

    return(.copula_gumbel_deriv(u1, u2, par, deriv = deriv, log = log))

  }

  if (.copula_backend() == "native" && family == "F" && deriv %in% c("u1", "u2", "par")) {

    return(.copula_frank_deriv(u1, u2, par, deriv = deriv, log = log))

  }

  if (.copula_backend() == "native" && family == "J" && deriv %in% c("u1", "u2", "par")) {

    return(.copula_joe_deriv(u1, u2, par, deriv = deriv, log = log))

  }

  if (.copula_backend() == "native" && family == "t" && deriv %in% c("u1", "u2", "par", "par2")) {

    return(.copula_t_deriv(u1, u2, par, par2, deriv = deriv, log = log))

  }

  .copula_require_vinecopula("the delegated VineCopula first-derivative backend")

  VineCopula::BiCopDeriv(

    u1,

    u2,

    family = .copula_family_number(family),

    par = par,

    par2 = par2,

    deriv = deriv,

    log = log

  )

}


.copula_deriv2 <- function(u1, u2, family, par, par2 = 0, deriv) {

  family <- .copula_family_code(family)

  if (.copula_backend() == "native" && family == "N" && deriv %in% c("u1", "u2", "par")) {

    return(.copula_gaussian_deriv2(u1, u2, par, deriv = deriv))

  }

  if (.copula_backend() == "native" && family == "C" && deriv %in% c("u1", "u2", "par")) {

    return(.copula_clayton_deriv2(u1, u2, par, deriv = deriv))

  }

  if (.copula_backend() == "native" && family == "G" && deriv %in% c("u1", "u2", "par")) {

    return(.copula_gumbel_deriv2(u1, u2, par, deriv = deriv))

  }

  if (.copula_backend() == "native" && family == "F" && deriv %in% c("u1", "u2", "par")) {

    return(.copula_frank_deriv2(u1, u2, par, deriv = deriv))

  }

  if (.copula_backend() == "native" && family == "J" && deriv %in% c("u1", "u2", "par")) {

    return(.copula_joe_deriv2(u1, u2, par, deriv = deriv))

  }

  if (.copula_backend() == "native" && family == "t" && deriv %in% c("u1", "u2", "par", "par2", "par1par2")) {

    return(.copula_t_deriv2(u1, u2, par, par2, deriv = deriv))

  }

  .copula_require_vinecopula("the delegated VineCopula second-derivative backend")

  VineCopula::BiCopDeriv2(

    u1,

    u2,

    family = .copula_family_number(family),

    par = par,

    par2 = par2,

    deriv = deriv

  )

}


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


.copula_dvine <- function(order, family, par, par2 = 0) {

  .copula_require_vinecopula("D-vine simulation setup via .copula_dvine()")

  VineCopula::D2RVine(order, .copula_family_numbers(family), par, par2)

}


.copula_rvine_sim <- function(n, rvm) {

  .copula_require_vinecopula("R-vine simulation via .copula_rvine_sim()")

  VineCopula::RVineSim(n, rvm)

}
