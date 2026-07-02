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
