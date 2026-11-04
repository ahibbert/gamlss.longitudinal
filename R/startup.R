.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion(pkgname)
  packageStartupMessage(sprintf(
    "%s %s - longitudinal GAMLSS models with copula dependence",
    pkgname,
    version
  ))
}
