.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion(pkgname)
  packageStartupMessage(
    sprintf(
      "This is %s %s. For overview type 'help(\"%s-package\")'.",
      pkgname,
      version,
      pkgname
    )
  )
  packageStartupMessage(
    "Warning: gamlss.longitudinal is in active development; please report issues at https://github.com/ahibbert/gamlss.longitudinal/issues"
  )
}
