#' Resolve the configured copula backend

#'

#' Reads `options(gamlss.longitudinal.copula_backend)` and normalises the value

#' to either `"native"` or `"vinecopula"`.

#'

#' @noRd

.copula_backend <- function() {

  backend <- getOption("gamlss.longitudinal.copula_backend", "native")

  backend <- tolower(as.character(backend)[1])

  if (!backend %in% c("vinecopula", "native")) {

    stop("Unknown copula backend: ", backend, call. = FALSE)

  }

  backend

}


#' Require the optional VineCopula backend

#'

#' Gives a contextual error when a requested operation needs `VineCopula` and it

#' is not installed.

#'

#' @noRd

.copula_require_vinecopula <- function(context = "this copula operation") {

  if (!requireNamespace("VineCopula", quietly = TRUE)) {

    stop(

      "VineCopula is required for ", context,

      ". Install VineCopula or use options(gamlss.longitudinal.copula_backend = 'native').",

      call. = FALSE

    )

  }

  invisible(TRUE)

}


#' Validate a public copula family code

#'

#' @noRd

.copula_family_code <- function(family) {

  if (!is.character(family) || length(family) != 1L || is.na(family)) {

    stop("Copula family must be a single character code.", call. = FALSE)

  }

  if (!family %in% c("N", "C", "F", "G", "J", "t")) {

    stop(

      "Unsupported copula family code '", family,

      "'. Use one of: N, C, F, G, J, t.",

      call. = FALSE

    )

  }

  family

}


#' Map a package copula code to a VineCopula family number

#'

#' @noRd

.copula_family_number <- function(family) {

  family <- .copula_family_code(family)

  switch(

    family,

    N = 1,

    C = 3,

    F = 5,

    G = 4,

    J = 6,

    t = 2

  )

}


#' Vectorised copula family number mapping

#'

#' Keeps independence code `0` available for diagnostic/reference calls.

#'

#' @noRd

.copula_family_numbers <- function(family) {

  vapply(family, function(x) {

    if (identical(x, 0) || identical(x, "0")) {

      return(0)

    }

    .copula_family_number(as.character(x))

  }, numeric(1), USE.NAMES = FALSE)

}


#' Clamp Gaussian copula correlation away from singular boundaries

#'

#' @noRd

.copula_gaussian_rho <- function(par) {

  pmin(pmax(as.numeric(par), -0.999999), 0.999999)

}


#' Clamp pseudo-observations to the open unit interval

#'

#' @noRd

.copula_clamp01 <- function(u) {

  pmin(pmax(as.numeric(u), 1e-12), 1 - 1e-12)

}


#' Recycle copula vector inputs to common length

#'

#' @noRd

.copula_recycle <- function(...) {

  args <- list(...)

  n <- max(vapply(args, length, integer(1)))

  lapply(args, rep, length.out = n)

}


#' Zero derivative helper for independence-limit copula paths

#'

#' @noRd

.copula_indep_deriv <- function(u1, deriv, log = FALSE) {

  rep(0, length(u1))

}


.copula_one_sided_par_deriv <- function(pdf_fun, u1, u2, par0, h, log = FALSE) {

  dens0 <- rep(1, length(u1))

  dens1 <- pdf_fun(u1, u2, par0 + h)

  if (isTRUE(log)) {

    (log(dens1) - log(dens0)) / h

  } else {

    (dens1 - dens0) / h

  }

}


.copula_central_par_deriv <- function(pdf_fun, u1, u2, par0, h, log = FALSE) {

  dens_plus <- pdf_fun(u1, u2, par0 + h)

  dens_minus <- pdf_fun(u1, u2, par0 - h)

  if (isTRUE(log)) {

    (log(dens_plus) - log(dens_minus)) / (2 * h)

  } else {

    (dens_plus - dens_minus) / (2 * h)

  }

}


.copula_one_sided_par_deriv2 <- function(pdf_fun, u1, u2, par0, h) {

  (pdf_fun(u1, u2, par0 + 2 * h) - 2 * pdf_fun(u1, u2, par0 + h) + 1) / h^2

}


.copula_central_par_deriv2 <- function(pdf_fun, u1, u2, par0, h) {

  (pdf_fun(u1, u2, par0 + h) - 2 + pdf_fun(u1, u2, par0 - h)) / h^2

}


