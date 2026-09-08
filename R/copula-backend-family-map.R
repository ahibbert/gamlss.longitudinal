#' Validate a public copula family code
#'
#' @noRd
.copula_family_code <- function(family) {
  if (!is.character(family) || length(family) != 1L || is.na(family)) {
    .gl_capability_stop(
      "gamlss_longitudinal_unsupported_copula_error",
      "Copula family must be a single character code."
    )
  }
  if (is.null(.gl_capability_copula_spec(family))) {
    .gl_capability_stop(
      "gamlss_longitudinal_unsupported_copula_error",
      paste0(
        "Unsupported copula family code '", family,
        "'. Use one of: ", paste(.gl_capability_all_copulas(), collapse = ", "), "."
      ),
      copula = family
    )
  }
  family
}

#' Map a package copula code to a VineCopula family number
#'
#' @noRd
.copula_family_number <- function(family) {
  family <- .copula_family_code(family)
  switch(family,
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
