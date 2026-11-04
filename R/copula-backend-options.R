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
