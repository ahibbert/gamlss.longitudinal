#' Extract a numeric column from copula-to-margin derivative input
#'
#' @keywords internal
#' @noRd
.copula_margin_derivative_numeric_column <- function(input, nm) {
  val <- if (is.data.frame(input)) input[[nm]] else input[, nm]
  suppressWarnings(as.numeric(val))
}
