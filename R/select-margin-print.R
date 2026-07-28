#' @rdname select_margin

#' @export

screen_margin <- function(...) {
  select_margin(...)
}

#' @export

print.margin_screen <- function(x, ..., n = 10L) {
  cat("\nMarginal Distribution Screen\n")

  cat("----------------------------\n")

  if (nrow(x) == 0L) {
    cat("No candidate families were retained.\n")

    return(invisible(x))
  }

  cat("Selected:", attr(x, "selected"), "\n\n")

  display <- utils::head(as.data.frame(x), n = n)

  unsupported <- "supported_by_longitudinal" %in% names(display) &

    any(!display$supported_by_longitudinal, na.rm = TRUE)

  display$supported_by_longitudinal <- NULL

  print(display, row.names = FALSE)

  if ("converged" %in% names(display) && any(!display$converged, na.rm = TRUE)) {
    cat("\nWarning: one or more printed time-intercept marginal fits did not report convergence.\n")
  }

  if (isTRUE(unsupported)) {
    cat("\nNote: one or more printed families are not currently supported by gamlss_longitudinal().\n")
  }

  invisible(x)
}
