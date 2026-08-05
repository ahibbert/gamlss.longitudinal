#' @export
print.gamlss_longitudinal_missingness_check <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nResponse Missingness Check\n")
  cat("--------------------------\n")
  print(x$response, digits = digits, row.names = FALSE)
  cat("\nAssessment:", x$assessment, "\n")
  cat(x$message, "\n", sep = "")
  if (!is.null(x$terms) && nrow(x$terms) > 0L) {
    cat("\nPredictor tests:\n")
    display <- x$terms
    display$p_value <- .gl_format_summary_p_value(display$p_value, digits + 1)
    print(display, digits = digits, row.names = FALSE)
  }
  invisible(x)
}
