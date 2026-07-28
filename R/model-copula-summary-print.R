#' @export

print.copula_time_summary <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCopula Dependence Summary\n")

  cat("-------------------------\n")

  if (!is.null(x$time_summary)) {
    cat("\nFitted dependence by time:\n")

    print(x$time_summary, digits = digits, row.names = FALSE)
  }

  if (!is.null(x$pair_summary)) {
    cat("\nAdjacent-pair dependence:\n")

    print(x$pair_summary, digits = digits, row.names = FALSE)
  }

  invisible(x)
}
