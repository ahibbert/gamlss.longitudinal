#' @export
print.gamlss_longitudinal_check <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  fmt_num <- function(value) {
    if (length(value) == 0L || is.null(value) || !is.finite(value)) {
      return("n/a")
    }

    formatC(value, digits = digits, format = "fg")
  }

  section <- function(title) {
    cat("\n", title, "\n", paste(rep("-", nchar(title)), collapse = ""), "\n", sep = "")
  }

  cat("\n")

  cat("GAMLSS Longitudinal Model Check\n")

  cat("===============================\n")

  cat("Margin: ", x$model$margin_dist, "    Copula: ", x$model$copula_dist, "\n", sep = "")

  cat("LogLik: ", fmt_num(x$fit$logLik), "    Converged: ", if (isTRUE(x$convergence$converged)) "yes" else "no", "\n", sep = "")

  if (!is.null(x$basic_checks)) {
    section("Basic Checks")

    basic_display <- x$basic_checks[, c("area", "status"), drop = FALSE]

    names(basic_display) <- c("Area", "Status")

    print(basic_display, row.names = FALSE, right = FALSE)

    cat("\nSummary: ", x$diagnostic_summary, "\n", sep = "")

    cat("Note: calibration and dependence rows are descriptive; they are not package pass/fail verdicts.\n")
  }

  section("Scores")

  print(x$scores, digits = digits, row.names = FALSE)

  section("PIT")

  print(x$pit, digits = digits, row.names = FALSE)

  if (!is.null(x$tail)) {
    section("Tail Calibration")

    print(x$tail, digits = digits, row.names = FALSE)
  }

  section("Residual Dependence")

  print(x$residual_dependence, digits = digits, row.names = FALSE)

  cat("\nUse check$checks for interpretation and check$flags for convergence or user-threshold review items.\n")

  invisible(x)
}
