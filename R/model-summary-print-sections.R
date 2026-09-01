#' Print smooth-term summary section
#'
#' @param smooth_terms Smooth-term table from `.gl_summary_smooth_terms()`.
#' @param digits Number of decimal places used by the print method.
#' @return Invisibly returns `smooth_terms`.
#' @noRd
.gl_print_summary_smooth_terms <- function(smooth_terms, digits) {
  cat("\nSmooth terms:\n")
  cat("--------------------\n")
  if (!is.null(smooth_terms) && nrow(smooth_terms) > 0) {
    smooth_disp <- smooth_terms
    if ("edf" %in% names(smooth_disp)) {
      smooth_disp$edf <- round(smooth_disp$edf, digits)
    }
    print(smooth_disp, row.names = FALSE)
    cat("Use plot(object) to visualize smooth and fixed terms with confidence bands.\n")
  } else {
    cat("None\n")
  }

  invisible(smooth_terms)
}

#' Print model-selection summary section
#'
#' @param fit Summary fit list from `.gl_build_summary_object()`.
#' @param digits Number of decimal places used by the print method.
#' @return Invisibly returns `fit`.
#' @noRd
.gl_print_summary_model_selection <- function(fit, digits) {
  heading <- if (identical(fit$criteria_status, "provisional_nonconverged")) {
    "Model Selection Criteria (PROVISIONAL; nonconverged fit):"
  } else {
    "Model Selection Criteria:"
  }
  cat("\n", heading, "\n", sep = "")
  cat("--------------------\n")
  if (!is.null(fit$model_selection)) {
    print(round(fit$model_selection, digits))
  } else {
    fit_tbl <- data.frame(
      metric = c("logLik", "AIC", "BIC"),
      value = c(fit$logLik, fit$AIC, fit$BIC),
      stringsAsFactors = FALSE
    )
    print(fit_tbl, row.names = FALSE, digits = digits)
  }

  invisible(fit)
}
