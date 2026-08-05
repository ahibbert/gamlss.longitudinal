#' Format summary p-values for console display
#'
#' @param x Numeric p-values.
#' @param digits Number of decimal places used for p-values.
#' @return A character vector with small positive p-values shown as less than
#'   the display threshold.
#' @noRd
.gl_format_summary_p_value <- function(x, digits) {
  x_raw <- x
  x <- round(x, digits)
  fmt_num <- function(v, d) ifelse(is.na(v), "NA", formatC(v, format = "f", digits = d))
  ifelse(
    is.na(x_raw),
    "NA",
    ifelse(x_raw > 0 & x_raw < 10^(-digits), paste0("<", formatC(10^(-digits), format = "f", digits = digits)), fmt_num(x, digits))
  )
}

#' Format summary coefficients for console display
#'
#' @param coef_tbl Coefficient table from `.gl_summary_coefficient_table()`.
#' @param digits Number of decimal places used by the print method.
#' @return A character data frame ready for aligned console printing.
#' @noRd
.gl_summary_coefficient_display <- function(coef_tbl, digits) {
  coef_tbl <- coef_tbl
  coef_tbl$estimate <- round(coef_tbl$estimate, digits)
  coef_tbl$std_error <- round(coef_tbl$std_error, digits)

  fmt_num <- function(v, d) ifelse(is.na(v), "NA", formatC(v, format = "f", digits = d))
  data.frame(
    term = as.character(coef_tbl$term),
    estimate = fmt_num(coef_tbl$estimate, digits),
    std_error = fmt_num(coef_tbl$std_error, digits),
    p_value = .gl_format_summary_p_value(coef_tbl$p_value, digits + 1),
    signif = ifelse(is.na(coef_tbl$signif), "", as.character(coef_tbl$signif)),
    parameter = as.character(coef_tbl$parameter),
    stringsAsFactors = FALSE
  )
}

#' Print grouped summary coefficient blocks
#'
#' @param coef_disp Character display table from `.gl_summary_coefficient_display()`.
#' @param param_order Preferred parameter print order.
#' @param prefix Prefix inserted before each table row.
#' @return Invisibly returns `coef_disp`.
#' @noRd
.gl_print_summary_coefficient_blocks <- function(
    coef_disp,
    param_order = c("mu", "sigma", "nu", "tau", "theta", "zeta"),
    prefix = "    ") {
  w_term <- max(nchar("term"), nchar(coef_disp$term, type = "width"), na.rm = TRUE)
  w_est <- max(nchar("estimate"), nchar(coef_disp$estimate, type = "width"), na.rm = TRUE)
  w_se <- max(nchar("std_error"), nchar(coef_disp$std_error, type = "width"), na.rm = TRUE)
  w_p <- max(nchar("p_value"), nchar(coef_disp$p_value, type = "width"), na.rm = TRUE)
  w_sig <- max(nchar("signif"), nchar(coef_disp$signif, type = "width"), na.rm = TRUE)

  format_row <- function(term, estimate, std_error, p_value, signif) {
    sprintf(
      "%-*s  %*s  %*s  %*s  %-*s",
      w_term, term,
      w_est, estimate,
      w_se, std_error,
      w_p, p_value,
      w_sig, signif
    )
  }

  print_coef_block <- function(block, prefix = "    ") {
    hdr <- format_row("term", "estimate", "std_error", "p_value", "signif")
    cat(prefix, hdr, "\n", sep = "")
    for (ii in seq_len(nrow(block))) {
      row_txt <- format_row(
        block$term[ii],
        block$estimate[ii],
        block$std_error[ii],
        block$p_value[ii],
        block$signif[ii]
      )
      cat(prefix, row_txt, "\n", sep = "")
    }
  }

  params_present <- unique(coef_disp$parameter)
  params_print <- c(param_order[param_order %in% params_present], setdiff(params_present, param_order))

  for (k in seq_along(params_print)) {
    p <- params_print[k]
    block <- coef_disp[coef_disp$parameter == p, c("term", "estimate", "std_error", "p_value", "signif"), drop = FALSE]
    cat(sprintf("  [%s]\n", p))
    print_coef_block(block, prefix = prefix)
    if (k < length(params_print)) cat("  --------------------\n")
  }

  invisible(coef_disp)
}
