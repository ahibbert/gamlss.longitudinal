#' Summarize PIT calibration for model checks
#'
#' @param pit Numeric PIT values, already clipped to `[0, 1]`.
#' @return One-row data frame with PIT calibration summary statistics.
#' @noRd
.gl_pit_calibration_stats <- function(pit) {
  ks_p <- tryCatch(stats::ks.test(pit, "punif")$p.value, error = function(e) NA_real_)

  data.frame(
    n = length(pit),
    mean = mean(pit, na.rm = TRUE),
    sd = stats::sd(pit, na.rm = TRUE),
    expected_sd = sqrt(1 / 12),
    ks_p_value = as.numeric(ks_p),
    stringsAsFactors = FALSE
  )
}

#' Summarize tail calibration for model checks
#'
#' @param pit Numeric PIT values, already clipped to `[0, 1]`.
#' @return List with full `tail_summary` table and compact `tail_stats` row.
#' @noRd
.gl_tail_calibration_stats <- function(pit) {
  tail_summary <- .gl_pit_tail_summary(pit)

  list(
    tail_summary = tail_summary,
    tail_stats = data.frame(
      tail_ratio_max = attr(tail_summary, "tail_ratio_max"),
      stringsAsFactors = FALSE
    )
  )
}
