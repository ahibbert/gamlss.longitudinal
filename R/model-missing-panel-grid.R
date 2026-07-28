#' Build the internal time lookup used when expanding fit panels
#'
#' @param dataset Prepared fitting data with `time` and `time_covariate` columns.
#' @return Data frame mapping internal numeric time to preserved time covariates.
#' @noRd
.gl_fit_time_lookup <- function(dataset) {
  time_lookup <- dataset[!duplicated(dataset$time), c("time", "time_covariate"), drop = FALSE]
  time_lookup[order(time_lookup$time), , drop = FALSE]
}

#' Build the full subject-time grid for structural missingness expansion
#'
#' @param dataset Prepared fitting data with `subject` and `time` columns.
#' @return Data frame containing all observed subjects crossed with all observed time points.
#' @noRd
.gl_fit_full_panel_grid <- function(dataset) {
  expand.grid(
    subject = sort(unique(dataset$subject)),
    time = sort(unique(dataset$time)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}
