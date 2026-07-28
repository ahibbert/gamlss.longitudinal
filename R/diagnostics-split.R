.gl_check_by_time <- function(by_time) {
  if (!is.logical(by_time) || length(by_time) != 1L || is.na(by_time)) {
    stop("'by_time' must be TRUE or FALSE.", call. = FALSE)
  }
}

.gl_resolve_diag_split <- function(by, diag_data, data = NULL, plot_name = "diagnostic") {
  if (is.null(by) || (is.character(by) && length(by) == 1L && !nzchar(by))) {
    return(NULL)
  }

  if (!is.character(by) || length(by) != 1L) {
    stop("'by' must be NULL or a single column name as a character string.", call. = FALSE)
  }

  if (by %in% c("time", "response_margin")) {
    return(as.factor(diag_data$time))
  }

  if (by %in% c("subject", "response_subject")) {
    return(as.factor(diag_data$subject))
  }

  if (is.null(data)) {
    stop("To split ", plot_name, " by '", by, "', provide data= containing that column.", call. = FALSE)
  }

  if (!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  if (!by %in% names(data)) {
    stop("Column '", by, "' not found in provided data.", call. = FALSE)
  }

  n_diag <- length(diag_data$response)

  if (nrow(data) == n_diag) {
    return(as.factor(data[[by]]))
  }

  keep_mask <- diag_data$keep_mask

  if (!is.null(keep_mask) && length(keep_mask) == nrow(data)) {
    vec <- data[[by]][keep_mask]

    if (length(vec) == n_diag) {
      return(as.factor(vec))
    }
  }

  keep_index <- diag_data$keep_index

  if (!is.null(keep_index) && length(keep_index) == n_diag && (n_diag == 0L || max(keep_index, na.rm = TRUE) <= nrow(data))) {
    return(as.factor(data[[by]][keep_index]))
  }

  stop(

    "Could not align data rows with ", plot_name, " rows for by='", by, "'. ",
    "Provide data with row count equal to either the diagnostic row count (", n_diag, ") or the original fit data rows.",
    call. = FALSE
  )
}

.gl_diag_split_info <- function(by_time, by, diag_data, data = NULL, plot_name = "diagnostic") {
  .gl_check_by_time(by_time)

  if (isTRUE(by_time) && is.null(by)) {
    by <- "time"
  } else if (isTRUE(by_time) && !is.null(by)) {
    warning("Both by_time and by were provided; using by='", by, "'.", call. = FALSE)
  }

  split_group <- .gl_resolve_diag_split(by, diag_data, data = data, plot_name = plot_name)

  list(by = by, group = split_group, split_by = !is.null(split_group))
}
