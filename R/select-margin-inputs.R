#' Normalize inputs for marginal distribution screening
#'
#' @keywords internal
#' @noRd
.select_margin_inputs <- function(
    response = NULL,
    data = NULL,
    response_var = NULL,
    type = NULL,
    time_intercepts = FALSE,
    time_var = NULL) {
  if (is.null(data) && is.data.frame(response) && !is.null(response_var)) {
    data <- response
    response <- NULL
  }

  if (!is.null(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    if (is.null(response_var) || !is.character(response_var) || length(response_var) != 1L) {
      stop("'response_var' must be a single column name when 'data' is supplied.", call. = FALSE)
    }
    if (!response_var %in% names(data)) {
      stop("response_var='", response_var, "' not found in 'data'.", call. = FALSE)
    }
    response <- data[[response_var]]
  }

  if (is.null(response)) {
    stop("Supply either 'response' or both 'data' and 'response_var'.", call. = FALSE)
  }

  response_all <- as.numeric(response)

  if (isTRUE(time_intercepts)) {
    if (is.null(data)) {
      stop("'time_intercepts = TRUE' requires 'data' and 'time_var'.", call. = FALSE)
    }
    if (is.null(time_var) || !is.character(time_var) || length(time_var) != 1L) {
      stop("'time_var' must be a single column name when 'time_intercepts = TRUE'.", call. = FALSE)
    }
    if (!time_var %in% names(data)) {
      stop("time_var='", time_var, "' not found in 'data'.", call. = FALSE)
    }
  }

  response <- response_all[is.finite(response_all)]
  if (length(response) < 3L) {
    stop("Need at least three finite response values to screen margins.", call. = FALSE)
  }

  if (is.null(type)) {
    is_count <- all(response >= 0) && all(abs(response - round(response)) < .Machine$double.eps^0.5)
    type <- if (is_count) {
      "counts"
    } else if (all(response > 0)) {
      "realplus"
    } else {
      "realAll"
    }
  }

  list(
    response = response,
    response_all = response_all,
    data = data,
    response_var = response_var,
    type = type,
    time_intercepts = isTRUE(time_intercepts),
    time_var = time_var
  )
}
