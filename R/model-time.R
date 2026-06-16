#' Prepare internal time index and subject identifiers for fitting
#' 
#' Check/update ordering for time. Convert time to numeric index and preserve original time covariate for formulas. 
#' Convert subject to character if factor. Check for missing/non-finite values in time and subject.
#'
#' @noRd
.gl_prepare_time_subject_columns <- function(dataset, time_var) {
  # Preserve the user-facing time covariate (including factor type) for formulas, while keeping an internal numeric time index for time-label-agnostic optimisation logic.

  dataset$time_covariate <- dataset$time
  time_covariate_is_factor <- is.factor(dataset$time_covariate)
  time_covariate_levels <- if (time_covariate_is_factor) levels(dataset$time_covariate) else NULL
  time_covariate_ordered <- if (time_covariate_is_factor) is.ordered(dataset$time_covariate) else FALSE

  if (is.factor(dataset$time_covariate)) {
    time_chr <- as.character(dataset$time_covariate)
    dataset$time <- match(time_chr, time_covariate_levels)

    if (anyNA(dataset$time)) {
      stop("ERROR: Failed to map factor time levels to internal numeric time index.")
    }

    dataset$time_covariate <- factor(time_chr, levels = time_covariate_levels, ordered = time_covariate_ordered)

    if (time_covariate_ordered) {
      time_contr <- contr.treatment(length(time_covariate_levels))

      if (length(time_covariate_levels) > 1) {
        colnames(time_contr) <- time_covariate_levels[-1]
      }

      contrasts(dataset$time_covariate) <- time_contr
    }
  } else if (is.numeric(dataset$time_covariate) || is.integer(dataset$time_covariate)) {
    dataset$time <- as.numeric(dataset$time_covariate)
  } else if (is.character(dataset$time_covariate)) {
    time_numeric <- suppressWarnings(as.numeric(dataset$time_covariate))

    if (anyNA(time_numeric)) {
      stop(
        "ERROR: time must be numeric-like unless supplied as factor.\n",
        "If time is categorical for formulas/interactions, convert it to factor before fitting."
      )
    }
    warning(
      "Converted character time variable '", time_var,
      "' to numeric for fitting; convert it to factor before fitting if visits are categorical.",
      call. = FALSE
    )
    dataset$time <- time_numeric
    dataset$time_covariate <- time_numeric
  } else {
    stop(
      "ERROR: Unsupported time variable type: ", class(dataset$time_covariate)[1],
      ". Use numeric/integer, numeric-like character, or factor."
    )
  }

  if (is.factor(dataset$subject)) {
    dataset$subject <- as.character(dataset$subject)
  }

  if (any(is.na(dataset$time)) || any(is.na(dataset$subject)) || any(!is.finite(dataset$time))) {
    stop("ERROR: time and subject variables cannot contain missing or non-finite values.")
  }

  list(
    dataset = dataset,
    time_covariate_is_factor = time_covariate_is_factor,
    time_covariate_levels = time_covariate_levels,
    time_covariate_ordered = time_covariate_ordered
  )
}
