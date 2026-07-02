.joint_selection_margin_candidates <- function(
    data,
    response_var,
    time_var,
    type,
    margin_families,
    time_intercepts,
    try.gamlss,
    trace) {
  if (!is.null(margin_families)) {
    if (!is.character(margin_families) || any(is.na(margin_families)) || any(!nzchar(margin_families))) {
      stop("'margin_families' must be a character vector of gamlss.dist family names.", call. = FALSE)
    }
    margin_families <- unique(margin_families)
    if (is.null(type)) {
      response <- as.numeric(data[[response_var]])
      response <- response[is.finite(response)]
      is_count <- all(response >= 0) && all(abs(response - round(response)) < .Machine$double.eps^0.5)
      type <- if (is_count) {
        "counts"
      } else if (all(response > 0)) {
        "realplus"
      } else {
        "realAll"
      }
    }
    out <- data.frame(
      family = margin_families,
      AIC = NA_real_,
      type = type,
      stringsAsFactors = FALSE
    )
    out$supported_by_longitudinal <- vapply(out$family, .joint_selection_margin_supported, logical(1))
    out$rank <- seq_len(nrow(out))
    out$delta_AIC <- NA_real_
    attr(out, "selected") <- if (nrow(out) > 0L) out$family[[1L]] else NA_character_
    attr(out, "response_type") <- type
    attr(out, "time_intercepts") <- isTRUE(time_intercepts)
    attr(out, "time_var") <- if (isTRUE(time_intercepts)) time_var else NULL
    class(out) <- c("margin_selection", "margin_screen", "data.frame")
    return(out)
  }

  margin_call <- quote(select_margin(
    data = data,
    response_var = response_var,
    time_var = if (isTRUE(time_intercepts)) time_var else NULL,
    type = type,
    families = NULL,
    time_intercepts = time_intercepts,
    try.gamlss = try.gamlss,
    trace = trace
  ))
  if (isTRUE(trace)) {
    eval(margin_call)
  } else {
    warnings <- character(0)
    withCallingHandlers(
      {
        utils::capture.output(ans <- eval(margin_call))
        ans
      },
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  }
}

.joint_selection_margin_supported <- function(family) {
  !is.null(.margin_family_object(family)) &&
    all(vapply(
      paste0(c("d", "p", "q"), family),
      exists,
      logical(1),
      envir = asNamespace("gamlss.dist"),
      inherits = FALSE
    ))
}

.joint_selection_check_column <- function(data, column, arg) {
  if (is.null(column) || !is.character(column) || length(column) != 1L || is.na(column) || !nzchar(column)) {
    stop("'", arg, "' must be a single column name.", call. = FALSE)
  }
  if (!column %in% names(data)) {
    stop(arg, "='", column, "' not found in 'data'.", call. = FALSE)
  }
  invisible(TRUE)
}

.joint_selection_count_pairs <- function(data, response_var, time_var, subject_var) {
  ord <- order(data[[subject_var]], data[[time_var]])
  data <- data[ord, , drop = FALSE]
  counts <- vapply(split(data, data[[subject_var]], drop = TRUE), function(subject_data) {
    subject_data <- subject_data[order(subject_data[[time_var]]), , drop = FALSE]
    if (nrow(subject_data) < 2L) {
      return(0L)
    }
    left <- seq_len(nrow(subject_data) - 1L)
    right <- left + 1L
    sum(is.finite(subject_data[[response_var]][left]) & is.finite(subject_data[[response_var]][right]))
  }, integer(1))
  sum(counts)
}
