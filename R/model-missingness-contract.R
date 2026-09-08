#' Describe contiguous observed-response segments within subjects
#'
#' @noRd
.gl_observation_segments <- function(response, time, subject) {
  n <- length(response)
  if (length(time) != n || length(subject) != n) {
    stop("'response', 'time', and 'subject' must have the same length.", call. = FALSE)
  }

  segment_id <- rep(NA_integer_, n)
  subject_rows <- split(seq_len(n), subject, drop = TRUE)
  subject_summary <- vector("list", length(subject_rows))
  next_segment <- 0L

  for (i in seq_along(subject_rows)) {
    rows <- subject_rows[[i]]
    rows <- rows[order(time[rows])]
    observed <- !is.na(response[rows])
    observed_positions <- which(observed)

    if (length(observed_positions) == 0L) {
      local_segments <- 0L
      interior_missing <- 0L
      gap_runs <- 0L
      leading_missing <- 0L
    } else {
      local_segment <- cumsum(c(TRUE, diff(observed_positions) > 1L))
      local_segments <- max(local_segment)
      segment_id[rows[observed_positions]] <- next_segment + local_segment
      next_segment <- next_segment + local_segments

      interior <- seq.int(min(observed_positions), max(observed_positions))
      interior_missing_pattern <- !observed[interior]
      interior_missing <- sum(interior_missing_pattern)
      gap_runs <- if (interior_missing == 0L) {
        0L
      } else {
        sum(interior_missing_pattern & !c(FALSE, head(interior_missing_pattern, -1L)))
      }
      leading_missing <- min(observed_positions) - 1L
    }

    subject_summary[[i]] <- data.frame(
      subject = as.character(subject[rows[1L]]),
      n_scheduled = length(rows),
      n_observed = sum(observed),
      n_segments = as.integer(local_segments),
      n_interior_missing = as.integer(interior_missing),
      n_gap_runs = as.integer(gap_runs),
      n_leading_missing = as.integer(leading_missing),
      stringsAsFactors = FALSE
    )
  }

  subject_summary <- if (length(subject_summary)) {
    do.call(rbind, subject_summary)
  } else {
    data.frame(
      subject = character(), n_scheduled = integer(), n_observed = integer(),
      n_segments = integer(), n_interior_missing = integer(), n_gap_runs = integer(),
      n_leading_missing = integer()
    )
  }
  rownames(subject_summary) <- NULL

  list(segment_id = segment_id, subject_summary = subject_summary)
}

#' Build and enforce the intermittent-gap likelihood contract
#'
#' @noRd
.gl_missingness_contract <- function(dataset, missingness = c("error", "segment")) {
  missingness <- match.arg(missingness)
  segments <- .gl_observation_segments(
    response = dataset$response,
    time = dataset$time,
    subject = dataset$subject
  )
  subject_summary <- segments$subject_summary
  gap_affected <- subject_summary$n_gap_runs > 0L
  delayed_entry_affected <- subject_summary$n_leading_missing > 0L
  affected <- gap_affected | delayed_entry_affected
  has_gaps <- any(gap_affected)
  has_delayed_entry <- any(delayed_entry_affected)
  requires_segment <- any(affected)

  out <- list(
    requested = missingness,
    objective = if (requires_segment) "segmented" else "ordinary",
    has_intermittent_gaps = has_gaps,
    has_delayed_entry = has_delayed_entry,
    n_subjects = nrow(subject_summary),
    n_subjects_with_gaps = sum(gap_affected),
    n_subjects_with_delayed_entry = sum(delayed_entry_affected),
    n_observed = sum(subject_summary$n_observed),
    n_missing_response = sum(subject_summary$n_scheduled - subject_summary$n_observed),
    n_interior_missing = sum(subject_summary$n_interior_missing),
    n_gap_runs = sum(subject_summary$n_gap_runs),
    n_leading_missing = sum(subject_summary$n_leading_missing),
    n_segments = sum(subject_summary$n_segments),
    n_within_segment_pairs = sum(pmax(subject_summary$n_observed - subject_summary$n_segments, 0L)),
    n_between_segment_transitions_omitted = sum(pmax(subject_summary$n_segments - 1L, 0L)),
    affected_subjects = subject_summary$subject[affected],
    subject_summary = subject_summary,
    segment_id = segments$segment_id,
    between_segment_assumption = if (requires_segment) "independent" else "not_applicable",
    criteria_status = if (requires_segment) "available_under_segment_independence" else "available",
    statement = if (requires_segment) {
      paste0(
        "The likelihood begins at each explicitly admitted observed segment; ",
        "segments separated by intermittent gaps are treated as independent."
      )
    } else {
      "The ordinary first-order likelihood is evaluated over the observed response sequence."
    },
    future_support = if (requires_segment) {
      "Numerical integration of dependence across intermittent gaps and a distinct delayed-entry likelihood are not implemented."
    } else {
      NA_character_
    }
  )
  class(out) <- "gamlss_longitudinal_missingness_contract"

  if (requires_segment && identical(missingness, "error")) {
    stop(structure(
      c(
        list(
          message = paste0(
            "Intermittent observation gaps or leading unobserved scheduled visits were detected for ",
            length(out$affected_subjects),
            " subject(s). The default likelihood is restricted to complete observed panels ",
            "and observed prefixes followed only by dropout. ",
            "Refit with missingness = \"segment\" to treat contiguous observed segments ",
            "as independent and explicitly admit the first observed visit after leading missing visits."
          ),
          call = NULL
        ),
        out
      ),
      class = c("gamlss.longitudinal_gap_error", "error", "condition")
    ))
  }

  out
}

#' Warn when a segmented likelihood has been explicitly requested
#'
#' @noRd
.gl_warn_segmented_missingness <- function(contract) {
  if (!identical(contract$objective, "segmented")) return(invisible(FALSE))

  warning(structure(
    c(
      list(
        message = paste0(
          "Intermittent observation gaps or leading unobserved scheduled visits were detected for ",
          length(contract$affected_subjects),
          " subject(s). Using the segmented likelihood requested by ",
          "missingness = \"segment\": observations in different contiguous segments ",
          "are treated as independent. AIC and BIC remain available under this explicit assumption. ",
          "Ordinary model-Hessian inference is unavailable for this composite objective; ",
          "use cluster-sandwich or bootstrap sensitivity analyses. ",
          "Numerical integration of dependence across gaps and a distinct delayed-entry likelihood are not implemented; ",
          "a future integration engine may add this support."
        ),
        call = NULL
      ),
      contract
    ),
    class = c("gamlss.longitudinal_segment_warning", "warning", "condition")
  ))
  invisible(TRUE)
}

#' Enforce inference methods supported by the missingness objective
#'
#' @noRd
.gl_require_supported_missingness_inference <- function(object, method) {
  objective <- object$likelihood_contract$objective %||%
    object$missingness$objective %||% "ordinary"
  if (!identical(objective, "segmented") || identical(method, "sandwich")) {
    return(invisible(TRUE))
  }

  stop(structure(
    list(
      message = paste0(
        "Model-Hessian inference is unavailable for a segmented likelihood. ",
        "Different observed segments were treated as independent. Use ",
        "vcov(object, method = \"sandwich\") or bootstrap sensitivity analysis."
      ),
      call = NULL,
      method = method,
      objective = objective
    ),
    class = c(
      "gamlss_longitudinal_segmented_inference_error",
      "error",
      "condition"
    )
  ))
}
