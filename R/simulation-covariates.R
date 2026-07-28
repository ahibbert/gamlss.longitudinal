#' Generate covariates for longitudinal simulations
#'
#' This helper is designed to be used inside the `covariates` argument of
#' [simulate_longitudinal_dataset()]. Subject-level variables are evaluated once
#' per subject and repeated across time; observation-level variables are
#' evaluated once per row of the long-format design.
#'
#' @param data Base long-format design data supplied by
#'   [simulate_longitudinal_dataset()] to a covariate callback.
#' @param subject Named list of subject-level covariate specifications.
#' @param observation Named list of observation-level covariate specifications.
#'
#' @return A data frame of covariates aligned row-for-row with `data`.
#' @export
simulate_longitudinal_covariates <- function(data, subject = list(), observation = list()) {
  if (!is.data.frame(data) || !".sim_subject_index" %in% names(data)) {
    stop(
      "data must be the base design supplied by simulate_longitudinal_dataset().",
      call. = FALSE
    )
  }
  out <- data.frame(row.names = seq_len(nrow(data)))
  subject_rows <- !duplicated(data$.sim_subject_index)
  subject_data <- data[subject_rows, , drop = FALSE]

  for (nm in names(subject)) {
    values <- .sim_eval_covariate_spec(subject[[nm]], subject_data, nm)
    if (length(values) == 1L) {
      values <- rep(values, nrow(subject_data))
    }
    if (length(values) != nrow(subject_data)) {
      stop("Subject covariate '", nm, "' must return one value per subject.", call. = FALSE)
    }
    out[[nm]] <- values[match(data$.sim_subject_index, subject_data$.sim_subject_index)]
  }

  for (nm in names(observation)) {
    values <- .sim_eval_covariate_spec(observation[[nm]], data, nm)
    if (length(values) == 1L) {
      values <- rep(values, nrow(data))
    }
    if (length(values) != nrow(data)) {
      stop("Observation covariate '", nm, "' must return one value per row.", call. = FALSE)
    }
    out[[nm]] <- values
  }

  out
}

.sim_eval_covariate_spec <- function(spec, data, label) {
  if (is.function(spec)) {
    return(spec(data))
  }
  spec
}

.sim_add_covariates <- function(long, covariates, n, n_time, subject_var) {
  if (is.null(covariates)) {
    return(long)
  }
  if (is.function(covariates)) {
    covariates <- covariates(long)
  }
  covariates <- as.data.frame(covariates, stringsAsFactors = FALSE)
  if (nrow(covariates) == n) {
    covariates <- covariates[long$.sim_subject_index, , drop = FALSE]
  } else if (nrow(covariates) != n * n_time) {
    stop(
      "covariates must have either n rows or n * length(times) rows.",
      call. = FALSE
    )
  }
  overlap <- intersect(names(long), names(covariates))
  if (length(overlap) > 0L) {
    covariates <- covariates[, setdiff(names(covariates), overlap), drop = FALSE]
  }
  if (ncol(covariates) == 0L) {
    return(long)
  }
  cbind(long, covariates)
}
