.gl_simulation_time_grid <- function(nd, object) {
  time_for_grid <- if ("time_covariate" %in% names(nd) && is.factor(nd$time_covariate)) {
    nd$time_covariate
  } else {
    nd$time
  }

  time_levels <- .gl_simulation_time_levels(time_for_grid, reference_time = object$response_margin)

  time_lookup <- stats::setNames(seq_along(time_levels), as.character(time_levels))

  time_idx <- unname(time_lookup[as.character(time_for_grid)])

  if (any(!is.finite(time_idx))) {
    stop("Could not map newdata time values to an ordered simulation grid.", call. = FALSE)
  }

  if (anyDuplicated(paste(nd$subject, as.character(time_for_grid), sep = "\r"))) {
    stop("newdata must contain at most one row per subject/time combination.", call. = FALSE)
  }

  list(
    time_for_grid = time_for_grid,
    time_levels = time_levels,
    time_idx = time_idx
  )
}

.gl_simulation_align_dependence_parameter <- function(param_vec, common_n, nd_eval, time_levels) {
  if (length(param_vec) == 0L) {
    return(rep(NA_real_, common_n))
  }

  if (length(param_vec) == common_n) {
    return(as.numeric(param_vec))
  }

  left_rows <- which(nd_eval$.gl_sim_time_idx < length(time_levels))

  if (length(param_vec) == length(left_rows)) {
    out <- rep(NA_real_, common_n)

    out[left_rows] <- as.numeric(param_vec)

    return(out)
  }

  rep(as.numeric(param_vec), length.out = common_n)
}
