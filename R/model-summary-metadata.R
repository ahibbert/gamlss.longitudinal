#' Build high-level summary metadata for a fitted model
#'
#' @noRd
.gl_summary_model_info <- function(object) {
  n_obs <- length(object$response)
  n_subjects <- length(unique(object$response_subject))
  n_timepoints <- length(unique(object$response_margin))

  n_fixed <- length(object$par)
  n_smooth_terms <- 0
  if (!is.null(object$par_s)) {
    n_smooth_terms <- sum(vapply(object$par_s, length, integer(1)))
  }

  edf_smooth <- NA_real_
  if (!is.null(object$df_s) && length(object$df_s) > 0) {
    df_vals <- suppressWarnings(as.numeric(unlist(object$df_s, use.names = FALSE)))
    df_vals <- df_vals[is.finite(df_vals)]
    if (length(df_vals) > 0) {
      edf_smooth <- sum(df_vals)
    }
  }

  list(
    margin_dist = if (!is.null(object$margin_dist$family[1])) as.character(object$margin_dist$family[1]) else NA_character_,
    copula_dist = object$copula_dist,
    n_obs = n_obs,
    n_subjects = n_subjects,
    n_timepoints = n_timepoints,
    n_fixed = n_fixed,
    n_smooth_terms = n_smooth_terms,
    edf_smooth = edf_smooth
  )
}
