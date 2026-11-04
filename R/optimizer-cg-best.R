#' Update the best CG raw log-likelihood tracker
#'
#' @noRd
.gl_update_cg_best_loglik <- function(
    candidate_loglik,
    best_raw_loglik,
    best_iteration,
    current_iteration) {
  if (is.finite(candidate_loglik) && candidate_loglik > best_raw_loglik) {
    best_raw_loglik <- candidate_loglik
    best_iteration <- current_iteration
  }

  list(
    best_raw_loglik = best_raw_loglik,
    best_iteration = best_iteration
  )
}

#' Compute CG raw log-likelihood drop from best seen value
#'
#' @noRd
.gl_cg_raw_loglik_drop <- function(
    best_raw_loglik,
    current_loglik,
    prevented_deterioration = FALSE,
    prevented_raw_loglik_drop = NA_real_) {
  raw_loglik_drop <- best_raw_loglik - current_loglik

  if (isTRUE(prevented_deterioration) && is.finite(prevented_raw_loglik_drop)) {
    raw_loglik_drop <- max(raw_loglik_drop, prevented_raw_loglik_drop, na.rm = TRUE)
  }

  raw_loglik_drop
}
