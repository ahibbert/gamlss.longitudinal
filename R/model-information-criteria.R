#' Count independent sampling units for longitudinal BIC
#'
#' BIC for a longitudinal likelihood uses subjects, not repeated responses, as
#' the independent sampling units. Only subjects contributing at least one
#' finite response are counted.
#'
#' @param response Response vector.
#' @param subject Subject identifier aligned with `response`.
#' @return Integer number of observed subjects.
#' @noRd
.gl_bic_subject_count <- function(response, subject) {
  if (is.null(response) || is.null(subject) || length(subject) != length(response)) {
    stop(
      "Subject-level BIC requires a subject identifier for every response.",
      call. = FALSE
    )
  }

  observed <- is.finite(response)
  if (!any(observed)) {
    stop("BIC requires at least one finite observed response.", call. = FALSE)
  }
  observed_subject <- subject[observed]
  if (anyNA(observed_subject)) {
    stop("Observed responses must have non-missing subject identifiers for BIC.", call. = FALSE)
  }

  as.integer(length(unique(observed_subject)))
}
