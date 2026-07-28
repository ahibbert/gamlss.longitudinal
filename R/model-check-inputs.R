#' Validate residual-dependence cutoff for model checks
#'
#' @param dependence_cor_cutoff Absolute correlation threshold.
#' @return Numeric scalar cutoff.
#' @noRd
.gl_validate_dependence_cor_cutoff <- function(dependence_cor_cutoff) {
  dependence_cor_cutoff <- as.numeric(dependence_cor_cutoff)

  if (length(dependence_cor_cutoff) != 1L || !is.finite(dependence_cor_cutoff) ||
    dependence_cor_cutoff <= 0 || dependence_cor_cutoff >= 1) {
    stop("'dependence_cor_cutoff' must be a single number between 0 and 1.", call. = FALSE)
  }

  dependence_cor_cutoff
}
