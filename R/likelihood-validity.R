#' Combine fixed likelihood contributions and record numerical failures
#'
#' Inclusion is determined only by the response missingness masks. Numerical
#' validity is recorded separately so that a failed contribution can never be
#' mistaken for missing data.
#'
#' @noRd
.gl_likelihood_combine_components <- function(
    margin_log_d,
    margin_included,
    copula_log_d,
    pair_included,
    pair_input_valid = rep(TRUE, length(pair_included))) {
  margin_included <- as.logical(margin_included)
  pair_included <- as.logical(pair_included)

  margin_contribution_valid <- is.finite(margin_log_d)
  pair_input_valid <- as.logical(pair_input_valid)
  pair_contribution_valid <- pair_input_valid & is.finite(copula_log_d)

  bad_margin <- which(margin_included & !margin_contribution_valid)
  bad_pair_input <- which(pair_included & !pair_input_valid)
  bad_pair_contribution <- which(pair_included & pair_input_valid & !is.finite(copula_log_d))

  marginal_loglik <- if (length(bad_margin) == 0L) {
    sum(margin_log_d[margin_included])
  } else {
    -Inf
  }
  copula_loglik <- if (length(bad_pair_input) == 0L && length(bad_pair_contribution) == 0L) {
    sum(copula_log_d[pair_included])
  } else {
    -Inf
  }
  joint_loglik <- if (is.finite(marginal_loglik) && is.finite(copula_loglik)) {
    marginal_loglik + copula_loglik
  } else {
    -Inf
  }

  codes <- character(0)
  if (length(bad_margin) > 0L) codes <- c(codes, "invalid_margin_contribution")
  if (length(bad_pair_input) > 0L) codes <- c(codes, "invalid_pair_input")
  if (length(bad_pair_contribution) > 0L) codes <- c(codes, "invalid_pair_contribution")

  valid <- length(codes) == 0L
  failure <- if (valid) {
    NULL
  } else {
    list(
      codes = codes,
      margin_rows = as.integer(bad_margin),
      pair_rows = as.integer(sort(unique(c(bad_pair_input, bad_pair_contribution)))),
      pair_input_rows = as.integer(bad_pair_input),
      pair_contribution_rows = as.integer(bad_pair_contribution),
      message = paste(codes, collapse = ", ")
    )
  }

  list(
    log_lik = c(
      marginal = marginal_loglik,
      copula = copula_loglik,
      joint = joint_loglik
    ),
    valid = valid,
    failure = failure,
    contribution_counts = c(
      marginal_included = sum(margin_included),
      pair_included = sum(pair_included),
      marginal_valid = sum(margin_included & margin_contribution_valid),
      pair_valid = sum(pair_included & pair_contribution_valid)
    ),
    margin_included = margin_included,
    pair_included = pair_included,
    margin_contribution_valid = margin_contribution_valid,
    pair_input_valid = pair_input_valid,
    pair_contribution_valid = pair_contribution_valid
  )
}
