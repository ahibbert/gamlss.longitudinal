#' Select margin and copula families by joint longitudinal fit
#'
#' `select_joint_distribution()` screens marginal distribution and copula
#' combinations by fitting intercept-only [gamlss_longitudinal()] models and
#' ranking their joint likelihood, AIC, or BIC. Unlike [select_margin()] and
#' [select_copula()], this selector evaluates each candidate as a full joint
#' longitudinal model, so it can take substantially longer on large candidate
#' sets.
#'
#' @param data Long-format data frame.
#' @param response_var Response column name in `data`.
#' @param time_var,subject_var Time and subject identifier column names.
#' @param type Optional `gamlss::fitDist()` type passed to [select_margin()].
#' @param margin_families Optional marginal family names. When supplied, these
#'   are used directly as the marginal candidate set; otherwise candidates come
#'   from [select_margin()].
#' @param copula_families Candidate copula family codes. Supported values are
#'   `"N"`, `"C"`, `"F"`, `"G"`, `"J"`, and `"t"`.
#' @param criterion Ranking criterion, one of `"AIC"`, `"BIC"`, or `"logLik"`.
#' @param min_pairs Minimum number of complete adjacent response pairs required
#'   before fitting candidates.
#' @param time_intercepts Logical; if `TRUE`, pass through to
#'   [select_margin()] and use time-specific intercepts for each marginal
#'   distribution parameter in the joint screening fits.
#' @param copula_time_intercepts Logical; if `TRUE`, use time-specific
#'   intercepts for the copula dependence parameter in the joint screening fits.
#'   Time is treated as a factor, not as a linear trend.
#' @param try.gamlss Passed to [select_margin()].
#' @param trace Logical; passed to [select_margin()] and used to decide whether
#'   candidate [gamlss_longitudinal()] fit output is shown.
#' @param progress Logical; if `TRUE`, print candidate-level progress messages
#'   with elapsed time and an estimated remaining runtime.
#' @param fit_args Optional named list of arguments overriding the default
#'   [gamlss_longitudinal()] screening fit settings.
#' @param keep_fits Logical; if `TRUE`, attach retained fitted models in the
#'   `"fits"` attribute.
#'
#' @return A data frame with one row per margin-copula combination and class
#'   `joint_distribution_selection`. Failed fits are retained with an `error`
#'   message and missing fit metrics. Nonconverged fits retain their provisional
#'   metrics and stop reason for audit, but have `selection_eligible = FALSE`,
#'   receive no rank, and can never be returned by [best_fit()].
#' @export
select_joint_distribution <- function(
    data,
    response_var,
    time_var = "time",
    subject_var = "subject",
    type = NULL,
    margin_families = NULL,
    copula_families = c("N", "C", "F", "G", "J", "t"),
    criterion = c("AIC", "BIC", "logLik"),
    min_pairs = 10,
    time_intercepts = FALSE,
    copula_time_intercepts = FALSE,
    try.gamlss = FALSE,
    trace = FALSE,
    progress = TRUE,
    fit_args = list(),
    keep_fits = FALSE) {
  criterion <- match.arg(criterion)
  if (!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  } else {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  .joint_selection_check_column(data, response_var, "response_var")
  .joint_selection_check_column(data, time_var, "time_var")
  .joint_selection_check_column(data, subject_var, "subject_var")
  if (!is.list(fit_args) ||
    (length(fit_args) > 0L && (is.null(names(fit_args)) || any(!nzchar(names(fit_args)))))) {
    stop("'fit_args' must be a named list.", call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("'progress' must be TRUE or FALSE.", call. = FALSE)
  }
  min_pairs <- as.integer(min_pairs)
  if (length(min_pairs) != 1L || !is.finite(min_pairs) || min_pairs < 1L) {
    stop("'min_pairs' must be a positive integer.", call. = FALSE)
  }
  n_pairs <- .joint_selection_count_pairs(
    data = data,
    response_var = response_var,
    time_var = time_var,
    subject_var = subject_var
  )
  if (n_pairs < min_pairs) {
    stop("At least ", min_pairs, " complete adjacent response pairs are required.", call. = FALSE)
  }

  copula_families <- unique(vapply(
    copula_families,
    .copula_family_code,
    character(1),
    USE.NAMES = FALSE
  ))
  if (length(copula_families) < 1L) {
    stop("'copula_families' must contain at least one supported family code.", call. = FALSE)
  }

  margin_selection <- .joint_selection_margin_candidates(
    data = data,
    response_var = response_var,
    time_var = time_var,
    type = type,
    margin_families = margin_families,
    time_intercepts = time_intercepts,
    try.gamlss = try.gamlss,
    trace = trace
  )
  margin_candidates <- as.data.frame(margin_selection)
  margin_candidates <- margin_candidates[
    isTRUE(length(margin_candidates$supported_by_longitudinal) > 0L) &
      margin_candidates$supported_by_longitudinal %in% TRUE, ,
    drop = FALSE
  ]
  if (nrow(margin_candidates) < 1L) {
    stop("No supported marginal families were retained for joint screening.", call. = FALSE)
  }

  combinations <- expand.grid(
    margin_family = margin_candidates$family,
    copula_family = copula_families,
    stringsAsFactors = FALSE
  )
  route_supported <- mapply(
    .gl_capability_route_supported,
    combinations$margin_family,
    combinations$copula_family,
    USE.NAMES = FALSE
  )
  excluded_routes <- combinations[!route_supported, , drop = FALSE]
  combinations <- combinations[route_supported, , drop = FALSE]
  if (nrow(combinations) < 1L) {
    .gl_capability_stop(
      "gamlss_longitudinal_unsupported_route_error",
      "No requested margin-copula combinations are present in the tested capability registry."
    )
  }

  candidate_fits <- .joint_selection_fit_candidates(
    combinations = combinations,
    data = data,
    response_var = response_var,
    time_var = time_var,
    subject_var = subject_var,
    fit_args = fit_args,
    time_intercepts = time_intercepts,
    copula_time_intercepts = copula_time_intercepts,
    trace = trace,
    progress = progress,
    keep_fits = keep_fits
  )

  out <- .joint_selection_finalize_result(
    rows = candidate_fits$rows,
    fit_store = candidate_fits$fit_store,
    n_pairs = n_pairs,
    criterion = criterion,
    margin_selection = margin_selection,
    time_intercepts = time_intercepts,
    time_var = time_var,
    copula_time_intercepts = copula_time_intercepts,
    keep_fits = keep_fits
  )
  attr(out, "excluded_by_capability_registry") <- excluded_routes
  out
}
