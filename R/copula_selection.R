#' Select a bivariate copula family from pseudo-observation pairs
#'
#' Screens supported copula families by maximising their native copula
#' log-likelihood on adjacent-time pseudo-observation pairs. Pseudo-observations
#' can be supplied directly as `u1`/`u2`, supplied as a row-aligned uniform
#' vector or column in `data`, or computed from a fitted
#' `gamlss.longitudinal` object.
#'
#' This helper is intended as a lightweight family-screening step. By default
#' it estimates a constant dependence parameter for each candidate family. With
#' `copula_time_intercepts = TRUE`, it estimates separate dependence parameters
#' for each observed adjacent time pair, treating time as a factor rather than
#' a linear trend. Richer covariate or smooth dependence structures should be
#' fitted afterwards with [gamlss_longitudinal()].
#'
#' @param data Optional long-format data frame.
#' @param object Optional fitted `gamlss.longitudinal` object.
#' @param u1,u2 Optional vectors of paired pseudo-observations.
#' @param u Optional row-aligned vector of pseudo-observations for `data`.
#' @param u_var Optional name of a pseudo-observation column in `data`.
#' @param response_var Optional response column name in `data`. Used to create
#'   pseudo-observations automatically when `u1`/`u2`, `u`, `u_var`, and
#'   `object` are not supplied.
#' @param margin_dist Optional `gamlss.dist` family object, or a
#'   [select_margin()] result, used to fit the temporary marginal model for
#'   automatic pseudo-observations. If omitted, [select_margin()] is called and
#'   a warning is issued.
#' @param mu.formula,sigma.formula,nu.formula,tau.formula Optional temporary
#'   marginal model formulas used when `select_copula()` creates
#'   pseudo-observations from `response_var`. If omitted, intercept-only
#'   formulas are used unless `time_intercepts = TRUE` or `margin_dist` is a
#'   time-intercept [select_margin()] result.
#' @param subject_var,time_var Subject and time column names used when building
#'   adjacent-time pairs from `data`.
#' @param families Candidate copula family codes. Supported values are `"N"`,
#'   `"C"`, `"F"`, `"G"`, `"J"`, and `"t"`.
#' @param lags Positive integer lag(s) used when forming adjacent pairs.
#' @param criterion Ranking criterion, one of `"AIC"`, `"BIC"`, or `"logLik"`.
#' @param t_df_grid Degrees-of-freedom grid used for the t-copula screen.
#' @param min_pairs Minimum number of complete pairs required.
#' @param time_intercepts Logical; when `select_copula()` creates
#'   pseudo-observations from `response_var`, use time-specific intercepts in
#'   the temporary marginal model and pass this setting through to
#'   [select_margin()] if the marginal family is auto-selected.
#' @param copula_time_intercepts Logical; if `TRUE`, screen each copula family
#'   with separate dependence parameters for each adjacent time-pair factor
#'   level. This does not impose a linear time trend.
#' @param tail_thresholds Numeric vector of lower-tail probabilities used for
#'   empirical-versus-fitted tail co-occurrence and conditional tail exceedance
#'   summaries.
#'
#' @return A data frame with one row per family and class
#'   `copula_selection`. The selected family is stored in the `selected`
#'   attribute. The main table includes fitted lower/upper tail co-occurrence
#'   and conditional exceedance probabilities at the smallest requested
#'   `tail_thresholds` value. Full empirical-versus-fitted tail diagnostic
#'   tables are stored in the `tail_cooccurrence` and
#'   `conditional_tail_exceedance` attributes.
#' @export
select_copula <- function(
    data = NULL,
    object = NULL,
    u1 = NULL,
    u2 = NULL,
    u = NULL,
    u_var = NULL,
    response_var = NULL,
    margin_dist = NULL,
    mu.formula = NULL,
    sigma.formula = NULL,
    nu.formula = NULL,
    tau.formula = NULL,
    subject_var = "subject",
    time_var = "time",
    families = c("N", "C", "F", "G", "J", "t"),
    lags = 1,
    criterion = c("AIC", "BIC", "logLik"),
    t_df_grid = c(3, 4, 6, 8, 12, 20, 30),
    min_pairs = 10,
    time_intercepts = FALSE,
    copula_time_intercepts = FALSE,
    tail_thresholds = c(0.05, 0.10, 0.20)) {
  criterion <- match.arg(criterion)
  families <- vapply(families, .copula_family_code, character(1), USE.NAMES = FALSE)
  lags <- as.integer(lags)
  if (length(lags) < 1L || any(!is.finite(lags)) || any(lags < 1L)) {
    stop("lags must contain positive integers.", call. = FALSE)
  }
  tail_thresholds <- .select_copula_tail_thresholds(tail_thresholds)

  pairs <- .select_copula_pairs(
    data = data,
    object = object,
    u1 = u1,
    u2 = u2,
    u = u,
    u_var = u_var,
    response_var = response_var,
    margin_dist = margin_dist,
    mu.formula = mu.formula,
    sigma.formula = sigma.formula,
    nu.formula = nu.formula,
    tau.formula = tau.formula,
    subject_var = subject_var,
    time_var = time_var,
    lags = lags,
    time_intercepts = time_intercepts
  )
  keep <- is.finite(pairs$u1) & is.finite(pairs$u2)
  pairs <- pairs[keep, , drop = FALSE]
  if (nrow(pairs) < min_pairs) {
    stop("At least ", min_pairs, " complete pseudo-observation pairs are required.", call. = FALSE)
  }
  if (isTRUE(copula_time_intercepts) && !"copula_time" %in% names(pairs)) {
    stop(
      "'copula_time_intercepts = TRUE' requires data, u/u_var, response_var, or object input with time information.",
      call. = FALSE
    )
  }

  fits <- lapply(families, function(family) {
    .select_copula_fit_family(
      u1 = pairs$u1,
      u2 = pairs$u2,
      family = family,
      t_df_grid = t_df_grid,
      copula_time = if (isTRUE(copula_time_intercepts)) pairs$copula_time else NULL
    )
  })
  tail_diagnostics <- .select_copula_tail_diagnostics(
    pairs = pairs,
    fits = fits,
    thresholds = tail_thresholds,
    copula_time_intercepts = copula_time_intercepts
  )
  out <- do.call(rbind, fits)
  out <- .select_copula_add_tail_summary(out, tail_diagnostics)
  out$n_pairs <- nrow(pairs)
  out <- out[order(out[[criterion]], decreasing = identical(criterion, "logLik")), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "selected") <- out$family[1]
  attr(out, "criterion") <- criterion
  attr(out, "margin_selection") <- attr(pairs, "margin_selection")
  attr(out, "pseudo_observation_source") <- attr(pairs, "pseudo_observation_source")
  attr(out, "copula_time_intercepts") <- isTRUE(copula_time_intercepts)
  attr(out, "copula_time_levels") <- if (isTRUE(copula_time_intercepts)) unique(pairs$copula_time) else NULL
  attr(out, "tail_thresholds") <- tail_thresholds
  attr(out, "tail_cooccurrence") <- tail_diagnostics$cooccurrence
  attr(out, "conditional_tail_exceedance") <- tail_diagnostics$conditional
  class(out) <- c("copula_selection", "data.frame")
  out
}
