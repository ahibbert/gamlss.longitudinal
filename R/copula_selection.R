#' Select a bivariate copula family from pseudo-observation pairs
#'
#' Screens supported copula families on adjacent-time pairs. Continuous
#' margins use the copula-density log-likelihood. Automatically constructed
#' discrete margins use the exact copula rectangle likelihood, matching the
#' likelihood used by [gamlss_longitudinal()]. Pseudo-observations
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
#'   `"C"`, `"F"`, `"G"`, `"J"`, and `"t"`. All six are available for every
#'   supported homogeneous marginal family.
#' @param lags Positive integer lag(s) used when forming adjacent pairs.
#' @param criterion Ranking criterion, one of `"AIC"`, `"BIC"`, or `"logLik"`.
#' @param t_df_grid Degrees-of-freedom grid used for the t-copula screen.
#' @param min_pairs Minimum number of complete pairs required.
#' @param time_intercepts Logical; when `select_copula()` creates
#'   pseudo-observations from `response_var`, use time-specific intercepts in
#'   the temporary marginal model and pass this setting through to
#'   [select_margin()] if the marginal family is auto-selected.
#' @param margin_n_cycles Positive integer maximum number of cycles used for
#'   the temporary marginal GAMLSS fit that creates pseudo-observations.
#' @param margin_args Named fixed marginal arguments passed to the temporary
#'   marginal fit. Binomial families require `list(bd = ...)` as described in
#'   [gamlss_longitudinal()].
#' @param pit_seed Finite numeric seed used for randomized PIT values when
#'   constructing descriptive diagnostics for a discrete margin. The seed does
#'   not affect the exact rectangle likelihood used to rank copulas.
#' @param copula_time_intercepts Logical; if `TRUE`, screen each copula family
#'   with separate dependence parameters for each adjacent time-pair factor
#'   level. This does not impose a linear time trend.
#' @param tail_thresholds Numeric vector of lower-tail probabilities used for
#'   empirical-versus-fitted tail co-occurrence and conditional tail exceedance
#'   summaries.
#' @param verbose Logical or integer verbosity level. The default, `1`, reports
#'   pair construction and candidate-level screening progress. Use `0` or
#'   `FALSE` for no progress output.
#'
#' @return A data frame with one row per retained family and class
#'   `copula_selection`. The selected family is stored in the `selected`
#'   attribute. Families removed by the route registry are stored in the
#'   `excluded_by_capability_registry` attribute. The main table includes fitted lower/upper tail co-occurrence
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
    margin_n_cycles = 100L,
    pit_seed = 1L,
    copula_time_intercepts = FALSE,
    tail_thresholds = c(0.05, 0.10, 0.20),
    margin_args = list(),
    verbose = 1) {
  verbose <- .gl_normalize_verbose(verbose)
  criterion <- match.arg(criterion)
  families <- vapply(families, .copula_family_code, character(1), USE.NAMES = FALSE)
  lags <- as.integer(lags)
  if (length(lags) < 1L || any(!is.finite(lags)) || any(lags < 1L)) {
    stop("lags must contain positive integers.", call. = FALSE)
  }
  tail_thresholds <- .select_copula_tail_thresholds(tail_thresholds)
  if (length(margin_n_cycles) != 1L || !is.finite(margin_n_cycles) || margin_n_cycles < 1) {
    stop("'margin_n_cycles' must be a positive integer.", call. = FALSE)
  }
  margin_n_cycles <- as.integer(margin_n_cycles)
  if (length(pit_seed) != 1L || !is.finite(pit_seed)) {
    stop("'pit_seed' must be one finite number.", call. = FALSE)
  }
  pit_seed <- as.integer(pit_seed)

  .gl_progress(verbose, "Preparing pseudo-observation pairs.", prefix = "select_copula")

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
    time_intercepts = time_intercepts,
    margin_n_cycles = margin_n_cycles,
    pit_seed = pit_seed,
    margin_args = margin_args,
    verbose = verbose
  )
  margin_family <- .select_copula_margin_family(
    margin_dist = margin_dist,
    object = object,
    margin_selection = attr(pairs, "margin_selection")
  )
  excluded_families <- character()
  if (!is.na(margin_family)) {
    supported <- vapply(
      families,
      function(family) .gl_capability_route_supported(margin_family, family),
      logical(1)
    )
    excluded_families <- families[!supported]
    families <- families[supported]
    if (!length(families)) {
      .gl_capability_stop(
        "gamlss_longitudinal_unsupported_route_error",
        paste0(
          "None of the requested copulas is a registered route for margin family '",
          margin_family,
          "'."
        ),
        margin_family = margin_family,
        copulas = excluded_families,
        context = "select_copula"
      )
    }
  }
  screen_likelihood <- attr(pairs, "screen_likelihood") %||% "copula_density"
  randomized_pit <- isTRUE(attr(pairs, "randomized_pit"))
  used_pit_seed <- attr(pairs, "pit_seed")
  exact_discrete <- identical(screen_likelihood, "exact_discrete_rectangle")
  if (exact_discrete && !all(c("upper1", "upper2", "lower1", "lower2") %in% names(pairs))) {
    stop("Exact discrete copula selection requires upper and lower CDF bounds.", call. = FALSE)
  }
  keep <- is.finite(pairs$u1) & is.finite(pairs$u2)
  if (exact_discrete) {
    keep <- keep & is.finite(pairs$upper1) & is.finite(pairs$upper2) &
      is.finite(pairs$lower1) & is.finite(pairs$lower2) &
      pairs$upper1 >= pairs$lower1 & pairs$upper2 >= pairs$lower2
    if (all(c("log_mass1", "log_mass2") %in% names(pairs))) {
      keep <- keep & is.finite(pairs$log_mass1) & is.finite(pairs$log_mass2)
    }
  }
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

  .gl_progress(
    verbose,
    "Prepared ", nrow(pairs), " complete pair(s); screening ",
    length(families), " copula family/families.",
    prefix = "select_copula"
  )

  fits <- lapply(seq_along(families), function(index) {
    family <- families[[index]]
    .gl_progress(
      verbose,
      "[", index, "/", length(families), "] Fitting ",
      .gl_copula_display_label(family), ".",
      prefix = "select_copula"
    )
    .select_copula_fit_family(
      u1 = pairs$u1,
      u2 = pairs$u2,
      family = family,
      t_df_grid = t_df_grid,
      copula_time = if (isTRUE(copula_time_intercepts)) pairs$copula_time else NULL,
      upper1 = if (exact_discrete) pairs$upper1 else NULL,
      upper2 = if (exact_discrete) pairs$upper2 else NULL,
      lower1 = if (exact_discrete) pairs$lower1 else NULL,
      lower2 = if (exact_discrete) pairs$lower2 else NULL,
      log_mass1 = if (exact_discrete && "log_mass1" %in% names(pairs)) pairs$log_mass1 else NULL,
      log_mass2 = if (exact_discrete && "log_mass2" %in% names(pairs)) pairs$log_mass2 else NULL,
      log_survival1 = if (exact_discrete && "log_survival1" %in% names(pairs)) pairs$log_survival1 else NULL,
      log_survival2 = if (exact_discrete && "log_survival2" %in% names(pairs)) pairs$log_survival2 else NULL
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
  out$screen_likelihood <- screen_likelihood
  out <- out[order(out[[criterion]], decreasing = identical(criterion, "logLik")), , drop = FALSE]
  rownames(out) <- NULL
  out$rank <- seq_len(nrow(out))
  best_value <- if (nrow(out)) out[[criterion]][[1L]] else NA_real_
  out$delta <- if (identical(criterion, "logLik")) best_value - out[[criterion]] else out[[criterion]] - best_value
  attr(out, "selected") <- out$family[1]
  attr(out, "criterion") <- criterion
  attr(out, "margin_selection") <- attr(pairs, "margin_selection")
  attr(out, "pseudo_observation_source") <- attr(pairs, "pseudo_observation_source")
  attr(out, "copula_time_intercepts") <- isTRUE(copula_time_intercepts)
  attr(out, "copula_time_levels") <- if (isTRUE(copula_time_intercepts)) unique(pairs$copula_time) else NULL
  attr(out, "tail_thresholds") <- tail_thresholds
  attr(out, "margin_family") <- margin_family
  attr(out, "excluded_by_capability_registry") <- excluded_families
  attr(out, "screen_likelihood") <- screen_likelihood
  attr(out, "screen_version") <- "2026.1"
  attr(out, "randomized_pit") <- randomized_pit
  attr(out, "pit_seed") <- if (randomized_pit) as.integer(used_pit_seed) else NA_integer_
  attr(out, "tail_cooccurrence") <- tail_diagnostics$cooccurrence
  attr(out, "conditional_tail_exceedance") <- tail_diagnostics$conditional
  class(out) <- c("copula_selection", "data.frame")
  .gl_progress(
    verbose,
    "Complete: selected ", .gl_copula_display_label(out$family[[1L]]), ".",
    prefix = "select_copula"
  )
  out
}
