#' Check a fitted longitudinal GAMLSS-copula model
#'
#' `check_model()` returns a compact descriptive diagnostic summary. Apart from
#' the fitted object's explicit convergence contract, reported calibration and
#' dependence statistics are not converted into package-defined pass/fail
#' decisions.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param include_vcov Logical; include variance-covariance inference metadata
#'   via [summary.gamlss.longitudinal()].
#' @param include_plots Logical; include standard diagnostic plot objects in
#'   `check$plots`.
#' @param dependence_cor_cutoff Optional user-supplied absolute lag-1
#'   Rosenblatt normal-score residual-correlation threshold. The default `NULL`
#'   reports the value descriptively without flagging it.
#' @param pit_seed Seed for randomized PIT values for discrete margins. The
#'   caller's random-number state is preserved.
#' @param residual_lags Positive integer lags to summarize. Lags with no usable
#'   pairs are retained and marked unavailable.
#' @param ... Passed to [summary.gamlss.longitudinal()] when `include_vcov` is
#'   `TRUE`.
#'
#' @return An object of class `gamlss_longitudinal_check`, including descriptive
#'   `basic_checks` and `checks` tables, any user-threshold `flags`, PIT method
#'   provenance, and residual-dependence scope by lag.
#' @export
check_model <- function(
    object,
    include_vcov = FALSE,
    include_plots = FALSE,
    dependence_cor_cutoff = NULL,
    pit_seed = 1L,
    residual_lags = 1:3,
    ...) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }

  dependence_cor_cutoff <- .gl_validate_dependence_cor_cutoff(dependence_cor_cutoff)
  residual_lags <- unique(as.integer(residual_lags))
  if (!length(residual_lags) || anyNA(residual_lags) || any(residual_lags < 1L)) {
    stop("'residual_lags' must contain positive integers.", call. = FALSE)
  }

  s <- summary(object, include_vcov = include_vcov, ...)
  discrete_pit <- identical(
    .gl_capability_likelihood_route(object$margin_dist),
    "exact_discrete_rectangle"
  )
  pit_out <- .gl_pit(object, randomize = discrete_pit, seed = pit_seed)
  pit <- pmin(pmax(pit_out$pit, 0), 1)
  pit_stats <- .gl_pit_calibration_stats(pit)
  tail_calibration <- .gl_tail_calibration_stats(pit)
  tail_summary <- tail_calibration$tail_summary
  tail_stats <- tail_calibration$tail_stats
  scores <- as.data.frame(
    as.list(proscore(object, type = c("logs", "mae", "mse", "dss"))),
    stringsAsFactors = FALSE
  )

  residual_dependence <- .gl_residual_dependence_summary(
    object,
    residual_lags = residual_lags
  )
  residual_dependence$scope_status <- ifelse(
    residual_dependence$n_pairs > 0L & is.finite(residual_dependence$normal_score_cor),
    "available",
    "unavailable_no_usable_pairs"
  )
  residual_dependence$cutoff <- if (is.null(dependence_cor_cutoff)) {
    NA_real_
  } else {
    dependence_cor_cutoff
  }
  residual_dependence$threshold_source <- if (is.null(dependence_cor_cutoff)) {
    "none_descriptive"
  } else {
    "user_supplied"
  }

  lag1_cor <- residual_dependence$normal_score_cor[match(1L, residual_dependence$lag)]
  if (length(lag1_cor) == 0L) lag1_cor <- NA_real_

  copula_summary <- tryCatch(copula_time_summary(object), error = function(e) NULL)
  vcov_method <- s$vcov$method %||% s$fit$vcov_method %||% NA_character_
  checks <- .gl_check_table(
    summary_obj = list(fit = s$fit, convergence = object$convergence),
    scores = scores,
    pit_stats = pit_stats,
    tail_stats = tail_stats,
    lag1_cor = lag1_cor,
    dependence_cor_cutoff = dependence_cor_cutoff,
    vcov_method = vcov_method
  )
  flags <- checks[
    checks$status %in% c("not_converged", "flagged", "review", "unavailable"),
    ,
    drop = FALSE
  ]

  out <- list(
    model = s$model,
    fit = s$fit,
    convergence = object$convergence,
    scores = scores,
    pit = pit_stats,
    pit_method = list(
      randomized = pit_out$randomized,
      seed = pit_out$seed,
      rng_state_preserved = TRUE
    ),
    tail = tail_summary,
    residual_dependence = residual_dependence,
    copula = copula_summary,
    basic_checks = .gl_basic_checks(checks),
    diagnostic_summary = .gl_basic_checks_result(checks),
    checks = checks,
    flags = flags,
    plots = if (isTRUE(include_plots)) {
      list(
        pithist = pithist(object, randomize = discrete_pit, seed = pit_seed, plot = TRUE),
        qqrplot = qqrplot(object, randomize = discrete_pit, seed = pit_seed, plot = TRUE),
        wormplot = wormplot(object, randomize = discrete_pit, seed = pit_seed, plot = TRUE),
        rootogram = rootogram(object, plot = TRUE)
      )
    } else {
      NULL
    }
  )
  class(out) <- "gamlss_longitudinal_check"
  out
}
