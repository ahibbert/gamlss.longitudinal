#' Check a fitted longitudinal GAMLSS-copula model

#'

#' `check_model()` turns diagnostics into a compact set of basic automated

#' checks for broad applied use. It does not replace visual inspection, but it

#' provides a stable first pass over convergence, marginal calibration, residual

#' dependence, and scoring summaries.

#'

#' @param object A fitted `gamlss.longitudinal` object.

#' @param include_vcov Logical; include variance-covariance inference metadata

#'   via [summary.gamlss.longitudinal()].

#' @param include_plots Logical; include standard diagnostic plot objects in

#'   `check$plots`. Visual review should usually use `plot(object)` or the

#'   explicit diagnostic helpers instead.

#' @param dependence_cor_cutoff Absolute lag-1 Rosenblatt normal-score residual

#'   correlation above which the dependence check is flagged. The default is a

#'   review threshold rather than a formal hypothesis test.

#' @param ... Passed to [summary.gamlss.longitudinal()] when `include_vcov` is

#'   `TRUE`.

#'

#' @return An object of class `gamlss_longitudinal_check`, including a compact

#'   `basic_checks` table, a full `checks` table, a `warnings` table containing

#'   failed checks, and overall `basic_checks_passed` and `basic_checks_result`

#'   fields.

#' @export

check_model <- function(
    object,
    include_vcov = FALSE,
    include_plots = FALSE,
    dependence_cor_cutoff = 0.25,
    ...) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }

  dependence_cor_cutoff <- .gl_validate_dependence_cor_cutoff(dependence_cor_cutoff)

  s <- summary(object, include_vcov = include_vcov, ...)

  pit_out <- .gl_pit(object, randomize = FALSE)

  pit <- pmin(pmax(pit_out$pit, 0), 1)

  pit_stats <- .gl_pit_calibration_stats(pit)

  tail_calibration <- .gl_tail_calibration_stats(pit)

  tail_summary <- tail_calibration$tail_summary

  tail_stats <- tail_calibration$tail_stats

  scores <- as.data.frame(as.list(proscore(object, type = c("logs", "mae", "mse", "dss"))), stringsAsFactors = FALSE)

  residual_dependence <- .gl_residual_dependence_summary(object, residual_lags = 1L)

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

  warnings <- checks[checks$status == "FAIL", , drop = FALSE]

  basic_checks_result <- .gl_basic_checks_result(checks)

  if (nrow(warnings) > 0L) {
    warning(
      paste0(
        "Basic model checks failed: ",
        paste(warnings$area, collapse = ", "),
        ". Review check$warnings and broader diagnostics."
      ),
      call. = FALSE
    )
  }

  out <- list(
    model = s$model,
    fit = s$fit,
    convergence = object$convergence,
    scores = scores,
    pit = pit_stats,
    tail = tail_summary,
    residual_dependence = transform(residual_dependence, cutoff = dependence_cor_cutoff),
    copula = copula_summary,
    basic_checks = .gl_basic_checks(checks),
    basic_checks_passed = !any(checks$status == "FAIL", na.rm = TRUE),
    basic_checks_result = basic_checks_result,
    checks = checks,
    warnings = warnings,
    plots = if (isTRUE(include_plots)) {
      list(
        pithist = pithist(object, plot = TRUE),
        qqrplot = qqrplot(object, plot = TRUE),
        wormplot = wormplot(object, plot = TRUE),
        rootogram = rootogram(object, plot = TRUE)
      )
    } else {
      NULL
    }
  )

  class(out) <- "gamlss_longitudinal_check"

  out
}
