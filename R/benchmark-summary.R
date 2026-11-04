# Consolidated support module for reviewer navigation.
# Source blocks are grouped mechanically from the files named below.

# ---- benchmark-summary.R ----

#' Summarise benchmark simulation results into win/tie/loss tables

#'

#' `summarise_benchmark_results()` turns rows from [run_coverage_simulations()]

#' or `benchmark_standard_models()`-style outputs into compact evidence tables.

#' Each metric is scored according to its estimand: lower error is better,

#' coverage targets 0.95, calibration/tail errors target 0, and runtime is

#' lower-is-better.

#'

#' @param results Data frame returned by [run_coverage_simulations()].

#' @param metrics Optional metric columns to summarise. Defaults to all known

#'   benchmark metrics present in `results`.

#' @param group_cols Columns defining one simulation case before method-level

#'   comparison.

#' @param tie_tolerance Relative gap from the case-best score treated as a tie.

#' @param absolute_tolerance Absolute score gap treated as a tie/win.

#'

#' @return An object of class `gamlss_longitudinal_benchmark_summary` with

#'   `summary`, `case_results`, and `metric_catalog` components.

#' @export

summarise_benchmark_results <- function(
    results,
    metrics = NULL,
    group_cols = c("family", "copula", "design", "n_subject", "n_time", "dependence", "missingness", "start_mode"),
    tie_tolerance = 0.05,
    absolute_tolerance = 1e-8) {
  results <- as.data.frame(results, stringsAsFactors = FALSE)

  if (!"method" %in% names(results)) {
    stop("'results' must contain a 'method' column.", call. = FALSE)
  }

  catalog <- .benchmark_metric_catalog()

  if (is.null(metrics)) {
    metrics <- intersect(catalog$metric, names(results))
  } else {
    metrics <- intersect(as.character(metrics), names(results))
  }

  metrics <- metrics[metrics %in% catalog$metric]

  if (length(metrics) == 0L) {
    stop("No supported benchmark metric columns were found in 'results'.", call. = FALSE)
  }

  group_cols <- intersect(group_cols, names(results))

  if (length(group_cols) == 0L) {
    results$.benchmark_case <- 1L

    group_cols <- ".benchmark_case"
  }

  group_key <- interaction(results[group_cols], drop = TRUE, lex.order = TRUE)

  case_rows <- list()

  for (metric in metrics) {
    meta <- catalog[match(metric, catalog$metric), , drop = FALSE]

    value <- suppressWarnings(as.numeric(results[[metric]]))

    score <- .benchmark_metric_score(value, meta$score_rule, meta$target)

    for (key in unique(group_key)) {
      idx <- which(group_key == key)

      finite_idx <- idx[is.finite(score[idx])]

      if (length(finite_idx) == 0L) next

      best <- min(score[finite_idx], na.rm = TRUE)

      labels <- vapply(score[idx], .benchmark_case_label, character(1),
        best_score = best,
        tie_tolerance = tie_tolerance,
        absolute_tolerance = absolute_tolerance
      )

      case_rows[[length(case_rows) + 1L]] <- data.frame(

        results[idx, group_cols, drop = FALSE],
        method = results$method[idx],
        metric = metric,
        estimand = meta$estimand,
        value = value[idx],
        score = score[idx],
        best_score = best,
        result = labels,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }

  case_results <- if (length(case_rows) == 0L) {
    data.frame()
  } else {
    do.call(rbind, case_rows)
  }

  if (nrow(case_results) == 0L) {
    summary <- data.frame()
  } else {
    split_key <- interaction(case_results$metric, case_results$method, drop = TRUE, lex.order = TRUE)

    summary_rows <- lapply(split(seq_len(nrow(case_results)), split_key), function(idx) {
      x <- case_results[idx, , drop = FALSE]

      finite <- is.finite(x$value)

      data.frame(
        metric = x$metric[[1L]],
        estimand = x$estimand[[1L]],
        method = x$method[[1L]],
        n = nrow(x),
        n_finite = sum(finite),
        wins = sum(x$result == "win", na.rm = TRUE),
        ties = sum(x$result == "tie", na.rm = TRUE),
        losses = sum(x$result == "loss", na.rm = TRUE),
        missing = sum(x$result == "missing", na.rm = TRUE),
        win_rate = mean(x$result == "win", na.rm = TRUE),
        win_or_tie_rate = mean(x$result %in% c("win", "tie"), na.rm = TRUE),
        median_value = if (any(finite)) stats::median(x$value[finite], na.rm = TRUE) else NA_real_,
        median_score = if (any(is.finite(x$score))) stats::median(x$score[is.finite(x$score)], na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    })

    summary <- do.call(rbind, summary_rows)

    summary <- summary[order(summary$estimand, summary$metric, -summary$win_or_tie_rate, summary$median_score), , drop = FALSE]

    rownames(summary) <- NULL
  }

  out <- list(
    summary = summary,
    case_results = case_results,
    metric_catalog = catalog[catalog$metric %in% metrics, , drop = FALSE]
  )

  class(out) <- "gamlss_longitudinal_benchmark_summary"

  out
}

# ---- benchmark-summary-metrics.R ----

.benchmark_metric_catalog <- function() {
  data.frame(
    metric = c(
      "benchmark_mean_rmse",
      "benchmark_mean_mae",
      "benchmark_mean_bias",
      "benchmark_rmse",
      "benchmark_mae",
      "benchmark_q90_mae",
      "benchmark_neg_log_score",
      "benchmark_upper_tail_error_90",
      "benchmark_theta_time_abs_error",
      "benchmark_interval_coverage_95",
      "benchmark_interval_width_95",
      "benchmark_pit_mean_abs_error",
      "benchmark_pit_ks_p_value",
      "benchmark_tail_error_lower_05",
      "benchmark_tail_error_upper_05",
      "smooth_eta_rmse",
      "smooth_eta_max_abs_error",
      "elapsed_sec"
    ),
    estimand = c(
      "mean",
      "mean",
      "mean",
      "observed_response",
      "observed_response",
      "quantile",
      "density",
      "tail",
      "dependence",
      "interval",
      "interval_width",
      "calibration",
      "calibration",
      "tail",
      "tail",
      "smooth",
      "smooth",
      "runtime"
    ),
    score_rule = c(
      "lower",
      "lower",
      "absolute",
      "lower",
      "lower",
      "lower",
      "lower",
      "absolute",
      "lower",
      "target",
      "lower",
      "lower",
      "higher",
      "absolute",
      "absolute",
      "lower",
      "lower",
      "lower"
    ),
    target = c(NA, NA, 0, NA, NA, NA, NA, 0, NA, 0.95, NA, NA, 1, 0, 0, NA, NA, NA),
    stringsAsFactors = FALSE
  )
}

.benchmark_metric_score <- function(value, rule, target) {
  if (identical(rule, "lower")) {
    return(value)
  }

  if (identical(rule, "higher")) {
    return(1 - value)
  }

  if (identical(rule, "absolute")) {
    return(abs(value))
  }

  abs(value - target)
}

.benchmark_case_label <- function(score, best_score, tie_tolerance, absolute_tolerance) {
  if (!is.finite(score) || !is.finite(best_score)) {
    return("missing")
  }

  gap <- score - best_score

  if (gap <= absolute_tolerance) {
    return("win")
  }

  scale <- max(abs(best_score), absolute_tolerance)

  if ((gap / scale) <= tie_tolerance) {
    return("tie")
  }

  "loss"
}

# ---- benchmark-summary-print.R ----

#' @export

print.gamlss_longitudinal_benchmark_summary <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nBenchmark Win/Tie/Loss Summary\n")

  cat("------------------------------\n")

  if (nrow(x$summary) == 0L) {
    cat("No finite benchmark metrics to summarise.\n")
  } else {
    print(.benchmark_print_metric_labels(x$summary), digits = digits, row.names = FALSE)
  }

  invisible(x)
}
