source(file.path("paper", "R", "08-simulation-sensitivity-correlation-misspecification", "standard-model-benchmarking", "00-benchmark-setup.R"))

run_dir <- bmk_env("GAMLSS_LONGITUDINAL_BENCHMARK_RUN_DIR", "")
if (!nzchar(run_dir)) {
  latest <- file.path(bmk_output_root, "latest_run_dir.txt")
  if (!file.exists(latest)) {
    stop("No run directory supplied and latest_run_dir.txt does not exist.", call. = FALSE)
  }
  run_dir <- readLines(latest, warn = FALSE)[1L]
}
run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)

results_path <- file.path(run_dir, "benchmark_results_by_rep.csv")
coef_path <- file.path(run_dir, "coefficient_results_by_rep.csv")
status_path <- file.path(run_dir, "primary_status_by_rep.csv")
complexity_path <- file.path(run_dir, "fit_complexity_by_rep.csv")
dependence_path <- file.path(run_dir, "dependence_recovery_by_rep.csv")

if (!file.exists(results_path)) stop("Missing ", results_path, call. = FALSE)
if (!file.exists(coef_path)) stop("Missing ", coef_path, call. = FALSE)

results <- utils::read.csv(results_path, stringsAsFactors = FALSE, check.names = FALSE)
coef <- utils::read.csv(coef_path, stringsAsFactors = FALSE, check.names = FALSE)
status <- if (file.exists(status_path)) utils::read.csv(status_path, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
complexity <- if (file.exists(complexity_path)) utils::read.csv(complexity_path, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
dependence <- if (file.exists(dependence_path)) utils::read.csv(dependence_path, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()

add_n_time <- function(x) {
  if (nrow(x) > 0L && !"n_time" %in% names(x)) x$n_time <- 4L
  x
}
results <- add_n_time(results)
coef <- add_n_time(coef)
status <- add_n_time(status)
complexity <- add_n_time(complexity)
dependence <- add_n_time(dependence)

metric_defs <- data.frame(
  metric = c(
    "benchmark_mean_rmse",
    "benchmark_q90_mae",
    "benchmark_upper_tail_error_90",
    "benchmark_neg_log_score",
    "benchmark_interval_coverage_95",
    "benchmark_interval_width_95",
    "elapsed_sec"
  ),
  score_rule = c("lower", "lower", "absolute", "lower", "target", "lower", "lower"),
  target = c(NA, NA, 0, NA, 0.95, NA, NA),
  stringsAsFactors = FALSE
)
metric_defs <- metric_defs[metric_defs$metric %in% names(results), , drop = FALSE]

score_metric <- function(values, rule, target) {
  values <- suppressWarnings(as.numeric(values))
  if (identical(rule, "absolute")) return(abs(values - target))
  if (identical(rule, "target")) return(abs(values - target))
  values
}

summarise_metric <- function(metric_row) {
  metric <- metric_row$metric
  rule <- metric_row$score_rule
  target <- metric_row$target
  group_cols <- c("generator", "scenario", "correlation", "correlation_level", "n_time", "family", "method")
  groups <- unique(results[group_cols])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    idx <- rep(TRUE, nrow(results))
    for (col in group_cols) idx <- idx & results[[col]] == groups[[col]][i]
    vals <- suppressWarnings(as.numeric(results[[metric]][idx]))
    scores <- score_metric(vals, rule, target)
    rows[[i]] <- cbind(
      groups[i, , drop = FALSE],
      data.frame(
        metric = metric,
        n = length(vals),
        n_finite = sum(is.finite(vals)),
        finite_rate = mean(is.finite(vals)),
        mean_value = mean(vals, na.rm = TRUE),
        median_value = stats::median(vals, na.rm = TRUE),
        mean_score = mean(scores, na.rm = TRUE),
        median_score = stats::median(scores, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    )
  }
  do.call(rbind, rows)
}

benchmark_summary <- if (nrow(metric_defs) > 0L) {
  do.call(rbind, lapply(seq_len(nrow(metric_defs)), function(i) summarise_metric(metric_defs[i, , drop = FALSE])))
} else {
  data.frame()
}

coef_numeric <- c("estimate", "std_error", "conf.low", "conf.high", "truth", "bias", "ci_width", "p_value")
for (nm in intersect(coef_numeric, names(coef))) {
  coef[[nm]] <- suppressWarnings(as.numeric(coef[[nm]]))
}
for (nm in c("ci_covers_truth", "false_positive")) {
  if (nm %in% names(coef)) {
    coef[[nm]] <- coef[[nm]] %in% c(TRUE, "TRUE", "True", "true", "1")
  }
}

coef_group_cols <- c("generator", "scenario", "correlation", "correlation_level", "n_time", "family", "method", "term")
coef_groups <- unique(coef[coef_group_cols])
coef_summary_rows <- vector("list", nrow(coef_groups))
for (i in seq_len(nrow(coef_groups))) {
  idx <- rep(TRUE, nrow(coef))
  for (col in coef_group_cols) idx <- idx & coef[[col]] == coef_groups[[col]][i]
  sub <- coef[idx, , drop = FALSE]
  est <- sub$estimate
  se <- sub$std_error
  truth <- sub$truth[is.finite(sub$truth)][1L]
  empirical_sd <- stats::sd(est, na.rm = TRUE)
  mean_reported_se <- mean(se, na.rm = TRUE)
  coef_summary_rows[[i]] <- cbind(
    coef_groups[i, , drop = FALSE],
    data.frame(
      truth = truth,
      n = nrow(sub),
      n_estimate = sum(is.finite(est)),
      mean_estimate = mean(est, na.rm = TRUE),
      bias = mean(est - sub$truth, na.rm = TRUE),
      rmse = sqrt(mean((est - sub$truth)^2, na.rm = TRUE)),
      empirical_sd = empirical_sd,
      mean_reported_se = mean_reported_se,
      se_calibration_ratio = mean_reported_se / empirical_sd,
      ci_coverage = mean(sub$ci_covers_truth, na.rm = TRUE),
      median_ci_width = stats::median(sub$ci_width, na.rm = TRUE),
      false_positive_rate = mean(sub$false_positive, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  )
}
coef_summary <- if (length(coef_summary_rows) > 0L) do.call(rbind, coef_summary_rows) else data.frame()

status_summary <- data.frame()
if (nrow(status) > 0L) {
  group_cols <- c("generator", "scenario", "correlation", "correlation_level", "n_time", "family", "method")
  groups <- unique(status[group_cols])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    idx <- rep(TRUE, nrow(status))
    for (col in group_cols) idx <- idx & status[[col]] == groups[[col]][i]
    sub <- status[idx, , drop = FALSE]
    success <- sub$success %in% c(TRUE, "TRUE", "True", "true", "1")
    converged <- if ("converged" %in% names(sub)) {
      sub$converged %in% c(TRUE, "TRUE", "True", "true", "1")
    } else {
      rep(NA, nrow(sub))
    }
    rows[[i]] <- cbind(
      groups[i, , drop = FALSE],
      data.frame(
        n = nrow(sub),
        success_rate = mean(success),
        convergence_rate = mean(converged[success], na.rm = TRUE),
        median_elapsed_sec = stats::median(suppressWarnings(as.numeric(sub$elapsed_sec)), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    )
  }
  status_summary <- do.call(rbind, rows)
}

complexity_summary <- data.frame()
if (nrow(complexity) > 0L) {
  for (nm in intersect(c("mean_df", "dependence_df", "total_df"), names(complexity))) {
    complexity[[nm]] <- suppressWarnings(as.numeric(complexity[[nm]]))
  }
  group_cols <- c("generator", "scenario", "correlation", "correlation_level", "n_time", "family", "method")
  groups <- unique(complexity[group_cols])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    idx <- rep(TRUE, nrow(complexity))
    for (col in group_cols) idx <- idx & complexity[[col]] == groups[[col]][i]
    sub <- complexity[idx, , drop = FALSE]
    rows[[i]] <- cbind(
      groups[i, , drop = FALSE],
      data.frame(
        n = nrow(sub),
        mean_df_median = stats::median(sub$mean_df, na.rm = TRUE),
        dependence_df_median = stats::median(sub$dependence_df, na.rm = TRUE),
        total_df_median = stats::median(sub$total_df, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    )
  }
  complexity_summary <- do.call(rbind, rows)
}

dependence_summary <- data.frame()
if (nrow(dependence) > 0L) {
  for (nm in intersect(c("dependence_n", "theta_mean_error", "theta_mae", "theta_rmse", "tau_mean_error", "tau_mae", "tau_rmse"), names(dependence))) {
    dependence[[nm]] <- suppressWarnings(as.numeric(dependence[[nm]]))
  }
  group_cols <- c("generator", "scenario", "correlation", "correlation_level", "n_time", "family", "method", "dependence_estimand")
  if ("dependence_scope" %in% names(dependence)) {
    group_cols <- c(group_cols, "dependence_scope")
  }
  groups <- unique(dependence[group_cols])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    idx <- rep(TRUE, nrow(dependence))
    for (col in group_cols) idx <- idx & dependence[[col]] == groups[[col]][i]
    sub <- dependence[idx, , drop = FALSE]
    rows[[i]] <- cbind(
      groups[i, , drop = FALSE],
      data.frame(
        n = nrow(sub),
        dependence_n_median = stats::median(sub$dependence_n, na.rm = TRUE),
        theta_mae_median = stats::median(sub$theta_mae, na.rm = TRUE),
        theta_rmse_median = stats::median(sub$theta_rmse, na.rm = TRUE),
        tau_mae_median = stats::median(sub$tau_mae, na.rm = TRUE),
        tau_rmse_median = stats::median(sub$tau_rmse, na.rm = TRUE),
        tau_mean_error_mean = mean(sub$tau_mean_error, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    )
  }
  dependence_summary <- do.call(rbind, rows)
}

bmk_write_csv(benchmark_summary, file.path(run_dir, "benchmark_summary.csv"))
bmk_write_csv(coef_summary, file.path(run_dir, "se_calibration_summary.csv"))
bmk_write_csv(status_summary, file.path(run_dir, "primary_status_summary.csv"))
bmk_write_csv(complexity_summary, file.path(run_dir, "fit_complexity_summary.csv"))
bmk_write_csv(dependence_summary, file.path(run_dir, "dependence_recovery_summary.csv"))

if (requireNamespace("ggplot2", quietly = TRUE) && nrow(coef_summary) > 0L) {
  p <- ggplot2::ggplot(
    coef_summary[coef_summary$term %in% c("x", "z_null"), , drop = FALSE],
    ggplot2::aes(x = method, y = se_calibration_ratio)
  ) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(family ~ n_time + scenario, scales = "free_x") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Mean reported SE / empirical SD") +
    ggplot2::theme_bw(base_size = 9)
  ggplot2::ggsave(file.path(run_dir, "se_calibration_ratio.png"), p, width = 12, height = 8, dpi = 150)

  p2 <- ggplot2::ggplot(
    coef_summary[coef_summary$term %in% c("x", "z_null"), , drop = FALSE],
    ggplot2::aes(x = method, y = ci_coverage)
  ) +
    ggplot2::geom_hline(yintercept = 0.95, linetype = 2, linewidth = 0.3) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(family ~ n_time + scenario, scales = "free_x") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Coefficient 95% CI coverage") +
    ggplot2::theme_bw(base_size = 9)
  ggplot2::ggsave(file.path(run_dir, "ci_coverage.png"), p2, width = 12, height = 8, dpi = 150)
}

report <- c(
  "# RS-Joint Standard-Model Benchmark Summary",
  "",
  paste0("- Run directory: `", run_dir, "`"),
  paste0("- Result rows: ", nrow(results)),
  paste0("- Coefficient rows: ", nrow(coef)),
  paste0("- Fit complexity rows: ", nrow(complexity)),
  paste0("- Dependence recovery rows: ", nrow(dependence)),
  "",
  "## Primary Fit Status",
  "",
  if (nrow(status_summary) > 0L) {
    utils::capture.output(print(status_summary, row.names = FALSE))
  } else {
    "No primary status rows found."
  },
  "",
  "## Benchmark Metrics",
  "",
  if (nrow(benchmark_summary) > 0L) {
    utils::capture.output(print(benchmark_summary, row.names = FALSE))
  } else {
    "No benchmark metric rows found."
  },
  "",
  "## SE Calibration",
  "",
  if (nrow(coef_summary) > 0L) {
    utils::capture.output(print(coef_summary, row.names = FALSE))
  } else {
    "No coefficient summary rows found."
  },
  "",
  "## Fit Complexity",
  "",
  if (nrow(complexity_summary) > 0L) {
    utils::capture.output(print(complexity_summary, row.names = FALSE))
  } else {
    "No fit complexity summary rows found."
  },
  "",
  "## Dependence Recovery",
  "",
  if (nrow(dependence_summary) > 0L) {
    utils::capture.output(print(dependence_summary, row.names = FALSE))
  } else {
    "No dependence recovery summary rows found."
  }
)
writeLines(report, file.path(run_dir, "benchmark_report.md"), useBytes = TRUE)

message("Summary complete: ", run_dir)
