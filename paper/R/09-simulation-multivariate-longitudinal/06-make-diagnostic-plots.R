source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))

mvt_load_package()
run_dir <- mvt_read_run_dir()
fig_dir <- file.path(run_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required for diagnostic plots.", call. = FALSE)
}

coef_path <- file.path(run_dir, "coefficient_summary.csv")
bench_path <- file.path(run_dir, "benchmark_summary.csv")
dep_path <- file.path(run_dir, "dependence_recovery_summary.csv")
if (!file.exists(coef_path) || !file.exists(bench_path) || !file.exists(dep_path)) {
  mvt_summarise_results(run_dir)
}

coef_summary <- mvt_read_optional_csv(coef_path)
bench_summary <- mvt_read_optional_csv(bench_path)
dep_summary <- mvt_read_optional_csv(dep_path)

prepare_plot_summary <- function(x) {
  if (!"scenario" %in% names(x)) {
    x$scenario <- "all"
  }
  scenario_label <- sub("_t[0-9]+$", "", x$scenario)
  scenario_label <- sub("^external_", "", scenario_label)
  scenario_label <- sub("^native_", "", scenario_label)
  scenario_label <- gsub("_", " ", scenario_label)
  x$scenario_label <- factor(scenario_label, levels = unique(scenario_label))
  x
}

scenario_points <- function() {
  ggplot2::geom_point(
    ggplot2::aes(
      colour = scenario_label,
      shape = scenario_label,
      group = scenario_label
    ),
    position = ggplot2::position_dodge(width = 0.65),
    size = 1.7
  )
}

diagnostic_theme <- function() {
  ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank()
    )
}

scenario_caption <- paste(
  "Each point is the scenario-level mean across completed replicates;",
  "colour and shape identify the dependence scenario."
)

if (nrow(coef_summary) > 0L && all(c("term", "metric", "mean", "method") %in% names(coef_summary))) {
  bias <- coef_summary[coef_summary$metric == "bias" & coef_summary$term %in% c("time", "x", "z"), , drop = FALSE]
  if (nrow(bias) > 0L) {
    bias <- prepare_plot_summary(bias)
    p <- ggplot2::ggplot(bias, ggplot2::aes(x = method, y = mean)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
      scenario_points() +
      ggplot2::facet_grid(term ~ n_time + family, scales = "free_x") +
      ggplot2::coord_flip() +
      diagnostic_theme() +
      ggplot2::labs(x = NULL, y = "Mean coefficient bias", caption = scenario_caption)
    ggplot2::ggsave(file.path(fig_dir, "coefficient_bias.png"), p, width = 12, height = 8, dpi = 150)
  }
}

if (nrow(coef_summary) > 0L && all(c("term", "metric", "mean", "method") %in% names(coef_summary))) {
  coverage <- coef_summary[
    coef_summary$metric == "ci_covers_truth" & coef_summary$term %in% c("time", "x", "z"),
    ,
    drop = FALSE
  ]
  if (nrow(coverage) > 0L) {
    coverage <- prepare_plot_summary(coverage)
    p <- ggplot2::ggplot(coverage, ggplot2::aes(x = method, y = mean)) +
      ggplot2::geom_hline(yintercept = 0.95, linetype = 2, linewidth = 0.3) +
      scenario_points() +
      ggplot2::facet_grid(term ~ n_time + family, scales = "free_x") +
      ggplot2::coord_flip() +
      diagnostic_theme() +
      ggplot2::labs(x = NULL, y = "Empirical 95% CI coverage", caption = scenario_caption)
    ggplot2::ggsave(file.path(fig_dir, "coefficient_coverage.png"), p, width = 12, height = 8, dpi = 150)
  }
}

if (nrow(bench_summary) > 0L && all(c("metric", "mean", "method") %in% names(bench_summary))) {
  elapsed <- bench_summary[bench_summary$metric == "elapsed_sec", , drop = FALSE]
  if (nrow(elapsed) > 0L) {
    elapsed <- prepare_plot_summary(elapsed)
    p <- ggplot2::ggplot(elapsed, ggplot2::aes(x = method, y = mean)) +
      scenario_points() +
      ggplot2::facet_grid(family ~ n_time, scales = "free_x") +
      ggplot2::coord_flip() +
      diagnostic_theme() +
      ggplot2::labs(x = NULL, y = "Mean elapsed seconds", caption = scenario_caption)
    ggplot2::ggsave(file.path(fig_dir, "runtime.png"), p, width = 10, height = 7, dpi = 150)

    p_log <- ggplot2::ggplot(elapsed, ggplot2::aes(x = mean, y = method)) +
      ggplot2::geom_point(
        ggplot2::aes(colour = scenario_label, shape = scenario_label, group = scenario_label),
        position = ggplot2::position_dodge(width = 0.65),
        size = 1.7
      ) +
      ggplot2::scale_x_log10() +
      ggplot2::facet_grid(family ~ n_time, scales = "free_x") +
      diagnostic_theme() +
      ggplot2::labs(x = "Mean elapsed seconds (log10 scale)", y = NULL, caption = scenario_caption)
    ggplot2::ggsave(file.path(fig_dir, "runtime_log.png"), p_log, width = 10, height = 7, dpi = 150)
  }

  pit <- bench_summary[bench_summary$metric == "benchmark_pit_mean_abs_error", , drop = FALSE]
  if (nrow(pit) > 0L) {
    pit <- prepare_plot_summary(pit)
    p <- ggplot2::ggplot(pit, ggplot2::aes(x = method, y = mean)) +
      scenario_points() +
      ggplot2::facet_grid(family ~ n_time, scales = "free_x") +
      ggplot2::coord_flip() +
      diagnostic_theme() +
      ggplot2::labs(x = NULL, y = "PIT mean absolute error from 0.5", caption = scenario_caption)
    ggplot2::ggsave(file.path(fig_dir, "pit_calibration.png"), p, width = 10, height = 7, dpi = 150)
  }
}

if (nrow(dep_summary) > 0L && all(c("metric", "mean", "method") %in% names(dep_summary))) {
  dep_rmse <- dep_summary[dep_summary$metric %in% c("theta_rmse", "tau_rmse"), , drop = FALSE]
  if (nrow(dep_rmse) > 0L) {
    dep_rmse <- prepare_plot_summary(dep_rmse)
    p <- ggplot2::ggplot(dep_rmse, ggplot2::aes(x = method, y = mean)) +
      scenario_points() +
      ggplot2::facet_grid(metric + family ~ n_time, scales = "free_x") +
      ggplot2::coord_flip() +
      diagnostic_theme() +
      ggplot2::labs(x = NULL, y = "Dependence recovery RMSE", caption = scenario_caption)
    ggplot2::ggsave(file.path(fig_dir, "dependence_recovery.png"), p, width = 10, height = 9, dpi = 150)
  }
}

tables_dir <- file.path(run_dir, "paper_tables")
read_primary_table <- function(stem) {
  mvt_read_optional_csv(file.path(tables_dir, paste0(stem, ".csv")))
}

primary_longer <- function(dat, metrics) {
  metrics <- intersect(metrics, names(dat))
  if (nrow(dat) == 0L || length(metrics) == 0L) return(data.frame())
  rows <- lapply(metrics, function(metric) {
    out <- dat
    out$metric <- metric
    out$value <- suppressWarnings(as.numeric(out[[metric]]))
    out
  })
  out <- mvt_bind_rows_fill(rows)
  out[is.finite(out$value), , drop = FALSE]
}

dist_primary <- read_primary_table("primary_distributional_fit")
if (nrow(dist_primary) > 0L) {
  likelihood_long <- primary_longer(
    dist_primary,
    c("logLik", "AIC", "BIC")
  )
  if (nrow(likelihood_long) > 0L) {
    likelihood_long <- prepare_plot_summary(likelihood_long)
    p <- ggplot2::ggplot(likelihood_long, ggplot2::aes(x = value, y = method)) +
      ggplot2::geom_point(
        ggplot2::aes(colour = scenario_label, shape = scenario_label, group = scenario_label),
        position = ggplot2::position_dodge(width = 0.65),
        size = 1.7
      ) +
      ggplot2::facet_wrap(ggplot2::vars(metric, family), scales = "free_x", ncol = 4) +
      diagnostic_theme() +
      ggplot2::labs(
        x = "Scenario-level mean",
        y = NULL,
        caption = paste(
          scenario_caption,
          "AIC and BIC are shown only for methods returning likelihood-based criteria."
        )
      )
    ggplot2::ggsave(file.path(fig_dir, "primary_likelihood_criteria.png"), p, width = 11, height = 8, dpi = 150)
  }

  variogram_long <- primary_longer(dist_primary, c("variogram_score_p05", "variogram_score_p2"))
  if (nrow(variogram_long) > 0L) {
    variogram_long <- prepare_plot_summary(variogram_long)
    p_vario <- ggplot2::ggplot(variogram_long, ggplot2::aes(x = value, y = method)) +
      ggplot2::geom_point(
        ggplot2::aes(colour = scenario_label, shape = scenario_label, group = scenario_label),
        position = ggplot2::position_dodge(width = 0.65),
        size = 1.7
      ) +
      ggplot2::facet_wrap(ggplot2::vars(metric, family), scales = "free_x", ncol = 4) +
      diagnostic_theme() +
      ggplot2::labs(
        x = "Scenario-level mean variogram score",
        y = NULL,
        caption = paste(scenario_caption, "Lower variogram scores indicate better multivariate distributional fit.")
      )
    ggplot2::ggsave(file.path(fig_dir, "primary_variogram_scores.png"), p_vario, width = 11, height = 7, dpi = 150)
  }

  score_long <- primary_longer(dist_primary, c("benchmark_neg_log_score", "logLik"))
  if (nrow(score_long) > 0L) {
    score_long <- prepare_plot_summary(score_long)
    p_score <- ggplot2::ggplot(score_long, ggplot2::aes(x = value, y = method)) +
      ggplot2::geom_point(
        ggplot2::aes(colour = scenario_label, shape = scenario_label, group = scenario_label),
        position = ggplot2::position_dodge(width = 0.65),
        size = 1.7
      ) +
      ggplot2::facet_wrap(ggplot2::vars(metric, family), scales = "free_x", ncol = 4) +
      diagnostic_theme() +
      ggplot2::labs(
        x = "Scenario-level mean",
        y = NULL,
        caption = scenario_caption
      )
    ggplot2::ggsave(file.path(fig_dir, "primary_distributional_fit.png"), p_score, width = 11, height = 6, dpi = 150)
  }
}

coef_primary <- read_primary_table("primary_marginal_covariate_recovery")
if (nrow(coef_primary) > 0L && all(c("method", "term", "family", "marginal_rmse") %in% names(coef_primary))) {
  coef_primary <- prepare_plot_summary(coef_primary)
  p <- ggplot2::ggplot(coef_primary, ggplot2::aes(x = method, y = marginal_rmse)) +
    scenario_points() +
    ggplot2::facet_grid(term ~ family, scales = "free_x") +
    ggplot2::coord_flip() +
    diagnostic_theme() +
    ggplot2::labs(x = NULL, y = "Marginal coefficient RMSE", caption = scenario_caption)
  ggplot2::ggsave(file.path(fig_dir, "primary_marginal_covariate_rmse.png"), p, width = 11, height = 8, dpi = 150)
}

uncertainty_primary <- read_primary_table("primary_marginal_uncertainty")
if (nrow(uncertainty_primary) > 0L && all(c("method", "term", "family", "ci_coverage", "se_ratio") %in% names(uncertainty_primary))) {
  uncertainty_primary <- prepare_plot_summary(uncertainty_primary)
  p_cov <- ggplot2::ggplot(uncertainty_primary, ggplot2::aes(x = method, y = ci_coverage)) +
    ggplot2::geom_hline(yintercept = 0.95, linetype = 2, linewidth = 0.3) +
    scenario_points() +
    ggplot2::facet_grid(term ~ family, scales = "free_x") +
    ggplot2::coord_flip() +
    diagnostic_theme() +
    ggplot2::labs(x = NULL, y = "Empirical 95% CI coverage", caption = scenario_caption)
  ggplot2::ggsave(file.path(fig_dir, "primary_uncertainty_coverage.png"), p_cov, width = 11, height = 8, dpi = 150)

  se_ratio_dat <- uncertainty_primary[
    is.finite(uncertainty_primary$se_ratio) &
      (!"n_reps" %in% names(uncertainty_primary) | uncertainty_primary$n_reps >= 2),
    ,
    drop = FALSE
  ]
  if (nrow(se_ratio_dat) > 0L) {
    se_ratio_dat <- prepare_plot_summary(se_ratio_dat)
    p_se <- ggplot2::ggplot(se_ratio_dat, ggplot2::aes(x = method, y = se_ratio)) +
      ggplot2::geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3) +
      scenario_points() +
      ggplot2::facet_grid(term ~ family, scales = "free_x") +
      ggplot2::coord_flip() +
      diagnostic_theme() +
      ggplot2::labs(
        x = NULL,
        y = "Mean reported SE / empirical SD",
        caption = paste(
          scenario_caption,
          "SE-ratio points require at least two completed replicates; one-replicate scenarios are omitted."
        )
      )
    ggplot2::ggsave(file.path(fig_dir, "primary_uncertainty_se_ratio.png"), p_se, width = 11, height = 8, dpi = 150)
  }
}

dep_primary <- read_primary_table("primary_dependence_recovery")
if (nrow(dep_primary) > 0L) {
  dep_long <- primary_longer(dep_primary, c("theta_rmse", "tau_rmse"))
  if (nrow(dep_long) > 0L) {
    dep_long <- prepare_plot_summary(dep_long)
    p <- ggplot2::ggplot(dep_long, ggplot2::aes(x = method, y = value)) +
      scenario_points() +
      ggplot2::facet_grid(metric ~ family, scales = "free_x") +
      ggplot2::coord_flip() +
      diagnostic_theme() +
      ggplot2::labs(x = NULL, y = "Dependence recovery RMSE", caption = scenario_caption)
    ggplot2::ggsave(file.path(fig_dir, "primary_dependence_recovery.png"), p, width = 11, height = 7, dpi = 150)
  }
}

message("Diagnostic plots written to: ", fig_dir)
