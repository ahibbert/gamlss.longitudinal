library(ggplot2)

input_path <- file.path(
  "results",
  "jss-replication",
  "expanded",
  "data",
  "03-joint-vs-separate-optimization-results.csv"
)

out_dir <- file.path(
  "results",
  "jss-replication",
  "expanded",
  "figures",
  "03-joint-vs-separate-100rep-diagnostics"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

results <- utils::read.csv(input_path, stringsAsFactors = FALSE)

as_num <- function(x) suppressWarnings(as.numeric(x))

numeric_cols <- c(
  "joint_loglik", "marginal_loglik", "copula_loglik", "elapsed_sec",
  "max_abs_param_error", "max_rel_param_error", "fitted_copula_tau",
  "true_copula_tau", "benchmark_interval_coverage_95",
  "benchmark_theta_time_abs_error"
)

for (col in intersect(numeric_cols, names(results))) {
  results[[col]] <- as_num(results[[col]])
}

results$method <- factor(results$method, levels = c("rs_separate", "rs_joint"))
results$family <- factor(results$family, levels = sort(unique(results$family)))
results$copula <- factor(results$copula, levels = c("N", "C", "G"))
results$design <- factor(results$design, levels = c("intercept", "covariate", "time_dependence"))
results$copula_tau_abs_error <- abs(results$fitted_copula_tau - results$true_copula_tau)
results$interval_coverage_abs_error <- abs(results$benchmark_interval_coverage_95 - 0.95)

plot_box_by_design <- function(data, metric, y_label, title, file_stem, hline = NULL) {
  data <- data[is.finite(data[[metric]]), , drop = FALSE]
  if (nrow(data) == 0L) {
    return(character())
  }

  paths <- character()
  pdf_path <- file.path(out_dir, paste0(file_stem, "-by-design.pdf"))
  grDevices::pdf(pdf_path, width = 11, height = 7.5)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (design_value in levels(results$design)) {
    plot_data <- data[data$design == design_value, , drop = FALSE]
    if (nrow(plot_data) == 0L) {
      next
    }

    p <- ggplot(plot_data, aes(x = method, y = .data[[metric]], fill = method)) +
      geom_boxplot(width = 0.65, outlier.alpha = 0.35, na.rm = TRUE) +
      facet_grid(family ~ copula, scales = "free_y") +
      labs(
        title = paste(title, "-", design_value),
        x = NULL,
        y = y_label
      ) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 25, hjust = 1)
      )

    if (!is.null(hline)) {
      p <- p + geom_hline(yintercept = hline, linetype = "dashed", colour = "grey35")
    }

    print(p)

    png_path <- file.path(out_dir, paste0(file_stem, "-", design_value, ".png"))
    ggsave(png_path, p, width = 11, height = 7.5, dpi = 320, bg = "white")
    paths <- c(paths, png_path)
  }

  c(pdf_path, paths)
}

plot_tau_scatter_by_design <- function(data) {
  data <- data[is.finite(data$fitted_copula_tau) & is.finite(data$true_copula_tau), , drop = FALSE]
  if (nrow(data) == 0L) {
    return(character())
  }

  paths <- character()
  pdf_path <- file.path(out_dir, "copula-tau-fitted-vs-true-by-design.pdf")
  grDevices::pdf(pdf_path, width = 11, height = 7.5)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (design_value in levels(results$design)) {
    plot_data <- data[data$design == design_value, , drop = FALSE]
    if (nrow(plot_data) == 0L) {
      next
    }

    p <- ggplot(plot_data, aes(x = true_copula_tau, y = fitted_copula_tau, colour = method)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey35") +
      geom_point(alpha = 0.55, size = 1.4) +
      facet_grid(family ~ copula, scales = "free") +
      labs(
        title = paste("Fitted versus true copula tau -", design_value),
        x = "True copula tau",
        y = "Fitted copula tau"
      ) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")

    print(p)

    png_path <- file.path(out_dir, paste0("copula-tau-fitted-vs-true-", design_value, ".png"))
    ggsave(png_path, p, width = 11, height = 7.5, dpi = 320, bg = "white")
    paths <- c(paths, png_path)
  }

  c(pdf_path, paths)
}

generated <- c(
  plot_box_by_design(
    results,
    "joint_loglik",
    "Joint log-likelihood",
    "Joint likelihood by method",
    "joint-loglik"
  ),
  plot_box_by_design(
    results,
    "marginal_loglik",
    "Marginal log-likelihood",
    "Marginal likelihood by method",
    "marginal-loglik"
  ),
  plot_box_by_design(
    results,
    "copula_loglik",
    "Copula log-likelihood",
    "Copula likelihood by method",
    "copula-loglik"
  ),
  plot_box_by_design(
    results,
    "max_abs_param_error",
    "Maximum absolute eta-scale parameter error",
    "Aggregate parameter error by method",
    "max-abs-param-error"
  ),
  plot_box_by_design(
    results,
    "copula_tau_abs_error",
    "Absolute copula tau error",
    "Copula tau error by method",
    "copula-tau-abs-error"
  ),
  plot_box_by_design(
    results,
    "benchmark_theta_time_abs_error",
    "Absolute time coefficient error for dependence",
    "Time-varying dependence coefficient error by method",
    "theta-time-abs-error"
  ),
  plot_box_by_design(
    results,
    "benchmark_interval_coverage_95",
    "Empirical 95% interval coverage",
    "Interval coverage by method",
    "interval-coverage-95",
    hline = 0.95
  ),
  plot_box_by_design(
    results,
    "interval_coverage_abs_error",
    "Absolute deviation from 0.95 coverage",
    "Interval coverage error by method",
    "interval-coverage-abs-error"
  ),
  plot_tau_scatter_by_design(results)
)

generated <- unique(generated[file.exists(generated)])
manifest <- data.frame(
  artifact = basename(generated),
  path = normalizePath(generated, winslash = "/", mustWork = FALSE),
  stringsAsFactors = FALSE
)

manifest_path <- file.path(out_dir, "plot-manifest.csv")
utils::write.csv(manifest, manifest_path, row.names = FALSE)

cat("Wrote", nrow(manifest), "plot artifacts to", out_dir, "\n")
cat("Manifest:", manifest_path, "\n")
