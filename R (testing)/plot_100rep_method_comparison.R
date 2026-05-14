#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

method_dirs <- data.frame(
  method = c(
    "RS separate",
    "RS dlcopdpar",
    "CG dlcopdpar",
    "gamlss2"
  ),
  dir = c(
    "results/bcpe_t_gamlss_comparison_rs_separate_absorbcons_100rep",
    "results/bcpe_t_gamlss_comparison_rs_dlcopdpar_absorbcons_100rep",
    "results/bcpe_t_gamlss_comparison_cg_dlcopdpar_absorbcons_delta025_lambda1_100rep",
    "results/bcpe_t_gamlss_comparison_rs_separate_absorbcons_100rep"
  ),
  source_model = c(
    "gamlss.longitudinal",
    "gamlss.longitudinal",
    "gamlss.longitudinal",
    "gamlss2"
  ),
  stringsAsFactors = FALSE
)

output_dir <- "results/bcpe_t_gamlss_comparison_method_comparison_100rep"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

smooth_ribbon_lower <- as.numeric(Sys.getenv("SMOOTH_RIBBON_LOWER", unset = "0.10"))
smooth_ribbon_upper <- as.numeric(Sys.getenv("SMOOTH_RIBBON_UPPER", unset = "0.90"))

read_fixed_summary <- function(method, dir, source_model) {
  path <- file.path(dir, "fixed_effects_bias_rmse_table.csv")
  fixed <- read.csv(path, stringsAsFactors = FALSE)
  fixed <- fixed[fixed$model == source_model, , drop = FALSE]
  fixed$method <- method
  fixed
}

read_smooth_summary <- function(method, dir, source_model) {
  rds_path <- file.path(dir, "all_results.rds")
  all_results <- readRDS(rds_path)
  smooth <- all_results$smooth
  smooth <- smooth[smooth$model == source_model, , drop = FALSE]

  smooth_stats <- aggregate(
    smooth_hat ~ scenario + n + d + parameter + s1,
    data = smooth,
    FUN = function(x) {
      c(
        median = stats::median(x, na.rm = TRUE),
        lower = unname(stats::quantile(x, smooth_ribbon_lower, na.rm = TRUE)),
        upper = unname(stats::quantile(x, smooth_ribbon_upper, na.rm = TRUE))
      )
    }
  )
  smooth_mat <- smooth_stats$smooth_hat
  if (!is.matrix(smooth_mat)) {
    smooth_mat <- do.call(rbind, smooth_mat)
  }
  smooth_stats$smooth_median <- smooth_mat[, "median"]
  smooth_stats$smooth_q_lower <- smooth_mat[, "lower"]
  smooth_stats$smooth_q_upper <- smooth_mat[, "upper"]
  smooth_stats$smooth_hat <- NULL

  smooth_truth <- aggregate(
    smooth_true ~ scenario + n + d + parameter + s1,
    data = smooth,
    FUN = mean,
    na.rm = TRUE
  )

  out <- merge(
    smooth_truth,
    smooth_stats,
    by = c("scenario", "n", "d", "parameter", "s1"),
    all.x = TRUE
  )
  out$method <- method
  out
}

read_method_metrics <- function(method, dir, source_model) {
  runs <- read.csv(file.path(dir, "fit_run_log.csv"), stringsAsFactors = FALSE)
  joint <- read.csv(
    file.path(dir, "joint_distribution_metrics_summary.csv"),
    stringsAsFactors = FALSE
  )
  fixed <- read.csv(
    file.path(dir, "fixed_effects_bias_rmse_table.csv"),
    stringsAsFactors = FALSE
  )
  smooth <- read.csv(
    file.path(dir, "smooth_integrated_metrics.csv"),
    stringsAsFactors = FALSE
  )
  se <- read.csv(
    file.path(dir, "fixed_effects_se_calibration.csv"),
    stringsAsFactors = FALSE
  )

  marginal_parameters <- c("mu", "sigma", "nu", "tau")
  fixed_i <- fixed[
    fixed$model == source_model &
      fixed$parameter %in% marginal_parameters,
    ,
    drop = FALSE
  ]
  fixed_nonintercept <- fixed_i[fixed_i$term != "intercept", , drop = FALSE]
  fixed_intercept <- fixed_i[fixed_i$term == "intercept", , drop = FALSE]

  se_i <- se[
    se$model == source_model &
      se$parameter %in% marginal_parameters &
      se$term != "intercept" &
      is.finite(se$se_to_empirical_sd),
    ,
    drop = FALSE
  ]

  get_smooth_irmse <- function(parameter) {
    value <- smooth$irmse[
      smooth$model == source_model &
        smooth$parameter == parameter
    ]
    if (length(value) == 0L) NA_real_ else value[1]
  }

  joint_i <- joint[joint$model == source_model, , drop = FALSE]
  runs_i <- runs[runs$model == source_model, , drop = FALSE]

  data.frame(
    method = method,
    success = paste0(sum(runs_i$success), "/", nrow(runs_i)),
    mean_logLik = joint_i$mean_logLik,
    sd_logLik = joint_i$sd_logLik,
    mean_df = mean(runs_i$df, na.rm = TRUE),
    mean_seconds = mean(runs_i$elapsed_sec, na.rm = TRUE),
    rosenblatt_z_lag1 = joint_i$mean_abs_rosenblatt_normal_lag1_cor,
    rosenblatt_cvm = joint_i$mean_rosenblatt_cvm,
    marginal_intercept_rmse = mean(fixed_intercept$rmse, na.rm = TRUE),
    marginal_nonintercept_rmse = mean(fixed_nonintercept$rmse, na.rm = TRUE),
    mu_smooth_irmse = get_smooth_irmse("mu"),
    sigma_smooth_irmse = get_smooth_irmse("sigma"),
    theta_smooth_irmse = get_smooth_irmse("theta"),
    se_to_sd = mean(se_i$se_to_empirical_sd, na.rm = TRUE),
    coverage_95 = mean(se_i$coverage_95, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

fmt_num <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))
}

tex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x)
  x
}

write_metric_table <- function(metrics, output_dir) {
  csv_path <- file.path(output_dir, "method_comparison_summary.csv")
  write.csv(metrics, csv_path, row.names = FALSE)

  lines <- c(
    "% Auto-generated method comparison summary.",
    "% Requires \\usepackage{booktabs}.",
    "",
    "\\begin{table}[ht]",
    "\\centering",
    "\\caption{Comparison of model performance across the 100-replication BCPE-$t$ simulation ($n = 500$, $d = 4$). Larger log-likelihood is better. Smaller Rosenblatt residual correlation, RMSE, IRMSE, and CvM values are better. SE/SD values closer to one indicate better standard-error calibration.}",
    "\\label{tab:bcpe-t-method-comparison-100rep}",
    "\\begin{tabular}{lrrrrrrrr}",
    "\\toprule",
    "Model & logLik & Rosen. $|r_1|$ & Fixed RMSE & $\\mu$ IRMSE & $\\sigma$ IRMSE & $\\theta$ IRMSE & SE/SD & Cover. \\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(metrics))) {
    lines <- c(
      lines,
      paste(
        tex_escape(metrics$method[i]),
        fmt_num(metrics$mean_logLik[i], 1),
        fmt_num(metrics$rosenblatt_z_lag1[i], 3),
        fmt_num(metrics$marginal_nonintercept_rmse[i], 3),
        fmt_num(metrics$mu_smooth_irmse[i], 3),
        fmt_num(metrics$sigma_smooth_irmse[i], 3),
        fmt_num(metrics$theta_smooth_irmse[i], 3),
        fmt_num(metrics$se_to_sd[i], 3),
        fmt_num(metrics$coverage_95[i], 3),
        sep = " & "
      ) |>
        paste0(" \\\\")
    )
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\begin{table}[ht]",
    "\\centering",
    "\\caption{Additional fit information for the same 100-replication comparison. Mean elapsed time is per fitted model replication. Mean degrees of freedom are available for gamlss2 only in these outputs.}",
    "\\label{tab:bcpe-t-method-comparison-100rep-fit-info}",
    "\\begin{tabular}{lrrrrr}",
    "\\toprule",
    "Model & Success & SD logLik & Mean df & Mean time (s) & Rosen. CvM \\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(metrics))) {
    lines <- c(
      lines,
      paste(
        tex_escape(metrics$method[i]),
        tex_escape(metrics$success[i]),
        fmt_num(metrics$sd_logLik[i], 1),
        fmt_num(metrics$mean_df[i], 1),
        fmt_num(metrics$mean_seconds[i], 1),
        fmt_num(metrics$rosenblatt_cvm[i], 3),
        sep = " & "
      ) |>
        paste0(" \\\\")
    )
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )

  tex_path <- file.path(output_dir, "method_comparison_summary.tex")
  writeLines(lines, tex_path, useBytes = TRUE)
}

fixed_summary <- do.call(
  rbind,
  Map(
    read_fixed_summary,
    method_dirs$method,
    method_dirs$dir,
    method_dirs$source_model
  )
)

metric_summary <- do.call(
  rbind,
  Map(
    read_method_metrics,
    method_dirs$method,
    method_dirs$dir,
    method_dirs$source_model
  )
)
write_metric_table(metric_summary, output_dir)

smooth_summary <- do.call(
  rbind,
  Map(
    read_smooth_summary,
    method_dirs$method,
    method_dirs$dir,
    method_dirs$source_model
  )
)

method_levels <- method_dirs$method
fixed_summary$method <- factor(fixed_summary$method, levels = method_levels)
smooth_summary$method <- factor(smooth_summary$method, levels = method_levels)

fixed_summary$parameter <- factor(
  fixed_summary$parameter,
  levels = c("mu", "sigma", "nu", "tau", "theta", "zeta")
)
fixed_summary$term <- factor(
  fixed_summary$term,
  levels = c("intercept", "x1", "x2", "t")
)
fixed_summary <- fixed_summary[!is.na(fixed_summary$parameter), , drop = FALSE]

smooth_summary$parameter <- factor(
  smooth_summary$parameter,
  levels = c("mu", "sigma", "theta")
)
smooth_summary <- smooth_summary[!is.na(smooth_summary$parameter), , drop = FALSE]

add_smooth_blanks <- function(smooth_summary) {
  template <- expand.grid(
    method = levels(smooth_summary$method),
    parameter = levels(smooth_summary$parameter),
    scenario = unique(smooth_summary$scenario),
    stringsAsFactors = FALSE
  )
  present <- unique(smooth_summary[c("method", "parameter", "scenario")])
  template$key <- paste(template$method, template$parameter, template$scenario, sep = "\r")
  present$key <- paste(present$method, present$parameter, present$scenario, sep = "\r")
  missing <- template[!template$key %in% present$key, , drop = FALSE]
  if (nrow(missing) == 0L) {
    return(smooth_summary)
  }

  missing_rows <- data.frame(
    scenario = missing$scenario,
    n = unique(smooth_summary$n)[1],
    d = unique(smooth_summary$d)[1],
    parameter = missing$parameter,
    s1 = NA_real_,
    smooth_true = NA_real_,
    smooth_median = NA_real_,
    smooth_q_lower = NA_real_,
    smooth_q_upper = NA_real_,
    method = missing$method,
    stringsAsFactors = FALSE
  )
  missing_rows$method <- factor(missing_rows$method, levels = levels(smooth_summary$method))
  missing_rows$parameter <- factor(
    missing_rows$parameter,
    levels = levels(smooth_summary$parameter)
  )

  rbind(smooth_summary, missing_rows)
}

smooth_summary <- add_smooth_blanks(smooth_summary)

draw_column_plot <- function(plot_list, path, width, height) {
  grDevices::png(path, width = width, height = height, units = "in", res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(1, length(plot_list))
  ))
  for (i in seq_along(plot_list)) {
    print(plot_list[[i]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
  }
  grid::popViewport()
}

fixed_parameters <- levels(droplevels(fixed_summary$parameter))
fixed_plots <- lapply(seq_along(fixed_parameters), function(i) {
  par_i <- fixed_parameters[i]
  dat_i <- fixed_summary[fixed_summary$parameter == par_i, , drop = FALSE]
  y_values <- unlist(dat_i[c("q05", "q95", "mean_estimate", "true_value")])
  y_values <- y_values[is.finite(y_values)]
  y_range <- range(y_values)
  y_pad <- diff(y_range) * 0.08
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.05

  ggplot(dat_i, aes(x = term, y = mean_estimate)) +
    geom_errorbar(
      aes(ymin = q05, ymax = q95),
      width = 0.18,
      color = "gray55",
      na.rm = TRUE
    ) +
    geom_point(size = 1.4, color = "black", na.rm = TRUE) +
    geom_point(
      aes(y = true_value),
      shape = 4,
      size = 2.3,
      stroke = 0.9,
      color = "#B3262E",
      na.rm = TRUE
    ) +
    facet_grid(method + scenario ~ ., scales = "fixed") +
    coord_cartesian(ylim = c(y_range[1] - y_pad, y_range[2] + y_pad)) +
    labs(
      x = NULL,
      y = if (i == 1L) "Estimate on eta scale" else NULL,
      title = as.character(par_i)
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.spacing = grid::unit(0.55, "lines"),
      plot.title = element_text(hjust = 0.5, size = 10),
      strip.text.y = if (i == length(fixed_parameters)) element_text(size = 8) else element_blank(),
      strip.background.y = if (i == length(fixed_parameters)) {
        element_rect(fill = "grey85", color = "grey40")
      } else {
        element_blank()
      }
    )
})

draw_column_plot(
  fixed_plots,
  file.path(output_dir, "fixed_effect_method_comparison_plot.png"),
  width = 15,
  height = 10
)

smooth_parameters <- levels(droplevels(smooth_summary$parameter))
smooth_plots <- lapply(seq_along(smooth_parameters), function(i) {
  par_i <- smooth_parameters[i]
  dat_i <- smooth_summary[smooth_summary$parameter == par_i, , drop = FALSE]
  y_values <- unlist(dat_i[c(
    "smooth_q_lower",
    "smooth_q_upper",
    "smooth_median",
    "smooth_true"
  )])
  y_values <- y_values[is.finite(y_values)]
  y_range <- range(y_values)
  y_pad <- diff(y_range) * 0.08
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.05

  ggplot(dat_i, aes(x = s1)) +
    geom_ribbon(
      aes(ymin = smooth_q_lower, ymax = smooth_q_upper),
      fill = "gray82",
      alpha = 0.85
    ) +
    geom_line(
      aes(y = smooth_median),
      color = "gray30",
      linewidth = 0.6,
      linetype = "dashed",
      na.rm = TRUE
    ) +
    geom_line(
      aes(y = smooth_true),
      color = "#B3262E",
      linewidth = 0.8,
      na.rm = TRUE
    ) +
    facet_grid(method + scenario ~ ., scales = "fixed") +
    coord_cartesian(ylim = c(y_range[1] - y_pad, y_range[2] + y_pad)) +
    labs(
      x = "s1",
      y = if (i == 1L) "Centered smooth contribution" else NULL,
      title = as.character(par_i)
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.spacing = grid::unit(0.55, "lines"),
      plot.title = element_text(hjust = 0.5, size = 10),
      strip.text.y = if (i == length(smooth_parameters)) element_text(size = 8) else element_blank(),
      strip.background.y = if (i == length(smooth_parameters)) {
        element_rect(fill = "grey85", color = "grey40")
      } else {
        element_blank()
      }
    )
})

draw_column_plot(
  smooth_plots,
  file.path(output_dir, "smooth_method_comparison_plot.png"),
  width = 12,
  height = 10
)

message("Wrote method comparison plots and tables to: ", output_dir)
