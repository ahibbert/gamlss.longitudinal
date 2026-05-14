#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

default_dirs <- c(
  "results/bcpe_t_gamlss_comparison_rs_separate_absorbcons_100rep",
  "results/bcpe_t_gamlss_comparison_rs_dlcopdpar_absorbcons_100rep"
)

args <- commandArgs(trailingOnly = TRUE)
out_dirs <- if (length(args) > 0L) args else default_dirs

model_labels <- c(
  "gamlss.longitudinal" = "gamlss.longitudinal",
  "gamlss2" = "gamlss2"
)

smooth_ribbon_lower <- as.numeric(Sys.getenv("SMOOTH_RIBBON_LOWER", unset = "0.10"))
smooth_ribbon_upper <- as.numeric(Sys.getenv("SMOOTH_RIBBON_UPPER", unset = "0.90"))
if (
  !is.finite(smooth_ribbon_lower) ||
    !is.finite(smooth_ribbon_upper) ||
    smooth_ribbon_lower < 0 ||
    smooth_ribbon_upper > 1 ||
    smooth_ribbon_lower >= smooth_ribbon_upper
) {
  stop("SMOOTH_RIBBON_LOWER and SMOOTH_RIBBON_UPPER must be valid probabilities.")
}

summarise_smooth_from_rds <- function(out_dir) {
  rds_path <- file.path(out_dir, "all_results.rds")
  if (!file.exists(rds_path)) {
    return(NULL)
  }

  all_results <- readRDS(rds_path)
  if (!is.list(all_results) || is.null(all_results$smooth)) {
    return(NULL)
  }

  smooth_raw <- all_results$smooth
  key_cols <- c("scenario", "model", "n", "d", "parameter", "s1")
  smooth_stats <- aggregate(
    smooth_hat ~ scenario + model + n + d + parameter + s1,
    data = smooth_raw,
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
    smooth_true ~ scenario + model + n + d + parameter + s1,
    data = smooth_raw,
    FUN = mean,
    na.rm = TRUE
  )

  merge(smooth_truth, smooth_stats, by = key_cols, all.x = TRUE)
}

regenerate_plots <- function(out_dir) {
  smooth_path <- file.path(out_dir, "smooth_pointwise_summary.csv")
  fixed_path <- file.path(out_dir, "fixed_effects_bias_rmse_table.csv")

  if (!file.exists(smooth_path)) {
    stop("Missing smooth summary: ", smooth_path)
  }
  if (!file.exists(fixed_path)) {
    stop("Missing fixed-effect summary: ", fixed_path)
  }

  smooth_summary <- summarise_smooth_from_rds(out_dir)
  if (is.null(smooth_summary)) {
    smooth_summary <- read.csv(smooth_path, stringsAsFactors = FALSE)
    smooth_summary$smooth_q_lower <- smooth_summary$smooth_q05
    smooth_summary$smooth_q_upper <- smooth_summary$smooth_q95
  }
  fixed_summary <- read.csv(fixed_path, stringsAsFactors = FALSE)

  smooth_summary$model <- factor(
    smooth_summary$model,
    levels = names(model_labels),
    labels = unname(model_labels)
  )
  fixed_summary$model <- factor(
    fixed_summary$model,
    levels = names(model_labels),
    labels = unname(model_labels)
  )

  smooth_summary$parameter <- factor(
    smooth_summary$parameter,
    levels = c("mu", "sigma", "theta")
  )
  smooth_summary <- smooth_summary[!is.na(smooth_summary$parameter), , drop = FALSE]
  fixed_summary$parameter <- factor(
    fixed_summary$parameter,
    levels = c("mu", "sigma", "nu", "tau", "theta", "zeta")
  )
  fixed_summary$term <- factor(
    fixed_summary$term,
    levels = c("intercept", "x1", "x2", "t")
  )

  smooth_parameters <- levels(droplevels(smooth_summary$parameter))
  smooth_panel_template <- unique(smooth_summary[c("scenario", "model")])
  smooth_plots <- lapply(seq_along(smooth_parameters), function(i) {
    par_i <- smooth_parameters[i]
    dat_i <- smooth_summary[smooth_summary$parameter == par_i, , drop = FALSE]
    if (nrow(dat_i) > 0L) {
      present_panels <- unique(dat_i[c("scenario", "model")])
      template_key <- paste(smooth_panel_template$scenario, smooth_panel_template$model, sep = "\r")
      present_key <- paste(present_panels$scenario, present_panels$model, sep = "\r")
      missing_panels <- smooth_panel_template[!template_key %in% present_key, , drop = FALSE]
      if (nrow(missing_panels) > 0L) {
        missing_rows <- data.frame(
          scenario = missing_panels$scenario,
          model = missing_panels$model,
          n = dat_i$n[1],
          d = dat_i$d[1],
          parameter = par_i,
          s1 = NA_real_,
          smooth_true = NA_real_,
          smooth_median = NA_real_,
          smooth_q_lower = NA_real_,
          smooth_q_upper = NA_real_,
          stringsAsFactors = FALSE
        )
        missing_rows$model <- factor(missing_rows$model, levels = levels(smooth_summary$model))
        missing_rows$parameter <- factor(
          missing_rows$parameter,
          levels = levels(smooth_summary$parameter)
        )
        dat_i <- rbind(dat_i, missing_rows)
      }
    }
    y_values <- unlist(dat_i[c(
      "smooth_q_lower",
      "smooth_q_upper",
      "smooth_median",
      "smooth_true"
    )])
    y_values <- y_values[is.finite(y_values)]
    y_range <- range(y_values)
    if (length(y_values) == 0L || !all(is.finite(y_range))) {
      y_range <- c(-1, 1)
    }
    y_pad <- diff(y_range) * 0.08
    if (!is.finite(y_pad) || y_pad == 0) {
      y_pad <- max(0.05, abs(y_range[1]) * 0.08)
    }

    ggplot(dat_i, aes(x = s1)) +
      geom_ribbon(
        aes(ymin = smooth_q_lower, ymax = smooth_q_upper),
        fill = "gray82",
        alpha = 0.85
      ) +
      geom_line(
        aes(y = smooth_median),
        color = "gray30",
        linewidth = 0.65,
        linetype = "dashed",
        na.rm = TRUE
      ) +
      geom_line(
        aes(y = smooth_true),
        color = "#B3262E",
        linewidth = 0.85,
        na.rm = TRUE
      ) +
      facet_grid(model + scenario ~ ., scales = "fixed") +
      coord_cartesian(ylim = c(y_range[1] - y_pad, y_range[2] + y_pad)) +
      labs(
        x = "s1",
        y = if (i == 1L) "Centered smooth contribution" else NULL,
        title = as.character(par_i)
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.spacing = grid::unit(0.7, "lines"),
        plot.title = element_text(hjust = 0.5, size = 10),
        strip.text.y = if (i == length(smooth_parameters)) {
          element_text()
        } else {
          element_blank()
        },
        strip.background.y = if (i == length(smooth_parameters)) {
          element_rect(fill = "grey85", color = "grey40")
        } else {
          element_blank()
        }
      )
  })

  smooth_plot_path <- file.path(out_dir, "smooth_recovery_plot.png")
  grDevices::png(
    smooth_plot_path,
    width = 11,
    height = 6,
    units = "in",
    res = 180
  )
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(1, length(smooth_plots))
  ))
  for (i in seq_along(smooth_plots)) {
    print(
      smooth_plots[[i]],
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i)
    )
  }
  grid::popViewport()
  grDevices::dev.off()

  fixed_parameters <- levels(droplevels(fixed_summary$parameter))
  fixed_plots <- lapply(seq_along(fixed_parameters), function(i) {
    par_i <- fixed_parameters[i]
    dat_i <- fixed_summary[fixed_summary$parameter == par_i, , drop = FALSE]
    y_values <- unlist(dat_i[c("q05", "q95", "mean_estimate", "true_value")])
    y_values <- y_values[is.finite(y_values)]
    y_range <- range(y_values)
    if (length(y_values) == 0L || !all(is.finite(y_range))) {
      y_range <- c(-1, 1)
    }
    y_pad <- diff(y_range) * 0.08
    if (!is.finite(y_pad) || y_pad == 0) {
      y_pad <- max(0.05, abs(y_range[1]) * 0.08)
    }

    ggplot(dat_i, aes(x = term, y = mean_estimate)) +
      geom_errorbar(
        aes(ymin = q05, ymax = q95),
        width = 0.18,
        color = "gray55",
        na.rm = TRUE
      ) +
      geom_point(size = 1.6, color = "black", na.rm = TRUE) +
      geom_point(
        aes(y = true_value),
        shape = 4,
        size = 2.7,
        stroke = 1.0,
        color = "#B3262E",
        na.rm = TRUE
      ) +
      facet_grid(model + scenario ~ ., scales = "fixed") +
      coord_cartesian(ylim = c(y_range[1] - y_pad, y_range[2] + y_pad)) +
      labs(
        x = NULL,
        y = if (i == 1L) "Estimate on eta scale" else NULL,
        title = as.character(par_i)
      ) +
      theme_bw(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1),
        panel.spacing = grid::unit(0.7, "lines"),
        plot.title = element_text(hjust = 0.5, size = 10),
        strip.text.y = if (i == length(fixed_parameters)) {
          element_text()
        } else {
          element_blank()
        },
        strip.background.y = if (i == length(fixed_parameters)) {
          element_rect(fill = "grey85", color = "grey40")
        } else {
          element_blank()
        }
      )
  })

  fixed_plot_path <- file.path(out_dir, "fixed_effect_recovery_plot.png")
  grDevices::png(
    fixed_plot_path,
    width = 14,
    height = 7,
    units = "in",
    res = 180
  )
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(1, length(fixed_plots))
  ))
  for (i in seq_along(fixed_plots)) {
    print(
      fixed_plots[[i]],
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i)
    )
  }
  grid::popViewport()
  grDevices::dev.off()

  message(
    "Regenerated shared-y recovery plots in: ",
    out_dir,
    " (smooth ribbon ",
    100 * smooth_ribbon_lower,
    "-",
    100 * smooth_ribbon_upper,
    "%)"
  )
}

invisible(lapply(out_dirs, regenerate_plots))
