.copula_v2_empirical_overlay_panel <- function(pair_data_plot,
                                               density_grid,
                                               plot1_style,
                                               contour_bins,
                                               is_grouped,
                                               copula_family_name,
                                               x_label,
                                               y_label) {
  title <- "Empirical Copula with Fitted Overlay"
  subtitle <- paste0("Red contours show the average fitted copula across paired observations; family = ", copula_family_name)

  if (all(!is.finite(density_grid$density))) {
    return(.copula_v2_message_plot(
      title = title,
      subtitle = subtitle,
      message = "No finite fitted copula density"
    ))
  }

  p <- ggplot2::ggplot(pair_data_plot, ggplot2::aes(x = .data$u1, y = .data$u2))

  if (plot1_style == "scatter") {
    p <- p +
      ggplot2::geom_point(color = "#4d4d4d", alpha = 0.45, size = 1.2)
  } else {
    p <- p +
      ggplot2::geom_bin2d(bins = 25, alpha = 0.8) +
      ggplot2::scale_fill_gradient(low = "white", high = "black", name = "Count")
  }

  p <- p +
    ggplot2::geom_contour(
      data = density_grid,
      ggplot2::aes(x = u1, y = u2, z = density),
      inherit.aes = FALSE,
      color = "#e41a1c",
      linewidth = 1.2,
      bins = contour_bins
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label
    ) +
    ggplot2::theme_minimal()

  if (is_grouped) {
    p <- p + ggplot2::facet_wrap(~split_group)
  }

  p
}

.copula_v2_cut_summary_panel <- function(quartile_df, plot2_cuts, is_grouped, tau_ylim) {
  if (nrow(quartile_df) == 0 || all(!is.finite(quartile_df$tau_emp)) || all(!is.finite(quartile_df$tau_fit))) {
    return(.copula_v2_message_plot(
      title = "Observed vs Fitted Correlation by Quantile Bin",
      subtitle = "Bins are formed from fitted copula strength",
      message = "No finite cut summaries"
    ))
  }

  cut_levels <- paste0("C", sort(unique(as.integer(sub("^C", "", quartile_df$cut_group)))))
  quartile_df$cut_group <- factor(quartile_df$cut_group, levels = cut_levels)

  p <- ggplot2::ggplot(quartile_df, ggplot2::aes(x = cut_group)) +
    ggplot2::geom_point(ggplot2::aes(y = tau_emp), color = "#4d4d4d", size = 2.8) +
    ggplot2::geom_line(ggplot2::aes(y = tau_emp, group = 1), color = "#4d4d4d", linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = tau_fit), color = "#e41a1c", size = 2.8, shape = 4, stroke = 1.1) +
    ggplot2::geom_line(ggplot2::aes(y = tau_fit, group = 1), color = "#e41a1c", linewidth = 0.8, linetype = "dashed") +
    ggplot2::labs(
      title = "Observed vs Fitted Correlation by Quantile Bin",
      subtitle = paste0("", plot2_cuts, " cuts formed from fitted copula strength"),
      x = "Cut",
      y = "Kendall's tau"
    ) +
    ggplot2::theme_minimal()

  if (is_grouped) {
    p <- p + ggplot2::facet_wrap(~split_group, scales = if (is.null(tau_ylim)) "free_y" else "fixed")
  }

  if (!is.null(tau_ylim)) {
    p <- p + ggplot2::coord_cartesian(ylim = tau_ylim)
  }

  p
}
