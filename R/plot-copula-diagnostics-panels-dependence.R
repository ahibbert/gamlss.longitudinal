.copula_v2_kendall_panel <- function(kendall_df) {
  if (nrow(kendall_df) == 0) {
    return(.copula_v2_message_plot(
      title = "Kendall Function Diagnostic",
      subtitle = "Empirical copula values compared with fitted copula values at observed pairs",
      message = "No finite Kendall diagnostic values"
    ))
  }

  ggplot2::ggplot(kendall_df, ggplot2::aes(x = fitted, y = empirical)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "#666666", linetype = "dashed") +
    ggplot2::geom_point(color = "#4d4d4d", alpha = 0.55, size = 1.2) +
    ggplot2::labs(
      title = "Kendall Function Diagnostic",
      subtitle = "Sorted empirical copula probabilities should track sorted fitted probabilities",
      x = "Fitted copula probability",
      y = "Empirical copula probability"
    ) +
    ggplot2::theme_minimal()
}

.copula_v2_tail_panel <- function(tail_long,
                                  title,
                                  message_subtitle,
                                  empty_message,
                                  plot_subtitle,
                                  y_label,
                                  ylim = NULL) {
  if (nrow(tail_long) == 0 || all(!is.finite(tail_long$probability))) {
    return(.copula_v2_message_plot(
      title = title,
      subtitle = message_subtitle,
      message = empty_message
    ))
  }

  p <- ggplot2::ggplot(tail_long, ggplot2::aes(x = threshold, y = probability, color = source, group = source)) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~tail) +
    ggplot2::scale_color_manual(values = c(Empirical = "#4d4d4d", Fitted = "#e41a1c")) +
    ggplot2::labs(
      title = title,
      subtitle = plot_subtitle,
      x = "Tail probability a",
      y = y_label,
      color = NULL
    ) +
    ggplot2::theme_minimal()

  if (!is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(ylim = ylim)
  }

  p
}

.copula_v2_residual_lag_panel <- function(lag_summary_df) {
  if (nrow(lag_summary_df) == 0 || all(!is.finite(lag_summary_df$cor_z))) {
    return(.copula_v2_message_plot(
      title = "Residual Dependence by Lag",
      subtitle = "Correlation of Rosenblatt normal scores within subject",
      message = "No finite residual lag correlations"
    ))
  }

  ggplot2::ggplot(lag_summary_df, ggplot2::aes(x = factor(lag), y = cor_z)) +
    ggplot2::geom_hline(yintercept = 0, color = "#666666", linetype = "dashed") +
    ggplot2::geom_col(fill = "#4d4d4d", alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = paste0("n=", n_pairs)), vjust = -0.35, size = 3) +
    ggplot2::coord_cartesian(ylim = c(-1, 1)) +
    ggplot2::labs(
      title = "Residual Dependence by Lag",
      subtitle = "Correlations should be close to zero after the Rosenblatt transform",
      x = "Lag",
      y = "Correlation"
    ) +
    ggplot2::theme_minimal()
}
