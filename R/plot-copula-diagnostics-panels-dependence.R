.copula_v2_kendall_panel <- function(kendall_df) {
  if (nrow(kendall_df) == 0) {
    return(.copula_v2_message_plot(
      title = "Kendall-function calibration",
      subtitle = "Empirical against fitted copula probabilities",
      message = "No finite Kendall diagnostic values"
    ))
  }

  cols <- .gl_plot_colours()
  ggplot2::ggplot(kendall_df, ggplot2::aes(x = fitted, y = empirical)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = cols[["reference"]], linetype = "dashed") +
    ggplot2::geom_point(color = cols[["observed"]], alpha = 0.55, size = 1.2) +
    ggplot2::labs(
      title = "Kendall-function calibration",
      subtitle = "Empirical against fitted copula probabilities",
      x = "Fitted copula probability",
      y = "Empirical copula probability"
    ) +
    theme_gamlss_longitudinal()
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

  cols <- .gl_plot_colours()
  p <- ggplot2::ggplot(tail_long, ggplot2::aes(x = threshold, y = probability, color = source, group = source)) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~tail) +
    ggplot2::scale_color_manual(values = c(Empirical = cols[["observed"]], Fitted = cols[["fitted"]])) +
    ggplot2::labs(
      title = title,
      subtitle = plot_subtitle,
      x = "Tail probability a",
      y = y_label,
      color = NULL
    ) +
    theme_gamlss_longitudinal()

  if (!is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(ylim = ylim)
  }

  p
}

.copula_v2_residual_lag_panel <- function(lag_summary_df) {
  if (nrow(lag_summary_df) == 0 || all(!is.finite(lag_summary_df$cor_z))) {
    return(.copula_v2_message_plot(
      title = "Residual dependence by lag",
      subtitle = "Rosenblatt-score correlations within subject",
      message = "No finite residual lag correlations"
    ))
  }

  cols <- .gl_plot_colours()
  ggplot2::ggplot(lag_summary_df, ggplot2::aes(x = factor(lag), y = cor_z)) +
    ggplot2::geom_hline(yintercept = 0, color = cols[["reference"]], linetype = "dashed") +
    ggplot2::geom_col(fill = cols[["observed"]], alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = paste0("n=", n_pairs)), vjust = -0.35, size = 3) +
    ggplot2::coord_cartesian(ylim = c(-1, 1)) +
    ggplot2::labs(
      title = "Residual dependence by lag",
      subtitle = "Rosenblatt-score correlations should be near zero",
      x = "Lag",
      y = "Correlation"
    ) +
    theme_gamlss_longitudinal()
}
