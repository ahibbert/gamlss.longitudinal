.copula_v2_rosenblatt_time_panel <- function(rosenblatt_df) {
  if (nrow(rosenblatt_df) == 0 || all(!is.finite(rosenblatt_df$z))) {
    return(.copula_v2_message_plot(
      title = "Rosenblatt scores by time",
      subtitle = "Conditional normal scores",
      message = "No finite Rosenblatt residuals"
    ))
  }

  cols <- .gl_plot_colours()
  ggplot2::ggplot(rosenblatt_df, ggplot2::aes(x = factor(time), y = z)) +
    ggplot2::geom_hline(yintercept = 0, color = cols[["reference"]], linetype = "dashed") +
    ggplot2::geom_boxplot(fill = cols[["interval"]], color = cols[["observed"]], outlier.alpha = 0.35) +
    ggplot2::labs(
      title = "Rosenblatt scores by time",
      subtitle = "Target: centred at 0 with similar spread",
      x = "Time",
      y = "Normal score"
    ) +
    theme_gamlss_longitudinal()
}

.copula_v2_rosenblatt_qq_panel <- function(rosenblatt_qq_df) {
  if (nrow(rosenblatt_qq_df) == 0) {
    return(.copula_v2_message_plot(
      title = "Rosenblatt Q-Q",
      subtitle = "Target: N(0, 1) conditional scores",
      message = "No finite Rosenblatt residuals"
    ))
  }

  cols <- .gl_plot_colours()
  ggplot2::ggplot(rosenblatt_qq_df, ggplot2::aes(x = theoretical, y = observed)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = cols[["reference"]], linetype = "dashed") +
    ggplot2::geom_point(color = cols[["observed"]], alpha = 0.55, size = 1.2) +
    ggplot2::labs(
      title = "Rosenblatt Q-Q",
      subtitle = "Target: N(0, 1) conditional scores",
      x = "Theoretical normal quantile",
      y = "Observed normal score"
    ) +
    theme_gamlss_longitudinal()
}

.copula_v2_rosenblatt_lag_panel <- function(rosenblatt_pair_df) {
  if (nrow(rosenblatt_pair_df) == 0 || all(!is.finite(rosenblatt_pair_df$z_prev)) || all(!is.finite(rosenblatt_pair_df$z_curr))) {
    return(.copula_v2_message_plot(
      title = "Rosenblatt Lag Plot",
      subtitle = "Current conditional vs previous score",
      message = "No finite Rosenblatt lag pairs"
    ))
  }

  cols <- .gl_plot_colours()
  conditional_pairs <- isTRUE(attr(
    rosenblatt_pair_df, "conditional_residual_pairs"
  ))
  x_label <- if (conditional_pairs) {
    expression(Phi^-1 * (R[t]))
  } else {
    expression(Phi^-1 * (U[t]))
  }
  y_label <- if (conditional_pairs) {
    expression(Phi^-1 * (R[t + 1]))
  } else {
    expression(Phi^-1 * (R[t + 1] ~ "|" ~ U[t]))
  }
  ggplot2::ggplot(rosenblatt_pair_df, ggplot2::aes(x = z_prev, y = z_curr)) +
    ggplot2::geom_hline(yintercept = 0, color = cols[["light"]]) +
    ggplot2::geom_vline(xintercept = 0, color = cols[["light"]]) +
    ggplot2::geom_point(color = cols[["observed"]], alpha = 0.35, size = 1.1) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = cols[["fitted"]], linewidth = 0.7) +
    ggplot2::labs(
      title = "Rosenblatt Lag Plot",
      subtitle = "Target: smooth is flat at 0",
      x = x_label,
      y = y_label
    ) +
    theme_gamlss_longitudinal()
}
