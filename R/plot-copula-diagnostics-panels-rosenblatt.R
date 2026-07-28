.copula_v2_rosenblatt_time_panel <- function(rosenblatt_df) {
  if (nrow(rosenblatt_df) == 0 || all(!is.finite(rosenblatt_df$z))) {
    return(.copula_v2_message_plot(
      title = "Rosenblatt Normal Scores by Time",
      subtitle = "Scores are qnorm of pairwise conditional Rosenblatt residuals",
      message = "No finite Rosenblatt residuals"
    ))
  }

  ggplot2::ggplot(rosenblatt_df, ggplot2::aes(x = factor(time), y = z)) +
    ggplot2::geom_hline(yintercept = 0, color = "#666666", linetype = "dashed") +
    ggplot2::geom_boxplot(fill = "#9ecae1", color = "#4d4d4d", outlier.alpha = 0.35) +
    ggplot2::labs(
      title = "Rosenblatt Normal Scores by Time",
      subtitle = "Each time point should be centered near zero with similar spread",
      x = "Time",
      y = "Normal score"
    ) +
    ggplot2::theme_minimal()
}

.copula_v2_rosenblatt_qq_panel <- function(rosenblatt_qq_df) {
  if (nrow(rosenblatt_qq_df) == 0) {
    return(.copula_v2_message_plot(
      title = "Rosenblatt Normal QQ",
      subtitle = "Conditional copula scores should follow N(0, 1)",
      message = "No finite Rosenblatt residuals"
    ))
  }

  ggplot2::ggplot(rosenblatt_qq_df, ggplot2::aes(x = theoretical, y = observed)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "#666666", linetype = "dashed") +
    ggplot2::geom_point(color = "#4d4d4d", alpha = 0.55, size = 1.2) +
    ggplot2::labs(
      title = "Rosenblatt Normal QQ",
      subtitle = "Conditional copula scores should follow N(0, 1)",
      x = "Theoretical normal quantile",
      y = "Observed normal score"
    ) +
    ggplot2::theme_minimal()
}

.copula_v2_rosenblatt_lag_panel <- function(rosenblatt_pair_df) {
  if (nrow(rosenblatt_pair_df) == 0 || all(!is.finite(rosenblatt_pair_df$z_prev)) || all(!is.finite(rosenblatt_pair_df$z_curr))) {
    return(.copula_v2_message_plot(
      title = "Rosenblatt Lag Plot",
      subtitle = "Current conditional score against previous marginal score",
      message = "No finite Rosenblatt lag pairs"
    ))
  }

  ggplot2::ggplot(rosenblatt_pair_df, ggplot2::aes(x = z_prev, y = z_curr)) +
    ggplot2::geom_hline(yintercept = 0, color = "#d9d9d9") +
    ggplot2::geom_vline(xintercept = 0, color = "#d9d9d9") +
    ggplot2::geom_point(color = "#4d4d4d", alpha = 0.35, size = 1.1) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#e41a1c", linewidth = 0.7) +
    ggplot2::labs(
      title = "Rosenblatt Lag Plot",
      subtitle = "The smooth should be approximately flat at zero",
      x = expression(Phi^-1 * (U[t])),
      y = expression(Phi^-1 * (R[t + 1] ~ "|" ~ U[t]))
    ) +
    ggplot2::theme_minimal()
}
