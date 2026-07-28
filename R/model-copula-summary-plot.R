#' Plot fitted copula trends by time

#'

#' @param x A `copula_time_summary` object or a fitted `gamlss.longitudinal` object.

#' @param ... Additional arguments (currently unused).

#' @param lags Integer lags to assess, measured in ordered time steps.

#' @param stat Character summary statistic for fitted values, one of "mean" or "median".

#' @param plot Logical; if TRUE, print the plot.

#'

#' @return Invisibly returns a list with the summary data and plot objects.

#' @export

plot.copula_time_summary <- function(x, ..., lags = 1, stat = c("mean", "median"), plot = TRUE) {
  summary_out <- if (inherits(x, "copula_time_summary")) {
    x
  } else {
    copula_time_summary(object = x, lags = lags, stat = stat)
  }

  time_summary <- summary_out$time_summary

  if (nrow(time_summary) == 0) {
    stop("No fitted copula summaries are available for plotting.")
  }

  time_summary$time <- as.factor(time_summary$time)

  p_theta <- ggplot2::ggplot(time_summary, ggplot2::aes(x = time, y = theta_fit, group = 1)) +
    ggplot2::geom_line(color = "#1f4e79", linewidth = 0.8) +
    ggplot2::geom_point(color = "#1f4e79", size = 2.5) +
    ggplot2::labs(
      title = "Fitted Copula Theta by Time",
      x = "Time",
      y = "Theta"
    ) +
    ggplot2::theme_minimal()

  p_tau <- ggplot2::ggplot(time_summary, ggplot2::aes(x = time, y = tau_fit, group = 1)) +
    ggplot2::geom_line(color = "#e41a1c", linewidth = 0.8) +
    ggplot2::geom_point(color = "#e41a1c", size = 2.5) +
    ggplot2::labs(
      title = "Fitted Copula Kendall's Tau by Time",
      x = "Time",
      y = "Tau"
    ) +
    ggplot2::theme_minimal()

  if ("zeta_fit" %in% names(time_summary)) {
    p_zeta <- ggplot2::ggplot(time_summary, ggplot2::aes(x = time, y = zeta_fit, group = 1)) +
      ggplot2::geom_line(color = "#4d4d4d", linewidth = 0.8) +
      ggplot2::geom_point(color = "#4d4d4d", size = 2.5) +
      ggplot2::labs(
        title = "Fitted Copula Zeta by Time",
        x = "Time",
        y = "Zeta"
      ) +
      ggplot2::theme_minimal()
  } else {
    p_zeta <- NULL
  }

  dashboard <- if (is.null(p_zeta)) {
    ggpubr::ggarrange(p_theta, p_tau, ncol = 1, nrow = 2)
  } else {
    ggpubr::ggarrange(p_theta, p_tau, p_zeta, ncol = 1, nrow = 3)
  }

  if (isTRUE(plot)) {
    print(dashboard)
  }

  invisible(list(
    summary = summary_out,
    p_theta = p_theta,
    p_tau = p_tau,
    p_zeta = p_zeta,
    dashboard = dashboard
  ))
}
