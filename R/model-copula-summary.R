#' Summarise fitted copula parameters by time

#'

#' @param object A fitted `gamlss.longitudinal` object.

#' @param lags Integer lags to assess, measured in ordered time steps.

#' @param stat Character summary statistic for fitted values, one of "mean" or "median".

#'

#' @return A data frame with fitted theta and tau summaries by time.

#' @export

copula_time_summary <- function(object, lags = 1, stat = c("mean", "median")) {

  if (!inherits(object, "gamlss.longitudinal")) {

    stop("'object' must be a fitted 'gamlss.longitudinal' object.")

  }


  stat <- match.arg(stat)

  copula_info <- get_copula_dist(object$copula_dist)

  has_zeta <- "zeta" %in% copula_info$parameters


  fit_data <- .copula_v2_fit_data(object)

  pair_data <- .copula_v2_pair_data(fit_data, lags = lags)


  agg_fun <- if (stat == "median") stats::median else mean


  time_summary <- do.call(rbind, lapply(split(fit_data, fit_data$time), function(x) {

    out <- data.frame(

      time = x$time[1],

      n_obs = nrow(x),

      theta_fit = agg_fun(x$theta_fit, na.rm = TRUE),

      tau_fit = agg_fun(x$tau_fit, na.rm = TRUE),

      stringsAsFactors = FALSE

    )

    if (has_zeta) {

      out$zeta_fit <- agg_fun(x$zeta_fit, na.rm = TRUE)

    }

    out

  }))


  pair_summary <- do.call(rbind, lapply(split(pair_data, pair_data$time_pair), function(x) {

    out <- data.frame(

      time_pair = x$time_pair[1],

      n_pairs = nrow(x),

      theta_pair = agg_fun(x$theta_pair, na.rm = TRUE),

      tau_pair = agg_fun(x$tau_fit, na.rm = TRUE),

      stringsAsFactors = FALSE

    )

    if (has_zeta) {

      out$zeta_pair <- agg_fun(x$zeta_pair, na.rm = TRUE)

    }

    out

  }))


  time_summary <- time_summary[order(time_summary$time), , drop = FALSE]


  if (!has_zeta) {

    # Keep the returned data tidy for one-parameter copulas.

    if ("zeta_fit" %in% names(fit_data)) {

      fit_data$zeta_fit <- NULL

    }

    if ("zeta_pair" %in% names(pair_data)) {

      pair_data$zeta_pair <- NULL

    }

  }


  out <- list(

    time_summary = time_summary,

    pair_summary = pair_summary,

    fit_data = fit_data,

    pair_data = pair_data

  )

  class(out) <- "copula_time_summary"

  out

}


#' @export

print.copula_time_summary <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  cat("\nCopula Dependence Summary\n")

  cat("-------------------------\n")

  if (!is.null(x$time_summary)) {

    cat("\nFitted dependence by time:\n")

    print(x$time_summary, digits = digits, row.names = FALSE)

  }

  if (!is.null(x$pair_summary)) {

    cat("\nAdjacent-pair dependence:\n")

    print(x$pair_summary, digits = digits, row.names = FALSE)

  }

  invisible(x)

}


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
