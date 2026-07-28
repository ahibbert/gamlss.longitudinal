.plot_method_fitted_quantile_panel <- function(x, quantiles, q_cols) {
  fc_fit_q <- procast(x, type = "quantile", at = quantiles)
  fc_fit_q$idx <- seq_len(nrow(fc_fit_q))

  p_fit_quant <- ggplot2::ggplot(fc_fit_q, ggplot2::aes(x = idx)) +
    ggplot2::geom_ribbon(ggplot2::aes_string(ymin = q_cols$q_low, ymax = q_cols$q_high), fill = "#4e79a7", alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes_string(y = q_cols$q_mid), color = "#1f4e79", linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = response), color = "black", alpha = 0.35, size = 0.9) +
    ggplot2::labs(
      title = "Fitted Forecast Quantiles",
      x = "Observation Index",
      y = "Response"
    ) +
    ggplot2::theme_minimal()

  list(plot = p_fit_quant, data = fc_fit_q)
}

.plot_method_newdata_quantile_input <- function(x, data = NULL, newdata = NULL, newdata_n = 8) {
  nd_use <- newdata
  if (is.null(nd_use) && !is.null(data) && is.data.frame(data)) {
    nd_use <- utils::head(data, newdata_n)
    # ensure quantile-only mode (response optional)
    if (is.null(x$var_map) || !"response" %in% x$var_map) {
      nd_use$response <- NA_real_
    } else {
      response_orig <- names(x$var_map)[x$var_map == "response"][1]
      if (!is.na(response_orig) && !response_orig %in% names(nd_use) && !"response" %in% names(nd_use)) {
        nd_use[[response_orig]] <- NA_real_
      }
    }
  }
  nd_use
}

.plot_method_newdata_quantile_panel <- function(x, nd_use, quantiles, q_cols) {
  p_new_quant <- NULL
  fc_new_q <- NULL

  if (!is.null(nd_use) && is.data.frame(nd_use) && nrow(nd_use) > 0) {
    fc_new_q <- tryCatch(
      procast.gamlss.longitudinal(x, type = "quantile", at = quantiles, newdata = nd_use),
      error = function(e) NULL
    )

    if (!is.null(fc_new_q)) {
      time_candidates <- c("time", if (!is.null(x$var_map)) names(x$var_map)[x$var_map == "time"] else character(0))
      person_candidates <- c("subject", if (!is.null(x$var_map)) names(x$var_map)[x$var_map == "subject"] else character(0))
      time_col <- time_candidates[time_candidates %in% names(nd_use)][1]
      person_col <- person_candidates[person_candidates %in% names(nd_use)][1]

      if (is.na(time_col) || is.null(time_col) || nchar(time_col) == 0) {
        fc_new_q$time_plot <- seq_len(nrow(fc_new_q))
      } else {
        fc_new_q$time_plot <- nd_use[[time_col]]
      }
      if (is.na(person_col) || is.null(person_col) || nchar(person_col) == 0) {
        fc_new_q$person_plot <- factor(seq_len(nrow(fc_new_q)))
      } else {
        fc_new_q$person_plot <- as.factor(nd_use[[person_col]])
      }

      p_new_quant <- ggplot2::ggplot(fc_new_q, ggplot2::aes(x = time_plot, color = person_plot, group = person_plot)) +
        ggplot2::geom_ribbon(ggplot2::aes_string(ymin = q_cols$q_low, ymax = q_cols$q_high, fill = "person_plot"), alpha = 0.14, color = NA, show.legend = FALSE) +
        ggplot2::geom_line(ggplot2::aes_string(y = q_cols$q_mid), linewidth = 0.8) +
        ggplot2::geom_point(ggplot2::aes_string(y = q_cols$q_mid), size = 1.7) +
        ggplot2::labs(
          title = "Newdata Forecast Quantiles",
          x = "Time",
          y = paste0("Predicted ", q_cols$q_mid)
        ) +
        ggplot2::theme_minimal()
    }
  }

  list(plot = p_new_quant, data = fc_new_q)
}
