#' Build a diagonal marginal panel for plot_dist()
#'
#' @noRd
.plot_dist_diagonal_panel <- function(
    i,
    margin_data,
    time_values,
    overlay,
    margin_dist,
    fit,
    fit_diag_data,
    grid_n) {
  input_data <- data.frame(X1 = margin_data[[i]]$response)
  x_lab <- latex2exp::TeX(paste("$Y_", i, "$"))

  p <- ggplot2::ggplot(input_data, ggplot2::aes(x = X1))
  if (overlay %in% c("margin", "model")) {
    p <- p + ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 30, na.rm = TRUE)
  } else {
    p <- p + ggplot2::geom_histogram(bins = 30, na.rm = TRUE)
  }
  p <- p + ggplot2::labs(x = x_lab)

  if (overlay == "margin") {
    params <- .plot_margin_constant_params(input_data$X1, margin_dist)
    density_grid <- .plot_margin_density_grid(input_data$X1, margin_dist, params, grid_n = grid_n)
    p <- p +
      ggplot2::geom_line(
        data = density_grid,
        ggplot2::aes(x = .data$response, y = .data$density),
        inherit.aes = FALSE,
        color = "#e41a1c",
        linewidth = 0.8
      ) +
      ggplot2::labs(y = "Density")
  }

  if (overlay == "model") {
    keep_time <- as.character(fit_diag_data$time) == as.character(time_values[i])
    density_grid <- .plot_margin_density_grid(
      fit_diag_data$response[keep_time],
      fit$margin_dist,
      lapply(fit_diag_data$params, function(x) x[keep_time]),
      grid_n = grid_n
    )
    p <- p +
      ggplot2::geom_line(
        data = density_grid,
        ggplot2::aes(x = .data$response, y = .data$density),
        inherit.aes = FALSE,
        color = "#e41a1c",
        linewidth = 0.8
      ) +
      ggplot2::labs(y = "Density")
  }

  p
}

#' Build an off-diagonal pairwise panel for plot_dist()
#'
#' @noRd
.plot_dist_offdiag_panel <- function(
    i,
    j,
    margin_data,
    margin_pseudo,
    offdiag_scale,
    transform,
    show_cor_stats,
    overlay,
    copula_spec,
    fit_pair_data,
    fit_copula_spec,
    time_values,
    grid_n,
    contour_bins) {
  offdiag_data <- .plot_dist_offdiag_data(
    i = i,
    j = j,
    margin_data = margin_data,
    margin_pseudo = margin_pseudo,
    offdiag_scale = offdiag_scale,
    transform = transform
  )
  input_data <- offdiag_data$input_data
  x_lab <- offdiag_data$x_lab
  y_lab <- offdiag_data$y_lab

  p <- ggplot2::ggplot(data = input_data, ggplot2::aes(x = X1, y = X2)) +
    ggplot2::geom_point(size = 0.4, alpha = 0.25, color = "black", na.rm = TRUE) +
    ggplot2::geom_density_2d(contour_var = "density", bins = 10, color = "black") +
    ggplot2::labs(x = x_lab, y = y_lab)

  if (show_cor_stats) {
    p <- .plot_dist_add_cor_stats(p, input_data)
  }

  if (offdiag_scale == "pseudo" && overlay == "copula") {
    p <- .plot_dist_add_selected_copula_overlay(
      p,
      input_data = input_data,
      copula_spec = copula_spec,
      transform = transform,
      grid_n = grid_n,
      contour_bins = contour_bins
    )
  }

  if (offdiag_scale == "pseudo" && overlay == "model") {
    p <- .plot_dist_add_model_copula_overlay(
      p,
      i = i,
      j = j,
      time_values = time_values,
      fit_pair_data = fit_pair_data,
      fit_copula_spec = fit_copula_spec,
      transform = transform,
      grid_n = grid_n,
      contour_bins = contour_bins
    )
  }

  p
}
