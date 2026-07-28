#' Build fitted overlay data for plot_dist()
#'
#' @noRd
.plot_dist_fit_overlay_data <- function(fit, dataset, num_margins) {
  fit_diag_data <- tryCatch(
    .gl_fitted_distribution(fit, newdata = dataset, require_response = TRUE),
    error = function(e) .gl_fitted_distribution(fit, newdata = NULL, require_response = TRUE)
  )
  fit_data <- tryCatch(
    .copula_v2_fit_data(fit, data = dataset),
    error = function(e) .copula_v2_fit_data(fit)
  )
  fit_pair_data <- .copula_v2_pair_data(fit_data, lags = seq_len(max(1L, num_margins - 1L)))
  fit_copula_spec <- get_copula_dist(fit$copula_dist)
  fit_copula_spec <- list(
    family = .copula_family_code(fit_copula_spec$copula_dist),
    par = NA_real_,
    par2 = 0,
    tau = NA_real_
  )

  list(
    diag_data = fit_diag_data,
    pair_data = fit_pair_data,
    copula_spec = fit_copula_spec
  )
}

#' Add correlation summary text to a plot_dist() off-diagonal panel
#'
#' @noRd
.plot_dist_add_cor_stats <- function(p, input_data) {
  if (nrow(input_data) >= 3) {
    pearson_r <- suppressWarnings(cor(input_data$X1, input_data$X2, method = "pearson", use = "complete.obs"))
    kendall_tau <- suppressWarnings(cor(input_data$X1, input_data$X2, method = "kendall", use = "complete.obs"))
    stats_lab <- sprintf("Pearson r = %.3f | Kendall tau = %.3f", pearson_r, kendall_tau)
  } else {
    stats_lab <- "Pearson r = NA | Kendall tau = NA"
  }

  p + ggplot2::labs(subtitle = stats_lab)
}

#' Add a selected copula overlay to a plot_dist() off-diagonal panel
#'
#' @noRd
.plot_dist_add_selected_copula_overlay <- function(p, input_data, copula_spec, transform, grid_n, contour_bins) {
  contour_grid <- .plot_copula_density_for_spec(
    data.frame(
      u1 = input_data$X1,
      u2 = input_data$X2,
      theta_pair = rep(copula_spec$par, nrow(input_data)),
      zeta_pair = rep(copula_spec$par2, nrow(input_data))
    ),
    copula_spec,
    grid_n = grid_n,
    max_pairs_overlay = 300
  )
  contour_grid <- .plot_copula_transform_grid(contour_grid, transform)
  contour_grid$X1 <- contour_grid$u1
  contour_grid$X2 <- contour_grid$u2

  p + ggplot2::geom_contour(
    data = contour_grid,
    ggplot2::aes(x = .data$X1, y = .data$X2, z = .data$density),
    inherit.aes = FALSE,
    color = "#e41a1c",
    linewidth = 0.8,
    bins = contour_bins
  )
}

#' Add a fitted-model copula overlay to a plot_dist() off-diagonal panel
#'
#' @noRd
.plot_dist_add_model_copula_overlay <- function(
    p,
    i,
    j,
    time_values,
    fit_pair_data,
    fit_copula_spec,
    transform,
    grid_n,
    contour_bins) {
  left_idx <- min(i, j)
  right_idx <- max(i, j)
  pd <- fit_pair_data[
    as.character(fit_pair_data$time_left) == as.character(time_values[left_idx]) &
      as.character(fit_pair_data$time_right) == as.character(time_values[right_idx]), ,
    drop = FALSE
  ]
  if (nrow(pd) <= 2L) {
    return(p)
  }

  contour_grid <- .plot_copula_density_for_spec(
    pd,
    fit_copula_spec,
    grid_n = grid_n,
    max_pairs_overlay = 300
  )
  contour_grid <- .plot_copula_transform_grid(contour_grid, transform)
  if (i < j) {
    contour_grid$X1 <- contour_grid$u1
    contour_grid$X2 <- contour_grid$u2
  } else {
    contour_grid$X1 <- contour_grid$u2
    contour_grid$X2 <- contour_grid$u1
  }

  p + ggplot2::geom_contour(
    data = contour_grid,
    ggplot2::aes(x = .data$X1, y = .data$X2, z = .data$density),
    inherit.aes = FALSE,
    color = "#e41a1c",
    linewidth = 0.8,
    bins = contour_bins
  )
}
