#' Compare fitted and empirical copula contour surfaces

#'

#' @param x A fitted `gamlss.longitudinal` object.

#' @param lags Integer lags to assess, measured in ordered time steps.

#' @param grid_n Grid size used for density surfaces.

#' @param max_pairs_overlay Maximum number of paired observations used for fitted surface averaging.

#' @param contour_bins Number of contour levels to draw in the surface panels.

#' @param transform Character; "uniform" compares surfaces on copula scale, "normal" compares them on z-scale.

#' @param diff_scale_limit Positive numeric; fixed symmetric color scale limit for the difference panel.

#' @param time_stratified Logical; if TRUE, compare surfaces by time pair.

#' @param plot Logical; if TRUE, print the dashboard.

#' @param ... Additional arguments reserved for future methods.

#'

#' @return Invisibly returns plots, grid-level surfaces, and numeric similarity metrics.

#' @export

plot.copula_contour_compare <- function(x, lags = 1, grid_n = 45, max_pairs_overlay = 300, contour_bins = 10, transform = "uniform", diff_scale_limit = 0.05, time_stratified = FALSE, plot = TRUE, ...) {
  if (!inherits(x, "gamlss.longitudinal")) {
    stop("'x' must be a fitted 'gamlss.longitudinal' object.")
  }

  object <- x

  controls <- .copula_contour_compare_controls(transform = transform, diff_scale_limit = diff_scale_limit)
  transform <- controls$transform
  diff_scale_limit <- controls$diff_scale_limit

  fit_data <- .copula_v2_fit_data(object)

  pair_data <- .copula_v2_pair_data(fit_data, lags = lags)

  copula_spec <- get_copula_dist(object$copula_dist)

  copula_family_name <- .copula_family_code(copula_spec$copula_dist)

  family_num <- tryCatch(
    {
      .copula_family_code(copula_family_name)
    },
    error = function(e) NA_character_
  )

  split_data <- if (isTRUE(time_stratified)) split(pair_data, pair_data$time_pair) else list(All = pair_data)

  grid_list <- lapply(names(split_data), function(nm) {
    pd <- split_data[[nm]]

    fit_grid <- .copula_v2_average_density_grid(
      family_num = family_num,
      pair_data = pd,
      grid_n = grid_n,
      max_pairs_overlay = max_pairs_overlay
    )

    # Build empirical surface on the same copula grid as fit_grid, then transform both

    # together if requested. This avoids grid mismatch artifacts in contouring.

    emp_grid <- .copula_v2_empirical_density_grid(pd, grid_n = grid_n, lims = c(0.02, 0.98, 0.02, 0.98))

    if (transform == "normal") {
      z1 <- stats::qnorm(.copula_v2_clamp01(fit_grid$u1))

      z2 <- stats::qnorm(.copula_v2_clamp01(fit_grid$u2))

      jacobian <- stats::dnorm(z1) * stats::dnorm(z2)

      fit_grid$u1 <- z1

      fit_grid$u2 <- z2

      fit_grid$density <- fit_grid$density * jacobian

      emp_grid$u1 <- z1

      emp_grid$u2 <- z2

      emp_grid$density <- emp_grid$density * jacobian
    } else {
      emp_grid <- emp_grid
    }

    # Merge on grid coordinates to ensure pointwise comparisons.

    g <- merge(

      emp_grid,
      fit_grid,
      by = c("u1", "u2"),
      suffixes = c("_emp", "_fit"),
      all = FALSE
    )

    g$density_diff <- g$density_fit - g$density_emp

    g$time_pair <- nm

    g
  })

  grid_df <- do.call(rbind, grid_list)

  metric_tables <- .copula_contour_compare_metric_tables(grid_df)
  metric_summary <- metric_tables$summary
  metric_overlap <- metric_tables$overlap

  axis_labels <- .copula_contour_compare_axis_labels(transform)
  x_label <- axis_labels$x
  y_label <- axis_labels$y

  p_emp <- ggplot2::ggplot(grid_df, ggplot2::aes(x = u1, y = u2, z = density_emp)) +
    ggplot2::geom_contour(color = "#4d4d4d", bins = contour_bins, linewidth = 0.9) +
    ggplot2::labs(title = "Empirical Copula Contours", x = x_label, y = y_label) +
    ggplot2::theme_minimal()

  p_fit <- ggplot2::ggplot(grid_df, ggplot2::aes(x = u1, y = u2, z = density_fit)) +
    ggplot2::geom_contour(color = "#e41a1c", bins = contour_bins, linewidth = 0.9) +
    ggplot2::labs(title = "Fitted Copula Contours", x = x_label, y = y_label) +
    ggplot2::theme_minimal()

  p_diff <- ggplot2::ggplot(grid_df, ggplot2::aes(x = u1, y = u2, fill = density_diff)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-diff_scale_limit, diff_scale_limit),
      oob = scales::squish,
      name = "Fit - Emp"
    ) +
    ggplot2::labs(title = "Contour Difference Surface", x = x_label, y = y_label) +
    ggplot2::theme_minimal()

  if (isTRUE(time_stratified)) {
    p_emp <- p_emp + ggplot2::facet_wrap(~time_pair)

    p_fit <- p_fit + ggplot2::facet_wrap(~time_pair)

    p_diff <- p_diff + ggplot2::facet_wrap(~time_pair)
  }

  dashboard <- ggpubr::ggarrange(p_emp, p_fit, p_diff, ncol = 1, nrow = 3)

  if (isTRUE(plot)) {
    print(dashboard)
  }

  invisible(list(
    plots = list(empirical_contours = p_emp, fitted_contours = p_fit, difference_surface = p_diff),
    dashboard = dashboard,
    grid = grid_df,
    metrics = list(summary = metric_summary, overlap = metric_overlap)
  ))
}
