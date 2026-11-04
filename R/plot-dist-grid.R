#' Build the full plot_dist() panel grid
#'
#' @noRd
.plot_dist_panel_grid <- function(num_margins,
                                  margin_data,
                                  margin_pseudo,
                                  time_values,
                                  overlay,
                                  margin_dist,
                                  fit,
                                  fit_diag_data,
                                  offdiag_scale,
                                  transform,
                                  show_cor_stats,
                                  copula_spec,
                                  fit_pair_data,
                                  fit_copula_spec,
                                  grid_n,
                                  contour_bins) {
  plots <- list()

  z <- 1
  for (i in seq_len(num_margins)) {
    for (j in seq_len(num_margins)) {
      if (i == j) {
        p <- .plot_dist_diagonal_panel(
          i = i,
          margin_data = margin_data,
          time_values = time_values,
          overlay = overlay,
          margin_dist = margin_dist,
          fit = fit,
          fit_diag_data = fit_diag_data,
          grid_n = grid_n
        )
      }
      if (i != j) {
        p <- .plot_dist_offdiag_panel(
          i = i,
          j = j,
          margin_data = margin_data,
          margin_pseudo = margin_pseudo,
          offdiag_scale = offdiag_scale,
          transform = transform,
          show_cor_stats = show_cor_stats,
          overlay = overlay,
          copula_spec = copula_spec,
          fit_pair_data = fit_pair_data,
          fit_copula_spec = fit_copula_spec,
          time_values = time_values,
          grid_n = grid_n,
          contour_bins = contour_bins
        )
      }

      plots[[z]] <- p
      z <- z + 1
    }
  }

  .plot_dist_share_panel_axes(plots, num_margins = num_margins)
}

#' Give plot_dist() panels comparable coordinate systems
#'
#' Diagonal histograms share their x and y ranges. Off-diagonal panels share a
#' single range across both axes so that the same pseudo-observation or response
#' value occupies the same position in every copula panel.
#'
#' @noRd
.plot_dist_share_panel_axes <- function(plots, num_margins) {
  panel_matrix <- matrix(seq_along(plots), nrow = num_margins, byrow = TRUE)
  diagonal <- diag(panel_matrix)
  off_diagonal <- setdiff(seq_along(plots), diagonal)

  plots <- .plot_dist_apply_shared_axes(
    plots,
    indices = diagonal,
    same_xy_range = FALSE
  )
  .plot_dist_apply_shared_axes(
    plots,
    indices = off_diagonal,
    same_xy_range = TRUE
  )
}

#' Apply common coordinate limits to selected plot_dist() panels
#'
#' @noRd
.plot_dist_apply_shared_axes <- function(plots, indices, same_xy_range) {
  if (length(indices) < 1L) {
    return(plots)
  }

  ranges <- lapply(plots[indices], .plot_dist_trained_axis_ranges)
  x_limits <- range(unlist(lapply(ranges, `[[`, "x")), finite = TRUE)
  y_limits <- range(unlist(lapply(ranges, `[[`, "y")), finite = TRUE)

  if (isTRUE(same_xy_range)) {
    shared_limits <- range(c(x_limits, y_limits), finite = TRUE)
    x_limits <- shared_limits
    y_limits <- shared_limits
  }

  x_limits <- .plot_dist_expand_degenerate_range(x_limits)
  y_limits <- .plot_dist_expand_degenerate_range(y_limits)

  for (index in indices) {
    plots[[index]] <- plots[[index]] + ggplot2::coord_cartesian(
      xlim = x_limits,
      ylim = y_limits
    )
  }

  plots
}

#' Extract the data-trained ranges of a ggplot axis
#'
#' @noRd
.plot_dist_trained_axis_ranges <- function(plot) {
  built <- suppressWarnings(ggplot2::ggplot_build(plot))
  list(
    x = built$layout$panel_scales_x[[1L]]$range$range,
    y = built$layout$panel_scales_y[[1L]]$range$range
  )
}

#' Expand a zero-width coordinate range
#'
#' @noRd
.plot_dist_expand_degenerate_range <- function(x) {
  if (length(x) != 2L || any(!is.finite(x))) {
    return(c(0, 1))
  }
  if (x[[1L]] < x[[2L]]) {
    return(x)
  }

  padding <- max(abs(x[[1L]]) * 0.05, 0.5)
  x + c(-padding, padding)
}
