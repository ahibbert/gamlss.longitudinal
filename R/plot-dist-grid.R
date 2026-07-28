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

  plots
}
