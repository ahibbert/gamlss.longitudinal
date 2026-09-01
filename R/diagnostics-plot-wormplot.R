#' @export
wormplot.gamlss.longitudinal <- function(object, randomize = NULL, seed = 1L, plot = TRUE, smooth = TRUE, by_time = FALSE, by = NULL, data = NULL, ...) {
  pit_out <- .gl_pit(object, randomize = randomize, seed = seed)

  split_info <- .gl_diag_split_info(by_time, by, pit_out$diag, data = data, plot_name = "wormplot")

  worm_frames <- .gl_worm_plot_frames(pit_out$pit, split_info)
  worm_df <- worm_frames$worm_df
  worm_band <- worm_frames$worm_band

  if (!plot) {
    return(worm_df)
  }

  theoretical <- detrended <- band_lower <- band_upper <- NULL

  p <- ggplot2::ggplot(worm_df, ggplot2::aes(x = theoretical, y = detrended)) +
    ggplot2::geom_ribbon(
      data = if (split_info$split_by) worm_df else worm_band,
      ggplot2::aes(x = theoretical, ymin = band_lower, ymax = band_upper),
      inherit.aes = FALSE,
      fill = "#9ecae1",
      alpha = 0.25
    ) +
    ggplot2::geom_point(size = 1.2, alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#666666") +
    ggplot2::labs(
      x = "Theoretical Normal Quantiles",
      y = "Detrended Quantiles",
      title = if (split_info$split_by) paste0("Worm Plot by ", split_info$by) else "Worm Plot"
    ) +
    ggplot2::theme_minimal()

  if (split_info$split_by) {
    p <- p + ggplot2::facet_wrap(~split_group, scales = "free")
  }

  if (smooth && nrow(worm_df) > 5) {
    p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#d95f0e", linewidth = 0.7)
  }

  p
}
