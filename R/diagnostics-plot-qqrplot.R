#' @export
qqrplot.gamlss.longitudinal <- function(object, randomize = FALSE, plot = TRUE, by_time = FALSE, by = NULL, data = NULL, ...) {
  pit_out <- .gl_pit(object, randomize = randomize)

  z <- stats::qnorm(pmin(pmax(pit_out$pit, .Machine$double.eps), 1 - .Machine$double.eps))

  split_info <- .gl_diag_split_info(by_time, by, pit_out$diag, data = data, plot_name = "QQ plot")

  qq_df <- .gl_qq_plot_frame(z, split_info)

  if (!plot) {
    return(qq_df)
  }

  theoretical <- observed <- NULL

  p <- ggplot2::ggplot(qq_df, ggplot2::aes(x = theoretical, y = observed)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.7) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#666666") +
    ggplot2::labs(
      x = "Theoretical Normal Quantiles",
      y = "Randomized Quantiles",
      title = if (split_info$split_by) paste0("QQ Plot by ", split_info$by) else "QQ Plot"
    ) +
    ggplot2::theme_minimal()

  if (split_info$split_by) {
    p <- p + ggplot2::facet_wrap(~split_group, scales = "free")
  }

  p
}
