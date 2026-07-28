#' @export
rootogram.gamlss.longitudinal <- function(object, bins = 20, plot = TRUE, by_time = FALSE, by = NULL, data = NULL, ...) {
  diag_data <- .gl_fitted_distribution(object, newdata = NULL, require_response = TRUE)

  y <- diag_data$response

  params <- diag_data$params

  split_info <- .gl_diag_split_info(by_time, by, diag_data, data = data, plot_name = "rootogram")

  if (length(y) == 0) {
    stop("No finite observations available for the rootogram.")
  }

  breaks <- pretty(range(y, na.rm = TRUE), n = bins)

  if (length(unique(breaks)) < 2) {
    breaks <- seq(min(y) - 0.5, max(y) + 0.5, length.out = bins + 1)
  }

  root_df <- .gl_rootogram_plot_frame(y, params, breaks, diag_data$family, split_info)

  if (!plot) {
    return(root_df)
  }

  midpoint <- root_diff <- NULL

  p <- ggplot2::ggplot(root_df, ggplot2::aes(x = midpoint, y = root_diff)) +
    ggplot2::geom_hline(yintercept = 0, color = "#666666") +
    ggplot2::geom_col(fill = "#2c7fb8", alpha = 0.8, width = diff(range(root_df$midpoint)) / max(length(root_df$midpoint), 1)) +
    ggplot2::labs(
      x = "Response",
      y = expression(sqrt(O) - sqrt(E)),
      title = if (split_info$split_by) paste0("Rootogram by ", split_info$by) else "Rootogram"
    ) +
    ggplot2::theme_minimal()

  if (split_info$split_by) {
    p <- p + ggplot2::facet_wrap(~split_group, scales = "free_y")
  }

  p
}
