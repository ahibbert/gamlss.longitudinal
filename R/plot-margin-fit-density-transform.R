#' Transform observed marginal responses for plotting
#'
#' @noRd
.plot_margin_transform_observed <- function(plot_data, response_scale) {
  if (response_scale == "response") {
    return(plot_data)
  }
  if (any(plot_data$response <= 0, na.rm = TRUE)) {
    stop("'response_scale = \"log\"' requires positive responses.", call. = FALSE)
  }
  plot_data$response <- log(plot_data$response)
  plot_data
}

#' Transform marginal density grids for plotting
#'
#' @noRd
.plot_margin_transform_density <- function(plot_data, response_scale) {
  if (response_scale == "response") {
    return(plot_data)
  }
  plot_data <- plot_data[is.finite(plot_data$response) & plot_data$response > 0, , drop = FALSE]
  plot_data <- .plot_margin_transform_observed(plot_data, response_scale)
  plot_data$density <- plot_data$density * exp(plot_data$response)
  plot_data
}
