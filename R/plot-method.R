#' Plot diagnostics dashboard for fitted `gamlss.longitudinal` objects
#'
#' Displays four ggplot-based diagnostic panels by default:
#' 1) PIT histogram
#' 2) QQ residual plot
#' 3) Worm plot
#' 4) Rootogram
#'
#' Fitted-data and newdata quantile forecast panels can be requested
#' explicitly with `include_fitted_quantiles` and `include_newdata_quantiles`.
#'
#' @param x A fitted `gamlss.longitudinal` object.
#' @param y Unused; included for S3 generic compatibility.
#' @param data Optional data frame used as fallback for `newdata` plotting.
#' @param newdata Optional data frame for the newdata forecast panel.
#' @param newdata_n Number of rows to use from `data` when `newdata` is NULL.
#' @param quantiles Quantiles for forecast panels.
#' @param include_fitted_quantiles Logical; if TRUE, include fitted-data
#'   forecast quantiles in the dashboard.
#' @param include_newdata_quantiles Logical; if TRUE, include newdata forecast
#'   quantiles in the dashboard. Uses `newdata`, or the first `newdata_n` rows
#'   of `data` if `newdata` is NULL.
#' @param randomize Logical; randomized PIT/residual diagnostics.
#' @param time_stratified Logical; if TRUE, show time-stratified diagnostic plots.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a list of generated plot/data objects.
#' @export
plot.gamlss.longitudinal <- function(
    x,
    y,
    data = NULL,
    newdata = NULL,
    newdata_n = 8,
    quantiles = c(0.1, 0.5, 0.9),
    include_fitted_quantiles = FALSE,
    include_newdata_quantiles = FALSE,
    randomize = TRUE,
    time_stratified = FALSE,
    ...) {
  if (!inherits(x, "gamlss.longitudinal")) {
    stop("'x' must be of class 'gamlss.longitudinal'.")
  }

  diagnostic_plots <- .plot_method_diagnostic_panels(x, randomize = randomize, time_stratified = time_stratified)

  include_any_quantiles <- isTRUE(include_fitted_quantiles) || isTRUE(include_newdata_quantiles)
  if (include_any_quantiles && length(quantiles) == 0L) {
    stop("'quantiles' must contain at least one probability when forecast quantile panels are requested.")
  }
  if (include_any_quantiles) {
    q_cols <- .plot_method_quantile_columns(quantiles)
  }

  p_fit_quant <- NULL
  fc_fit_q <- NULL
  if (isTRUE(include_fitted_quantiles)) {
    fitted_quantiles <- .plot_method_fitted_quantile_panel(x, quantiles = quantiles, q_cols = q_cols)
    p_fit_quant <- fitted_quantiles$plot
    fc_fit_q <- fitted_quantiles$data
  }

  # Optional newdata forecast quantiles
  nd_use <- NULL
  p_new_quant <- NULL
  fc_new_q <- NULL
  if (isTRUE(include_newdata_quantiles)) {
    nd_use <- .plot_method_newdata_quantile_input(x, data = data, newdata = newdata, newdata_n = newdata_n)
    newdata_quantiles <- .plot_method_newdata_quantile_panel(x, nd_use = nd_use, quantiles = quantiles, q_cols = q_cols)
    p_new_quant <- newdata_quantiles$plot
    fc_new_q <- newdata_quantiles$data
  }

  if (isTRUE(include_newdata_quantiles) && is.null(p_new_quant)) {
    p_new_quant <- .plot_method_empty_plot("Newdata Forecast Quantiles", "Provide 'newdata' or 'data' for this panel")
  }

  dashboard_plots <- unname(diagnostic_plots)
  if (!is.null(p_fit_quant)) {
    dashboard_plots <- c(dashboard_plots, list(p_fit_quant))
  }
  if (!is.null(p_new_quant)) {
    dashboard_plots <- c(dashboard_plots, list(p_new_quant))
  }

  dashboard <- .plot_method_arrange_dashboard(dashboard_plots)
  print(dashboard)

  invisible(list(
    diagnostics = diagnostic_plots,
    forecasts = list(fitted_quantiles = p_fit_quant, newdata_quantiles = p_new_quant),
    fitted_data = fc_fit_q,
    newdata_data = fc_new_q,
    dashboard = dashboard
  ))
}
