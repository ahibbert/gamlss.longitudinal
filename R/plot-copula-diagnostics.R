#' Plot copula diagnostics for a fitted gamlss.longitudinal object

#'

#' @param x A fitted `gamlss.longitudinal` object.

#' @param lags Integer lags to assess, measured in ordered time steps.

#' @param grid_n Grid size used for contour averaging.

#' @param max_pairs_overlay Maximum number of paired observations used for the fitted overlay.

#' @param transform Character; "uniform" (default) shows empirical copula on

#'   the unit interval, "normal" transforms to standard normal scale.

#' @param plot1_style Character; "bins" (default) draws a binned empirical layer, "scatter" draws points.

#' @param contour_bins Integer number of contour levels for the fitted copula overlay in plot 1.

#' @param time_stratified Logical; if TRUE, facet both plots by time pair.

#' @param by Optional grouping variable name for stratified plots. Defaults to

#'   time-pair grouping when NULL. Use `data` for covariates not stored on the

#'   fitted pair object (for example gender).

#' @param data Optional data frame used when grouping by a covariate via `by`.

#' @param tau_ylim Optional numeric vector of length 2 specifying y-axis limits

#'   for Kendall's tau chart(s). If `NULL` (default), y-axis scales are automatic.

#' @param plot2_cuts Integer number of quantile-based cuts used in plot 2 (default 10).

#' @param tail_thresholds Numeric vector of lower-tail probabilities used for

#'   tail co-occurrence and conditional exceedance diagnostics.

#' @param residual_lags Integer lags used for Rosenblatt normal-score

#'   autocorrelation diagnostics.

#' @param dashboard_ncol Number of columns in the combined diagnostic dashboard.

#' @param plot Logical; if TRUE, print the dashboard.

#' @param ... Additional arguments reserved for future methods.

#'

#' @return Invisibly returns a list with plot objects and summaries.

#' @export

plot_copula_diagnostics <- function(
    x,
    lags = 1,
    grid_n = 35,
    max_pairs_overlay = 300,
    transform = "normal",
    plot1_style = "bins",
    contour_bins = 8,
    time_stratified = FALSE,
    by = NULL,
    data = NULL,
    tau_ylim = NULL,
    plot2_cuts = 10,
    tail_thresholds = c(0.05, 0.10, 0.20),
    residual_lags = 1:3,
    dashboard_ncol = 3,
    plot = TRUE,
    ...) {
  plot.copula(
    x = x,
    lags = lags,
    grid_n = grid_n,
    max_pairs_overlay = max_pairs_overlay,
    transform = transform,
    plot1_style = plot1_style,
    contour_bins = contour_bins,
    time_stratified = time_stratified,
    by = by,
    data = data,
    tau_ylim = tau_ylim,
    plot2_cuts = plot2_cuts,
    tail_thresholds = tail_thresholds,
    residual_lags = residual_lags,
    dashboard_ncol = dashboard_ncol,
    plot = plot,
    ...
  )
}
