#' Normalize copula diagnostics plot controls
#'
#' Centralizes the argument validation and coercion for the copula
#' diagnostics plotting entry point. Kept separate so reviewers can inspect
#' plotting controls without stepping through the full plotting assembly.
#'
#' @noRd
.copula_v2_normalize_plot_controls <- function(transform,
                                               plot1_style,
                                               contour_bins,
                                               time_stratified,
                                               plot2_cuts,
                                               tau_ylim,
                                               tail_thresholds,
                                               residual_lags,
                                               dashboard_ncol) {
  if (!transform %in% c("uniform", "normal")) {
    stop("'transform' must be either 'uniform' or 'normal'.")
  }

  if (!plot1_style %in% c("bins", "scatter")) {
    stop("'plot1_style' must be either 'bins' or 'scatter'.")
  }

  if (!is.numeric(contour_bins) || length(contour_bins) != 1 || !is.finite(contour_bins) || contour_bins < 1) {
    stop("'contour_bins' must be a single finite number >= 1.")
  }
  contour_bins <- as.integer(round(contour_bins))

  if (!is.logical(time_stratified) || length(time_stratified) != 1 || is.na(time_stratified)) {
    stop("'time_stratified' must be TRUE or FALSE.")
  }

  if (!is.numeric(plot2_cuts) || length(plot2_cuts) != 1 || !is.finite(plot2_cuts) || plot2_cuts < 2) {
    stop("'plot2_cuts' must be a single finite number >= 2.")
  }
  plot2_cuts <- as.integer(round(plot2_cuts))

  if (!is.null(tau_ylim)) {
    if (!is.numeric(tau_ylim) || length(tau_ylim) != 2 || any(!is.finite(tau_ylim)) || tau_ylim[1] >= tau_ylim[2]) {
      stop("'tau_ylim' must be NULL or a numeric vector of length 2 with tau_ylim[1] < tau_ylim[2].")
    }
    tau_ylim <- as.numeric(tau_ylim)
  }

  tail_thresholds <- sort(unique(as.numeric(tail_thresholds)))
  tail_thresholds <- tail_thresholds[is.finite(tail_thresholds) & tail_thresholds > 0 & tail_thresholds < 0.5]
  if (length(tail_thresholds) == 0) {
    tail_thresholds <- c(0.05, 0.10, 0.20)
  }

  residual_lags <- sort(unique(as.integer(residual_lags)))
  residual_lags <- residual_lags[residual_lags > 0]
  if (length(residual_lags) == 0) {
    residual_lags <- 1:3
  }

  if (!is.numeric(dashboard_ncol) || length(dashboard_ncol) != 1 || !is.finite(dashboard_ncol) || dashboard_ncol < 1) {
    stop("'dashboard_ncol' must be a single positive integer.")
  }
  dashboard_ncol <- as.integer(round(dashboard_ncol))

  list(
    transform = transform,
    plot1_style = plot1_style,
    contour_bins = contour_bins,
    time_stratified = time_stratified,
    plot2_cuts = plot2_cuts,
    tau_ylim = tau_ylim,
    tail_thresholds = tail_thresholds,
    residual_lags = residual_lags,
    dashboard_ncol = dashboard_ncol
  )
}
