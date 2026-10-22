#' Plot copula diagnostics for a fitted gamlss.longitudinal object

#'

#' @param x A fitted `gamlss.longitudinal` object.

#' @param type Diagnostic set to display. `"overview"` (the default) shows four
#'   complementary checks; use `"calibration"`, `"residuals"`, `"tails"`, or
#'   `"all"` for focused or complete dashboards.

#' @param empirical Empirical layer for the copula overlay: `"bins"` or
#'   `"points"`.

#' @param n_bins Number of bins used for the fitted-versus-observed dependence
#'   calibration panel.

#' @param ncol Number of dashboard columns.

#' @param lags Integer lags to assess, measured in ordered time steps.

#' @param grid_n Grid size used for contour averaging.

#' @param max_pairs_overlay Maximum number of paired observations used for the fitted overlay.

#' @param transform Character; `"normal"` (the default) shows standard-normal

#'   scores, while `"uniform"` shows the empirical copula on the unit interval.

#' @param contour_bins Integer number of contour levels for the fitted copula overlay in plot 1.

#' @param time_stratified Logical; if TRUE, facet both plots by time pair.

#' @param by Optional grouping variable name for stratified plots. Defaults to

#'   time-pair grouping when NULL. Use `data` for covariates not stored on the

#'   fitted pair object (for example gender).

#' @param data Optional data frame used when grouping by a covariate via `by`.

#' @param tau_ylim Optional numeric vector of length 2 specifying y-axis limits

#'   for Kendall's tau chart(s). If `NULL` (default), y-axis scales are automatic.

#' @param tail_thresholds Numeric vector of lower-tail probabilities used for

#'   tail co-occurrence and conditional exceedance diagnostics.

#' @param residual_lags Integer lags used for Rosenblatt normal-score

#'   autocorrelation diagnostics.

#' @param randomized_residuals Logical; for discrete margins, use seeded
#'   copula-conditional latent completions and exact randomized conditional PIT
#'   residuals. The default is `TRUE`; set to `FALSE` to withhold diagnostics
#'   that require latent continuous values.
#' @param randomized_seed Finite numeric seed used for the discrete
#'   randomization. The caller's random-number state is preserved.
#' @param randomized_draws Positive integer number of copula-conditional latent
#'   completions per observed adjacent pair for discrete copula, Kendall, and
#'   tail panels. Multiple completions reduce Monte Carlo sensitivity while
#'   retaining the observed probability rectangles.

#' @param plot Logical; if TRUE, print the dashboard.

#' @param ... Deprecated aliases `plot1_style`, `plot2_cuts`, and
#'   `dashboard_ncol` are accepted temporarily.
#' @param plot1_style Deprecated alias for `empirical`, retained by the S3
#'   method for compatibility.
#' @param plot2_cuts Deprecated alias for `n_bins`, retained by the S3 method
#'   for compatibility.
#' @param dashboard_ncol Deprecated alias for `ncol`, retained by the S3 method
#'   for compatibility.

#'

#' @return Invisibly returns a `gamlss_longitudinal_plot` object containing the
#'   displayed dashboard, all component plots, and their diagnostic data.

#' @export

plot_copula_diagnostics <- function(
    x,
    lags = 1,
    grid_n = 35,
    max_pairs_overlay = 300,
    transform = "normal",
    plot1_style = NULL,
    contour_bins = 8,
    time_stratified = FALSE,
    by = NULL,
    data = NULL,
    tau_ylim = NULL,
    plot2_cuts = NULL,
    tail_thresholds = c(0.05, 0.10, 0.20),
    residual_lags = 1:3,
    dashboard_ncol = NULL,
    plot = TRUE,
    type = c("overview", "calibration", "residuals", "tails", "all"),
    empirical = c("bins", "points"),
    n_bins = 10,
    ncol = 2,
    randomized_residuals = TRUE,
    randomized_seed = 1L,
    randomized_draws = 5L,
    ...) {
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || anyNA(names(dots)) || any(!nzchar(names(dots))))) {
    stop("All arguments supplied through `...` must be named.", call. = FALSE)
  }
  if (anyDuplicated(names(dots))) {
    stop("Arguments supplied through `...` must not be duplicated.", call. = FALSE)
  }
  unknown <- names(dots)
  if (length(unknown)) {
    stop("Unused argument", if (length(unknown) > 1L) "s" else "", ": ",
         paste0("`", unknown, "`", collapse = ", "), ".", call. = FALSE)
  }
  if (!is.null(plot1_style)) {
    if (!missing(empirical)) stop("Supply only `empirical`; do not also supply deprecated `plot1_style`.", call. = FALSE)
    .gl_warn_interface_deprecation(
      "plot1_style_argument",
      "`plot1_style` is deprecated; use `empirical = \"bins\"` or `empirical = \"points\"`.",
      "gamlss.longitudinal_deprecated_plot_argument"
    )
    empirical <- if (identical(plot1_style, "scatter")) "points" else plot1_style
    plot1_style <- NULL
  }
  if (!is.null(plot2_cuts)) {
    if (!missing(n_bins)) stop("Supply only `n_bins`; do not also supply deprecated `plot2_cuts`.", call. = FALSE)
    .gl_warn_interface_deprecation(
      "plot2_cuts_argument", "`plot2_cuts` is deprecated; use `n_bins`.",
      "gamlss.longitudinal_deprecated_plot_argument"
    )
    n_bins <- plot2_cuts
    plot2_cuts <- NULL
  }
  if (!is.null(dashboard_ncol)) {
    if (!missing(ncol)) stop("Supply only `ncol`; do not also supply deprecated `dashboard_ncol`.", call. = FALSE)
    .gl_warn_interface_deprecation(
      "dashboard_ncol_argument", "`dashboard_ncol` is deprecated; use `ncol`.",
      "gamlss.longitudinal_deprecated_plot_argument"
    )
    ncol <- dashboard_ncol
    dashboard_ncol <- NULL
  }
  plot.copula(
    x = x,
    type = type,
    empirical = empirical,
    n_bins = n_bins,
    ncol = ncol,
    lags = lags,
    grid_n = grid_n,
    max_pairs_overlay = max_pairs_overlay,
    transform = transform,
    contour_bins = contour_bins,
    time_stratified = time_stratified,
    by = by,
    data = data,
    tau_ylim = tau_ylim,
    tail_thresholds = tail_thresholds,
    residual_lags = residual_lags,
    randomized_residuals = randomized_residuals,
    randomized_seed = randomized_seed,
    randomized_draws = randomized_draws,
    plot = plot
  )
}
