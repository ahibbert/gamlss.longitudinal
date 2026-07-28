#' Plot marginal and pairwise distribution diagnostics
#'
#' @param dataset Optional long-format data frame. Required unless `fit`
#'   contains stored data.
#' @param margin_dist Optional `gamlss.dist` family object used for diagonal
#'   marginal overlays.
#' @param subject_var,time_var,response_var Column names identifying subjects,
#'   time points, and responses. Required for raw data and inferred from `fit`
#'   when possible.
#' @param offdiag_scale Character; show off-diagonal panels on response or
#'   pseudo-observation scale.
#' @param transform Character; transform pseudo-observation off-diagonal panels
#'   to normal-score or uniform scale.
#' @param show_cor_stats Logical; include Pearson and Kendall correlations in
#'   off-diagonal panel subtitles.
#' @param fit Optional fitted `gamlss.longitudinal` object used for final-model
#'   margin and copula overlays.
#' @param overlay Character; `"none"` preserves the historical plot,
#'   `"margin"` overlays fitted marginal densities on diagonal panels,
#'   `"copula"` overlays a selected copula on pseudo-observation off-diagonal
#'   panels, and `"model"` overlays the final fitted model using `fit`.
#' @param copula_dist Optional `copula_selection` result or one-row selection table
#'   used when `overlay = "copula"`.
#' @param grid_n Grid size used for density overlays.
#' @param contour_bins Number of copula contour levels.
#' @param ... Compatibility arguments passed from `plotDist()` to
#'   `plot_dist()`.
#'
#' @return A `ggpubr` arranged plot object.
#' @export
plot_dist <- function(
    dataset = NULL,
    margin_dist = NULL,
    subject_var = NULL,
    time_var = NULL,
    response_var = NULL,
    offdiag_scale = c("pseudo", "response"),
    transform = c("normal", "uniform"),
    show_cor_stats = TRUE,
    fit = NULL,
    overlay = NULL,
    copula_dist = NULL,
    grid_n = 80,
    contour_bins = 8,
    ...) {
  .plot_reject_old_call_args(sys.call(), old_args = c("dist", "family", "copula"))
  .plot_reject_old_args(list(...), old_args = c("dist", "family", "copula"))
  controls <- .plot_dist_resolve_controls(
    fit = fit,
    offdiag_scale = offdiag_scale,
    transform = transform,
    overlay = overlay,
    copula_dist = copula_dist,
    margin_dist = margin_dist
  )
  offdiag_scale <- controls$offdiag_scale
  transform <- controls$transform
  overlay <- controls$overlay

  plot_data <- .plot_dist_normalise_data(
    dataset = dataset,
    fit = fit,
    subject_var = subject_var,
    time_var = time_var,
    response_var = response_var
  )
  time_values <- .plot_dist_time_values(plot_data$time)
  num_margins <- length(time_values)
  if (num_margins < 1L) {
    stop("No time points are available to plot.", call. = FALSE)
  }
  if (!is.null(margin_dist)) {
    margin_dist <- .plot_margin_resolve_family(margin_dist)
  } else if (!is.null(fit)) {
    margin_dist <- fit$margin_dist
  }

  overlay_state <- .plot_dist_overlay_state(
    fit = fit,
    dataset = dataset,
    num_margins = num_margins,
    copula_dist = copula_dist
  )
  fit_diag_data <- overlay_state$fit_diag_data
  fit_pair_data <- overlay_state$fit_pair_data
  fit_copula_spec <- overlay_state$fit_copula_spec
  copula_spec <- overlay_state$copula_spec

  margin_slices <- .plot_dist_margin_slices(plot_data, time_values)
  margin_data <- margin_slices$margin_data
  margin_pseudo <- margin_slices$margin_pseudo

  ## plot.new()
  # par(mfrow=c(1,num_margins))

  # Historical base graphics margin inspection is now covered by plot_margin_fit().
  # invisible(readline(prompt="Press [enter] to continue"))

  plots <- .plot_dist_panel_grid(
    num_margins = num_margins,
    margin_data = margin_data,
    margin_pseudo = margin_pseudo,
    time_values = time_values,
    overlay = overlay,
    margin_dist = margin_dist,
    fit = fit,
    fit_diag_data = fit_diag_data,
    offdiag_scale = offdiag_scale,
    transform = transform,
    show_cor_stats = show_cor_stats,
    copula_spec = copula_spec,
    fit_pair_data = fit_pair_data,
    fit_copula_spec = fit_copula_spec,
    grid_n = grid_n,
    contour_bins = contour_bins
  )
  ggpubr::ggarrange(plotlist = plots, ncol = num_margins, nrow = num_margins)
}
