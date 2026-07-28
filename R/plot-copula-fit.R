#' Plot empirical pseudo-observation pairs with a fitted copula overlay
#'
#' `plot_copula_fit()` and `plot_copula_overlay()` standardise the
#' copula-screen visual check used in examples and vignettes. They accept raw
#' pseudo-observations plus a `select_copula()` result, or a final
#' `gamlss.longitudinal` fit whose row-specific fitted copula density is
#' averaged over paired observations.
#'
#' @param data Optional long-format data frame. If a fitted
#'   `gamlss.longitudinal` object is supplied here, it is treated as `fit`.
#' @param copula_dist A `copula_selection` result, one-row selection data frame,
#'   family code(s), or `"best"`/`"auto"` to screen all supported families.
#'   Character values are fitted to the supplied pseudo-observation pairs before
#'   the overlay is drawn. Required unless `fit` is supplied.
#' @param fit Optional fitted `gamlss.longitudinal` object.
#' @param object Optional alias for `fit`.
#' @param u1,u2 Optional direct paired pseudo-observations.
#' @param u Optional row-aligned pseudo-observation vector for `data`.
#' @param u_var Optional pseudo-observation column name in `data`.
#' @param response_var Optional response column name used to create temporary
#'   pseudo-observations when `u`, `u_var`, and `u1`/`u2` are absent.
#' @param margin_dist Optional margin family used when creating temporary
#'   pseudo-observations from `response_var`.
#' @param subject_var,time_var Subject and time column names.
#' @param lags Positive integer lag(s) used when forming pairs.
#' @param by_time Logical; if `TRUE`, facet the copula overlay by adjacent
#'   time pair.
#' @param transform Character; either `"normal"` or `"uniform"`.
#' @param grid_n Grid size for fitted density contours.
#' @param bins Number of empirical two-dimensional bins.
#' @param contour_bins Number of fitted contour levels.
#' @param max_pairs_overlay Maximum paired observations used when averaging
#'   fitted copula densities.
#' @param plot Logical; if TRUE, print the plot.
#' @param ... Additional arguments reserved for future methods.
#'
#' @return Invisibly returns a list containing the plot, pair data, and fitted
#'   density grid.
#' @export
plot_copula_fit <- function(
    data = NULL,
    copula_dist = NULL,
    fit = NULL,
    object = NULL,
    u1 = NULL,
    u2 = NULL,
    u = NULL,
    u_var = NULL,
    response_var = NULL,
    margin_dist = NULL,
    subject_var = "subject",
    time_var = "time",
    lags = 1,
    by_time = FALSE,
    transform = c("normal", "uniform"),
    grid_n = 80,
    bins = 28,
    contour_bins = 8,
    max_pairs_overlay = 300,
    plot = TRUE,
    ...) {
  .plot_reject_old_call_args(sys.call(), old_args = c("copula", "selected_fit"))
  .plot_reject_old_args(list(...), old_args = c("copula", "selected_fit"))
  transform <- match.arg(transform)
  if (!is.logical(by_time) || length(by_time) != 1L || is.na(by_time)) {
    stop("'by_time' must be TRUE or FALSE.", call. = FALSE)
  }

  plot_data <- .plot_copula_fit_inputs(
    data = data,
    copula_dist = copula_dist,
    fit = fit,
    object = object,
    u1 = u1,
    u2 = u2,
    u = u,
    u_var = u_var,
    response_var = response_var,
    margin_dist = margin_dist,
    subject_var = subject_var,
    time_var = time_var,
    lags = lags,
    by_time = by_time
  )

  pair_data <- plot_data$pair_data
  spec <- plot_data$spec
  selection <- plot_data$selection
  density_grid <- .plot_copula_fit_density_grid(
    pair_data = pair_data,
    spec = spec,
    by_time = by_time,
    grid_n = grid_n,
    max_pairs_overlay = max_pairs_overlay,
    refit_by_group = is.null(plot_data$fit)
  )
  density_grid <- .plot_copula_transform_grid(density_grid, transform)
  pair_plot <- .copula_v2_transform_data(pair_data, transform = transform)

  p <- .plot_copula_fit_plot(
    pair_plot = pair_plot,
    density_grid = density_grid,
    spec = spec,
    transform = transform,
    bins = bins,
    contour_bins = contour_bins,
    by_time = by_time
  )

  if (isTRUE(plot)) {
    print(p)
  }

  invisible(list(
    plot = p,
    pair_data = pair_data,
    density = density_grid,
    copula = spec,
    selection = selection
  ))
}
