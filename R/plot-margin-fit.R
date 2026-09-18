#' Plot an observed marginal response with a fitted GAMLSS density overlay
#'
#' `plot_margin_fit()` replaces the common exploratory use of
#' `gamlss::histDist()` with a ggplot-based helper that also understands final
#' `gamlss.longitudinal` fits. With raw data it fits the supplied family as an
#' intercept-only marginal model. With a fitted longitudinal model it overlays
#' the average row-specific fitted marginal density.
#'
#' @param x Numeric response vector, data frame, or fitted
#'   `gamlss.longitudinal` object.
#' @param margin_dist A `gamlss.dist` family object, family name, or
#'   `margin_selection` result. Required for raw data unless `fit` is supplied.
#' @param data Optional data frame containing the response.
#' @param fit Optional fitted `gamlss.longitudinal` object. When supplied, the
#'   fitted marginal density is overlaid using the model distribution and
#'   row-specific fitted parameters.
#' @param response_var Response column name when `x` or `data` is a data frame.
#' @param bins Number of histogram bins for continuous margins and the maximum
#'   number of displayed bins for discrete margins. Discrete responses use
#'   unit-width bins while the observed integer support fits within this limit;
#'   wider supports are grouped into shared, equal-width integer ranges.
#' @param separate_zero Logical; for discrete margins, keep zero in its own
#'   unit-width bin when a wide response range is grouped. This is useful for
#'   displaying the excess mass at zero in zero-inflated distributions.
#' @param grid_n Number of grid points used for a fitted continuous density.
#'   Discrete probability masses are evaluated at every integer in the observed
#'   response range.
#' @param response_scale Plot the density on the original response scale or,
#'   for positive responses, on the log-response scale.
#' @param time_intercepts Logical; for raw data, fit time-specific intercepts
#'   for each marginal distribution parameter before drawing the overlay.
#' @param by_time Logical; facet the plot by time. For fitted longitudinal
#'   models, fitted densities are averaged within each time point.
#' @param by Optional name of a factor or other grouping column used to facet
#'   the plot. For fitted models, supply `data` when grouping by a column other
#'   than time. When both `by_time = TRUE` and `by` are supplied, `by` takes
#'   precedence.
#' @param time_var Time column name used when `time_intercepts = TRUE` or
#'   `by_time = TRUE` for raw data.
#' @param fit_control Control object passed to the internal `gamlss()` overlay
#'   fit for raw-data plots. Defaults to `gamlss::gamlss.control(n.cyc = 50)`.
#' @param plot Logical; if TRUE, print the plot.
#' @param ... Additional arguments reserved for future methods.
#'
#' @return Invisibly returns a list containing the plot, observed data, and
#'   fitted density grid.
#' @export
plot_margin_fit <- function(
    x = NULL,
    margin_dist = NULL,
    data = NULL,
    fit = NULL,
    response_var = "response",
    bins = 30,
    separate_zero = FALSE,
    grid_n = 200,
    response_scale = c("response", "log"),
    time_intercepts = FALSE,
    by_time = FALSE,
    by = NULL,
    time_var = "time",
    fit_control = gamlss::gamlss.control(n.cyc = 50),
    plot = TRUE,
    ...) {
  .plot_reject_old_call_args(sys.call(), old_args = c("family"))
  .plot_reject_old_args(list(...), old_args = c("family"))
  response_scale <- match.arg(response_scale)
  if (!is.logical(separate_zero) || length(separate_zero) != 1L || is.na(separate_zero)) {
    stop("'separate_zero' must be TRUE or FALSE.", call. = FALSE)
  }

  if (inherits(x, "gamlss.longitudinal") && is.null(fit)) {
    fit <- x
    x <- NULL
  }

  plot_data <- if (!is.null(fit)) {
    if (!inherits(fit, "gamlss.longitudinal")) {
      stop("'fit' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
    }
    .plot_margin_fit_from_object(
      fit = fit,
      x = x,
      data = data,
      grid_n = grid_n,
      by_time = by_time,
      by = by,
      response_scale = response_scale
    )
  } else {
    .plot_margin_fit_from_raw(
      x = x,
      margin_dist = margin_dist,
      data = data,
      response_var = response_var,
      grid_n = grid_n,
      response_scale = response_scale,
      time_intercepts = time_intercepts,
      by_time = by_time,
      by = by,
      time_var = time_var,
      fit_control = fit_control
    )
  }

  p <- .plot_margin_fit_plot(
    obs = plot_data$obs,
    density_grid = plot_data$density_grid,
    margin_dist = plot_data$margin_dist,
    bins = bins,
    separate_zero = separate_zero,
    by_time = plot_data$split_by
  )

  if (isTRUE(plot)) {
    print(p)
  }

  invisible(.gl_new_plot_result(
    list(
      plot = p,
      data = plot_data$obs,
      density = plot_data$density_grid,
      by = plot_data$by
    ),
    primary = "plot"
  ))
}
