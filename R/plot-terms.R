#' Plot term effects for a fitted longitudinal model

#'

#' @param x A fitted `gamlss.longitudinal` object.

#' @param y Unused; included for compatibility with older calls.

#' @param data Optional data frame used to recover factor levels, transformed

#'   covariate scales, and interaction plotting metadata. Factor fixed terms

#'   are grouped into one point-and-interval plot per model term and parameter.

#' @param ci_level Confidence level for pointwise intervals.

#' @param ncol Number of columns in the combined dashboard.

#' @param include_intercept Logical; include intercept terms in fixed-effect

#'   plots.

#' @param plot_interactions Logical; include fixed-effect interaction terms.

#' @param ci_col,fit_col Colours for interval and fitted-term layers.

#' @param show_legend Logical; include plot captions/legends where available.

#' @param smooth_even_grid Logical; draw smooth terms on an evenly spaced grid.

#' @param smooth_grid_n Number of grid points for smooth-term plots.

#' @param paginate Logical; print one chart at a time for large dashboards.

#' @param ... Additional arguments reserved for future use.

#'

#' @return Invisibly returns a list with smooth-term, fixed-term, and dashboard

#'   plot objects.

#' @export

plot_terms <- function(
    x,
    y = NULL,
    data = NULL,
    ci_level = 0.95,
    ncol = 4,
    include_intercept = FALSE,
    plot_interactions = FALSE,
    ci_col = "red",
    fit_col = "black",
    show_legend = TRUE,
    smooth_even_grid = TRUE,
    smooth_grid_n = 200,
    paginate = FALSE,
    ...) {
  .plot_terms_gamlss_longitudinal(
    x = x,
    y = y,
    data = data,
    ci_level = ci_level,
    ncol = ncol,
    include_intercept = include_intercept,
    plot_interactions = plot_interactions,
    ci_col = ci_col,
    fit_col = fit_col,
    show_legend = show_legend,
    smooth_even_grid = smooth_even_grid,
    smooth_grid_n = smooth_grid_n,
    paginate = paginate,
    ...
  )
}

#' @rdname plot_terms

#' @usage NULL

#' @rawNamespace export(plot.terms)

plot.terms <- function(x, ...) {
  .Deprecated("plot_terms", package = "gamlss.longitudinal")

  plot_terms(x, ...)
}

.plot_terms_gamlss_longitudinal <- function(
    x,
    y = NULL,
    data = NULL,
    ci_level = 0.95,
    ncol = 4,
    include_intercept = FALSE,
    plot_interactions = FALSE,
    ci_col = "red",
    fit_col = "black",
    show_legend = TRUE,
    smooth_even_grid = TRUE,
    smooth_grid_n = 200,
    paginate = FALSE,
    ...) {
  if (!inherits(x, "gamlss.longitudinal")) {
    stop("'x' must be a fitted 'gamlss.longitudinal' object.")
  }

  cat("\n=== Plotting term effects for gamlss.longitudinal object ===\n")

  counts <- .plot_terms_count(
    x,
    data = data,
    include_intercept = include_intercept,
    plot_interactions = plot_interactions
  )

  cat(sprintf(
    "Found %d smooth terms and %d fixed terms (total: %d plots).\n\n",
    counts$smooth, counts$fixed, counts$total
  ))

  if (counts$total == 0) {
    warning("No term plots to display.")

    return(invisible(list(smooth_terms = list(), fixed_terms = list())))
  }

  vcov_obj <- .resolve_vcov(x, numderiv = FALSE, extra_args = list(method = "analytical"))

  smooth_results <- list()

  fixed_results <- list()

  if (counts$smooth > 0) {
    smooth_results <- plot_smooth_terms(
      object = x,
      vcov_obj = vcov_obj,
      data = data,
      ci_level = ci_level,
      ncol = ncol,
      ci_col = ci_col,
      fit_col = fit_col,
      even_grid = smooth_even_grid,
      grid_n = smooth_grid_n,
      setup_mfrow = FALSE,
      show_legend = show_legend
    )
  }

  if (counts$fixed > 0) {
    fixed_results <- plot_fixed_terms(
      object = x,
      vcov_obj = vcov_obj,
      ci_level = ci_level,
      ncol = ncol,
      include_intercept = include_intercept,
      plot_interactions = plot_interactions,
      ci_col = ci_col,
      fit_col = fit_col,
      setup_mfrow = FALSE,
      data = data,
      show_legend = show_legend
    )
  }

  plot_objects <- .plot_terms_collect_plot_objects(smooth_results, fixed_results)

  dashboard <- .plot_terms_render_dashboard(plot_objects, ncol = ncol, paginate = paginate)

  out <- list(
    smooth_terms = smooth_results,
    fixed_terms = fixed_results,
    dashboard = dashboard
  )
  attr(out, "inference_contract") <- list(
    fixed = attr(fixed_results, "inference_contract"),
    smooth = attr(smooth_results, "inference_contract")
  )
  invisible(out)
}
