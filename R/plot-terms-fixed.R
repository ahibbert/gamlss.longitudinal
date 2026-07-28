#' Plot all fixed terms with confidence bands

#'

#' This utility plots fixed-effect term contributions for a fitted

#' `gamlss.longitudinal` object using coefficient uncertainty from

#' `vcov.gamlss.longitudinal()`.

#'

#' For each fixed-effect design-matrix column \eqn{x_j}, it plots

#' \eqn{x_j \hat{\beta}_j} with pointwise confidence bands

#' \eqn{x_j \hat{\beta}_j \pm z_{\alpha/2}\sqrt{x_j^2 \mathrm{Var}(\hat{\beta}_j)}}.

#'

#' @param object A fitted `gamlss.longitudinal` object.

#' @param vcov_obj Optional output from `vcov(object, ...)`. If `NULL`, this is

#' computed internally with the analytical vcov path.

#' @param ci_level Confidence level for pointwise intervals.

#' @param ncol Number of plot columns (defaults to 2 or fewer if needed).

#' @param include_intercept Logical; include intercept columns in plots.

#' @param plot_interactions Logical; include interaction columns in plots.

#' @param ci_col Color for confidence bands.

#' @param fit_col Color for fitted fixed-term line.

#' @param ci_lty Line type for confidence bands.

#' @param fit_lwd Line width for fitted fixed-term line.

#' @param sort_x Logical; sort x-values before drawing lines.

#' @param fallback_to_index Logical; if x has one unique value, use index on x-axis.

#' @param setup_mfrow Logical; if TRUE (default), configure par(mfrow) inside

#' this function. Set FALSE when caller configures layout.

#' @param data Optional data frame used to detect factor columns and show

#' factor levels on x-axis for categorical fixed terms. Factor terms are grouped

#' into one point-and-interval plot per model term and parameter.

#' @param factor_pch Point symbol for factor-level estimates.

#' @param factor_cex Point size for factor-level estimates.

#' @param show_legend Logical; if TRUE, draw a small legend in each panel.

#'

#' @return Invisibly returns a nested list with x, fitted values, standard

#' errors, and confidence limits for each fixed term.

#' @export

plot_fixed_terms <- function(
    object,
    vcov_obj = NULL,
    ci_level = 0.95,
    ncol = NULL,
    include_intercept = FALSE,
    plot_interactions = FALSE,
    ci_col = "red",
    fit_col = "black",
    ci_lty = 2,
    fit_lwd = 2,
    sort_x = TRUE,
    fallback_to_index = TRUE,
    setup_mfrow = TRUE,
    data = NULL,
    factor_pch = 16,
    factor_cex = 1.2,
    show_legend = TRUE) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be of class 'gamlss.longitudinal'.")
  }

  if (is.null(vcov_obj)) {
    vcov_obj <- .resolve_vcov(object, numderiv = FALSE, extra_args = list(method = "analytical"))
  }

  if (!is.list(vcov_obj) || is.null(vcov_obj$vcov) || is.null(vcov_obj$vcov$overall)) {
    stop("'vcov_obj' must be the output of vcov.gamlss.longitudinal() with vcov$overall present.")
  }

  V <- vcov_obj$vcov$overall

  if (is.null(rownames(V)) || is.null(colnames(V))) {
    stop("vcov$overall must have row and column names matching fixed coefficients.")
  }

  z <- qnorm((1 + ci_level) / 2)

  gg_add <- .plot_fixed_terms_gg_add

  data_for_terms <- data

  if ((is.null(data_for_terms) || !is.data.frame(data_for_terms)) && !is.null(object$dataset)) {
    data_for_terms <- object$dataset
  }

  plot_specs <- .plot_fixed_terms_plot_specs(
    object = object,
    V = V,
    data_for_terms = data_for_terms,
    include_intercept = include_intercept,
    plot_interactions = plot_interactions
  )

  n_plots <- length(plot_specs)

  if (n_plots == 0) {
    warning("No fixed terms found to plot with matching vcov entries.")

    return(invisible(list()))
  }

  render_result <- .plot_fixed_terms_render_specs(
    plot_specs = plot_specs,
    object = object,
    V = V,
    z = z,
    ci_col = ci_col,
    fit_col = fit_col,
    ci_level = ci_level,
    factor_cex = factor_cex,
    show_legend = show_legend,
    gg_add = gg_add,
    sort_x = sort_x,
    fallback_to_index = fallback_to_index,
    fit_lwd = fit_lwd
  )

  out <- render_result$out

  plot_objects <- render_result$plot_objects

  if (length(plot_objects) > 0) {
    dashboard <- .plot_fixed_terms_dashboard(plot_objects, ncol = ncol)

    if (setup_mfrow) {
      grid::grid.newpage()

      grid::pushViewport(grid::viewport(layout = grid::grid.layout(dashboard$nrow, dashboard$ncol)))

      for (i_plot in seq_along(plot_objects)) {
        r <- ((i_plot - 1) %/% dashboard$ncol) + 1

        c <- ((i_plot - 1) %% dashboard$ncol) + 1

        print(plot_objects[[i_plot]], vp = grid::viewport(layout.pos.row = r, layout.pos.col = c))
      }

      grid::popViewport()
    }

    out$plots <- plot_objects

    out$dashboard <- dashboard
  }

  invisible(out)
}
