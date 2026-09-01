#' Plot all smooth terms with confidence bands
#'
#' This utility plots every smooth term in a fitted `gamlss.longitudinal` object
#' and computes pointwise confidence bands using the smooth coefficient
#' covariance matrices returned by `vcov.gamlss.longitudinal()`.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param vcov_obj Optional output from `vcov(object, ...)`. If `NULL`, this is
#' computed internally with the analytical vcov path.
#' @param data Optional data frame containing original covariates used for the
#' x-axis variable of each smooth.
#' @param ci_level Confidence level for pointwise intervals.
#' @param ncol Number of plot columns (defaults to 2 or fewer if needed).
#' @param ci_col Color for confidence bands.
#' @param fit_col Color for fitted smooth line.
#' @param ci_lty Line type for confidence bands.
#' @param fit_lwd Line width for fitted smooth.
#' @param sort_x Logical; sort points by x before plotting lines.
#' @param even_grid Logical; if TRUE, plot smooths on an evenly spaced x-grid
#' built over observed x-range.
#' @param grid_n Number of grid points when `even_grid = TRUE`.
#' @param fallback_to_index Logical; if x variable cannot be inferred, plot
#' against row index.
#' @param setup_mfrow Logical; if TRUE (default), configure par(mfrow) inside
#' this function. Set FALSE when caller configures layout.
#' @param show_legend Logical; if TRUE, draw a small legend in each panel.
#'
#' @return Invisibly returns a nested list with x, fitted values, standard
#' errors, and approximate conditional pointwise confidence limits for each
#' smooth term. The bands omit fixed-smooth, between-smooth, and
#' smoothing-parameter uncertainty; the returned object records this contract.
#' @export
plot_smooth_terms <- function(
    object,
    vcov_obj = NULL,
    data = NULL,
    ci_level = 0.95,
    ncol = NULL,
    ci_col = "red",
    fit_col = "black",
    ci_lty = 2,
    fit_lwd = 2,
    sort_x = TRUE,
    even_grid = TRUE,
    grid_n = 200,
    fallback_to_index = TRUE,
    setup_mfrow = TRUE,
    show_legend = TRUE) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be of class 'gamlss.longitudinal'.")
  }

  if (is.null(vcov_obj)) {
    vcov_obj <- .resolve_vcov(object, numderiv = FALSE, extra_args = list(method = "analytical"))
  }

  if (!is.list(vcov_obj) || is.null(vcov_obj$vcov)) {
    stop("'vcov_obj' must be the output of vcov.gamlss.longitudinal().")
  }

  smooth_vcov_list <- vcov_obj$vcov$smooth_vcov

  smooth_se_list <- vcov_obj$vcov$smooth_se

  z <- qnorm((1 + ci_level) / 2)

  smooth_index <- .plot_smooth_terms_index(object)

  n_plots <- length(smooth_index)

  if (n_plots == 0) {
    warning("No smooth terms found to plot.")

    return(invisible(.gl_attach_inference_contract(
      list(),
      .gl_inference_contract(
        "smooth_term_pointwise", coefficient_names = character(),
        validity_status = "not_applicable"
      )
    )))
  }

  out <- list()

  plot_objects <- list()

  for (i in seq_len(n_plots)) {
    par_name <- smooth_index[[i]]$par_name

    s_name <- smooth_index[[i]]$s_name

    B <- object$model_matrix$s[[par_name]][[s_name]]

    beta_s <- object$par_s[[par_name]][[s_name]]

    x_info <- .plot_smooth_terms_x_info(
      par_name = par_name,
      s_name = s_name,
      B = B,
      object = object,
      data = data,
      fallback_to_index = fallback_to_index
    )

    x <- x_info$x

    fitted_smooth <- as.numeric(B %*% beta_s)

    smooth_vcov <- NULL

    smooth_se <- NULL

    if (!is.null(smooth_vcov_list) && !is.null(smooth_vcov_list[[par_name]])) {
      smooth_vcov <- smooth_vcov_list[[par_name]][[s_name]]
    }

    if (!is.null(smooth_se_list) && !is.null(smooth_se_list[[par_name]])) {
      smooth_se <- smooth_se_list[[par_name]][[s_name]]
    }

    smooth_fit_se <- .plot_smooth_terms_fit_se(B, smooth_vcov, smooth_se)

    ci_lower <- fitted_smooth - z * smooth_fit_se

    ci_upper <- fitted_smooth + z * smooth_fit_se

    main_title <- paste(par_name, s_name, sep = ": ")

    ylab_text <- paste("smooth(", x_info$x_var, ")", sep = "")

    plot_df <- .plot_smooth_terms_plot_df(
      x = x,
      fitted_smooth = fitted_smooth,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      sort_x = sort_x,
      even_grid = even_grid,
      grid_n = grid_n
    )

    y_lim <- .plot_fixed_terms_y_limits(plot_df$fitted, plot_df$ci_lower, plot_df$ci_upper)

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = fitted)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), fill = ci_col, alpha = 0.16) +
      ggplot2::geom_line(color = fit_col, linewidth = fit_lwd) +
      ggplot2::labs(title = main_title, x = x_info$x_var, y = ylab_text)

    if (!is.null(y_lim)) {
      p <- p + ggplot2::coord_cartesian(ylim = y_lim)
    }

    if (show_legend) {
      p <- p + ggplot2::labs(caption = paste("fit /", round(ci_level * 100), "% CI"))
    }

    p <- p + ggplot2::theme_minimal()

    plot_objects[[length(plot_objects) + 1]] <- p

    if (is.null(out[[par_name]])) out[[par_name]] <- list()

    out[[par_name]][[s_name]] <- list(
      x = x,
      fitted = fitted_smooth,
      se = smooth_fit_se,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      plot = p
    )
  }

  if (length(plot_objects) > 0) {
    if (is.null(ncol)) {
      ncol <- min(2, n_plots)
    }

    nrow <- ceiling(length(plot_objects) / ncol)

    dashboard <- list(plotlist = plot_objects, ncol = ncol, nrow = nrow)

    if (setup_mfrow) {
      grid::grid.newpage()

      grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol)))

      for (i_plot in seq_along(plot_objects)) {
        r <- ((i_plot - 1) %/% ncol) + 1

        c <- ((i_plot - 1) %% ncol) + 1

        print(plot_objects[[i_plot]], vp = grid::viewport(layout.pos.row = r, layout.pos.col = c))
      }

      grid::popViewport()
    }

    out$plots <- plot_objects

    out$dashboard <- dashboard
  }

  invisible(.gl_attach_inference_contract(
    out,
    .gl_inference_contract(
      "smooth_term_pointwise",
      coefficient_names = .gl_smooth_coefficient_names(object),
      validity_status = "approximate"
    )
  ))
}
