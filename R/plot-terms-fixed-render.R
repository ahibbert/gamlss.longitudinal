.plot_fixed_terms_render_specs <- function(
    plot_specs,
    object,
    V,
    z,
    ci_col,
    fit_col,
    ci_level,
    factor_cex,
    show_legend,
    gg_add,
    sort_x,
    fallback_to_index,
    fit_lwd) {
  out <- list()

  plot_objects <- list()

  for (i in seq_along(plot_specs)) {
    spec <- plot_specs[[i]]

    par_name <- spec$par_name

    if (identical(spec$type, "factor")) {
      fg <- spec$group

      levs <- fg$levels

      term_data <- .plot_fixed_terms_factor_data(fg, par_name, object, V, z)

      x_plot <- term_data$x_plot

      fitted_term <- term_data$fitted

      term_se <- term_data$se

      ci_lower <- term_data$ci_lower

      ci_upper <- term_data$ci_upper

      p <- .plot_fixed_terms_factor_plot(
        term_data = term_data,
        fg = fg,
        par_name = par_name,
        ci_col = ci_col,
        fit_col = fit_col,
        factor_cex = factor_cex,
        ci_level = ci_level,
        show_legend = show_legend,
        gg_add = gg_add
      )

      if (is.null(out[[par_name]])) out[[par_name]] <- list()

      out[[par_name]][[fg$var_name]] <- list(
        coefficient = paste(par_name, fg$var_name, sep = "."),
        x = x_plot,
        levels = levs,
        fitted = fitted_term,
        se = term_se,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        plot = p
      )

      plot_objects[[length(plot_objects) + 1]] <- p
    } else if (identical(spec$type, "interaction_factor_factor")) {
      ig <- spec$group

      term_data <- .plot_fixed_terms_factor_factor_data(ig, par_name, object, V, z)

      panel_levels <- term_data$panel_levels

      other_levels <- term_data$other_levels

      x_plot <- term_data$x_plot

      plot_df <- term_data$plot_df

      p <- .plot_fixed_terms_factor_factor_plot(
        term_data = term_data,
        ig = ig,
        par_name = par_name,
        factor_cex = factor_cex,
        ci_level = ci_level,
        show_legend = show_legend,
        gg_add = gg_add
      )

      if (is.null(out[[par_name]])) out[[par_name]] <- list()

      out[[par_name]][[ig$interaction_name]] <- list(
        coefficient = paste(par_name, ig$interaction_name, sep = "."),
        x = x_plot,
        levels = panel_levels,
        series = other_levels,
        fitted = plot_df$fitted,
        se = plot_df$se,
        ci_lower = plot_df$ci_lower,
        ci_upper = plot_df$ci_upper,
        plot_data = plot_df,
        plot = p
      )

      plot_objects[[length(plot_objects) + 1]] <- p
    } else if (identical(spec$type, "interaction_factor")) {
      term_data <- .plot_fixed_terms_interaction_factor_data(spec, par_name, object, V, z)

      if (is.null(term_data)) next

      col_name <- term_data$col_name

      coef_name <- term_data$coef_name

      x_levels <- term_data$x_levels

      x_labels <- term_data$x_labels

      x_plot <- term_data$x_plot

      fitted_term <- term_data$fitted

      term_se <- term_data$se

      ci_lower <- term_data$ci_lower

      ci_upper <- term_data$ci_upper

      p <- .plot_fixed_terms_interaction_factor_plot(
        term_data = term_data,
        par_name = par_name,
        ci_col = ci_col,
        fit_col = fit_col,
        factor_cex = factor_cex,
        ci_level = ci_level,
        show_legend = show_legend,
        gg_add = gg_add
      )

      if (is.null(out[[par_name]])) out[[par_name]] <- list()

      out[[par_name]][[col_name]] <- list(
        coefficient = coef_name,
        x = x_levels,
        levels = x_labels,
        fitted = fitted_term,
        se = term_se,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        plot = p
      )

      plot_objects[[length(plot_objects) + 1]] <- p
    } else {
      term_data <- .plot_fixed_terms_continuous_data(
        spec,
        par_name,
        object,
        V,
        z,
        sort_x = sort_x,
        fallback_to_index = fallback_to_index
      )

      col_name <- term_data$col_name

      coef_name <- term_data$coef_name

      x_plot <- term_data$x_plot

      fitted_term <- term_data$fitted

      term_se <- term_data$se

      ci_lower <- term_data$ci_lower

      ci_upper <- term_data$ci_upper

      p <- .plot_fixed_terms_continuous_plot(
        term_data = term_data,
        par_name = par_name,
        ci_col = ci_col,
        fit_col = fit_col,
        fit_lwd = fit_lwd,
        ci_level = ci_level,
        show_legend = show_legend,
        gg_add = gg_add
      )

      if (is.null(out[[par_name]])) out[[par_name]] <- list()

      out[[par_name]][[col_name]] <- list(
        coefficient = coef_name,
        x = x_plot,
        fitted = fitted_term,
        se = term_se,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        plot = p
      )

      plot_objects[[length(plot_objects) + 1]] <- p
    }
  }

  list(out = out, plot_objects = plot_objects)
}
