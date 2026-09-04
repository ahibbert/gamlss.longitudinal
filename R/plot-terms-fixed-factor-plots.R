.plot_fixed_terms_factor_plot <- function(term_data,
                                          fg,
                                          par_name,
                                          ci_col,
                                          fit_col,
                                          factor_cex,
                                          ci_level,
                                          show_legend,
                                          gg_add = .plot_fixed_terms_gg_add) {
  p <- ggplot2::ggplot(
    term_data$plot_df[term_data$plot_df$keep, , drop = FALSE],
    ggplot2::aes(x = x, y = fitted)
  )
  p <- gg_add(p, ggplot2::geom_hline(yintercept = 0, color = "grey70", linetype = 3), "geom_hline")
  p <- gg_add(p, ggplot2::geom_point(color = fit_col, size = factor_cex), "geom_point")
  p <- gg_add(p, ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), color = ci_col, width = 0.15), "geom_errorbar")
  p <- gg_add(p, ggplot2::scale_x_continuous(breaks = term_data$x_plot, labels = fg$levels), "scale_x_continuous")
  p <- gg_add(
    p,
    ggplot2::labs(
      title = paste(par_name, fg$var_name, sep = ": "),
      x = fg$var_name,
      y = paste("fixed contribution:", paste(par_name, fg$var_name, sep = "."))
    ),
    "labs"
  )
  p <- gg_add(p, theme_gamlss_longitudinal(), "theme_gamlss_longitudinal")

  if (!is.null(term_data$y_lim)) {
    p <- gg_add(p, ggplot2::coord_cartesian(ylim = term_data$y_lim), "coord_cartesian")
  }

  if (show_legend) {
    p <- gg_add(p, ggplot2::labs(caption = .plot_fixed_terms_ci_caption("estimate", ci_level)), "labs")
  }

  p
}

.plot_fixed_terms_factor_factor_plot <- function(term_data,
                                                 ig,
                                                 par_name,
                                                 factor_cex,
                                                 ci_level,
                                                 show_legend,
                                                 gg_add = .plot_fixed_terms_gg_add) {
  pg <- ig$panel_group
  og <- ig$other_group
  plot_df <- term_data$plot_df

  p <- ggplot2::ggplot(
    plot_df[plot_df$keep, , drop = FALSE],
    ggplot2::aes(x = x, y = fitted, color = group, group = group)
  )
  p <- gg_add(p, ggplot2::geom_hline(yintercept = 0, color = "grey70", linetype = 3), "geom_hline")
  p <- gg_add(p, ggplot2::geom_line(linewidth = 0.8), "geom_line")
  p <- gg_add(p, ggplot2::geom_point(size = factor_cex), "geom_point")
  p <- gg_add(p, ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), width = 0.15), "geom_errorbar")
  p <- gg_add(p, ggplot2::scale_x_continuous(breaks = term_data$x_plot, labels = term_data$x_labels), "scale_x_continuous")
  group_levels <- unique(as.character(plot_df$group[plot_df$keep]))
  p <- gg_add(
    p,
    ggplot2::scale_color_manual(
      name = og$var_name,
      values = stats::setNames(.gl_group_colours(length(group_levels)), group_levels)
    ),
    "scale_color_manual"
  )
  p <- gg_add(
    p,
    ggplot2::labs(
      title = paste(par_name, ig$interaction_name, sep = ": "),
      x = pg$var_name,
      y = paste("fixed contribution:", paste(par_name, ig$interaction_name, sep = "."))
    ),
    "labs"
  )
  p <- gg_add(p, theme_gamlss_longitudinal(), "theme_gamlss_longitudinal")

  if (!is.null(term_data$y_lim)) {
    p <- gg_add(p, ggplot2::coord_cartesian(ylim = term_data$y_lim), "coord_cartesian")
  }

  if (show_legend) {
    p <- gg_add(p, ggplot2::labs(caption = .plot_fixed_terms_ci_caption("estimate", ci_level)), "labs")
  }

  p
}

.plot_fixed_terms_interaction_factor_plot <- function(term_data,
                                                      par_name,
                                                      ci_col,
                                                      fit_col,
                                                      factor_cex,
                                                      ci_level,
                                                      show_legend,
                                                      gg_add = .plot_fixed_terms_gg_add) {
  p <- ggplot2::ggplot(term_data$plot_df, ggplot2::aes(x = x, y = fitted))
  p <- gg_add(p, ggplot2::geom_hline(yintercept = 0, color = "grey70", linetype = 3), "geom_hline")
  p <- gg_add(p, ggplot2::geom_point(color = fit_col, size = factor_cex), "geom_point")
  p <- gg_add(p, ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), color = ci_col, width = 0.15), "geom_errorbar")
  p <- gg_add(p, ggplot2::scale_x_continuous(breaks = term_data$x_plot, labels = term_data$x_labels), "scale_x_continuous")
  p <- gg_add(
    p,
    ggplot2::labs(
      title = paste(par_name, term_data$col_name, sep = ": "),
      x = paste(term_data$col_name, "(interaction level)"),
      y = paste("fixed contribution:", term_data$coef_name)
    ),
    "labs"
  )
  p <- gg_add(p, theme_gamlss_longitudinal(), "theme_gamlss_longitudinal")

  if (!is.null(term_data$y_lim)) {
    p <- gg_add(p, ggplot2::coord_cartesian(ylim = term_data$y_lim), "coord_cartesian")
  }

  if (show_legend) {
    p <- gg_add(p, ggplot2::labs(caption = .plot_fixed_terms_ci_caption("estimate", ci_level)), "labs")
  }

  p
}
