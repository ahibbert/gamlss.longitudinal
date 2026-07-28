.plot_fixed_terms_continuous_plot <- function(term_data,
                                              par_name,
                                              ci_col,
                                              fit_col,
                                              fit_lwd,
                                              ci_level,
                                              show_legend,
                                              gg_add = .plot_fixed_terms_gg_add) {
  main_title <- paste(par_name, term_data$col_name, sep = ": ")
  ylab_text <- paste("fixed contribution:", term_data$coef_name)

  p <- ggplot2::ggplot(term_data$plot_df, ggplot2::aes(x = x, y = fitted))
  p <- gg_add(p, ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), fill = ci_col, alpha = 0.16), "geom_ribbon")
  p <- gg_add(p, ggplot2::geom_line(color = fit_col, linewidth = fit_lwd), "geom_line")
  p <- gg_add(p, ggplot2::labs(title = main_title, x = term_data$xlab_text, y = ylab_text), "labs")
  p <- gg_add(p, ggplot2::theme_minimal(), "theme_minimal")

  if (!is.null(term_data$y_lim)) {
    p <- gg_add(p, ggplot2::coord_cartesian(ylim = term_data$y_lim), "coord_cartesian")
  }

  if (show_legend) {
    p <- gg_add(p, ggplot2::labs(caption = .plot_fixed_terms_ci_caption("fit", ci_level)), "labs")
  }

  p
}
