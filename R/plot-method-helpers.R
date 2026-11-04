.plot_method_quantile_col_name <- function(prob) {
  paste0("q", gsub("^0\\.", "", format(prob, trim = TRUE)))
}

.plot_method_empty_plot <- function(title_txt, msg_txt) {
  ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_text(label = msg_txt) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title_txt) +
    ggplot2::theme_void()
}

.plot_method_diagnostic_panels <- function(x, randomize, time_stratified) {
  list(
    pithist = pithist(x, bins = 20, randomize = randomize, plot = TRUE, by_time = time_stratified),
    qqrplot = qqrplot(x, randomize = randomize, plot = TRUE, by_time = time_stratified),
    wormplot = wormplot(x, randomize = randomize, plot = TRUE, by_time = time_stratified),
    rootogram = rootogram(x, bins = 20, plot = TRUE, by_time = time_stratified)
  )
}

.plot_method_quantile_columns <- function(quantiles) {
  list(
    q_low = .plot_method_quantile_col_name(min(quantiles)),
    q_high = .plot_method_quantile_col_name(max(quantiles)),
    q_mid = .plot_method_quantile_col_name(if (0.5 %in% quantiles) 0.5 else quantiles[ceiling(length(quantiles) / 2)])
  )
}

.plot_method_arrange_dashboard <- function(dashboard_plots) {
  ggpubr::ggarrange(
    plotlist = dashboard_plots,
    ncol = 2,
    nrow = ceiling(length(dashboard_plots) / 2)
  )
}
