.copula_v2_dashboard_plot_list <- function(empirical_overlay,
                                           quartile_correlation,
                                           rosenblatt_by_time,
                                           rosenblatt_qq,
                                           rosenblatt_lag,
                                           kendall_function,
                                           tail_cooccurrence,
                                           conditional_tail_exceedance,
                                           residual_lag_correlation) {
  list(
    empirical_overlay = empirical_overlay,
    quartile_correlation = quartile_correlation,
    kendall_function = kendall_function,
    rosenblatt_by_time = rosenblatt_by_time,
    rosenblatt_qq = rosenblatt_qq,
    rosenblatt_lag = rosenblatt_lag,
    tail_cooccurrence = tail_cooccurrence,
    conditional_tail_exceedance = conditional_tail_exceedance,
    residual_lag_correlation = residual_lag_correlation
  )
}

.copula_v2_arrange_dashboard <- function(dashboard_plots, dashboard_ncol) {
  do.call(
    ggpubr::ggarrange,
    c(
      dashboard_plots,
      list(
        ncol = min(dashboard_ncol, length(dashboard_plots)),
        nrow = ceiling(length(dashboard_plots) / dashboard_ncol)
      )
    )
  )
}

.copula_v2_plot_result <- function(plots,
                                   dashboard,
                                   fit_data,
                                   pair_data,
                                   pair_data_uniform,
                                   rosenblatt,
                                   rosenblatt_pairs,
                                   quartile_summary,
                                   kendall_summary,
                                   tail_summary,
                                   conditional_tail_summary,
                                   residual_lag_summary) {
  list(
    plots = plots,
    dashboard = dashboard,
    fit_data = fit_data,
    pair_data = pair_data,
    pair_data_uniform = pair_data_uniform,
    rosenblatt = rosenblatt,
    rosenblatt_pairs = rosenblatt_pairs,
    quartile_summary = quartile_summary,
    kendall_summary = kendall_summary,
    tail_summary = tail_summary,
    conditional_tail_summary = conditional_tail_summary,
    residual_lag_summary = residual_lag_summary
  )
}
