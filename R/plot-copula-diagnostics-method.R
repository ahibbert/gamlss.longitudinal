#' @rdname plot_copula_diagnostics

#' @method plot copula

#' @export

plot.copula <- function(x, lags = 1, grid_n = 35, max_pairs_overlay = 300, transform = "normal", plot1_style = "bins", contour_bins = 8, time_stratified = FALSE, by = NULL, data = NULL, tau_ylim = NULL, plot2_cuts = 10, tail_thresholds = c(0.05, 0.10, 0.20), residual_lags = 1:3, dashboard_ncol = 3, plot = TRUE, ...) {
  if (!inherits(x, "gamlss.longitudinal")) {
    stop("'x' must be a fitted 'gamlss.longitudinal' object.")
  }

  object <- x

  controls <- .copula_v2_normalize_plot_controls(
    transform = transform,
    plot1_style = plot1_style,
    contour_bins = contour_bins,
    time_stratified = time_stratified,
    plot2_cuts = plot2_cuts,
    tau_ylim = tau_ylim,
    tail_thresholds = tail_thresholds,
    residual_lags = residual_lags,
    dashboard_ncol = dashboard_ncol
  )
  transform <- controls$transform
  plot1_style <- controls$plot1_style
  contour_bins <- controls$contour_bins
  time_stratified <- controls$time_stratified
  plot2_cuts <- controls$plot2_cuts
  tau_ylim <- controls$tau_ylim
  tail_thresholds <- controls$tail_thresholds
  residual_lags <- controls$residual_lags
  dashboard_ncol <- controls$dashboard_ncol

  fit_data <- .copula_v2_fit_data(object)

  pair_data_uniform <- .copula_v2_pair_data(fit_data, lags = lags)

  if (isTRUE(time_stratified) && is.null(by)) {
    by <- "time_pair"
  } else if (isTRUE(time_stratified) && !is.null(by)) {
    warning("Both time_stratified and by were supplied; using by='", by, "'.", call. = FALSE)
  }

  pair_data_uniform <- .copula_v2_attach_group(pair_data_uniform, object = object, by = by, data = data)

  pair_data_plot <- pair_data_uniform

  # Apply transform to pair data if requested

  if (transform == "normal") {
    pair_data_plot <- .copula_v2_transform_data(pair_data_plot, transform = "normal")
  }

  copula_spec <- get_copula_dist(object$copula_dist)

  copula_family_name <- .copula_family_code(copula_spec$copula_dist)

  family_num <- tryCatch(
    {
      .copula_family_code(copula_family_name)
    },
    error = function(e) NA_character_
  )

  is_grouped <- !is.null(by) || isTRUE(time_stratified)

  density_grid <- .copula_v2_density_grid_for_plot(
    family_num = family_num,
    pair_data_uniform = pair_data_uniform,
    grid_n = grid_n,
    max_pairs_overlay = max_pairs_overlay,
    transform = transform,
    is_grouped = is_grouped
  )

  axis_labels <- .copula_v2_axis_labels(transform)
  x_label <- axis_labels$x
  y_label <- axis_labels$y

  p1 <- .copula_v2_empirical_overlay_panel(
    pair_data_plot = pair_data_plot,
    density_grid = density_grid,
    plot1_style = plot1_style,
    contour_bins = contour_bins,
    is_grouped = is_grouped,
    copula_family_name = copula_family_name,
    x_label = x_label,
    y_label = y_label
  )

  quartile_df <- .copula_v2_cut_summary_by_group(
    pair_data_plot = pair_data_plot,
    plot2_cuts = plot2_cuts,
    is_grouped = is_grouped
  )

  p2 <- .copula_v2_cut_summary_panel(
    quartile_df = quartile_df,
    plot2_cuts = plot2_cuts,
    is_grouped = is_grouped,
    tau_ylim = tau_ylim
  )

  rosenblatt_df <- tryCatch(

    .copula_v2_rosenblatt_series(fit_data, family_num),
    error = function(e) data.frame()
  )

  rosenblatt_pair_df <- tryCatch(

    .copula_v2_rosenblatt_pair_data(pair_data_uniform, family_num),
    error = function(e) data.frame()
  )

  p_ros_time <- .copula_v2_rosenblatt_time_panel(rosenblatt_df)

  rosenblatt_qq_df <- .copula_v2_rosenblatt_qq_data(rosenblatt_df)

  p_ros_qq <- .copula_v2_rosenblatt_qq_panel(rosenblatt_qq_df)

  p_ros_lag <- .copula_v2_rosenblatt_lag_panel(rosenblatt_pair_df)

  kendall_df <- tryCatch(

    .copula_v2_kendall_diagnostic(pair_data_uniform, family_num),
    error = function(e) data.frame()
  )

  p_kendall <- .copula_v2_kendall_panel(kendall_df)

  tail_df <- tryCatch(

    .copula_v2_tail_diagnostics(pair_data_uniform, family_num, thresholds = tail_thresholds),
    error = function(e) data.frame()
  )

  cond_tail_df <- .copula_v2_conditional_tail_diagnostics(tail_df)

  tail_long <- .copula_v2_tail_long_data(tail_df)

  p_tail <- .copula_v2_tail_panel(
    tail_long = tail_long,
    title = "Tail Co-occurrence",
    message_subtitle = "Observed joint tail probability against fitted copula probability",
    empty_message = "No finite tail diagnostics",
    plot_subtitle = "Lower: P(Ut <= a, Ut+1 <= a); Upper: P(Ut >= 1-a, Ut+1 >= 1-a)",
    y_label = "Joint probability"
  )

  cond_tail_long <- .copula_v2_tail_long_data(cond_tail_df)

  p_cond_tail <- .copula_v2_tail_panel(
    tail_long = cond_tail_long,
    title = "Conditional Tail Exceedance",
    message_subtitle = "Observed conditional tail probability against fitted copula probability",
    empty_message = "No finite conditional tail diagnostics",
    plot_subtitle = "Lower: P(Ut+1 <= a | Ut <= a); Upper: P(Ut+1 >= 1-a | Ut >= 1-a)",
    y_label = "Conditional probability",
    ylim = c(0, 1)
  )

  lag_summary_df <- tryCatch(

    .copula_v2_rosenblatt_lag_summary(rosenblatt_df, lag_values = residual_lags),
    error = function(e) data.frame()
  )

  p_lag_summary <- .copula_v2_residual_lag_panel(lag_summary_df)

  dashboard_plots <- .copula_v2_dashboard_plot_list(
    empirical_overlay = p1,
    quartile_correlation = p2,
    rosenblatt_by_time = p_ros_time,
    rosenblatt_qq = p_ros_qq,
    rosenblatt_lag = p_ros_lag,
    kendall_function = p_kendall,
    tail_cooccurrence = p_tail,
    conditional_tail_exceedance = p_cond_tail,
    residual_lag_correlation = p_lag_summary
  )

  dashboard <- .copula_v2_arrange_dashboard(dashboard_plots, dashboard_ncol = dashboard_ncol)

  if (isTRUE(plot)) {
    print(dashboard)
  }

  invisible(.copula_v2_plot_result(
    plots = dashboard_plots,
    dashboard = dashboard,
    fit_data = fit_data,
    pair_data = pair_data_plot,
    pair_data_uniform = pair_data_uniform,
    rosenblatt = rosenblatt_df,
    rosenblatt_pairs = rosenblatt_pair_df,
    quartile_summary = quartile_df,
    kendall_summary = kendall_df,
    tail_summary = tail_df,
    conditional_tail_summary = cond_tail_df,
    residual_lag_summary = lag_summary_df
  ))
}
