#' @rdname plot_copula_fit
#' @export
plot_copula_overlay <- function(
    data = NULL,
    copula_dist = NULL,
    fit = NULL,
    object = NULL,
    u1 = NULL,
    u2 = NULL,
    u = NULL,
    u_var = NULL,
    response_var = NULL,
    margin_dist = NULL,
    subject_var = "subject",
    time_var = "time",
    lags = 1,
    by_time = FALSE,
    transform = c("normal", "uniform"),
    grid_n = 80,
    bins = 28,
    contour_bins = 8,
    max_pairs_overlay = 300,
    plot = TRUE,
    ...) {
  .plot_reject_old_call_args(sys.call(), old_args = c("copula", "selected_fit"))
  plot_copula_fit(
    data = data,
    copula_dist = copula_dist,
    fit = fit,
    object = object,
    u1 = u1,
    u2 = u2,
    u = u,
    u_var = u_var,
    response_var = response_var,
    margin_dist = margin_dist,
    subject_var = subject_var,
    time_var = time_var,
    lags = lags,
    by_time = by_time,
    transform = transform,
    grid_n = grid_n,
    bins = bins,
    contour_bins = contour_bins,
    max_pairs_overlay = max_pairs_overlay,
    plot = plot,
    ...
  )
}
