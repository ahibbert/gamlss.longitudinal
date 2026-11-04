.copula_contour_compare_controls <- function(transform, diff_scale_limit) {
  if (!transform %in% c("uniform", "normal")) {
    stop("'transform' must be either 'uniform' or 'normal'.")
  }

  if (!is.numeric(diff_scale_limit) || length(diff_scale_limit) != 1 || !is.finite(diff_scale_limit) || diff_scale_limit <= 0) {
    stop("'diff_scale_limit' must be a single positive numeric value.")
  }

  list(transform = transform, diff_scale_limit = diff_scale_limit)
}

.copula_contour_compare_metric_tables <- function(grid_df) {
  metric_list <- lapply(split(grid_df, grid_df$time_pair), function(g) {
    m <- .copula_v2_surface_metrics(g$density_emp, g$density_fit)

    out <- cbind(data.frame(time_pair = unique(g$time_pair), stringsAsFactors = FALSE), m$summary)

    if (nrow(m$overlap) > 0) {
      overlap <- cbind(data.frame(time_pair = unique(g$time_pair), stringsAsFactors = FALSE), m$overlap)
    } else {
      overlap <- data.frame()
    }

    list(summary = out, overlap = overlap)
  })

  list(
    summary = do.call(rbind, lapply(metric_list, function(x) x$summary)),
    overlap = do.call(rbind, lapply(metric_list, function(x) x$overlap))
  )
}

.copula_contour_compare_axis_labels <- function(transform) {
  list(
    x = if (transform == "normal") expression(Phi^-1 * (U[t])) else expression(U[t]),
    y = if (transform == "normal") expression(Phi^-1 * (U[t + 1])) else expression(U[t + 1])
  )
}
