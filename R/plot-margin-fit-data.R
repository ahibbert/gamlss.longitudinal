#' @keywords internal
#' @noRd
.plot_margin_fit_from_object <- function(fit, x = NULL, data = NULL, grid_n = 200, by_time = FALSE, response_scale = "response") {
  newdata <- if (!is.null(data)) {
    data
  } else if (is.data.frame(x)) {
    x
  } else {
    NULL
  }
  diag_data <- .gl_fitted_distribution(fit, newdata = newdata, require_response = TRUE)
  margin_dist <- fit$margin_dist
  obs <- data.frame(
    response = diag_data$response,
    split_group = as.character(diag_data$time),
    stringsAsFactors = FALSE
  )
  density_grid <- .plot_margin_density_grid(
    y = diag_data$response,
    family = margin_dist,
    params = diag_data$params,
    grid_n = grid_n,
    group = if (isTRUE(by_time)) diag_data$time else NULL
  )

  list(
    margin_dist = margin_dist,
    obs = .plot_margin_transform_observed(obs, response_scale),
    density_grid = .plot_margin_transform_density(density_grid, response_scale)
  )
}

#' @keywords internal
#' @noRd
.plot_margin_fit_from_raw <- function(
    x = NULL,
    margin_dist = NULL,
    data = NULL,
    response_var = "response",
    grid_n = 200,
    response_scale = "response",
    time_intercepts = FALSE,
    by_time = FALSE,
    time_var = "time",
    fit_control = gamlss::gamlss.control(n.cyc = 50)) {
  margin_dist <- .plot_margin_resolve_family(margin_dist)
  y <- .plot_margin_response(x = x, data = data, response_var = response_var)
  raw_data <- if (!is.null(data)) as.data.frame(data, stringsAsFactors = FALSE) else if (is.data.frame(x)) as.data.frame(x, stringsAsFactors = FALSE) else NULL
  raw_time <- NULL
  if (isTRUE(time_intercepts) || isTRUE(by_time)) {
    if (is.null(raw_data) || !time_var %in% names(raw_data)) {
      stop("Raw-data time overlays require 'data' or data-frame 'x' containing 'time_var'.", call. = FALSE)
    }
    raw_time <- raw_data[[time_var]]
  }

  if (isTRUE(time_intercepts)) {
    obs <- data.frame(
      response = y,
      split_group = if (isTRUE(by_time)) as.character(raw_time) else "All",
      stringsAsFactors = FALSE
    )
    params <- .plot_margin_time_intercept_params(y, raw_time, margin_dist, fit_control = fit_control)
    density_grid <- if (is.null(params)) {
      .plot_margin_empty_density()
    } else {
      .plot_margin_density_grid(
        y,
        margin_dist,
        params,
        grid_n = grid_n,
        group = if (isTRUE(by_time)) raw_time else NULL
      )
    }
  } else if (isTRUE(by_time)) {
    obs <- data.frame(response = y, split_group = as.character(raw_time), stringsAsFactors = FALSE)
    overlay_failed <- FALSE
    obs_split <- split(obs, obs$split_group)
    density_grid <- do.call(rbind, lapply(names(obs_split), function(group_name) {
      df <- obs_split[[group_name]]
      params <- .plot_margin_constant_params(df$response, margin_dist, warn = FALSE, fit_control = fit_control)
      grid <- if (is.null(params)) {
        overlay_failed <<- TRUE
        .plot_margin_empty_density(group_name)
      } else {
        .plot_margin_density_grid(df$response, margin_dist, params, grid_n = grid_n)
      }
      if (nrow(grid) > 0L) {
        grid$split_group <- group_name
      }
      grid
    }))
    if (isTRUE(overlay_failed)) {
      warning(
        "The ", .plot_margin_family_name(margin_dist), " marginal overlay fit did not converge for at least one time group; no fitted density was drawn for failed group(s).",
        call. = FALSE
      )
    }
  } else {
    obs <- data.frame(response = y, split_group = "All", stringsAsFactors = FALSE)
    params <- .plot_margin_constant_params(y, margin_dist, fit_control = fit_control)
    density_grid <- if (is.null(params)) {
      .plot_margin_empty_density()
    } else {
      .plot_margin_density_grid(y, margin_dist, params, grid_n = grid_n)
    }
  }

  list(
    margin_dist = margin_dist,
    obs = .plot_margin_transform_observed(obs, response_scale),
    density_grid = .plot_margin_transform_density(density_grid, response_scale)
  )
}

#' @keywords internal
#' @noRd
.plot_margin_fit_plot <- function(obs, density_grid, margin_dist, bins = 30, by_time = FALSE) {
  family_name <- .plot_margin_family_name(margin_dist)
  p <- ggplot2::ggplot(obs, ggplot2::aes(x = .data$response)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = bins, fill = "#d9d9d9", color = "white", na.rm = TRUE) +
    ggplot2::geom_line(data = density_grid, ggplot2::aes(x = .data$response, y = .data$density), inherit.aes = FALSE, color = "#e41a1c", linewidth = 1) +
    ggplot2::labs(
      title = paste0("Marginal Fit: ", family_name),
      x = "Response",
      y = "Density"
    ) +
    ggplot2::theme_minimal()

  if (isTRUE(by_time)) {
    p <- p + ggplot2::facet_wrap(~split_group, scales = "free_y")
  }

  p
}
