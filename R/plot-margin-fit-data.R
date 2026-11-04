#' @keywords internal
#' @noRd
.plot_margin_fit_from_object <- function(
    fit,
    x = NULL,
    data = NULL,
    grid_n = 200,
    by_time = FALSE,
    by = NULL,
    response_scale = "response") {
  newdata <- if (!is.null(data)) {
    data
  } else if (is.data.frame(x)) {
    x
  } else {
    NULL
  }
  diag_data <- .gl_fitted_distribution(fit, newdata = newdata, require_response = TRUE)
  split_info <- .gl_diag_split_info(
    by_time = by_time,
    by = by,
    diag_data = diag_data,
    data = newdata,
    plot_name = "marginal fit plot"
  )
  margin_dist <- fit$margin_dist
  obs <- data.frame(
    response = diag_data$response,
    split_group = if (isTRUE(split_info$split_by)) split_info$group else "All",
    stringsAsFactors = FALSE
  )
  density_grid <- .plot_margin_density_grid(
    y = diag_data$response,
    family = margin_dist,
    params = diag_data$params,
    grid_n = grid_n,
    group = if (isTRUE(split_info$split_by)) split_info$group else NULL
  )

  list(
    margin_dist = margin_dist,
    obs = .plot_margin_transform_observed(obs, response_scale),
    density_grid = .plot_margin_transform_density(density_grid, response_scale),
    by = split_info$by,
    split_by = split_info$split_by
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
    by = NULL,
    time_var = "time",
    fit_control = gamlss::gamlss.control(n.cyc = 50)) {
  margin_dist <- .plot_margin_resolve_family(margin_dist)
  y <- .plot_margin_response(x = x, data = data, response_var = response_var)
  raw_data <- if (!is.null(data)) as.data.frame(data, stringsAsFactors = FALSE) else if (is.data.frame(x)) as.data.frame(x, stringsAsFactors = FALSE) else NULL
  raw_time <- NULL
  needs_time <- isTRUE(time_intercepts) || isTRUE(by_time) ||
    (!is.null(by) && is.character(by) && length(by) == 1L && by %in% c("time", "response_margin"))
  if (needs_time) {
    if (is.null(raw_data) || !time_var %in% names(raw_data)) {
      stop("Raw-data time overlays require 'data' or data-frame 'x' containing 'time_var'.", call. = FALSE)
    }
    raw_time <- raw_data[[time_var]]
  }
  split_info <- .plot_margin_raw_split_info(
    by_time = by_time,
    by = by,
    raw_data = raw_data,
    raw_time = raw_time
  )
  split_group <- split_info$group

  if (isTRUE(time_intercepts)) {
    obs <- data.frame(
      response = y,
      split_group = if (isTRUE(split_info$split_by)) split_group else "All",
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
        group = if (isTRUE(split_info$split_by)) split_group else NULL
      )
    }
  } else if (isTRUE(split_info$split_by)) {
    obs <- data.frame(response = y, split_group = split_group, stringsAsFactors = FALSE)
    overlay_failed <- FALSE
    obs_split <- split(obs, obs$split_group, drop = TRUE)
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
        "The ", .plot_margin_family_name(margin_dist), " marginal overlay fit did not converge for at least one facet group; no fitted density was drawn for failed group(s).",
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
    density_grid = .plot_margin_transform_density(density_grid, response_scale),
    by = split_info$by,
    split_by = split_info$split_by
  )
}

#' Resolve raw-data faceting for a marginal fit plot
#'
#' @keywords internal
#' @noRd
.plot_margin_raw_split_info <- function(by_time, by, raw_data, raw_time) {
  .gl_check_by_time(by_time)

  if (!is.null(by) && (!is.character(by) || length(by) != 1L || is.na(by))) {
    stop("'by' must be NULL or a single column name as a character string.", call. = FALSE)
  }
  if (is.character(by) && length(by) == 1L && !nzchar(by)) {
    by <- NULL
  }

  if (isTRUE(by_time) && is.null(by)) {
    by <- "time"
  } else if (isTRUE(by_time) && !is.null(by)) {
    warning("Both by_time and by were provided; using by='", by, "'.", call. = FALSE)
  }

  if (is.null(by)) {
    return(list(by = NULL, group = NULL, split_by = FALSE))
  }
  if (by %in% c("time", "response_margin")) {
    return(list(by = by, group = as.factor(raw_time), split_by = TRUE))
  }
  if (is.null(raw_data)) {
    stop("To facet the marginal fit plot by '", by, "', provide data= containing that column.", call. = FALSE)
  }
  if (!by %in% names(raw_data)) {
    stop("Column '", by, "' not found in provided data.", call. = FALSE)
  }

  list(by = by, group = as.factor(raw_data[[by]]), split_by = TRUE)
}

#' Aggregate a wide discrete support for display
#'
#' @keywords internal
#' @noRd
.plot_margin_discrete_display_data <- function(
    obs,
    density_grid,
    max_bins = 30L,
    separate_zero = FALSE) {
  if (length(max_bins) != 1L || !is.numeric(max_bins) || !is.finite(max_bins) ||
      max_bins < 1 || max_bins != as.integer(max_bins)) {
    stop("'bins' must be one positive integer.", call. = FALSE)
  }
  if (!is.logical(separate_zero) || length(separate_zero) != 1L || is.na(separate_zero)) {
    stop("'separate_zero' must be TRUE or FALSE.", call. = FALSE)
  }
  max_bins <- as.integer(max_bins)

  observed_values <- obs$response[is.finite(obs$response)]
  fitted_values <- density_grid$response[is.finite(density_grid$response)]
  support_values <- c(observed_values, fitted_values)
  if (!length(support_values)) {
    return(list(
      obs = data.frame(
        response = numeric(), probability = numeric(), split_group = character(),
        bin_start = numeric(), bin_end = numeric(), bin_width = numeric(),
        stringsAsFactors = FALSE
      ),
      density = data.frame(
        response = numeric(), density = numeric(), split_group = character(),
        bin_start = numeric(), bin_end = numeric(), bin_width = numeric(),
        stringsAsFactors = FALSE
      ),
      bin_width = 1L,
      support_size = 0L
    ))
  }

  support_min <- floor(min(support_values))
  support_max <- ceiling(max(support_values))
  support_size <- support_max - support_min + 1
  support <- seq.int(support_min, support_max)
  use_separate_zero <- isTRUE(separate_zero) && 0 %in% support

  if (use_separate_zero && support_size > 1L && max_bins < 2L) {
    stop("'bins' must be at least 2 when 'separate_zero = TRUE' and non-zero values are present.", call. = FALSE)
  }

  if (use_separate_zero) {
    negative_size <- max(0L, as.integer(-support_min))
    positive_size <- max(0L, as.integer(support_max))
    side_count <- as.integer(negative_size > 0L) + as.integer(positive_size > 0L)
    available_nonzero_bins <- max_bins - 1L
    if (available_nonzero_bins < side_count) {
      stop("'bins' is too small to keep zero separate from both negative and positive values.", call. = FALSE)
    }
    nonzero_bin_width <- if (negative_size + positive_size == 0L) {
      1L
    } else {
      candidate_widths <- seq_len(max(negative_size, positive_size))
      required_bins <- ceiling(negative_size / candidate_widths) +
        ceiling(positive_size / candidate_widths)
      as.integer(candidate_widths[which(required_bins <= available_nonzero_bins)[1L]])
    }

    bin_start <- support
    bin_end <- support
    positive <- support > 0
    if (any(positive)) {
      bin_start[positive] <- 1 +
        floor((support[positive] - 1) / nonzero_bin_width) * nonzero_bin_width
      bin_end[positive] <- pmin(
        support_max,
        bin_start[positive] + nonzero_bin_width - 1
      )
    }
    negative <- support < 0
    if (any(negative)) {
      bin_end[negative] <- -1 -
        floor((-1 - support[negative]) / nonzero_bin_width) * nonzero_bin_width
      bin_start[negative] <- pmax(
        support_min,
        bin_end[negative] - nonzero_bin_width + 1
      )
    }
    unique_starts <- sort(unique(bin_start))
    bin_id <- match(bin_start, unique_starts) - 1L
    bin_width <- max(bin_end - bin_start + 1L)
  } else {
    bin_width <- max(1L, as.integer(ceiling(support_size / max_bins)))
    bin_id <- floor((support - support_min) / bin_width)
    bin_start <- support_min + bin_id * bin_width
    bin_end <- bin_start + bin_width - 1L
  }
  bin_lookup <- unique(data.frame(
    bin_id = bin_id,
    bin_start = bin_start,
    bin_end = bin_end,
    stringsAsFactors = FALSE
  ))
  bin_lookup <- bin_lookup[order(bin_lookup$bin_id), , drop = FALSE]
  bin_lookup$response <- (bin_lookup$bin_start + bin_lookup$bin_end) / 2
  bin_lookup$bin_width <- bin_lookup$bin_end - bin_lookup$bin_start + 1

  add_bin_columns <- function(data) {
    lookup_rows <- match(data$bin_id, bin_lookup$bin_id)
    data$bin_start <- bin_lookup$bin_start[lookup_rows]
    data$bin_end <- bin_lookup$bin_end[lookup_rows]
    data$response <- bin_lookup$response[lookup_rows]
    data$bin_width <- bin_lookup$bin_width[lookup_rows]
    data
  }

  obs_keep <- is.finite(obs$response) & !is.na(obs$split_group)
  obs_work <- obs[obs_keep, c("response", "split_group"), drop = FALSE]
  obs_work$split_group <- as.character(obs_work$split_group)
  obs_work$bin_id <- bin_id[match(obs_work$response, support)]
  obs_binned <- stats::aggregate(
    x = list(count = rep.int(1L, nrow(obs_work))),
    by = list(split_group = obs_work$split_group, bin_id = obs_work$bin_id),
    FUN = sum
  )
  obs_binned$probability <- ave(
    obs_binned$count,
    obs_binned$split_group,
    FUN = function(value) value / sum(value)
  )
  obs_binned <- add_bin_columns(obs_binned)

  density_keep <- is.finite(density_grid$response) & !is.na(density_grid$split_group)
  density_work <- density_grid[density_keep, c("response", "density", "split_group"), drop = FALSE]
  density_work$split_group <- as.character(density_work$split_group)
  if (nrow(density_work)) {
    density_work$bin_id <- bin_id[match(density_work$response, support)]
    density_binned <- stats::aggregate(
      x = list(density = density_work$density),
      by = list(
        split_group = density_work$split_group,
        bin_id = density_work$bin_id
      ),
      FUN = function(value) {
        if (!any(is.finite(value))) NA_real_ else sum(value[is.finite(value)])
      }
    )
    density_binned <- add_bin_columns(density_binned)
  } else {
    density_binned <- data.frame(
      response = numeric(), density = numeric(), split_group = character(),
      bin_start = numeric(), bin_end = numeric(), bin_width = numeric(),
      stringsAsFactors = FALSE
    )
  }

  obs_binned <- obs_binned[order(obs_binned$split_group, obs_binned$bin_id), , drop = FALSE]
  density_binned <- density_binned[order(density_binned$split_group, density_binned$bin_id), , drop = FALSE]
  rownames(obs_binned) <- NULL
  rownames(density_binned) <- NULL

  list(
    obs = obs_binned,
    density = density_binned,
    bin_width = bin_width,
    support_size = as.integer(support_size)
  )
}

#' @keywords internal
#' @noRd
.plot_margin_discrete_step_data <- function(density_grid) {
  if (!nrow(density_grid)) {
    return(data.frame(
      response = numeric(), density = numeric(), split_group = character(),
      stringsAsFactors = FALSE
    ))
  }

  group_rows <- split(seq_len(nrow(density_grid)), density_grid$split_group, drop = TRUE)
  step_data <- lapply(group_rows, function(rows) {
    data <- density_grid[rows, , drop = FALSE]
    data <- data[order(data$bin_start), , drop = FALSE]
    data.frame(
      response = c(rbind(data$bin_start - 0.5, data$bin_end + 0.5)),
      density = rep(data$density, each = 2L),
      split_group = rep(data$split_group, each = 2L),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, step_data)
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
.plot_margin_fit_plot <- function(
    obs,
    density_grid,
    margin_dist,
    bins = 30,
    separate_zero = FALSE,
    by_time = FALSE) {
  family_name <- .plot_margin_family_name(margin_dist)
  is_discrete <- .plot_margin_is_discrete(margin_dist)
  cols <- .gl_plot_colours()
  panel_levels <- NULL
  if (isTRUE(by_time)) {
    panel_levels <- if (is.factor(obs$split_group)) {
      levels(obs$split_group)[levels(obs$split_group) %in% as.character(obs$split_group)]
    } else {
      unique(as.character(obs$split_group))
    }
  }
  discrete_display <- NULL
  if (is_discrete) {
    discrete_display <- .plot_margin_discrete_display_data(
      obs = obs,
      density_grid = density_grid,
      max_bins = bins,
      separate_zero = separate_zero
    )
    obs <- discrete_display$obs
    density_grid <- discrete_display$density
  }
  if (isTRUE(by_time)) {
    obs$split_group <- factor(as.character(obs$split_group), levels = panel_levels)
    density_grid$split_group <- factor(
      as.character(density_grid$split_group),
      levels = panel_levels
    )
  }
  histogram_layer <- if (is_discrete) {
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$bin_start - 0.5,
        xmax = .data$bin_end + 0.5,
        ymin = 0,
        ymax = .data$probability
      ),
      fill = cols[["neutral"]],
      color = NA,
      na.rm = TRUE
    )
  } else {
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = bins,
      fill = cols[["neutral"]],
      color = NA,
      na.rm = TRUE
    )
  }
  overlay_layer <- if (is_discrete) {
    ggplot2::geom_line(
      data = .plot_margin_discrete_step_data(density_grid),
      ggplot2::aes(
        x = .data$response,
        y = .data$density,
        group = .data$split_group
      ),
      inherit.aes = FALSE,
      color = cols[["fitted"]],
      linewidth = 1
    )
  } else {
    ggplot2::geom_line(
      data = density_grid,
      ggplot2::aes(x = .data$response, y = .data$density),
      inherit.aes = FALSE,
      color = cols[["fitted"]],
      linewidth = 1
    )
  }
  p <- ggplot2::ggplot(obs, ggplot2::aes(x = .data$response)) +
    histogram_layer +
    overlay_layer +
    ggplot2::labs(
      title = paste0("Marginal Fit: ", family_name),
      x = "Response",
      y = if (is_discrete) "Probability" else "Density"
    ) +
    theme_gamlss_longitudinal(base_size = 12) +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 8, 8, 14))

  if (isTRUE(by_time)) {
    p <- p + ggplot2::facet_wrap(~split_group, scales = "free_y")
  }

  p
}
