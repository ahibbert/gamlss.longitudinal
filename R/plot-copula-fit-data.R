#' @keywords internal
#' @noRd
.plot_copula_fit_inputs <- function(
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
    by = NULL) {
  split_info <- .plot_copula_split_info(by_time = by_time, by = by)
  if (inherits(data, "gamlss.longitudinal") && is.null(fit) && is.null(object)) {
    fit <- data
    data <- NULL
  }
  if (is.null(fit) && !is.null(object)) {
    fit <- object
  }

  if (!is.null(fit)) {
    if (!inherits(fit, "gamlss.longitudinal")) {
      stop("'fit' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
    }
    fit_data <- .copula_v2_fit_data(fit, data = data)
    pair_data <- .copula_v2_pair_data(fit_data, lags = lags)
    copula_spec <- get_copula_dist(fit$copula_dist)
    spec <- list(
      family = .copula_family_code(copula_spec$copula_dist),
      par = NA_real_,
      par2 = 0,
      tau = NA_real_
    )
    selection <- NULL
  } else {
    if (is.null(copula_dist)) {
      stop("Supply 'copula_dist' or a fitted 'fit'.", call. = FALSE)
    }
    pair_data <- .select_copula_pairs(
      data = data,
      object = NULL,
      u1 = u1,
      u2 = u2,
      u = u,
      u_var = u_var,
      response_var = response_var,
      margin_dist = margin_dist,
      mu.formula = NULL,
      sigma.formula = NULL,
      nu.formula = NULL,
      tau.formula = NULL,
      subject_var = subject_var,
      time_var = time_var,
      lags = lags
    )
  }

  pair_data <- pair_data[is.finite(pair_data$u1) & is.finite(pair_data$u2), , drop = FALSE]
  if (nrow(pair_data) < 3L) {
    stop("Need at least three finite pseudo-observation pairs to plot a copula overlay.", call. = FALSE)
  }
  if (is.null(fit)) {
    resolved <- .plot_copula_resolve_spec(copula_dist, pair_data = pair_data, min_pairs = 3L)
    spec <- resolved$spec
    selection <- resolved$selection
    pair_data$theta_pair <- rep(spec$par, nrow(pair_data))
    pair_data$zeta_pair <- rep(spec$par2, nrow(pair_data))
  }

  if (isTRUE(split_info$split_by)) {
    pair_data <- .plot_copula_attach_group(
      pair_data = pair_data,
      by = split_info$by,
      data = data,
      fit = fit,
      subject_var = subject_var,
      time_var = time_var
    )
  }

  list(
    fit = fit,
    pair_data = pair_data,
    spec = spec,
    selection = selection,
    by = split_info$by,
    split_by = split_info$split_by
  )
}

#' Resolve faceting for a copula fit plot
#'
#' @keywords internal
#' @noRd
.plot_copula_split_info <- function(by_time, by) {
  .gl_check_by_time(by_time)
  if (!is.null(by) && (!is.character(by) || length(by) != 1L || is.na(by))) {
    stop("'by' must be NULL or a single column name as a character string.", call. = FALSE)
  }
  if (is.character(by) && length(by) == 1L && !nzchar(by)) {
    by <- NULL
  }
  if (isTRUE(by_time) && is.null(by)) {
    by <- "time_pair"
  } else if (isTRUE(by_time) && !is.null(by)) {
    warning("Both by_time and by were provided; using by='", by, "'.", call. = FALSE)
  }
  list(by = by, split_by = !is.null(by))
}

#' Attach a requested facet variable to adjacent copula pairs
#'
#' @keywords internal
#' @noRd
.plot_copula_attach_group <- function(pair_data, by, data, fit, subject_var, time_var) {
  if (by %in% c("time", "time_pair", "copula_time")) {
    pair_data$split_group <- .plot_copula_pair_time_group(pair_data)
    return(pair_data)
  }
  if (by %in% names(pair_data)) {
    pair_data$split_group <- as.factor(pair_data[[by]])
    return(pair_data[!is.na(pair_data$split_group), , drop = FALSE])
  }
  if (is.null(data)) {
    stop("To facet the copula fit plot by '", by, "', provide data= containing that column.", call. = FALSE)
  }
  if (!all(c("subject", "time_left") %in% names(pair_data))) {
    stop("Faceting direct u1/u2 pairs requires a grouping column already present in the pair data.", call. = FALSE)
  }

  grouping_object <- fit
  if (is.null(grouping_object)) {
    var_map <- list()
    var_map[[subject_var]] <- "subject"
    var_map[[time_var]] <- "time"
    grouping_object <- list(var_map = var_map)
  }
  .copula_v2_attach_group(pair_data, object = grouping_object, by = by, data = data)
}

#' @keywords internal
#' @noRd
.plot_copula_fit_density_grid <- function(pair_data, spec, by_time = FALSE, grid_n = 80, max_pairs_overlay = 300, refit_by_group = FALSE) {
  if (isTRUE(by_time)) {
    return(do.call(rbind, lapply(split(pair_data, pair_data$split_group), function(group_data) {
      if (isTRUE(refit_by_group)) {
        group_fit <- .select_copula_fit_family(
          u1 = group_data$u1,
          u2 = group_data$u2,
          family = spec$family,
          t_df_grid = c(3, 4, 6, 8, 12, 20, 30)
        )
        group_data$theta_pair <- rep(group_fit$par[[1L]], nrow(group_data))
        group_data$zeta_pair <- rep(group_fit$par2[[1L]], nrow(group_data))
      }
      grid <- .plot_copula_density_for_spec(group_data, spec, grid_n = grid_n, max_pairs_overlay = max_pairs_overlay)
      grid$split_group <- as.character(group_data$split_group[[1L]])
      grid
    })))
  }

  .plot_copula_density_for_spec(pair_data, spec, grid_n = grid_n, max_pairs_overlay = max_pairs_overlay)
}

#' @keywords internal
#' @noRd
.plot_copula_fit_plot <- function(pair_plot, density_grid, spec, transform = "normal", bins = 28, contour_bins = 8, by_time = FALSE) {
  x_label <- if (transform == "normal") expression(Phi^-1 * (U[t])) else expression(U[t])
  y_label <- if (transform == "normal") expression(Phi^-1 * (U[t + 1])) else expression(U[t + 1])

  if (isTRUE(by_time)) {
    panel_levels <- if (is.factor(pair_plot$split_group)) {
      levels(pair_plot$split_group)[
        levels(pair_plot$split_group) %in% as.character(pair_plot$split_group)
      ]
    } else {
      unique(as.character(pair_plot$split_group))
    }
    pair_plot$split_group <- factor(
      as.character(pair_plot$split_group),
      levels = panel_levels
    )
    density_grid$split_group <- factor(
      as.character(density_grid$split_group),
      levels = panel_levels
    )
  }

  p <- ggplot2::ggplot(pair_plot, ggplot2::aes(x = .data$u1, y = .data$u2)) +
    ggplot2::geom_bin2d(ggplot2::aes(fill = ggplot2::after_stat(density)), bins = bins, alpha = 0.85) +
    ggplot2::geom_contour(
      data = density_grid,
      ggplot2::aes(x = .data$u1, y = .data$u2, z = .data$density),
      inherit.aes = FALSE,
      color = .gl_plot_colours()[["fitted"]],
      linewidth = 0.8,
      bins = contour_bins
    ) +
    ggplot2::scale_fill_gradient(
      low = "#F7F9FB",
      high = .gl_plot_colours()[["observed"]],
      guide = "none"
    ) +
    ggplot2::labs(
      title = paste0("Copula Overlay: ", spec$family),
      subtitle = "Darker bins show empirical density; blue contours show the fitted copula",
      x = x_label,
      y = y_label
    ) +
    theme_gamlss_longitudinal()

  if (isTRUE(by_time)) {
    p <- p + ggplot2::facet_wrap(~split_group, scales = "fixed")
  }

  p
}
