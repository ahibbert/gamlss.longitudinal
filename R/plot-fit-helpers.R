#' Plot diagnostics dashboard for fitted `gamlss.longitudinal` objects
#'
#' Displays four ggplot-based diagnostic panels by default:
#' 1) PIT histogram
#' 2) QQ residual plot
#' 3) Worm plot
#' 4) Rootogram
#'
#' Fitted-data and newdata quantile forecast panels can be requested
#' explicitly with `include_fitted_quantiles` and `include_newdata_quantiles`.
#'
#' @param x A fitted `gamlss.longitudinal` object.
#' @param y Unused; included for S3 generic compatibility.
#' @param data Optional data frame used as fallback for `newdata` plotting.
#' @param newdata Optional data frame for the newdata forecast panel.
#' @param newdata_n Number of rows to use from `data` when `newdata` is NULL.
#' @param quantiles Quantiles for forecast panels.
#' @param include_fitted_quantiles Logical; if TRUE, include fitted-data
#'   forecast quantiles in the dashboard.
#' @param include_newdata_quantiles Logical; if TRUE, include newdata forecast
#'   quantiles in the dashboard. Uses `newdata`, or the first `newdata_n` rows
#'   of `data` if `newdata` is NULL.
#' @param randomize Logical; randomized PIT/residual diagnostics.
#' @param time_stratified Logical; if TRUE, show time-stratified diagnostic plots.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a list of generated plot/data objects.
#' @export
plot.gamlss.longitudinal = function(
  x,
  y,
  data = NULL,
  newdata = NULL,
  newdata_n = 8,
  quantiles = c(0.1, 0.5, 0.9),
  include_fitted_quantiles = FALSE,
  include_newdata_quantiles = FALSE,
  randomize = TRUE,
  time_stratified = FALSE,
  ...
) {
  if(!inherits(x, "gamlss.longitudinal")) {
    stop("'x' must be of class 'gamlss.longitudinal'.")
  }

  q_col_name = function(prob) {
    paste0("q", gsub("^0\\.", "", format(prob, trim = TRUE)))
  }

  make_empty_plot = function(title_txt, msg_txt) {
    ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_text(label = msg_txt) +
      ggplot2::xlim(-1, 1) +
      ggplot2::ylim(-1, 1) +
      ggplot2::labs(title = title_txt) +
      ggplot2::theme_void()
  }

  # 1-4: Standard diagnostics
  p_diag1 = pithist(x, bins = 20, randomize = randomize, plot = TRUE, by_time = time_stratified)
  p_diag2 = qqrplot(x, randomize = randomize, plot = TRUE, by_time = time_stratified)
  p_diag3 = wormplot(x, randomize = randomize, plot = TRUE, by_time = time_stratified)
  p_diag4 = rootogram(x, bins = 20, plot = TRUE, by_time = time_stratified)

  include_any_quantiles = isTRUE(include_fitted_quantiles) || isTRUE(include_newdata_quantiles)
  if(include_any_quantiles && length(quantiles) == 0L) {
    stop("'quantiles' must contain at least one probability when forecast quantile panels are requested.")
  }
  if(include_any_quantiles) {
    q_low = q_col_name(min(quantiles))
    q_high = q_col_name(max(quantiles))
    q_mid = q_col_name(if(0.5 %in% quantiles) 0.5 else quantiles[ceiling(length(quantiles) / 2)])
  }

  p_fit_quant = NULL
  fc_fit_q = NULL
  if(isTRUE(include_fitted_quantiles)) {
    fc_fit_q = procast(x, type = "quantile", at = quantiles)
    fc_fit_q$idx = seq_len(nrow(fc_fit_q))

    p_fit_quant = ggplot2::ggplot(fc_fit_q, ggplot2::aes(x = idx)) +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = q_low, ymax = q_high), fill = "#4e79a7", alpha = 0.25) +
      ggplot2::geom_line(ggplot2::aes_string(y = q_mid), color = "#1f4e79", linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(y = response), color = "black", alpha = 0.35, size = 0.9) +
      ggplot2::labs(
        title = "Fitted Forecast Quantiles",
        x = "Observation Index",
        y = "Response"
      ) +
      ggplot2::theme_minimal()
  }

  # Optional newdata forecast quantiles
  nd_use = NULL
  if(isTRUE(include_newdata_quantiles)) {
    nd_use = newdata
    if(is.null(nd_use) && !is.null(data) && is.data.frame(data)) {
      nd_use = utils::head(data, newdata_n)
      # ensure quantile-only mode (response optional)
      if(is.null(x$var_map) || !"response" %in% x$var_map) {
        nd_use$response = NA_real_
      } else {
        response_orig = names(x$var_map)[x$var_map == "response"][1]
        if(!is.na(response_orig) && !response_orig %in% names(nd_use) && !"response" %in% names(nd_use)) {
          nd_use[[response_orig]] = NA_real_
        }
      }
    }
  }

  p_new_quant = NULL
  fc_new_q = NULL
  if(!is.null(nd_use) && is.data.frame(nd_use) && nrow(nd_use) > 0) {
    fc_new_q = tryCatch(
      procast.gamlss.longitudinal(x, type = "quantile", at = quantiles, newdata = nd_use),
      error = function(e) NULL
    )

    if(!is.null(fc_new_q)) {
      time_candidates = c("time", if(!is.null(x$var_map)) names(x$var_map)[x$var_map == "time"] else character(0))
      person_candidates = c("subject", if(!is.null(x$var_map)) names(x$var_map)[x$var_map == "subject"] else character(0))
      time_col = time_candidates[time_candidates %in% names(nd_use)][1]
      person_col = person_candidates[person_candidates %in% names(nd_use)][1]

      if(is.na(time_col) || is.null(time_col) || nchar(time_col) == 0) {
        fc_new_q$time_plot = seq_len(nrow(fc_new_q))
      } else {
        fc_new_q$time_plot = nd_use[[time_col]]
      }
      if(is.na(person_col) || is.null(person_col) || nchar(person_col) == 0) {
        fc_new_q$person_plot = factor(seq_len(nrow(fc_new_q)))
      } else {
        fc_new_q$person_plot = as.factor(nd_use[[person_col]])
      }

      p_new_quant = ggplot2::ggplot(fc_new_q, ggplot2::aes(x = time_plot, color = person_plot, group = person_plot)) +
        ggplot2::geom_ribbon(ggplot2::aes_string(ymin = q_low, ymax = q_high, fill = "person_plot"), alpha = 0.14, color = NA, show.legend = FALSE) +
        ggplot2::geom_line(ggplot2::aes_string(y = q_mid), linewidth = 0.8) +
        ggplot2::geom_point(ggplot2::aes_string(y = q_mid), size = 1.7) +
        ggplot2::labs(
          title = "Newdata Forecast Quantiles",
          x = "Time",
          y = paste0("Predicted ", q_mid)
        ) +
        ggplot2::theme_minimal()
    }
  }

  if(isTRUE(include_newdata_quantiles) && is.null(p_new_quant)) {
    p_new_quant = make_empty_plot("Newdata Forecast Quantiles", "Provide 'newdata' or 'data' for this panel")
  }

  dashboard_plots = list(p_diag1, p_diag2, p_diag3, p_diag4)
  if(!is.null(p_fit_quant)) {
    dashboard_plots = c(dashboard_plots, list(p_fit_quant))
  }
  if(!is.null(p_new_quant)) {
    dashboard_plots = c(dashboard_plots, list(p_new_quant))
  }

  dashboard = ggpubr::ggarrange(
    plotlist = dashboard_plots,
    ncol = 2,
    nrow = ceiling(length(dashboard_plots) / 2)
  )
  print(dashboard)

  invisible(list(
    diagnostics = list(pithist = p_diag1, qqrplot = p_diag2, wormplot = p_diag3, rootogram = p_diag4),
    forecasts = list(fitted_quantiles = p_fit_quant, newdata_quantiles = p_new_quant),
    fitted_data = fc_fit_q,
    newdata_data = fc_new_q,
    dashboard = dashboard
  ))
}

#' @keywords internal
#' @noRd
.plot_margin_resolve_family <- function(margin_dist) {
  if (inherits(margin_dist, "margin_selection") || inherits(margin_dist, "margin_screen")) {
    margin_dist <- best_fit_family(margin_dist)
  }
  if (is.character(margin_dist) && length(margin_dist) == 1L) {
    margin_dist <- .margin_family_object(margin_dist)
  }
  if (is.null(margin_dist) || is.null(margin_dist$family) || is.null(margin_dist$parameters)) {
    stop("'margin_dist' must be a gamlss.dist family object, family name, or margin_selection result.", call. = FALSE)
  }
  margin_dist
}

#' @keywords internal
#' @noRd
.plot_reject_old_args <- function(dots, old_args) {
  dot_names <- names(dots)
  old_used <- intersect(dot_names[!is.na(dot_names) & nzchar(dot_names)], old_args)
  if (length(old_used) > 0L) {
    replacements <- c(
      family = "margin_dist",
      dist = "margin_dist",
      copula = "copula_dist",
      selected_fit = "copula_dist"
    )
    msg <- vapply(old_used, function(arg) {
      paste0("'", arg, "' has been removed; use '", replacements[[arg]], "' instead")
    }, character(1), USE.NAMES = FALSE)
    stop(paste(msg, collapse = "; "), call. = FALSE)
  }
  if (length(dots) > 0L) {
    stop("Unused argument(s): ", paste(dot_names, collapse = ", "), call. = FALSE)
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.plot_reject_old_call_args <- function(call, old_args) {
  call_names <- names(as.list(call)[-1L])
  old_used <- intersect(call_names[!is.na(call_names) & nzchar(call_names)], old_args)
  if (length(old_used) == 0L) {
    return(invisible(NULL))
  }
  replacements <- c(
    family = "margin_dist",
    dist = "margin_dist",
    copula = "copula_dist",
    selected_fit = "copula_dist"
  )
  msg <- vapply(old_used, function(arg) {
    paste0("'", arg, "' has been removed; use '", replacements[[arg]], "' instead")
  }, character(1), USE.NAMES = FALSE)
  stop(paste(msg, collapse = "; "), call. = FALSE)
}

#' @keywords internal
#' @noRd
.plot_margin_family_name <- function(family) {
  if (is.character(family$family)) {
    family$family[1]
  } else {
    as.character(family$family[[1L]])
  }
}

#' @keywords internal
#' @noRd
.plot_margin_response <- function(x = NULL, data = NULL, response_var = "response") {
  if (!is.null(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    if (!is.character(response_var) || length(response_var) != 1L || !response_var %in% names(data)) {
      stop("'response_var' must name a response column in 'data'.", call. = FALSE)
    }
    return(as.numeric(data[[response_var]]))
  }
  if (is.data.frame(x)) {
    if (!is.character(response_var) || length(response_var) != 1L || !response_var %in% names(x)) {
      stop("'response_var' must name a response column in 'x'.", call. = FALSE)
    }
    return(as.numeric(x[[response_var]]))
  }
  as.numeric(x)
}

#' @keywords internal
#' @noRd
.plot_margin_constant_params <- function(y, family, warn = TRUE, fit_control = gamlss::gamlss.control(n.cyc = 50)) {
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  if (length(y) < 3L) {
    stop("Need at least three finite response values to fit a marginal overlay.", call. = FALSE)
  }

  fit <- NULL
  invisible(utils::capture.output({
    fit <- suppressWarnings(suppressMessages(
      gamlss::gamlss(y ~ 1, family = family, trace = FALSE, control = fit_control)
    ))
  }))
  if (!.plot_margin_check_fit(fit, family, warn = warn)) {
    return(NULL)
  }
  params <- lapply(names(family$parameters), function(par_name) {
    as.numeric(stats::fitted(fit, what = par_name))
  })
  names(params) <- names(family$parameters)
  params
}

#' @keywords internal
#' @noRd
.plot_margin_time_intercept_params <- function(y, time, family, warn = TRUE, fit_control = gamlss::gamlss.control(n.cyc = 50)) {
  y <- as.numeric(y)
  time <- as.character(time)
  keep <- is.finite(y) & !is.na(time)
  if (sum(keep) < 3L) {
    stop("Need at least three finite response values with non-missing time values to fit a time-intercept marginal overlay.", call. = FALSE)
  }

  fit_data <- data.frame(
    y = y[keep],
    time_intercept = factor(time[keep], levels = unique(time[keep])),
    stringsAsFactors = FALSE
  )
  formula_time <- stats::as.formula("y ~ time_intercept")
  fit_args <- list(formula = formula_time, family = family, data = fit_data, trace = FALSE, control = fit_control)
  for (par_name in setdiff(names(family$parameters), "mu")) {
    fit_args[[paste0(par_name, ".formula")]] <- stats::as.formula("~ time_intercept")
  }
  fit <- NULL
  invisible(utils::capture.output({
    fit <- suppressWarnings(suppressMessages(do.call(gamlss::gamlss, fit_args)))
  }))
  if (!.plot_margin_check_fit(fit, family, warn = warn)) {
    return(NULL)
  }

  params <- lapply(names(family$parameters), function(par_name) {
    out <- rep(NA_real_, length(y))
    out[keep] <- as.numeric(stats::fitted(fit, what = par_name))
    out
  })
  names(params) <- names(family$parameters)
  params
}

#' @keywords internal
#' @noRd
.plot_margin_check_fit <- function(fit, family, warn = TRUE) {
  if (is.null(fit)) {
    stop("The marginal overlay fit failed.", call. = FALSE)
  }
  if (identical(fit$converged, FALSE)) {
    if (isTRUE(warn)) {
      family_name <- .plot_margin_family_name(family)
      warning(
        "The ", family_name, " marginal overlay fit did not converge; no fitted density was drawn.",
        call. = FALSE
      )
    }
    return(FALSE)
  }
  TRUE
}

#' @keywords internal
#' @noRd
.plot_margin_empty_density <- function(group = "All") {
  data.frame(
    response = numeric(),
    density = numeric(),
    split_group = as.character(group)[0L],
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.plot_margin_support_bounds <- function(family, params = NULL) {
  family_name <- .plot_margin_family_name(family)
  qfun <- tryCatch(
    get(paste0("q", family_name), envir = asNamespace("gamlss.dist"), inherits = FALSE),
    error = function(e) NULL
  )
  if (!is.function(qfun)) {
    return(c(lower = -Inf, upper = Inf))
  }

  q_args <- list()
  for (par_name in names(family$parameters)) {
    par_value <- NULL
    if (!is.null(params) && par_name %in% names(params)) {
      par_vec <- as.numeric(params[[par_name]])
      par_vec <- par_vec[is.finite(par_vec)]
      if (length(par_vec) > 0L) {
        par_value <- stats::median(par_vec)
      }
    }
    if (is.null(par_value)) {
      par_value <- tryCatch(
        eval(formals(qfun)[[par_name]], envir = baseenv()),
        error = function(e) NA_real_
      )
      par_value <- as.numeric(par_value)[1L]
    }
    if (is.finite(par_value)) {
      q_args[[par_name]] <- par_value
    }
  }

  bounds <- tryCatch(
    do.call(qfun, c(list(p = c(0, 1)), q_args)),
    error = function(e) {
      tryCatch(
        do.call(qfun, c(list(p = c(.Machine$double.eps, 1 - .Machine$double.eps)), q_args)),
        error = function(e2) c(-Inf, Inf)
      )
    }
  )
  bounds <- as.numeric(bounds)
  if (length(bounds) < 2L) {
    return(c(lower = -Inf, upper = Inf))
  }
  c(
    lower = if (is.finite(bounds[1L])) bounds[1L] else -Inf,
    upper = if (is.finite(bounds[2L])) bounds[2L] else Inf
  )
}

#' @keywords internal
#' @noRd
.plot_margin_density_grid <- function(y, family, params, grid_n = 200, group = NULL) {
  y <- as.numeric(y)
  keep <- is.finite(y)
  y <- y[keep]
  params <- lapply(params, function(x) {
    x <- as.numeric(x)
    if (length(x) == length(keep)) {
      x[keep]
    } else {
      rep(x, length.out = length(y))
    }
  })
  if (!is.null(group)) {
    group <- as.character(group)[keep]
  }

  if (length(y) < 1L) {
    return(data.frame(response = numeric(), density = numeric(), split_group = character()))
  }

  family_name <- .plot_margin_family_name(family)
  support_bounds <- .plot_margin_support_bounds(family, params)
  response_range <- range(y, na.rm = TRUE)
  pad <- diff(response_range) * 0.04
  if (!is.finite(pad) || pad <= 0) {
    pad <- max(abs(response_range), 1) * 0.04
  }

  build_one <- function(idx, group_name = "All") {
    yy <- y[idx]
    x_min <- min(yy, na.rm = TRUE) - pad
    x_max <- max(yy, na.rm = TRUE) + pad
    if (is.finite(support_bounds["lower"])) {
      x_min <- max(x_min, support_bounds["lower"])
      if (x_min <= support_bounds["lower"]) {
        yy_inside <- yy[is.finite(yy) & yy > support_bounds["lower"]]
        if (length(yy_inside) > 0L) {
          range_floor <- support_bounds["lower"] + (x_max - support_bounds["lower"]) * 0.01
          x_min <- max(
            as.numeric(stats::quantile(yy_inside, probs = 0.025, type = 8, names = FALSE)),
            range_floor,
            support_bounds["lower"] + max(diff(range(yy_inside, na.rm = TRUE)), abs(support_bounds["lower"]), 1) * 1e-6
          )
        }
      }
    }
    if (is.finite(support_bounds["upper"])) {
      x_max <- min(x_max, support_bounds["upper"])
    }
    if (!is.finite(x_min) || !is.finite(x_max) || x_max <= x_min) {
      x_min <- min(yy, na.rm = TRUE)
      x_max <- max(yy, na.rm = TRUE)
    }
    x_grid <- seq(x_min, x_max, length.out = grid_n)
    density <- vapply(x_grid, function(x_value) {
      par_i <- lapply(params, function(p) p[idx])
      d_value <- suppressWarnings(.gl_call_family_fun("d", family_name, rep(x_value, length(idx)), par_i))
      d_value[!is.finite(d_value)] <- NA_real_
      mean(d_value, na.rm = TRUE)
    }, numeric(1), USE.NAMES = FALSE)
    density[!is.finite(density)] <- NA_real_
    data.frame(
      response = x_grid,
      density = density,
      split_group = group_name,
      stringsAsFactors = FALSE
    )
  }

  if (is.null(group)) {
    return(build_one(seq_along(y)))
  }

  groups <- unique(group[!is.na(group)])
  out <- lapply(groups, function(g) build_one(which(group == g), group_name = g))
  do.call(rbind, out)
}

#' Plot an observed marginal response with a fitted GAMLSS density overlay
#'
#' `plot_margin_fit()` replaces the common exploratory use of
#' `gamlss::histDist()` with a ggplot-based helper that also understands final
#' `gamlss.longitudinal` fits. With raw data it fits the supplied family as an
#' intercept-only marginal model. With a fitted longitudinal model it overlays
#' the average row-specific fitted marginal density.
#'
#' @param x Numeric response vector, data frame, or fitted
#'   `gamlss.longitudinal` object.
#' @param margin_dist A `gamlss.dist` family object, family name, or
#'   `margin_selection` result. Required for raw data unless `fit` is supplied.
#' @param data Optional data frame containing the response.
#' @param fit Optional fitted `gamlss.longitudinal` object. When supplied, the
#'   fitted marginal density is overlaid using the model distribution and
#'   row-specific fitted parameters.
#' @param response_var Response column name when `x` or `data` is a data frame.
#' @param bins Number of histogram bins.
#' @param grid_n Number of grid points used for the fitted density.
#' @param response_scale Plot the density on the original response scale or,
#'   for positive responses, on the log-response scale.
#' @param time_intercepts Logical; for raw data, fit time-specific intercepts
#'   for each marginal distribution parameter before drawing the overlay.
#' @param by_time Logical; facet the plot by time. For fitted longitudinal
#'   models, fitted densities are averaged within each time point.
#' @param time_var Time column name used when `time_intercepts = TRUE` or
#'   `by_time = TRUE` for raw data.
#' @param fit_control Control object passed to the internal `gamlss()` overlay
#'   fit for raw-data plots. Defaults to `gamlss::gamlss.control(n.cyc = 50)`.
#' @param plot Logical; if TRUE, print the plot.
#' @param ... Additional arguments reserved for future methods.
#'
#' @return Invisibly returns a list containing the plot, observed data, and
#'   fitted density grid.
#' @export
plot_margin_fit <- function(
  x = NULL,
  margin_dist = NULL,
  data = NULL,
  fit = NULL,
  response_var = "response",
  bins = 30,
  grid_n = 200,
  response_scale = c("response", "log"),
  time_intercepts = FALSE,
  by_time = FALSE,
  time_var = "time",
  fit_control = gamlss::gamlss.control(n.cyc = 50),
  plot = TRUE,
  ...
) {
  .plot_reject_old_call_args(sys.call(), old_args = c("family"))
  .plot_reject_old_args(list(...), old_args = c("family"))
  response_scale <- match.arg(response_scale)

  if (inherits(x, "gamlss.longitudinal") && is.null(fit)) {
    fit <- x
    x <- NULL
  }

  transform_margin_obs <- function(plot_data) {
    if (response_scale == "response") {
      return(plot_data)
    }
    if (any(plot_data$response <= 0, na.rm = TRUE)) {
      stop("'response_scale = \"log\"' requires positive responses.", call. = FALSE)
    }
    plot_data$response <- log(plot_data$response)
    plot_data
  }

  transform_margin_density <- function(plot_data) {
    if (response_scale == "response") {
      return(plot_data)
    }
    plot_data <- plot_data[is.finite(plot_data$response) & plot_data$response > 0, , drop = FALSE]
    plot_data <- transform_margin_obs(plot_data)
    plot_data$density <- plot_data$density * exp(plot_data$response)
    plot_data
  }

  if (!is.null(fit)) {
    if (!inherits(fit, "gamlss.longitudinal")) {
      stop("'fit' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
    }
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
    obs <- transform_margin_obs(obs)
    density_grid <- transform_margin_density(density_grid)
  } else {
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
    } else {
      if (isTRUE(by_time)) {
        obs <- data.frame(response = y, split_group = as.character(raw_time), stringsAsFactors = FALSE)
        overlay_failed <- FALSE
        density_grid <- do.call(rbind, lapply(names(split(obs, obs$split_group)), function(group_name) {
          df <- split(obs, obs$split_group)[[group_name]]
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
    }

    obs <- transform_margin_obs(obs)
    density_grid <- transform_margin_density(density_grid)
  }

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

  if (isTRUE(plot)) {
    print(p)
  }

  invisible(list(plot = p, data = obs, density = density_grid))
}

#' @keywords internal
#' @noRd
.plot_copula_selection_spec <- function(copula) {
  if (inherits(copula, "copula_selection")) {
    copula <- as.data.frame(copula)[1L, , drop = FALSE]
  }
  if (is.character(copula) && length(copula) == 1L) {
    family <- .copula_family_code(copula)
    par <- switch(family, N = 0, F = 0, t = 0, C = 1e-8, G = 1 + 1e-8, J = 1 + 1e-8, 0)
    par2 <- if (identical(family, "t")) 4 else 0
    return(list(family = family, par = par, par2 = par2, tau = 0))
  }
  if (is.data.frame(copula) && nrow(copula) >= 1L) {
    family <- .copula_family_code(as.character(copula$family[1L]))
    par <- if ("par" %in% names(copula)) as.numeric(copula$par[1L]) else NA_real_
    par2 <- if ("par2" %in% names(copula)) as.numeric(copula$par2[1L]) else 0
    tau <- if ("tau" %in% names(copula)) as.numeric(copula$tau[1L]) else NA_real_
    if (!is.finite(par) && is.finite(tau)) {
      par <- .copula_tau_to_par(family, tau)
    }
    return(list(family = family, par = par, par2 = par2, tau = tau))
  }
  stop("'copula' must be a copula_selection result, one-row selection data frame, or family code.", call. = FALSE)
}

#' @keywords internal
#' @noRd
.plot_copula_candidate_families <- function(copula) {
  if (!is.character(copula) || length(copula) < 1L || any(is.na(copula))) {
    stop("'copula_dist' must be a copula_selection result, one-row selection data frame, or family code.", call. = FALSE)
  }
  if (length(copula) == 1L && tolower(copula) %in% c("best", "auto")) {
    return(c("N", "C", "F", "G", "J", "t"))
  }
  vapply(copula, .copula_family_code, character(1), USE.NAMES = FALSE)
}

#' @keywords internal
#' @noRd
.plot_copula_resolve_spec <- function(copula, pair_data = NULL, min_pairs = 3L) {
  if (is.character(copula)) {
    if (is.null(pair_data)) {
      return(list(spec = .plot_copula_selection_spec(copula[1L]), selection = NULL))
    }
    pair_data <- pair_data[is.finite(pair_data$u1) & is.finite(pair_data$u2), , drop = FALSE]
    if (nrow(pair_data) < min_pairs) {
      stop("Need at least ", min_pairs, " finite pseudo-observation pairs to fit the copula overlay.", call. = FALSE)
    }
    selection <- select_copula(
      u1 = pair_data$u1,
      u2 = pair_data$u2,
      families = .plot_copula_candidate_families(copula),
      min_pairs = min_pairs
    )
    return(list(spec = .plot_copula_selection_spec(selection), selection = selection))
  }
  list(spec = .plot_copula_selection_spec(copula), selection = if (inherits(copula, "copula_selection")) copula else NULL)
}

#' @keywords internal
#' @noRd
.plot_copula_density_for_spec <- function(pair_data, spec, grid_n, max_pairs_overlay) {
  pair_data$theta_pair <- if ("theta_pair" %in% names(pair_data)) pair_data$theta_pair else rep(spec$par, nrow(pair_data))
  pair_data$zeta_pair <- if ("zeta_pair" %in% names(pair_data)) pair_data$zeta_pair else rep(spec$par2, nrow(pair_data))
  .copula_v2_average_density_grid(
    family_num = spec$family,
    pair_data = pair_data,
    grid_n = grid_n,
    max_pairs_overlay = max_pairs_overlay
  )
}

#' @keywords internal
#' @noRd
.plot_copula_transform_grid <- function(density_grid, transform) {
  if (identical(transform, "normal")) {
    z1 <- stats::qnorm(.copula_v2_clamp01(density_grid$u1))
    z2 <- stats::qnorm(.copula_v2_clamp01(density_grid$u2))
    density_grid$density <- density_grid$density * stats::dnorm(z1) * stats::dnorm(z2)
    density_grid$u1 <- z1
    density_grid$u2 <- z2
  }
  density_grid
}

#' Plot empirical pseudo-observation pairs with a fitted copula overlay
#'
#' `plot_copula_fit()` and `plot_copula_overlay()` standardise the
#' copula-screen visual check used in examples and vignettes. They accept raw
#' pseudo-observations plus a `select_copula()` result, or a final
#' `gamlss.longitudinal` fit whose row-specific fitted copula density is
#' averaged over paired observations.
#'
#' @param data Optional long-format data frame. If a fitted
#'   `gamlss.longitudinal` object is supplied here, it is treated as `fit`.
#' @param copula_dist A `copula_selection` result, one-row selection data frame,
#'   family code(s), or `"best"`/`"auto"` to screen all supported families.
#'   Character values are fitted to the supplied pseudo-observation pairs before
#'   the overlay is drawn. Required unless `fit` is supplied.
#' @param fit Optional fitted `gamlss.longitudinal` object.
#' @param object Optional alias for `fit`.
#' @param u1,u2 Optional direct paired pseudo-observations.
#' @param u Optional row-aligned pseudo-observation vector for `data`.
#' @param u_var Optional pseudo-observation column name in `data`.
#' @param response_var Optional response column name used to create temporary
#'   pseudo-observations when `u`, `u_var`, and `u1`/`u2` are absent.
#' @param margin_dist Optional margin family used when creating temporary
#'   pseudo-observations from `response_var`.
#' @param subject_var,time_var Subject and time column names.
#' @param lags Positive integer lag(s) used when forming pairs.
#' @param by_time Logical; if `TRUE`, facet the copula overlay by adjacent
#'   time pair.
#' @param transform Character; either `"normal"` or `"uniform"`.
#' @param grid_n Grid size for fitted density contours.
#' @param bins Number of empirical two-dimensional bins.
#' @param contour_bins Number of fitted contour levels.
#' @param max_pairs_overlay Maximum paired observations used when averaging
#'   fitted copula densities.
#' @param plot Logical; if TRUE, print the plot.
#' @param ... Additional arguments reserved for future methods.
#'
#' @return Invisibly returns a list containing the plot, pair data, and fitted
#'   density grid.
#' @export
plot_copula_fit <- function(
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
  ...
) {
  .plot_reject_old_call_args(sys.call(), old_args = c("copula", "selected_fit"))
  .plot_reject_old_args(list(...), old_args = c("copula", "selected_fit"))
  transform <- match.arg(transform)
  if (!is.logical(by_time) || length(by_time) != 1L || is.na(by_time)) {
    stop("'by_time' must be TRUE or FALSE.", call. = FALSE)
  }
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

  if (isTRUE(by_time)) {
    pair_data$split_group <- .plot_copula_pair_time_group(pair_data)
  }

  density_grid <- if (isTRUE(by_time)) {
    do.call(rbind, lapply(split(pair_data, pair_data$split_group), function(group_data) {
      if (is.null(fit)) {
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
    }))
  } else {
    .plot_copula_density_for_spec(pair_data, spec, grid_n = grid_n, max_pairs_overlay = max_pairs_overlay)
  }
  density_grid <- .plot_copula_transform_grid(density_grid, transform)
  pair_plot <- .copula_v2_transform_data(pair_data, transform = transform)

  x_label <- if (transform == "normal") expression(Phi^-1 * (U[t])) else expression(U[t])
  y_label <- if (transform == "normal") expression(Phi^-1 * (U[t + 1])) else expression(U[t + 1])

  p <- ggplot2::ggplot(pair_plot, ggplot2::aes(x = .data$u1, y = .data$u2)) +
    ggplot2::geom_bin2d(ggplot2::aes(fill = ggplot2::after_stat(density)), bins = bins, alpha = 0.85) +
    ggplot2::geom_contour(
      data = density_grid,
      ggplot2::aes(x = .data$u1, y = .data$u2, z = .data$density),
      inherit.aes = FALSE,
      color = "#e41a1c",
      linewidth = 0.8,
      bins = contour_bins
    ) +
    ggplot2::scale_fill_gradient(low = "#f7f7f7", high = "#4d4d4d", name = "Empirical") +
    ggplot2::labs(
      title = paste0("Copula Overlay: ", spec$family),
      x = x_label,
      y = y_label
    ) +
    ggplot2::theme_minimal()

  if (isTRUE(by_time)) {
    p <- p + ggplot2::facet_wrap(~split_group)
  }

  if (isTRUE(plot)) {
    print(p)
  }

  invisible(list(
    plot = p,
    pair_data = pair_data,
    density = density_grid,
    copula = spec,
    selection = selection
  ))
}

.plot_copula_pair_time_group <- function(pair_data) {
  if ("time_pair" %in% names(pair_data)) {
    return(factor(pair_data$time_pair, levels = unique(pair_data$time_pair)))
  }
  if ("copula_time" %in% names(pair_data)) {
    return(factor(pair_data$copula_time, levels = unique(pair_data$copula_time)))
  }
  stop("'by_time = TRUE' requires pair data with adjacent time-pair labels.", call. = FALSE)
}

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
  ...
) {
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

#' @keywords internal
#' @noRd
.plot_dist_fit_var <- function(fit, data_names, role, explicit = NULL) {
  if (!is.null(explicit)) {
    return(explicit)
  }
  candidates <- character(0)
  if (!is.null(fit)) {
    stored <- fit[[paste0(role, "_var")]]
    candidates <- c(candidates, stored)
    if (!is.null(fit$var_map)) {
      candidates <- c(candidates, names(fit$var_map)[fit$var_map == role])
    }
  }
  candidates <- c(candidates, role)
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  hit <- candidates[candidates %in% data_names]
  if (length(hit) > 0L) {
    return(hit[1L])
  }
  NULL
}

#' @keywords internal
#' @noRd
.plot_dist_time_values <- function(time) {
  if (is.factor(time)) {
    lev <- levels(time)
    lev[lev %in% as.character(unique(time))]
  } else {
    values <- unique(time)
    if (is.numeric(values) || is.integer(values)) sort(values) else sort(as.character(values))
  }
}

#' @keywords internal
#' @noRd
.plot_dist_normalise_data <- function(dataset, fit, subject_var, time_var, response_var) {
  if (is.null(dataset)) {
    if (!is.null(fit) && !is.null(fit$dataset)) {
      dataset <- fit$dataset
    } else {
      stop("Supply 'dataset', or supply a fitted 'fit' with stored data.", call. = FALSE)
    }
  }
  dataset <- as.data.frame(dataset, stringsAsFactors = FALSE)

  if (!is.null(fit)) {
    subject_var <- .plot_dist_fit_var(fit, names(dataset), "subject", subject_var)
    time_var <- .plot_dist_fit_var(fit, names(dataset), "time", time_var)
    response_var <- .plot_dist_fit_var(fit, names(dataset), "response", response_var)
  }

  missing_inputs <- c(
    subject_var = is.null(subject_var),
    time_var = is.null(time_var),
    response_var = is.null(response_var)
  )
  if (any(missing_inputs)) {
    stop(
      "Raw-data plotting requires 'subject_var', 'time_var', and 'response_var'.",
      call. = FALSE
    )
  }

  missing_cols <- setdiff(c(subject_var, time_var, response_var), names(dataset))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in 'dataset': ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  data.frame(
    subject = dataset[[subject_var]],
    time = dataset[[time_var]],
    response = dataset[[response_var]],
    stringsAsFactors = FALSE
  )
}

#' Plot marginal and pairwise distribution diagnostics
#'
#' @param dataset Optional long-format data frame. Required unless `fit`
#'   contains stored data.
#' @param margin_dist Optional `gamlss.dist` family object used for diagonal
#'   marginal overlays.
#' @param subject_var,time_var,response_var Column names identifying subjects,
#'   time points, and responses. Required for raw data and inferred from `fit`
#'   when possible.
#' @param offdiag_scale Character; show off-diagonal panels on response or
#'   pseudo-observation scale.
#' @param transform Character; transform pseudo-observation off-diagonal panels
#'   to normal-score or uniform scale.
#' @param show_cor_stats Logical; include Pearson and Kendall correlations in
#'   off-diagonal panel subtitles.
#' @param fit Optional fitted `gamlss.longitudinal` object used for final-model
#'   margin and copula overlays.
#' @param overlay Character; `"none"` preserves the historical plot,
#'   `"margin"` overlays fitted marginal densities on diagonal panels,
#'   `"copula"` overlays a selected copula on pseudo-observation off-diagonal
#'   panels, and `"model"` overlays the final fitted model using `fit`.
#' @param copula_dist Optional `copula_selection` result or one-row selection table
#'   used when `overlay = "copula"`.
#' @param grid_n Grid size used for density overlays.
#' @param contour_bins Number of copula contour levels.
#' @param ... Compatibility arguments passed from `plotDist()` to
#'   `plot_dist()`.
#'
#' @return A `ggpubr` arranged plot object.
#' @export
plot_dist <- function (
  dataset = NULL,
  margin_dist = NULL,
  subject_var = NULL,
  time_var = NULL,
  response_var = NULL,
  offdiag_scale = c("pseudo", "response"),
  transform = c("normal", "uniform"),
  show_cor_stats = TRUE,
  fit = NULL,
  overlay = NULL,
  copula_dist = NULL,
  grid_n = 80,
  contour_bins = 8,
  ...
) {
  .plot_reject_old_call_args(sys.call(), old_args = c("dist", "family", "copula"))
  .plot_reject_old_args(list(...), old_args = c("dist", "family", "copula"))
  if (!is.null(fit) && !inherits(fit, "gamlss.longitudinal")) {
    stop("'fit' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }

  offdiag_scale <- match.arg(offdiag_scale)
  transform <- match.arg(transform)
  overlay_choices <- c("none", "margin", "copula", "model")
  overlay <- if (is.null(overlay)) {
    if (!is.null(fit)) "model" else "none"
  } else {
    match.arg(overlay, overlay_choices)
  }
  if (overlay == "model" && is.null(fit)) {
    stop("'fit' is required when overlay = 'model'.", call. = FALSE)
  }
  if (overlay == "copula" && is.null(copula_dist)) {
    stop("'copula_dist' is required when overlay = 'copula'.", call. = FALSE)
  }
  if (overlay %in% c("copula", "model") && offdiag_scale != "pseudo") {
    warning("Copula overlays are only drawn when offdiag_scale = 'pseudo'.", call. = FALSE)
  }
  if (overlay == "margin" && is.null(margin_dist) && is.null(fit)) {
    stop("'margin_dist' is required when overlay = 'margin'.", call. = FALSE)
  }

  plot_data <- .plot_dist_normalise_data(
    dataset = dataset,
    fit = fit,
    subject_var = subject_var,
    time_var = time_var,
    response_var = response_var
  )
  time_values <- .plot_dist_time_values(plot_data$time)
  num_margins=length(time_values)
  if (num_margins < 1L) {
    stop("No time points are available to plot.", call. = FALSE)
  }
  if (!is.null(margin_dist)) {
    margin_dist <- .plot_margin_resolve_family(margin_dist)
  } else if (!is.null(fit)) {
    margin_dist <- fit$margin_dist
  }

  fit_diag_data <- NULL
  fit_pair_data <- NULL
  fit_copula_spec <- NULL
  if (!is.null(fit)) {
    fit_diag_data <- tryCatch(
      .gl_fitted_distribution(fit, newdata = dataset, require_response = TRUE),
      error = function(e) .gl_fitted_distribution(fit, newdata = NULL, require_response = TRUE)
    )
    fit_data <- tryCatch(
      .copula_v2_fit_data(fit, data = dataset),
      error = function(e) .copula_v2_fit_data(fit)
    )
    fit_pair_data <- .copula_v2_pair_data(fit_data, lags = seq_len(max(1L, num_margins - 1L)))
    fit_copula_spec <- get_copula_dist(fit$copula_dist)
    fit_copula_spec <- list(
      family = .copula_family_code(fit_copula_spec$copula_dist),
      par = NA_real_,
      par2 = 0,
      tau = NA_real_
    )
  }
  copula_spec <- if (!is.null(copula_dist)) .plot_copula_selection_spec(copula_dist) else NULL

  margin_data=list()
  margin_pseudo=list()
  for (i in seq_len(num_margins)) {
    margin_data[[i]] <- plot_data[as.character(plot_data$time) == as.character(time_values[i]), c("subject", "response")]

    r <- rank(margin_data[[i]]$response, ties.method = "average", na.last = "keep")
    n_obs <- sum(!is.na(margin_data[[i]]$response))
    u <- r / (n_obs + 1)
    margin_pseudo[[i]] <- data.frame(subject = margin_data[[i]]$subject, u = u)
  }

  ##plot.new()
  #par(mfrow=c(1,num_margins))

  # Historical base graphics margin inspection is now covered by plot_margin_fit().
  #invisible(readline(prompt="Press [enter] to continue"))

  plots=list()

  z=1
  for (i in seq_len(num_margins)) {
    for (j in seq_len(num_margins)) {
      if(i==j) {
        input_data=data.frame(X1 = margin_data[[i]]$response)
        x_lab <- latex2exp::TeX(paste("$Y_",i,"$"))

        p <- ggplot2::ggplot(input_data, ggplot2::aes(x=X1))
        if (overlay %in% c("margin", "model")) {
          p <- p + ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins=30, na.rm=TRUE)
        } else {
          p <- p + ggplot2::geom_histogram(bins=30, na.rm=TRUE)
        }
        p <- p + ggplot2::labs(x = x_lab)

        if (overlay == "margin") {
          params <- .plot_margin_constant_params(input_data$X1, margin_dist)
          density_grid <- .plot_margin_density_grid(input_data$X1, margin_dist, params, grid_n = grid_n)
          p <- p +
            ggplot2::geom_line(
              data = density_grid,
              ggplot2::aes(x = .data$response, y = .data$density),
              inherit.aes = FALSE,
              color = "#e41a1c",
              linewidth = 0.8
            ) +
            ggplot2::labs(y = "Density")
        }

        if (overlay == "model") {
          keep_time <- as.character(fit_diag_data$time) == as.character(time_values[i])
          density_grid <- .plot_margin_density_grid(
            fit_diag_data$response[keep_time],
            fit$margin_dist,
            lapply(fit_diag_data$params, function(x) x[keep_time]),
            grid_n = grid_n
          )
          p <- p +
            ggplot2::geom_line(
              data = density_grid,
              ggplot2::aes(x = .data$response, y = .data$density),
              inherit.aes = FALSE,
              color = "#e41a1c",
              linewidth = 0.8
            ) +
            ggplot2::labs(y = "Density")
        }
      }
      if(i!=j) {
        if (offdiag_scale == "pseudo") {
          input_data <- merge(
            margin_pseudo[[i]],
            margin_pseudo[[j]],
            by = "subject",
            suffixes = c(".i", ".j"),
            all = FALSE
          )
          input_data <- input_data[complete.cases(input_data$u.i, input_data$u.j), c("u.i", "u.j")]
          names(input_data) <- c("X1", "X2")
          x_lab <- latex2exp::TeX(paste("$U_",i,"$"))
          y_lab <- latex2exp::TeX(paste("$U_",j,"$"))
          if (identical(transform, "normal")) {
            input_data$X1 <- stats::qnorm(.copula_v2_clamp01(input_data$X1))
            input_data$X2 <- stats::qnorm(.copula_v2_clamp01(input_data$X2))
            x_lab <- latex2exp::TeX(paste0("$\\Phi^{-1}(U_{", i, "})$"))
            y_lab <- latex2exp::TeX(paste0("$\\Phi^{-1}(U_{", j, "})$"))
          }
        } else {
          input_data <- merge(
            margin_data[[i]],
            margin_data[[j]],
            by = "subject",
            suffixes = c(".i", ".j"),
            all = FALSE
          )
          input_data <- input_data[complete.cases(input_data$response.i, input_data$response.j), c("response.i", "response.j")]
          names(input_data) <- c("X1", "X2")
          x_lab <- latex2exp::TeX(paste("$Y_",i,"$"))
          y_lab <- latex2exp::TeX(paste("$Y_",j,"$"))
        }

        p=ggplot2::ggplot(data=input_data,ggplot2::aes(x=X1,y=X2)) +
          ggplot2::geom_point(size=0.4, alpha=0.25, color="black", na.rm=TRUE) +
          ggplot2::geom_density_2d(contour_var="density",bins=10,color="black") +
          ggplot2::labs(x = x_lab, y = y_lab)

        if (show_cor_stats) {
          if (nrow(input_data) >= 3) {
            pearson_r <- suppressWarnings(cor(input_data$X1, input_data$X2, method = "pearson", use = "complete.obs"))
            kendall_tau <- suppressWarnings(cor(input_data$X1, input_data$X2, method = "kendall", use = "complete.obs"))
            stats_lab <- sprintf("Pearson r = %.3f | Kendall tau = %.3f", pearson_r, kendall_tau)
          } else {
            stats_lab <- "Pearson r = NA | Kendall tau = NA"
          }

          p <- p + ggplot2::labs(subtitle = stats_lab)
        }

        if (offdiag_scale == "pseudo" && overlay == "copula") {
          contour_grid <- .plot_copula_density_for_spec(
            data.frame(
              u1 = input_data$X1,
              u2 = input_data$X2,
              theta_pair = rep(copula_spec$par, nrow(input_data)),
              zeta_pair = rep(copula_spec$par2, nrow(input_data))
            ),
            copula_spec,
            grid_n = grid_n,
            max_pairs_overlay = 300
          )
          contour_grid <- .plot_copula_transform_grid(contour_grid, transform)
          contour_grid$X1 <- contour_grid$u1
          contour_grid$X2 <- contour_grid$u2
          p <- p + ggplot2::geom_contour(
            data = contour_grid,
            ggplot2::aes(x = .data$X1, y = .data$X2, z = .data$density),
            inherit.aes = FALSE,
            color = "#e41a1c",
            linewidth = 0.8,
            bins = contour_bins
          )
        }

        if (offdiag_scale == "pseudo" && overlay == "model") {
          left_idx <- min(i, j)
          right_idx <- max(i, j)
          pd <- fit_pair_data[
            as.character(fit_pair_data$time_left) == as.character(time_values[left_idx]) &
              as.character(fit_pair_data$time_right) == as.character(time_values[right_idx]),
            ,
            drop = FALSE
          ]
          if (nrow(pd) > 2L) {
            contour_grid <- .plot_copula_density_for_spec(
              pd,
              fit_copula_spec,
              grid_n = grid_n,
              max_pairs_overlay = 300
            )
            contour_grid <- .plot_copula_transform_grid(contour_grid, transform)
            if (i < j) {
              contour_grid$X1 <- contour_grid$u1
              contour_grid$X2 <- contour_grid$u2
            } else {
              contour_grid$X1 <- contour_grid$u2
              contour_grid$X2 <- contour_grid$u1
            }
            p <- p + ggplot2::geom_contour(
              data = contour_grid,
              ggplot2::aes(x = .data$X1, y = .data$X2, z = .data$density),
              inherit.aes = FALSE,
              color = "#e41a1c",
              linewidth = 0.8,
              bins = contour_bins
            )
          }
        }
      }

      plots[[z]]=p
      z=z+1
    }
  }
  ggpubr::ggarrange(plotlist=plots,ncol=num_margins,nrow=num_margins)

}

#' @rdname plot_dist
#' @export
plotDist <- function(...) {
  warning("plotDist() is retained for compatibility; use plot_dist() instead.", call. = FALSE)
  plot_dist(...)
}
