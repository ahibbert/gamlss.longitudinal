#' Compare fitted count margins across time
#'
#' Overlay several fitted count distributions on the same empirical density
#' histogram in each time panel. The fitted curves are marginal probability
#' masses averaged over the rows in a panel and divided by the displayed bin
#' width, so fitted
#' `gamlss.longitudinal` models with covariates remain directly comparable.
#'
#' Supply either a named list of fitted `gamlss.longitudinal` models in `fits`,
#' or raw `data` and a vector/list of candidate `families`. In the raw-data
#' route, a separate time-intercept model is fitted for every parameter of each
#' candidate family.
#'
#' @param fits Optional list of fitted `gamlss.longitudinal` models. All models
#'   must have been fitted to the same observed response-time pairs.
#' @param data Optional data frame used for a raw-data comparison.
#' @param families Candidate `gamlss.dist` family names or family objects for
#'   the raw-data route. Defaults to the five count families in the example.
#' @param response_var Response column name in `data`.
#' @param time_var Time column name in `data`.
#' @param bins Maximum number of displayed count bins. Adjacent integer values
#'   are combined when the observed support is wider than this value.
#' @param kde Logical; add a boundary-corrected Gaussian kernel density
#'   reference line in each panel.
#' @param kde_adjust Positive bandwidth multiplier for the KDE reference.
#' @param kde_n Number of points used to draw each KDE line.
#' @param facet_scales Scales passed to [ggplot2::facet_wrap()]. Use `"fixed"`
#'   to compare density heights across time or `"free_y"` to emphasize
#'   within-time shape.
#' @param fit_control Control object used by the raw-data GAMLSS fits.
#' @param plot Logical; if `TRUE`, print the plot.
#'
#' @return Invisibly returns a `gamlss_longitudinal_plot` containing `plot`,
#'   `data`, `density`, and `kde` components.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_margin_comparison(
#'   data = dat,
#'   families = c("ZISICHEL", "ZIBNB", "ZIPIG", "PO", "NBI"),
#'   response_var = "response",
#'   time_var = "time"
#' )
#'
#' plot_margin_comparison(list(
#'   ZISICHEL = fit_zisichel,
#'   ZIBNB = fit_zibnb,
#'   ZIPIG = fit_zipig,
#'   PO = fit_po,
#'   NBI = fit_nbi
#' ))
#' }
plot_margin_comparison <- function(
    fits = NULL,
    data = NULL,
    families = c("ZISICHEL", "ZIBNB", "ZIPIG", "PO", "NBI"),
    response_var = "response",
    time_var = "time",
    bins = 30,
    kde = TRUE,
    kde_adjust = 1,
    kde_n = 512,
    facet_scales = c("free_y", "fixed"),
    fit_control = gamlss::gamlss.control(n.cyc = 50),
    plot = TRUE) {
  facet_scales <- match.arg(facet_scales)
  .plot_margin_comparison_validate_controls(
    fits = fits,
    data = data,
    bins = bins,
    kde = kde,
    kde_adjust = kde_adjust,
    kde_n = kde_n
  )

  comparison <- if (!is.null(fits)) {
    .plot_margin_comparison_from_fits(fits)
  } else {
    .plot_margin_comparison_from_raw(
      data = data,
      families = families,
      response_var = response_var,
      time_var = time_var,
      fit_control = fit_control
    )
  }

  display <- .plot_margin_comparison_display_data(
    obs = comparison$obs,
    density_grid = comparison$density,
    max_bins = bins
  )
  model_levels <- unique(as.character(comparison$density$model))
  kde_grid <- if (isTRUE(kde)) {
    .plot_margin_comparison_kde(
      obs = comparison$obs,
      bin_width = display$bin_width,
      support = display$support,
      adjust = kde_adjust,
      n = kde_n
    )
  } else {
    data.frame(response = numeric(), density = numeric(), split_group = character())
  }

  panel_levels <- unique(as.character(comparison$obs$split_group))
  display$obs$split_group <- factor(display$obs$split_group, levels = panel_levels)
  display$density$split_group <- factor(display$density$split_group, levels = panel_levels)
  display$density$model <- factor(display$density$model, levels = model_levels)
  kde_grid$split_group <- factor(kde_grid$split_group, levels = panel_levels)

  p <- .plot_margin_comparison_plot(
    obs = display$obs,
    density_grid = display$density,
    kde_grid = kde_grid,
    bin_width = display$bin_width,
    facet_scales = facet_scales
  )
  if (isTRUE(plot)) {
    print(p)
  }

  invisible(.gl_new_plot_result(
    list(
      plot = p,
      data = display$obs,
      density = display$density,
      kde = kde_grid
    ),
    primary = "plot"
  ))
}

#' @keywords internal
#' @noRd
.plot_margin_comparison_validate_controls <- function(fits, data, bins, kde, kde_adjust, kde_n) {
  if (is.null(fits) == is.null(data)) {
    stop("Supply exactly one of 'fits' or 'data'.", call. = FALSE)
  }
  if (length(bins) != 1L || !is.numeric(bins) || !is.finite(bins) ||
      bins < 1 || bins != as.integer(bins)) {
    stop("'bins' must be one positive integer.", call. = FALSE)
  }
  if (length(kde) != 1L || !is.logical(kde) || is.na(kde)) {
    stop("'kde' must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(kde_adjust) != 1L || !is.numeric(kde_adjust) ||
      !is.finite(kde_adjust) || kde_adjust <= 0) {
    stop("'kde_adjust' must be one positive number.", call. = FALSE)
  }
  if (length(kde_n) != 1L || !is.numeric(kde_n) || !is.finite(kde_n) ||
      kde_n < 2 || kde_n != as.integer(kde_n)) {
    stop("'kde_n' must be an integer of at least two.", call. = FALSE)
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.plot_margin_comparison_labels <- function(families, supplied_names = NULL) {
  labels <- vapply(families, .plot_margin_family_name, character(1))
  if (!is.null(supplied_names)) {
    use_name <- !is.na(supplied_names) & nzchar(supplied_names)
    labels[use_name] <- supplied_names[use_name]
  }
  make.unique(labels)
}

#' @keywords internal
#' @noRd
.plot_margin_comparison_signature <- function(response, time) {
  sort(paste(
    as.character(time),
    format(as.numeric(response), digits = 17, scientific = FALSE, trim = TRUE),
    sep = "\r"
  ))
}

#' @keywords internal
#' @noRd
.plot_margin_comparison_from_fits <- function(fits) {
  if (inherits(fits, "gamlss.longitudinal")) {
    fits <- list(fits)
  }
  if (!is.list(fits) || length(fits) < 2L ||
      !all(vapply(fits, inherits, logical(1), what = "gamlss.longitudinal"))) {
    stop("'fits' must contain at least two fitted 'gamlss.longitudinal' models.", call. = FALSE)
  }

  families <- lapply(fits, function(fit) fit$margin_dist)
  if (!all(vapply(families, .plot_margin_is_discrete, logical(1)))) {
    stop("All compared marginal families must be discrete.", call. = FALSE)
  }
  labels <- .plot_margin_comparison_labels(families, names(fits))
  diagnostics <- lapply(fits, .gl_fitted_distribution, require_response = TRUE)
  reference_signature <- .plot_margin_comparison_signature(
    diagnostics[[1L]]$response,
    diagnostics[[1L]]$time
  )
  same_data <- vapply(diagnostics[-1L], function(x) {
    identical(
      .plot_margin_comparison_signature(x$response, x$time),
      reference_signature
    )
  }, logical(1))
  if (!all(same_data)) {
    stop("All models in 'fits' must use the same observed response-time pairs.", call. = FALSE)
  }

  obs <- data.frame(
    response = diagnostics[[1L]]$response,
    split_group = as.character(diagnostics[[1L]]$time),
    stringsAsFactors = FALSE
  )
  density <- do.call(rbind, lapply(seq_along(fits), function(i) {
    grid <- .plot_margin_density_grid(
      y = diagnostics[[i]]$response,
      family = families[[i]],
      params = diagnostics[[i]]$params,
      group = diagnostics[[i]]$time
    )
    grid$model <- labels[[i]]
    grid
  }))
  rownames(density) <- NULL
  list(obs = obs, density = density)
}

#' @keywords internal
#' @noRd
.plot_margin_comparison_from_raw <- function(data, families, response_var, time_var, fit_control) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  missing_columns <- setdiff(c(response_var, time_var), names(data))
  if (length(missing_columns)) {
    stop("Column(s) not found in 'data': ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  if (is.character(families)) {
    family_names <- names(families)
    families <- lapply(families, .plot_margin_resolve_family)
  } else if (is.list(families)) {
    family_names <- names(families)
    families <- lapply(families, .plot_margin_resolve_family)
  } else {
    families <- list(.plot_margin_resolve_family(families))
    family_names <- NULL
  }
  if (length(families) < 2L) {
    stop("'families' must contain at least two candidate distributions.", call. = FALSE)
  }
  if (!all(vapply(families, .plot_margin_is_discrete, logical(1)))) {
    stop("All compared marginal families must be discrete.", call. = FALSE)
  }
  labels <- .plot_margin_comparison_labels(families, family_names)

  fitted <- lapply(seq_along(families), function(i) {
    tryCatch(
      .plot_margin_fit_from_raw(
        data = data,
        margin_dist = families[[i]],
        response_var = response_var,
        time_intercepts = TRUE,
        by_time = TRUE,
        time_var = time_var,
        fit_control = fit_control
      ),
      error = function(e) {
        warning(
          "Could not fit ", labels[[i]], " for the marginal comparison: ",
          conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
  })
  usable <- vapply(fitted, function(x) !is.null(x) && nrow(x$density_grid) > 0L, logical(1))
  if (sum(usable) < 2L) {
    stop("Fewer than two candidate marginal models produced usable fitted probabilities.", call. = FALSE)
  }
  fitted <- fitted[usable]
  labels <- labels[usable]
  density <- do.call(rbind, lapply(seq_along(fitted), function(i) {
    grid <- fitted[[i]]$density_grid
    grid$model <- labels[[i]]
    grid
  }))
  rownames(density) <- NULL
  list(obs = fitted[[1L]]$obs, density = density)
}

#' @keywords internal
#' @noRd
.plot_margin_comparison_display_data <- function(obs, density_grid, max_bins = 30L) {
  observed_values <- obs$response[is.finite(obs$response)]
  fitted_values <- density_grid$response[is.finite(density_grid$response)]
  support_values <- c(observed_values, fitted_values)
  if (!length(support_values)) {
    stop("No finite response values are available to plot.", call. = FALSE)
  }
  if (any(abs(observed_values - round(observed_values)) > sqrt(.Machine$double.eps))) {
    stop("Count-margin comparisons require integer-valued responses.", call. = FALSE)
  }

  support_min <- floor(min(support_values))
  support_max <- ceiling(max(support_values))
  support_size <- support_max - support_min + 1
  bin_width <- max(1L, as.integer(ceiling(support_size / as.integer(max_bins))))

  add_bins <- function(x) {
    x$bin_start <- support_min + x$bin_id * bin_width
    x$bin_end <- pmin(x$bin_start + bin_width - 1, support_max)
    x$response <- (x$bin_start + x$bin_end) / 2
    x
  }

  obs_keep <- is.finite(obs$response) & !is.na(obs$split_group)
  obs_work <- obs[obs_keep, c("response", "split_group"), drop = FALSE]
  obs_work$split_group <- as.character(obs_work$split_group)
  obs_work$bin_id <- floor((obs_work$response - support_min) / bin_width)
  obs_binned <- stats::aggregate(
    x = list(count = rep.int(1L, nrow(obs_work))),
    by = list(split_group = obs_work$split_group, bin_id = obs_work$bin_id),
    FUN = sum
  )
  obs_binned$probability <- ave(
    obs_binned$count,
    obs_binned$split_group,
    FUN = function(x) x / sum(x)
  )
  obs_binned <- add_bins(obs_binned)

  density_keep <- is.finite(density_grid$response) &
    is.finite(density_grid$density) &
    !is.na(density_grid$split_group) &
    !is.na(density_grid$model)
  density_work <- density_grid[density_keep, c("response", "density", "split_group", "model"), drop = FALSE]
  density_work$split_group <- as.character(density_work$split_group)
  density_work$model <- as.character(density_work$model)
  density_work$bin_id <- floor((density_work$response - support_min) / bin_width)
  density_binned <- stats::aggregate(
    x = list(density = density_work$density),
    by = list(
      split_group = density_work$split_group,
      model = density_work$model,
      bin_id = density_work$bin_id
    ),
    FUN = sum
  )
  density_binned <- add_bins(density_binned)

  list(
    obs = obs_binned[order(obs_binned$split_group, obs_binned$bin_id), , drop = FALSE],
    density = density_binned[order(density_binned$split_group, density_binned$model, density_binned$bin_id), , drop = FALSE],
    bin_width = bin_width,
    support = c(lower = support_min, upper = support_max)
  )
}

#' @keywords internal
#' @noRd
.plot_margin_comparison_kde <- function(obs, bin_width, support, adjust = 1, n = 512L) {
  groups <- unique(as.character(obs$split_group[!is.na(obs$split_group)]))
  out <- lapply(groups, function(group_name) {
    y <- obs$response[as.character(obs$split_group) == group_name & is.finite(obs$response)]
    if (!length(y)) {
      return(NULL)
    }
    bandwidth <- suppressWarnings(stats::bw.nrd0(y))
    if (!is.finite(bandwidth) || bandwidth <= 0) {
      bandwidth <- max(0.5, bin_width / 2)
    }
    # Reflection at zero prevents a count-data reference curve from leaking
    # probability into the negative half-line. Multiplication by two undoes
    # the doubled sample size used by stats::density().
    estimate <- stats::density(
      c(y, -y),
      bw = bandwidth * adjust,
      from = max(0, support[["lower"]] - bin_width / 2),
      to = support[["upper"]] + bin_width / 2,
      n = as.integer(n)
    )
    data.frame(
      response = estimate$x,
      density = 2 * estimate$y,
      split_group = group_name,
      stringsAsFactors = FALSE
    )
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) {
    return(data.frame(response = numeric(), density = numeric(), split_group = character()))
  }
  do.call(rbind, out)
}

#' @keywords internal
#' @noRd
.plot_margin_comparison_plot <- function(
    obs,
    density_grid,
    kde_grid,
    bin_width,
    facet_scales = "free_y") {
  models <- if (is.factor(density_grid$model)) {
    levels(density_grid$model)
  } else {
    unique(as.character(density_grid$model))
  }
  model_colours <- stats::setNames(.gl_group_colours(length(models)), models)
  cols <- .gl_plot_colours()
  show_kde <- nrow(kde_grid) > 0L

  p <- ggplot2::ggplot(obs, ggplot2::aes(x = .data$response)) +
    ggplot2::geom_col(
      ggplot2::aes(y = .data$probability / bin_width),
      width = bin_width,
      fill = cols[["neutral"]],
      colour = "white",
      linewidth = 0.15,
      alpha = 0.8,
      na.rm = TRUE
    ) +
    ggplot2::geom_step(
      data = density_grid,
      ggplot2::aes(
        x = .data$response,
        y = .data$density / bin_width,
        colour = .data$model,
        group = .data$model
      ),
      inherit.aes = FALSE,
      direction = "mid",
      linewidth = 0.9,
      na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(values = model_colours) +
    ggplot2::labs(
      title = "Marginal count-distribution fit across time",
      subtitle = if (show_kde) {
        "Bars: observed density; coloured steps: fitted marginal density; dashed line: KDE reference"
      } else {
        "Bars: observed density; coloured steps: fitted marginal density"
      },
      x = "Response",
      y = "Probability density",
      colour = "Distribution",
      linetype = NULL,
      caption = if (show_kde) {
        "The KDE is boundary-corrected; it is a descriptive guide and does not model a structural zero mass."
      } else {
        NULL
      }
    ) +
    ggplot2::facet_wrap(~split_group, scales = facet_scales) +
    theme_gamlss_longitudinal()

  if (show_kde) {
    p <- p +
      ggplot2::geom_line(
        data = kde_grid,
        ggplot2::aes(
          x = .data$response,
          y = .data$density,
          linetype = "KDE reference"
        ),
        inherit.aes = FALSE,
        colour = cols[["observed"]],
        linewidth = 0.8,
        na.rm = TRUE
      ) +
      ggplot2::scale_linetype_manual(values = c("KDE reference" = "22"))
  }
  p
}
