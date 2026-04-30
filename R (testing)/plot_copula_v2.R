# Copula diagnostics v2 for testing a gamlss.longitudinal fit.
#
# Two plots only:
# 1) empirical copula with fitted overlay
# 2) observed vs fitted Kendall correlation by quartile

utils::globalVariables(c("u1", "u2", "quartile", "tau_emp", "tau_fit", "density", "x_id", "time_pair", "split_group"))

.copula_v2_clamp01 <- function(x) {
  pmin(pmax(x, 0.001), 0.999)
}

.copula_v2_tau_from_par <- function(family_num, par, par2 = NA_real_) {
  # Handle NA inputs immediately
  if (!is.finite(par)) {
    return(NA_real_)
  }

  tau <- tryCatch({
    if (is.finite(par2)) {
      suppressWarnings(VineCopula::BiCopPar2Tau(family = family_num, par = par, par2 = par2))
    } else {
      suppressWarnings(VineCopula::BiCopPar2Tau(family = family_num, par = par, par2 = 0))
    }
  }, error = function(e) NA_real_)

  if (is.finite(tau)) {
    return(as.numeric(tau))
  }

  # Fallback for Gaussian/Clayton copulas using formula
  if (family_num %in% c(1, 2) && is.finite(par)) {
    return(2 / pi * asin(max(min(par, 0.999999), -0.999999)))
  }

  NA_real_
}

.copula_v2_message_plot <- function(title, subtitle, message) {
  ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_text(label = message, size = 4) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_void()
}

.copula_v2_fit_data <- function(object) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.")
  }

  copula_spec <- get_copula_dist(object$copula_dist)
  
  # Extract proper copula family name for BiCopName
  copula_family_name <- object$copula_dist
  if (!is.character(copula_family_name) || nchar(copula_family_name) == 0) {
    copula_family_name <- copula_spec$copula_dist
  }
  if (grepl("^[0-9]+$", copula_family_name)) {
    family_map <- c("1" = "N", "2" = "C", "3" = "G", "4" = "F", "5" = "J", "6" = "BB1", "7" = "BB6", "8" = "BB7", "9" = "BB8", "10" = "T")
    if (copula_family_name %in% names(family_map)) {
      copula_family_name <- family_map[[copula_family_name]]
    }
  }
  
  eta_out <- calc_eta(
    par_cov = object$par,
    mm = object$model_matrix,
    margin_dist = object$margin_dist,
    copula_link = copula_spec$copula_link,
    par_s = object$par_s
  )

  response <- object$response
  subject <- object$response_subject
  time <- object$response_margin

  # Extract only margin parameters that are actually in eta_out$eta_inv
  margin_param_names <- names(object$margin_dist$parameters)
  margin_params <- list()
  for (param_name in margin_param_names) {
    if (param_name %in% names(eta_out$eta_inv)) {
      margin_params[[param_name]] <- eta_out$eta_inv[[param_name]]
    }
  }

  theta_fit <- if ("theta" %in% names(eta_out$eta_inv)) eta_out$eta_inv$theta else numeric(0)
  zeta_fit <- if ("zeta" %in% names(eta_out$eta_inv)) eta_out$eta_inv$zeta else numeric(0)

  # Align response-side vectors to a common leading length.
  margin_min_n <- if (length(margin_params) > 0) {
    min(vapply(margin_params, length, integer(1)))
  } else {
    length(response)
  }

  common_n <- min(length(response), length(subject), length(time), margin_min_n)
  if (!is.finite(common_n) || common_n < 1) {
    stop("No finite fitted observations are available for copula diagnostics.")
  }

  response <- response[seq_len(common_n)]
  subject <- subject[seq_len(common_n)]
  time <- time[seq_len(common_n)]
  margin_params <- lapply(margin_params, function(x) x[seq_len(common_n)])

  align_copula_param <- function(param_vec) {
    n_resp <- common_n
    if (length(param_vec) == 0) {
      return(rep(NA_real_, n_resp))
    }

    # Full-row parameterization.
    if (length(param_vec) == n_resp) {
      return(param_vec)
    }

    # Pair-row parameterization: parameters correspond to times 1:(T-1) only.
    margin_names <- sort(unique(time))
    left_time_rows <- which(time %in% margin_names[seq_len(max(1, length(margin_names) - 1))])
    if (length(param_vec) == length(left_time_rows)) {
      out <- rep(NA_real_, n_resp)
      out[left_time_rows] <- param_vec
      return(out)
    }

    # Fallback for unexpected lengths.
    rep(param_vec, length.out = n_resp)
  }

  theta_fit <- align_copula_param(theta_fit)
  zeta_fit <- align_copula_param(zeta_fit)

  # Filter by finite values
  keep <- is.finite(response)
  for (param_name in names(margin_params)) {
    keep <- keep & is.finite(margin_params[[param_name]])
  }

  response <- response[keep]
  subject <- subject[keep]
  time <- time[keep]
  margin_params <- lapply(margin_params, function(x) x[keep])
  theta_fit <- theta_fit[keep]
  zeta_fit <- zeta_fit[keep]

  if (length(response) == 0) {
    stop("No finite fitted observations are available for copula diagnostics.")
  }

  # Convert margin_dist$family to family name if needed
  family_name <- object$margin_dist$family[1]
  if (!is.character(family_name)) {
    family_name <- object$margin_dist$family[1]$family
  }

  u <- .gl_call_family_fun("p", family_name, response, margin_params)
  u <- .copula_v2_clamp01(u)

  family_num <- tryCatch({
    as.numeric(VineCopula::BiCopName(copula_family_name))
  }, error = function(e) NA_numeric_)
  
  # Compute tau_fit, suppressing coercion warnings
  tau_fit <- suppressWarnings(
    vapply(seq_along(theta_fit), function(i) {
      .copula_v2_tau_from_par(family_num, theta_fit[i], zeta_fit[i])
    }, numeric(1), USE.NAMES = FALSE)
  )

  data.frame(
    subject = subject,
    time = time,
    response = response,
    u = u,
    theta_fit = theta_fit,
    zeta_fit = zeta_fit,
    tau_fit = tau_fit,
    stringsAsFactors = FALSE
  )
}

.copula_v2_pair_data <- function(fit_data, lags = 1) {
  time_vec <- fit_data$time
  time_levels <- if (is.factor(time_vec)) {
    lev <- levels(time_vec)
    lev[lev %in% as.character(unique(time_vec))]
  } else {
    u <- unique(time_vec)
    if (is.numeric(u) || is.integer(u)) sort(u) else sort(as.character(u))
  }
  if (length(time_levels) < 2) {
    stop("Need at least two time points to build copula pair diagnostics.")
  }

  time_lookup <- setNames(seq_along(time_levels), as.character(time_levels))
  fit_data$time_idx <- unname(time_lookup[as.character(fit_data$time)])
  if (any(!is.finite(fit_data$time_idx))) {
    stop("Could not map time values to an ordered index for copula pair diagnostics.")
  }

  lag_values <- sort(unique(as.integer(lags)))
  lag_values <- lag_values[lag_values > 0]
  if (length(lag_values) == 0) {
    lag_values <- 1L
  }

  pair_list <- list()
  idx <- 1L

  for (lag_value in lag_values) {
    for (subject_id in unique(fit_data$subject)) {
      subject_rows <- fit_data[fit_data$subject == subject_id, , drop = FALSE]
      subject_rows <- subject_rows[order(subject_rows$time_idx), , drop = FALSE]
      if (nrow(subject_rows) < 2) next

      for (j in seq_len(nrow(subject_rows) - lag_value)) {
        k <- j + lag_value
        if (k > nrow(subject_rows)) next

        t1 <- subject_rows$time[j]
        t2 <- subject_rows$time[k]
        t1_idx <- subject_rows$time_idx[j]
        t2_idx <- subject_rows$time_idx[k]
        if ((t2_idx - t1_idx) != lag_value) next

        row1 <- subject_rows[j, , drop = FALSE]
        row2 <- subject_rows[k, , drop = FALSE]

        # Match likelihood indexing: pair (t, t+lag) uses the left-row copula parameter.
        theta_pair <- as.numeric(row1$theta_fit)
        zeta_pair <- as.numeric(row1$zeta_fit)
        tau_pair <- as.numeric(row1$tau_fit)
        if (!is.finite(theta_pair)) theta_pair <- NA_real_
        if (!is.finite(zeta_pair)) zeta_pair <- NA_real_
        if (!is.finite(tau_pair)) tau_pair <- NA_real_

        pair_list[[idx]] <- data.frame(
          subject = subject_id,
          time_left = as.character(t1),
          time_right = as.character(t2),
          time_pair = paste0("T", as.character(t1), " vs T", as.character(t2)),
          lag = lag_value,
          u1 = row1$u,
          u2 = row2$u,
          theta_pair = theta_pair,
          zeta_pair = zeta_pair,
          tau_fit = tau_pair,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }

  if (length(pair_list) == 0) {
    stop("No complete subject-time pairs were found for copula diagnostics.")
  }

  do.call(rbind, pair_list)
}

.copula_v2_attach_group <- function(pair_data, object, by, data = NULL) {
  if (is.null(by) || (is.character(by) && length(by) == 1 && !nzchar(by))) {
    pair_data$split_group <- factor(pair_data$time_pair)
    return(pair_data)
  }

  if (!is.character(by) || length(by) != 1) {
    stop("'by' must be NULL or a single column name as a character string.")
  }

  if (by %in% c("time", "time_pair")) {
    pair_data$split_group <- factor(pair_data$time_pair)
    return(pair_data)
  }
  if (by %in% c("subject", "lag") && by %in% names(pair_data)) {
    pair_data$split_group <- factor(pair_data[[by]])
    return(pair_data)
  }

  if (is.null(data)) {
    stop("To split plot.copula by '", by, "', provide data= containing that column.")
  }

  df <- as.data.frame(data, stringsAsFactors = FALSE)
  if (!is.null(object$var_map)) {
    for (old_name in names(object$var_map)) {
      new_name <- object$var_map[[old_name]]
      if (old_name %in% names(df) && !new_name %in% names(df)) {
        names(df)[names(df) == old_name] <- new_name
      }
    }
  }

  by_col <- by
  if (!by_col %in% names(df) && !is.null(object$var_map) && by %in% names(object$var_map)) {
    mapped_col <- object$var_map[[by]]
    if (mapped_col %in% names(df)) {
      by_col <- mapped_col
    }
  }
  if (!by_col %in% names(df)) {
    stop("Column '", by, "' not found in provided data after internal name mapping.")
  }
  if (!all(c("subject", "time") %in% names(df))) {
    stop("Provided data must contain subject and time columns (or names mappable via object$var_map) to split by '", by, "'.")
  }

  key_df <- paste(df$subject, as.character(df$time), sep = "::")
  key_pair <- paste(pair_data$subject, as.character(pair_data$time_left), sep = "::")
  matched <- df[[by_col]][match(key_pair, key_df)]

  if (all(is.na(matched))) {
    by_subj <- tapply(df[[by_col]], as.character(df$subject), function(v) {
      vv <- unique(v[!is.na(v)])
      if (length(vv) == 1) vv else NA
    })
    matched <- by_subj[as.character(pair_data$subject)]
  }

  pair_data$split_group <- factor(matched)
  pair_data <- pair_data[!is.na(pair_data$split_group), , drop = FALSE]

  if (nrow(pair_data) == 0) {
    stop("No valid paired rows remained after grouping by '", by, "'.")
  }

  pair_data
}

.copula_v2_transform_data <- function(data, transform = "uniform") {
  # Transform uniform [0,1] data to normal scale or other scales
  if (transform == "normal") {
    # Clamp to avoid infinite values from qnorm at 0 or 1
    data$u1 <- stats::qnorm(.copula_v2_clamp01(data$u1))
    data$u2 <- stats::qnorm(.copula_v2_clamp01(data$u2))
  }
  data
}

.copula_v2_average_density_grid <- function(family_num, pair_data, grid_n = 35, max_pairs_overlay = 300) {
  grid <- seq(0.02, 0.98, length.out = grid_n)
  grid_df <- expand.grid(u1 = grid, u2 = grid)

  pair_data <- pair_data[is.finite(pair_data$theta_pair), , drop = FALSE]
  if (nrow(pair_data) == 0) {
    grid_df$density <- NA_real_
    return(grid_df)
  }

  if (nrow(pair_data) > max_pairs_overlay) {
    set.seed(1)
    pair_data <- pair_data[sample(seq_len(nrow(pair_data)), max_pairs_overlay), , drop = FALSE]
  }

  density_sum <- rep(0, nrow(grid_df))
  density_count <- 0L

  for (i in seq_len(nrow(pair_data))) {
    par <- pair_data$theta_pair[i]
    par2 <- pair_data$zeta_pair[i]
    density_i <- tryCatch({
      if (is.finite(par2)) {
        VineCopula::BiCopPDF(grid_df$u1, grid_df$u2, family = family_num, par = par, par2 = par2)
      } else {
        VineCopula::BiCopPDF(grid_df$u1, grid_df$u2, family = family_num, par = par, par2 = 0)
      }
    }, error = function(e) rep(NA_real_, nrow(grid_df)))

    if (all(!is.finite(density_i))) next
    density_i[!is.finite(density_i)] <- 0
    density_sum <- density_sum + density_i
    density_count <- density_count + 1L
  }

  if (density_count == 0L) {
    grid_df$density <- NA_real_
  } else {
    grid_df$density <- density_sum / density_count
  }

  grid_df
}

.copula_v2_empirical_density_grid <- function(pair_data, grid_n = 35, lims = NULL) {
  if (is.null(lims)) {
    x_rng <- range(pair_data$u1, na.rm = TRUE)
    y_rng <- range(pair_data$u2, na.rm = TRUE)
    x_pad <- max(0.001, 0.025 * diff(x_rng))
    y_pad <- max(0.001, 0.025 * diff(y_rng))
    lims <- c(x_rng[1] - x_pad, x_rng[2] + x_pad, y_rng[1] - y_pad, y_rng[2] + y_pad)
  }
  kde <- MASS::kde2d(pair_data$u1, pair_data$u2, n = grid_n, lims = lims)
  data.frame(
    u1 = rep(kde$x, each = length(kde$y)),
    u2 = rep(kde$y, times = length(kde$x)),
    density = as.vector(kde$z),
    stringsAsFactors = FALSE
  )
}

.copula_v2_surface_metrics <- function(emp_density, fit_density, overlap_probs = c(0.7, 0.85, 0.95)) {
  emp <- as.numeric(emp_density)
  fit <- as.numeric(fit_density)
  ok <- is.finite(emp) & is.finite(fit)
  emp <- emp[ok]
  fit <- fit[ok]

  if (length(emp) == 0) {
    return(list(summary = data.frame(), overlap = data.frame()))
  }

  # Scale both surfaces to unit mass before computing distance metrics.
  emp <- pmax(emp, 0)
  fit <- pmax(fit, 0)
  emp <- emp / max(sum(emp), .Machine$double.eps)
  fit <- fit / max(sum(fit), .Machine$double.eps)

  summary_df <- data.frame(
    rmse = sqrt(mean((fit - emp)^2)),
    mae = mean(abs(fit - emp)),
    surface_cor = suppressWarnings(stats::cor(emp, fit, use = "complete.obs")),
    stringsAsFactors = FALSE
  )

  overlap_df <- do.call(rbind, lapply(overlap_probs, function(p) {
    thr_emp <- stats::quantile(emp, probs = p, na.rm = TRUE, type = 7)
    thr_fit <- stats::quantile(fit, probs = p, na.rm = TRUE, type = 7)
    mask_emp <- emp >= thr_emp
    mask_fit <- fit >= thr_fit
    union_n <- sum(mask_emp | mask_fit)
    iou <- if (union_n == 0) NA_real_ else sum(mask_emp & mask_fit) / union_n
    data.frame(level_prob = p, contour_iou = iou, stringsAsFactors = FALSE)
  }))

  list(summary = summary_df, overlap = overlap_df)
}

#' Compare fitted and empirical copula contour surfaces
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param lags Integer lags to assess, measured in ordered time steps.
#' @param grid_n Grid size used for density surfaces.
#' @param max_pairs_overlay Maximum number of paired observations used for fitted surface averaging.
#' @param contour_bins Number of contour levels to draw in the surface panels.
#' @param transform Character; "uniform" compares surfaces on copula scale, "normal" compares them on z-scale.
#' @param diff_scale_limit Positive numeric; fixed symmetric color scale limit for the difference panel.
#' @param time_stratified Logical; if TRUE, compare surfaces by time pair.
#' @param plot Logical; if TRUE, print the dashboard.
#'
#' @return Invisibly returns plots, grid-level surfaces, and numeric similarity metrics.
plot.copula_contour_compare <- function(object, lags = 1, grid_n = 45, max_pairs_overlay = 300, contour_bins = 10, transform = "uniform", diff_scale_limit = 0.05, time_stratified = FALSE, plot = TRUE, ...) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.")
  }

  if (!transform %in% c("uniform", "normal")) {
    stop("'transform' must be either 'uniform' or 'normal'.")
  }

  if (!is.numeric(diff_scale_limit) || length(diff_scale_limit) != 1 || !is.finite(diff_scale_limit) || diff_scale_limit <= 0) {
    stop("'diff_scale_limit' must be a single positive numeric value.")
  }

  fit_data <- .copula_v2_fit_data(object)
  pair_data <- .copula_v2_pair_data(fit_data, lags = lags)

  copula_spec <- get_copula_dist(object$copula_dist)
  copula_family_name <- object$copula_dist
  if (!is.character(copula_family_name) || nchar(copula_family_name) == 0) {
    copula_family_name <- copula_spec$copula_dist
  }
  if (grepl("^[0-9]+$", copula_family_name)) {
    family_map <- c("1" = "N", "2" = "C", "3" = "G", "4" = "F", "5" = "J", "6" = "BB1", "7" = "BB6", "8" = "BB7", "9" = "BB8", "10" = "T")
    if (copula_family_name %in% names(family_map)) {
      copula_family_name <- family_map[[copula_family_name]]
    }
  }

  family_num <- tryCatch({
    as.numeric(VineCopula::BiCopName(copula_family_name))
  }, error = function(e) NA_numeric_)

  split_data <- if (isTRUE(time_stratified)) split(pair_data, pair_data$time_pair) else list(All = pair_data)

  grid_list <- lapply(names(split_data), function(nm) {
    pd <- split_data[[nm]]

    fit_grid <- .copula_v2_average_density_grid(
      family_num = family_num,
      pair_data = pd,
      grid_n = grid_n,
      max_pairs_overlay = max_pairs_overlay
    )

    # Build empirical surface on the same copula grid as fit_grid, then transform both
    # together if requested. This avoids grid mismatch artifacts in contouring.
    emp_grid <- .copula_v2_empirical_density_grid(pd, grid_n = grid_n, lims = c(0.02, 0.98, 0.02, 0.98))

    if (transform == "normal") {
      z1 <- stats::qnorm(.copula_v2_clamp01(fit_grid$u1))
      z2 <- stats::qnorm(.copula_v2_clamp01(fit_grid$u2))
      jacobian <- stats::dnorm(z1) * stats::dnorm(z2)
      fit_grid$u1 <- z1
      fit_grid$u2 <- z2
      fit_grid$density <- fit_grid$density * jacobian

      emp_grid$u1 <- z1
      emp_grid$u2 <- z2
      emp_grid$density <- emp_grid$density * jacobian
    } else {
      emp_grid <- emp_grid
    }

    # Merge on grid coordinates to ensure pointwise comparisons.
    g <- merge(
      emp_grid,
      fit_grid,
      by = c("u1", "u2"),
      suffixes = c("_emp", "_fit"),
      all = FALSE
    )
    g$density_diff <- g$density_fit - g$density_emp
    g$time_pair <- nm
    g
  })

  grid_df <- do.call(rbind, grid_list)

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

  metric_summary <- do.call(rbind, lapply(metric_list, function(x) x$summary))
  metric_overlap <- do.call(rbind, lapply(metric_list, function(x) x$overlap))

  x_label <- if (transform == "normal") expression(Phi^-1 * (U[t])) else expression(U[t])
  y_label <- if (transform == "normal") expression(Phi^-1 * (U[t + 1])) else expression(U[t + 1])

  p_emp <- ggplot2::ggplot(grid_df, ggplot2::aes(x = u1, y = u2, z = density_emp)) +
    ggplot2::geom_contour(color = "#4d4d4d", bins = contour_bins, linewidth = 0.9) +
    ggplot2::labs(title = "Empirical Copula Contours", x = x_label, y = y_label) +
    ggplot2::theme_minimal()

  p_fit <- ggplot2::ggplot(grid_df, ggplot2::aes(x = u1, y = u2, z = density_fit)) +
    ggplot2::geom_contour(color = "#e41a1c", bins = contour_bins, linewidth = 0.9) +
    ggplot2::labs(title = "Fitted Copula Contours", x = x_label, y = y_label) +
    ggplot2::theme_minimal()

  p_diff <- ggplot2::ggplot(grid_df, ggplot2::aes(x = u1, y = u2, fill = density_diff)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-diff_scale_limit, diff_scale_limit),
      oob = scales::squish,
      name = "Fit - Emp"
    ) +
    ggplot2::labs(title = "Contour Difference Surface", x = x_label, y = y_label) +
    ggplot2::theme_minimal()

  if (isTRUE(time_stratified)) {
    p_emp <- p_emp + ggplot2::facet_wrap(~time_pair)
    p_fit <- p_fit + ggplot2::facet_wrap(~time_pair)
    p_diff <- p_diff + ggplot2::facet_wrap(~time_pair)
  }

  dashboard <- ggpubr::ggarrange(p_emp, p_fit, p_diff, ncol = 1, nrow = 3)

  if (isTRUE(plot)) {
    print(dashboard)
  }

  invisible(list(
    plots = list(empirical_contours = p_emp, fitted_contours = p_fit, difference_surface = p_diff),
    dashboard = dashboard,
    grid = grid_df,
    metrics = list(summary = metric_summary, overlap = metric_overlap)
  ))
}

#' Plot copula diagnostics for a fitted gamlss.longitudinal object
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param lags Integer lags to assess, measured in ordered time steps.
#' @param grid_n Grid size used for contour averaging.
#' @param max_pairs_overlay Maximum number of paired observations used for the fitted overlay.
#' @param transform Character; "uniform" (default) shows empirical copula on [0,1], "normal" transforms to standard normal scale.
#' @param plot1_style Character; "bins" (default) draws a binned empirical layer, "scatter" draws points.
#' @param contour_bins Integer number of contour levels for the fitted copula overlay in plot 1.
#' @param time_stratified Logical; if TRUE, facet both plots by time pair.
#' @param by Optional grouping variable name for stratified plots. Defaults to
#'   time-pair grouping when NULL. Use `data` for covariates not stored on the
#'   fitted pair object (for example gender).
#' @param data Optional data frame used when grouping by a covariate via `by`.
#' @param tau_ylim Optional numeric vector of length 2 specifying y-axis limits
#'   for Kendall's tau chart(s). If `NULL` (default), y-axis scales are automatic.
#' @param plot2_cuts Integer number of quantile-based cuts used in plot 2 (default 10).
#' @param plot Logical; if TRUE, print the dashboard.
#'
#' @return Invisibly returns a list with plot objects and summaries.
plot.copula <- function(object, lags = 1, grid_n = 35, max_pairs_overlay = 300, transform = "normal", plot1_style = "bins", contour_bins = 8, time_stratified = FALSE, by = NULL, data = NULL, tau_ylim = NULL, plot2_cuts = 10, plot = TRUE, ...) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.")
  }

  # Validate and apply transformation
  if (!transform %in% c("uniform", "normal")) {
    stop("'transform' must be either 'uniform' or 'normal'.")
  }

  if (!plot1_style %in% c("bins", "scatter")) {
    stop("'plot1_style' must be either 'bins' or 'scatter'.")
  }

  if (!is.numeric(contour_bins) || length(contour_bins) != 1 || !is.finite(contour_bins) || contour_bins < 1) {
    stop("'contour_bins' must be a single finite number >= 1.")
  }
  contour_bins <- as.integer(round(contour_bins))

  if (!is.logical(time_stratified) || length(time_stratified) != 1 || is.na(time_stratified)) {
    stop("'time_stratified' must be TRUE or FALSE.")
  }

  if (!is.numeric(plot2_cuts) || length(plot2_cuts) != 1 || !is.finite(plot2_cuts) || plot2_cuts < 2) {
    stop("'plot2_cuts' must be a single finite number >= 2.")
  }
  plot2_cuts <- as.integer(round(plot2_cuts))

  if (!is.null(tau_ylim)) {
    if (!is.numeric(tau_ylim) || length(tau_ylim) != 2 || any(!is.finite(tau_ylim)) || tau_ylim[1] >= tau_ylim[2]) {
      stop("'tau_ylim' must be NULL or a numeric vector of length 2 with tau_ylim[1] < tau_ylim[2].")
    }
    tau_ylim <- as.numeric(tau_ylim)
  }

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
  
  # Ensure proper copula family name for BiCopName
  copula_family_name <- object$copula_dist
  if (!is.character(copula_family_name) || nchar(copula_family_name) == 0) {
    copula_family_name <- copula_spec$copula_dist
  }
  
  # If still a numeric code, convert to family name
  if (grepl("^[0-9]+$", copula_family_name)) {
    # Map numeric codes back to family names
    family_map <- c("1" = "N", "2" = "C", "3" = "G", "4" = "F", "5" = "J", "6" = "BB1", "7" = "BB6", "8" = "BB7", "9" = "BB8", "10" = "T")
    if (copula_family_name %in% names(family_map)) {
      copula_family_name <- family_map[[copula_family_name]]
    }
  }

  family_num <- tryCatch({
    as.numeric(VineCopula::BiCopName(copula_family_name))
  }, error = function(e) NA_numeric_)

  is_grouped <- !is.null(by) || isTRUE(time_stratified)

  if (is_grouped) {
    density_list <- lapply(split(pair_data_uniform, pair_data_uniform$split_group), function(x) {
      grid_i <- .copula_v2_average_density_grid(
        family_num = family_num,
        pair_data = x,
        grid_n = grid_n,
        max_pairs_overlay = max_pairs_overlay
      )
      grid_i$split_group <- as.character(x$split_group[1])
      grid_i
    })
    density_grid <- do.call(rbind, density_list)
  } else {
    density_grid <- .copula_v2_average_density_grid(
      family_num = family_num,
      pair_data = pair_data_uniform,
      grid_n = grid_n,
      max_pairs_overlay = max_pairs_overlay
    )
  }

  # Apply transform to density grid if requested
  if (transform == "normal") {
    # Transform coordinates to normal scale
    z1 <- stats::qnorm(.copula_v2_clamp01(density_grid$u1))
    z2 <- stats::qnorm(.copula_v2_clamp01(density_grid$u2))
    
    # Apply Jacobian correction: multiply by phi(z1) * phi(z2)
    # where phi is the standard normal PDF
    jacobian_correction <- stats::dnorm(z1) * stats::dnorm(z2)
    
    density_grid$u1 <- z1
    density_grid$u2 <- z2
    density_grid$density <- density_grid$density * jacobian_correction
  }

  # Set axis labels based on transform
  x_label <- if (transform == "normal") {
    expression(Phi^-1 * (U[t]))
  } else {
    expression(U[t])
  }
  
  y_label <- if (transform == "normal") {
    expression(Phi^-1 * (U[t + 1]))
  } else {
    expression(U[t + 1])
  }

  p1 <- ggplot2::ggplot(pair_data_plot, ggplot2::aes_string(x = "u1", y = "u2"))

  if (plot1_style == "scatter") {
    p1 <- p1 +
      ggplot2::geom_point(color = "#4d4d4d", alpha = 0.45, size = 1.2)
  } else {
    p1 <- p1 +
      ggplot2::geom_bin2d(bins = 25, alpha = 0.8) +
      ggplot2::scale_fill_gradient(low = "white", high = "black", name = "Count")
  }

  p1 <- p1 +
    ggplot2::geom_contour(
      data = density_grid,
      ggplot2::aes(x = u1, y = u2, z = density),
      inherit.aes = FALSE,
      color = "#e41a1c",
      linewidth = 1.2,
      bins = contour_bins
    ) +
    ggplot2::labs(
      title = "Empirical Copula with Fitted Overlay",
      subtitle = paste0("Red contours show the average fitted copula across paired observations; family = ", copula_family_name),
      x = x_label,
      y = y_label
    ) +
    ggplot2::theme_minimal()

  if (is_grouped) {
    p1 <- p1 + ggplot2::facet_wrap(~split_group)
  }

  if (all(!is.finite(density_grid$density))) {
    p1 <- .copula_v2_message_plot(
      title = "Empirical Copula with Fitted Overlay",
      subtitle = paste0("Red contours show the average fitted copula across paired observations; family = ", copula_family_name),
      message = "No finite fitted copula density"
    )
  }

  build_cut_summary <- function(df, split_name = NULL) {
    if (nrow(df) < 1) {
      return(data.frame())
    }

    # Use rank-based bins to avoid collapsed quantile cuts when many fitted tau values are tied.
    df <- df[is.finite(df$tau_fit), , drop = FALSE]
    if (nrow(df) < 1) {
      return(data.frame())
    }

    effective_cuts <- min(plot2_cuts, nrow(df))
    cut_labels <- paste0("C", seq_len(effective_cuts))
    tau_rank <- rank(df$tau_fit, ties.method = "first", na.last = "keep")
    df$cut_group <- cut(tau_rank, breaks = effective_cuts, include.lowest = TRUE, labels = cut_labels)

    out <- do.call(rbind, lapply(split(df, df$cut_group), function(x) {
      tau_emp <- suppressWarnings(stats::cor(x$u1, x$u2, method = "kendall", use = "complete.obs"))
      tau_fit <- mean(x$tau_fit, na.rm = TRUE)
      data.frame(
        cut_group = as.character(x$cut_group[1]),
        tau_emp = tau_emp,
        tau_fit = tau_fit,
        n_pairs = nrow(x),
        stringsAsFactors = FALSE
      )
    }))

    if (!is.null(split_name)) {
      out$split_group <- split_name
    }
    out
  }

  if (is_grouped) {
    quartile_list <- lapply(split(pair_data_plot, pair_data_plot$split_group), function(x) {
      build_cut_summary(x, split_name = as.character(x$split_group[1]))
    })
    quartile_df <- do.call(rbind, quartile_list)
  } else {
    quartile_df <- build_cut_summary(pair_data_plot)
  }

  if (nrow(quartile_df) == 0 || all(!is.finite(quartile_df$tau_emp)) || all(!is.finite(quartile_df$tau_fit))) {
    p2 <- .copula_v2_message_plot(
      title = "Observed vs Fitted Correlation by Quantile Bin",
      subtitle = "Bins are formed from fitted copula strength",
      message = "No finite cut summaries"
    )
  } else {
    cut_levels <- paste0("C", sort(unique(as.integer(sub("^C", "", quartile_df$cut_group)))))
    quartile_df$cut_group <- factor(quartile_df$cut_group, levels = cut_levels)
    p2 <- ggplot2::ggplot(quartile_df, ggplot2::aes(x = cut_group)) +
      ggplot2::geom_point(ggplot2::aes(y = tau_emp), color = "#4d4d4d", size = 2.8) +
      ggplot2::geom_line(ggplot2::aes(y = tau_emp, group = 1), color = "#4d4d4d", linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(y = tau_fit), color = "#e41a1c", size = 2.8, shape = 4, stroke = 1.1) +
      ggplot2::geom_line(ggplot2::aes(y = tau_fit, group = 1), color = "#e41a1c", linewidth = 0.8, linetype = "dashed") +
      ggplot2::labs(
        title = "Observed vs Fitted Correlation by Quantile Bin",
        subtitle = paste0("", plot2_cuts, " cuts formed from fitted copula strength"),
        x = "Cut",
        y = "Kendall's tau"
      ) +
      ggplot2::theme_minimal()

    if (is_grouped) {
      p2 <- p2 + ggplot2::facet_wrap(~split_group, scales = if (is.null(tau_ylim)) "free_y" else "fixed")
    }

    if (!is.null(tau_ylim)) {
      p2 <- p2 + ggplot2::coord_cartesian(ylim = tau_ylim)
    }
  }

  dashboard <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2)

  if (isTRUE(plot)) {
    print(dashboard)
  }

  invisible(list(
    plots = list(
      empirical_overlay = p1,
      quartile_correlation = p2
    ),
    dashboard = dashboard,
    fit_data = fit_data,
    pair_data = pair_data_plot,
    quartile_summary = quartile_df
  ))
}
