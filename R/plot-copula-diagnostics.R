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

      suppressWarnings(.copula_par_to_tau(family = family_num, par = par, par2 = par2))

    } else {

      suppressWarnings(.copula_par_to_tau(family = family_num, par = par, par2 = 0))

    }

  }, error = function(e) NA_real_)


  if (is.finite(tau)) {

    return(as.numeric(tau))

  }


  # Fallback for Gaussian copulas using formula.

  if (identical(family_num, "N") && is.finite(par)) {

    return(2 / pi * asin(max(min(par, 0.999999), -0.999999)))

  }


  NA_real_

}


.copula_v2_bicop_cdf <- function(u1, u2, family_num, par, par2 = NA_real_) {

  n <- max(length(u1), length(u2), length(par), length(par2))

  if (!is.character(family_num) || length(family_num) != 1L || n < 1) return(rep(NA_real_, length(u1)))

  u1 <- rep(u1, length.out = n)

  u2 <- rep(u2, length.out = n)

  par <- rep(par, length.out = n)

  par2 <- rep(par2, length.out = n)

  vapply(seq_len(n), function(i) {

    if (!is.finite(u1[i]) || !is.finite(u2[i]) || !is.finite(par[i])) return(NA_real_)

    tryCatch({

      .copula_cdf(

        u1[i],

        u2[i],

        family = family_num,

        par = par[i],

        par2 = if (is.finite(par2[i])) par2[i] else 0

      )

    }, error = function(e) NA_real_)

  }, numeric(1), USE.NAMES = FALSE)

}


.copula_v2_bicop_cond_u2_given_u1 <- function(u1, u2, family_num, par, par2 = NA_real_) {

  n <- max(length(u1), length(u2), length(par), length(par2))

  if (!is.character(family_num) || length(family_num) != 1L || n < 1) return(rep(NA_real_, length(u1)))

  u1 <- rep(u1, length.out = n)

  u2 <- rep(u2, length.out = n)

  par <- rep(par, length.out = n)

  par2 <- rep(par2, length.out = n)

  out <- vapply(seq_len(n), function(i) {

    if (!is.finite(u1[i]) || !is.finite(u2[i]) || !is.finite(par[i])) return(NA_real_)

    tryCatch({

      # BiCopHfunc1 gives dC(u1, u2) / du1, i.e. F(U2 <= u2 | U1 = u1).

      .copula_hfunc1(

        u1[i],

        u2[i],

        family = family_num,

        par = par[i],

        par2 = if (is.finite(par2[i])) par2[i] else 0

      )

    }, error = function(e) NA_real_)

  }, numeric(1), USE.NAMES = FALSE)

  .copula_v2_clamp01(as.numeric(out))

}


.copula_v2_rosenblatt_pair_data <- function(pair_data, family_num) {

  pair_data$rosenblatt <- .copula_v2_bicop_cond_u2_given_u1(

    pair_data$u1,

    pair_data$u2,

    family_num = family_num,

    par = pair_data$theta_pair,

    par2 = pair_data$zeta_pair

  )

  pair_data$z <- stats::qnorm(.copula_v2_clamp01(pair_data$rosenblatt))

  pair_data$z_prev <- stats::qnorm(.copula_v2_clamp01(pair_data$u1))

  pair_data$z_curr <- pair_data$z

  pair_data

}


.copula_v2_rosenblatt_series <- function(fit_data, family_num) {

  pair_data <- .copula_v2_pair_data(fit_data, lags = 1)

  pair_data <- .copula_v2_rosenblatt_pair_data(pair_data, family_num)


  out <- fit_data[, c("subject", "time", "u"), drop = FALSE]

  out$key <- paste(out$subject, as.character(out$time), sep = "::")

  out$rosenblatt <- NA_real_


  first_idx <- ave(seq_len(nrow(out)), out$subject, FUN = function(x) x == min(x))

  out$rosenblatt[as.logical(first_idx)] <- out$u[as.logical(first_idx)]


  pair_key <- paste(pair_data$subject, as.character(pair_data$time_right), sep = "::")

  out$rosenblatt <- ifelse(

    is.na(out$rosenblatt),

    pair_data$rosenblatt[match(out$key, pair_key)],

    out$rosenblatt

  )

  out$rosenblatt <- .copula_v2_clamp01(out$rosenblatt)

  out$z <- stats::qnorm(out$rosenblatt)

  out[is.finite(out$z), c("subject", "time", "rosenblatt", "z"), drop = FALSE]

}


.copula_v2_kendall_diagnostic <- function(pair_data, family_num) {

  if (nrow(pair_data) < 2) return(data.frame())


  emp_copula <- vapply(seq_len(nrow(pair_data)), function(i) {

    mean(pair_data$u1 <= pair_data$u1[i] & pair_data$u2 <= pair_data$u2[i], na.rm = TRUE)

  }, numeric(1))


  fit_copula <- .copula_v2_bicop_cdf(

    pair_data$u1,

    pair_data$u2,

    family_num = family_num,

    par = pair_data$theta_pair,

    par2 = pair_data$zeta_pair

  )


  ok <- is.finite(emp_copula) & is.finite(fit_copula)

  data.frame(

    empirical = sort(emp_copula[ok]),

    fitted = sort(fit_copula[ok]),

    stringsAsFactors = FALSE

  )

}


.copula_v2_tail_diagnostics <- function(pair_data, family_num, thresholds = c(0.05, 0.10, 0.20)) {

  rows <- lapply(thresholds, function(alpha) {

    c_lower <- .copula_v2_bicop_cdf(

      rep(alpha, nrow(pair_data)),

      rep(alpha, nrow(pair_data)),

      family_num = family_num,

      par = pair_data$theta_pair,

      par2 = pair_data$zeta_pair

    )

    upper_cut <- 1 - alpha

    c_upper_cut <- .copula_v2_bicop_cdf(

      rep(upper_cut, nrow(pair_data)),

      rep(upper_cut, nrow(pair_data)),

      family_num = family_num,

      par = pair_data$theta_pair,

      par2 = pair_data$zeta_pair

    )

    lower_fit <- mean(c_lower, na.rm = TRUE)

    upper_fit <- mean(1 - 2 * upper_cut + c_upper_cut, na.rm = TRUE)


    data.frame(

      threshold = alpha,

      tail = c("Lower", "Upper"),

      empirical = c(

        mean(pair_data$u1 <= alpha & pair_data$u2 <= alpha, na.rm = TRUE),

        mean(pair_data$u1 >= upper_cut & pair_data$u2 >= upper_cut, na.rm = TRUE)

      ),

      fitted = c(lower_fit, upper_fit),

      stringsAsFactors = FALSE

    )

  })

  do.call(rbind, rows)

}


.copula_v2_conditional_tail_diagnostics <- function(tail_df) {

  if (nrow(tail_df) == 0) return(tail_df)

  out <- tail_df

  out$empirical <- out$empirical / out$threshold

  out$fitted <- out$fitted / out$threshold

  out$empirical <- pmin(pmax(out$empirical, 0), 1)

  out$fitted <- pmin(pmax(out$fitted, 0), 1)

  out

}


.copula_v2_rosenblatt_lag_summary <- function(rosenblatt_df, lag_values = 1:3) {

  lag_values <- sort(unique(as.integer(lag_values)))

  lag_values <- lag_values[lag_values > 0]

  rows <- lapply(lag_values, function(lag_value) {

    pair_list <- lapply(split(rosenblatt_df, rosenblatt_df$subject), function(x) {

      x <- x[order(x$time), , drop = FALSE]

      if (nrow(x) <= lag_value) return(NULL)

      data.frame(

        z_prev = x$z[seq_len(nrow(x) - lag_value)],

        z_curr = x$z[(lag_value + 1):nrow(x)]

      )

    })

    pair_list <- pair_list[!vapply(pair_list, is.null, logical(1))]

    if (length(pair_list) == 0) {

      return(data.frame(lag = lag_value, cor_z = NA_real_, n_pairs = 0L))

    }

    pairs <- do.call(rbind, pair_list)

    data.frame(

      lag = lag_value,

      cor_z = suppressWarnings(stats::cor(pairs$z_prev, pairs$z_curr, use = "complete.obs")),

      n_pairs = nrow(pairs)

    )

  })

  do.call(rbind, rows)

}


.copula_v2_message_plot <- function(title, subtitle, message) {

  ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x = x, y = y)) +

    ggplot2::geom_text(label = message, size = 4) +

    ggplot2::xlim(-1, 1) +

    ggplot2::ylim(-1, 1) +

    ggplot2::labs(title = title, subtitle = subtitle) +

    ggplot2::theme_void()

}


.copula_v2_fit_data <- function(object, data = NULL) {

  if (!inherits(object, "gamlss.longitudinal")) {

    stop("'object' must be a fitted 'gamlss.longitudinal' object.")

  }


  copula_spec <- get_copula_dist(object$copula_dist)


  copula_family_name <- .copula_family_code(copula_spec$copula_dist)


  if (is.null(data)) {

    mm_use <- object$model_matrix

    response <- object$response

    subject <- object$response_subject

    time <- object$response_margin

  } else {

    nd <- .gl_prepare_newdata_internal(object, data, require_response = TRUE)

    mm_use <- do.call(

      create_model_matrices,

      list(

        mu.formula = object$formulas_int$mu,

        sigma.formula = object$formulas_int$sigma,

        nu.formula = object$formulas_int$nu,

        tau.formula = object$formulas_int$tau,

        theta.formula = object$formulas_int$theta,

        zeta.formula = object$formulas_int$zeta,

        margin.family = object$margin_dist,

        copula.family = object$copula_dist,

        copula.link = copula_spec$copula_link,

        dataset = nd,

        quiet_gamlss2 = TRUE,

        preserve_factor_levels = TRUE

      )

    )

    mm_use <- .gl_align_model_matrix_columns(mm_use, object$model_matrix)

    response <- nd$response

    subject <- nd$subject

    time <- nd$time

  }


  eta_out <- calc_eta(

    par_cov = object$par,

    mm = mm_use,

    margin_dist = object$margin_dist,

    copula_link = copula_spec$copula_link,

    par_s = object$par_s

  )


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

    .copula_family_code(copula_family_name)

  }, error = function(e) NA_character_)


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

        .copula_pdf(grid_df$u1, grid_df$u2, family = family_num, par = par, par2 = par2)

      } else {

        .copula_pdf(grid_df$u1, grid_df$u2, family = family_num, par = par, par2 = 0)

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

#' @param x A fitted `gamlss.longitudinal` object.

#' @param lags Integer lags to assess, measured in ordered time steps.

#' @param grid_n Grid size used for density surfaces.

#' @param max_pairs_overlay Maximum number of paired observations used for fitted surface averaging.

#' @param contour_bins Number of contour levels to draw in the surface panels.

#' @param transform Character; "uniform" compares surfaces on copula scale, "normal" compares them on z-scale.

#' @param diff_scale_limit Positive numeric; fixed symmetric color scale limit for the difference panel.

#' @param time_stratified Logical; if TRUE, compare surfaces by time pair.

#' @param plot Logical; if TRUE, print the dashboard.

#' @param ... Additional arguments reserved for future methods.

#'

#' @return Invisibly returns plots, grid-level surfaces, and numeric similarity metrics.

#' @export

plot.copula_contour_compare <- function(x, lags = 1, grid_n = 45, max_pairs_overlay = 300, contour_bins = 10, transform = "uniform", diff_scale_limit = 0.05, time_stratified = FALSE, plot = TRUE, ...) {

  if (!inherits(x, "gamlss.longitudinal")) {

    stop("'x' must be a fitted 'gamlss.longitudinal' object.")

  }

  object <- x


  if (!transform %in% c("uniform", "normal")) {

    stop("'transform' must be either 'uniform' or 'normal'.")

  }


  if (!is.numeric(diff_scale_limit) || length(diff_scale_limit) != 1 || !is.finite(diff_scale_limit) || diff_scale_limit <= 0) {

    stop("'diff_scale_limit' must be a single positive numeric value.")

  }


  fit_data <- .copula_v2_fit_data(object)

  pair_data <- .copula_v2_pair_data(fit_data, lags = lags)


  copula_spec <- get_copula_dist(object$copula_dist)

  copula_family_name <- .copula_family_code(copula_spec$copula_dist)


  family_num <- tryCatch({

    .copula_family_code(copula_family_name)

  }, error = function(e) NA_character_)


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

#' @param x A fitted `gamlss.longitudinal` object.

#' @param lags Integer lags to assess, measured in ordered time steps.

#' @param grid_n Grid size used for contour averaging.

#' @param max_pairs_overlay Maximum number of paired observations used for the fitted overlay.

#' @param transform Character; "uniform" (default) shows empirical copula on

#'   the unit interval, "normal" transforms to standard normal scale.

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

#' @param tail_thresholds Numeric vector of lower-tail probabilities used for

#'   tail co-occurrence and conditional exceedance diagnostics.

#' @param residual_lags Integer lags used for Rosenblatt normal-score

#'   autocorrelation diagnostics.

#' @param dashboard_ncol Number of columns in the combined diagnostic dashboard.

#' @param plot Logical; if TRUE, print the dashboard.

#' @param ... Additional arguments reserved for future methods.

#'

#' @return Invisibly returns a list with plot objects and summaries.

#' @export

plot_copula_diagnostics <- function(

  x,

  lags = 1,

  grid_n = 35,

  max_pairs_overlay = 300,

  transform = "normal",

  plot1_style = "bins",

  contour_bins = 8,

  time_stratified = FALSE,

  by = NULL,

  data = NULL,

  tau_ylim = NULL,

  plot2_cuts = 10,

  tail_thresholds = c(0.05, 0.10, 0.20),

  residual_lags = 1:3,

  dashboard_ncol = 3,

  plot = TRUE,

  ...

) {

  plot.copula(

    x = x,

    lags = lags,

    grid_n = grid_n,

    max_pairs_overlay = max_pairs_overlay,

    transform = transform,

    plot1_style = plot1_style,

    contour_bins = contour_bins,

    time_stratified = time_stratified,

    by = by,

    data = data,

    tau_ylim = tau_ylim,

    plot2_cuts = plot2_cuts,

    tail_thresholds = tail_thresholds,

    residual_lags = residual_lags,

    dashboard_ncol = dashboard_ncol,

    plot = plot,

    ...

  )

}


#' @rdname plot_copula_diagnostics

#' @method plot copula

#' @export

plot.copula <- function(x, lags = 1, grid_n = 35, max_pairs_overlay = 300, transform = "normal", plot1_style = "bins", contour_bins = 8, time_stratified = FALSE, by = NULL, data = NULL, tau_ylim = NULL, plot2_cuts = 10, tail_thresholds = c(0.05, 0.10, 0.20), residual_lags = 1:3, dashboard_ncol = 3, plot = TRUE, ...) {

  if (!inherits(x, "gamlss.longitudinal")) {

    stop("'x' must be a fitted 'gamlss.longitudinal' object.")

  }

  object <- x


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


  tail_thresholds <- sort(unique(as.numeric(tail_thresholds)))

  tail_thresholds <- tail_thresholds[is.finite(tail_thresholds) & tail_thresholds > 0 & tail_thresholds < 0.5]

  if (length(tail_thresholds) == 0) {

    tail_thresholds <- c(0.05, 0.10, 0.20)

  }


  residual_lags <- sort(unique(as.integer(residual_lags)))

  residual_lags <- residual_lags[residual_lags > 0]

  if (length(residual_lags) == 0) {

    residual_lags <- 1:3

  }


  if (!is.numeric(dashboard_ncol) || length(dashboard_ncol) != 1 || !is.finite(dashboard_ncol) || dashboard_ncol < 1) {

    stop("'dashboard_ncol' must be a single positive integer.")

  }

  dashboard_ncol <- as.integer(round(dashboard_ncol))


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


  family_num <- tryCatch({

    .copula_family_code(copula_family_name)

  }, error = function(e) NA_character_)


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


  rosenblatt_df <- tryCatch(

    .copula_v2_rosenblatt_series(fit_data, family_num),

    error = function(e) data.frame()

  )


  rosenblatt_pair_df <- tryCatch(

    .copula_v2_rosenblatt_pair_data(pair_data_uniform, family_num),

    error = function(e) data.frame()

  )


  if (nrow(rosenblatt_df) == 0 || all(!is.finite(rosenblatt_df$z))) {

    p_ros_time <- .copula_v2_message_plot(

      title = "Rosenblatt Normal Scores by Time",

      subtitle = "Scores are qnorm of pairwise conditional Rosenblatt residuals",

      message = "No finite Rosenblatt residuals"

    )

  } else {

    p_ros_time <- ggplot2::ggplot(rosenblatt_df, ggplot2::aes(x = factor(time), y = z)) +

      ggplot2::geom_hline(yintercept = 0, color = "#666666", linetype = "dashed") +

      ggplot2::geom_boxplot(fill = "#9ecae1", color = "#4d4d4d", outlier.alpha = 0.35) +

      ggplot2::labs(

        title = "Rosenblatt Normal Scores by Time",

        subtitle = "Each time point should be centered near zero with similar spread",

        x = "Time",

        y = "Normal score"

      ) +

      ggplot2::theme_minimal()

  }


  rosenblatt_z <- rosenblatt_df$z[is.finite(rosenblatt_df$z)]

  if (length(rosenblatt_z) == 0) {

    p_ros_qq <- .copula_v2_message_plot(

      title = "Rosenblatt Normal QQ",

      subtitle = "Conditional copula scores should follow N(0, 1)",

      message = "No finite Rosenblatt residuals"

    )

  } else {

    rosenblatt_qq_df <- data.frame(

      theoretical = stats::qnorm(stats::ppoints(length(rosenblatt_z))),

      observed = sort(rosenblatt_z),

      stringsAsFactors = FALSE

    )

    p_ros_qq <- ggplot2::ggplot(rosenblatt_qq_df, ggplot2::aes(x = theoretical, y = observed)) +

      ggplot2::geom_abline(intercept = 0, slope = 1, color = "#666666", linetype = "dashed") +

      ggplot2::geom_point(color = "#4d4d4d", alpha = 0.55, size = 1.2) +

      ggplot2::labs(

        title = "Rosenblatt Normal QQ",

        subtitle = "Conditional copula scores should follow N(0, 1)",

        x = "Theoretical normal quantile",

        y = "Observed normal score"

      ) +

      ggplot2::theme_minimal()

  }


  if (nrow(rosenblatt_pair_df) == 0 || all(!is.finite(rosenblatt_pair_df$z_prev)) || all(!is.finite(rosenblatt_pair_df$z_curr))) {

    p_ros_lag <- .copula_v2_message_plot(

      title = "Rosenblatt Lag Plot",

      subtitle = "Current conditional score against previous marginal score",

      message = "No finite Rosenblatt lag pairs"

    )

  } else {

    p_ros_lag <- ggplot2::ggplot(rosenblatt_pair_df, ggplot2::aes(x = z_prev, y = z_curr)) +

      ggplot2::geom_hline(yintercept = 0, color = "#d9d9d9") +

      ggplot2::geom_vline(xintercept = 0, color = "#d9d9d9") +

      ggplot2::geom_point(color = "#4d4d4d", alpha = 0.35, size = 1.1) +

      ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#e41a1c", linewidth = 0.7) +

      ggplot2::labs(

        title = "Rosenblatt Lag Plot",

        subtitle = "The smooth should be approximately flat at zero",

        x = expression(Phi^-1 * (U[t])),

        y = expression(Phi^-1 * (R[t + 1] ~ "|" ~ U[t]))

      ) +

      ggplot2::theme_minimal()

  }


  kendall_df <- tryCatch(

    .copula_v2_kendall_diagnostic(pair_data_uniform, family_num),

    error = function(e) data.frame()

  )


  if (nrow(kendall_df) == 0) {

    p_kendall <- .copula_v2_message_plot(

      title = "Kendall Function Diagnostic",

      subtitle = "Empirical copula values compared with fitted copula values at observed pairs",

      message = "No finite Kendall diagnostic values"

    )

  } else {

    p_kendall <- ggplot2::ggplot(kendall_df, ggplot2::aes(x = fitted, y = empirical)) +

      ggplot2::geom_abline(intercept = 0, slope = 1, color = "#666666", linetype = "dashed") +

      ggplot2::geom_point(color = "#4d4d4d", alpha = 0.55, size = 1.2) +

      ggplot2::labs(

        title = "Kendall Function Diagnostic",

        subtitle = "Sorted empirical copula probabilities should track sorted fitted probabilities",

        x = "Fitted copula probability",

        y = "Empirical copula probability"

      ) +

      ggplot2::theme_minimal()

  }


  tail_df <- tryCatch(

    .copula_v2_tail_diagnostics(pair_data_uniform, family_num, thresholds = tail_thresholds),

    error = function(e) data.frame()

  )

  cond_tail_df <- .copula_v2_conditional_tail_diagnostics(tail_df)


  tail_long <- if (nrow(tail_df) > 0) {

    rbind(

      data.frame(threshold = tail_df$threshold, tail = tail_df$tail, source = "Empirical", probability = tail_df$empirical),

      data.frame(threshold = tail_df$threshold, tail = tail_df$tail, source = "Fitted", probability = tail_df$fitted)

    )

  } else {

    data.frame()

  }


  if (nrow(tail_long) == 0 || all(!is.finite(tail_long$probability))) {

    p_tail <- .copula_v2_message_plot(

      title = "Tail Co-occurrence",

      subtitle = "Observed joint tail probability against fitted copula probability",

      message = "No finite tail diagnostics"

    )

  } else {

    p_tail <- ggplot2::ggplot(tail_long, ggplot2::aes(x = threshold, y = probability, color = source, group = source)) +

      ggplot2::geom_point(size = 2.4) +

      ggplot2::geom_line(linewidth = 0.8) +

      ggplot2::facet_wrap(~tail) +

      ggplot2::scale_color_manual(values = c(Empirical = "#4d4d4d", Fitted = "#e41a1c")) +

      ggplot2::labs(

        title = "Tail Co-occurrence",

        subtitle = "Lower: P(Ut <= a, Ut+1 <= a); Upper: P(Ut >= 1-a, Ut+1 >= 1-a)",

        x = "Tail probability a",

        y = "Joint probability",

        color = NULL

      ) +

      ggplot2::theme_minimal()

  }


  cond_tail_long <- if (nrow(cond_tail_df) > 0) {

    rbind(

      data.frame(threshold = cond_tail_df$threshold, tail = cond_tail_df$tail, source = "Empirical", probability = cond_tail_df$empirical),

      data.frame(threshold = cond_tail_df$threshold, tail = cond_tail_df$tail, source = "Fitted", probability = cond_tail_df$fitted)

    )

  } else {

    data.frame()

  }


  if (nrow(cond_tail_long) == 0 || all(!is.finite(cond_tail_long$probability))) {

    p_cond_tail <- .copula_v2_message_plot(

      title = "Conditional Tail Exceedance",

      subtitle = "Observed conditional tail probability against fitted copula probability",

      message = "No finite conditional tail diagnostics"

    )

  } else {

    p_cond_tail <- ggplot2::ggplot(cond_tail_long, ggplot2::aes(x = threshold, y = probability, color = source, group = source)) +

      ggplot2::geom_point(size = 2.4) +

      ggplot2::geom_line(linewidth = 0.8) +

      ggplot2::facet_wrap(~tail) +

      ggplot2::scale_color_manual(values = c(Empirical = "#4d4d4d", Fitted = "#e41a1c")) +

      ggplot2::coord_cartesian(ylim = c(0, 1)) +

      ggplot2::labs(

        title = "Conditional Tail Exceedance",

        subtitle = "Lower: P(Ut+1 <= a | Ut <= a); Upper: P(Ut+1 >= 1-a | Ut >= 1-a)",

        x = "Tail probability a",

        y = "Conditional probability",

        color = NULL

      ) +

      ggplot2::theme_minimal()

  }


  lag_summary_df <- tryCatch(

    .copula_v2_rosenblatt_lag_summary(rosenblatt_df, lag_values = residual_lags),

    error = function(e) data.frame()

  )


  if (nrow(lag_summary_df) == 0 || all(!is.finite(lag_summary_df$cor_z))) {

    p_lag_summary <- .copula_v2_message_plot(

      title = "Residual Dependence by Lag",

      subtitle = "Correlation of Rosenblatt normal scores within subject",

      message = "No finite residual lag correlations"

    )

  } else {

    p_lag_summary <- ggplot2::ggplot(lag_summary_df, ggplot2::aes(x = factor(lag), y = cor_z)) +

      ggplot2::geom_hline(yintercept = 0, color = "#666666", linetype = "dashed") +

      ggplot2::geom_col(fill = "#4d4d4d", alpha = 0.8) +

      ggplot2::geom_text(ggplot2::aes(label = paste0("n=", n_pairs)), vjust = -0.35, size = 3) +

      ggplot2::coord_cartesian(ylim = c(-1, 1)) +

      ggplot2::labs(

        title = "Residual Dependence by Lag",

        subtitle = "Correlations should be close to zero after the Rosenblatt transform",

        x = "Lag",

        y = "Correlation"

      ) +

      ggplot2::theme_minimal()

  }


  dashboard_plots <- list(p1, p2, p_ros_time, p_ros_qq, p_ros_lag, p_kendall, p_tail, p_cond_tail, p_lag_summary)

  dashboard <- do.call(

    ggpubr::ggarrange,

    c(

      dashboard_plots,

      list(

        ncol = min(dashboard_ncol, length(dashboard_plots)),

        nrow = ceiling(length(dashboard_plots) / dashboard_ncol)

      )

    )

  )


  if (isTRUE(plot)) {

    print(dashboard)

  }


  invisible(list(

    plots = list(

      empirical_overlay = p1,

      quartile_correlation = p2,

      rosenblatt_by_time = p_ros_time,

      rosenblatt_qq = p_ros_qq,

      rosenblatt_lag = p_ros_lag,

      kendall_function = p_kendall,

      tail_cooccurrence = p_tail,

      conditional_tail_exceedance = p_cond_tail,

      residual_lag_correlation = p_lag_summary

    ),

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


