# Silence NSE checks for ggplot mappings.
utils::globalVariables(c("theoretical", "observed", "detrended", "midpoint", "root_diff", "band_lower", "band_upper", "time", "split_group"))

#' Diagnostic generics for fitted longitudinal GAMLSS-copula models
#'
#' These generics dispatch to methods for `gamlss.longitudinal` objects and
#' return either diagnostic plot objects, scoring summaries, or forecast data.
#'
#' @name topmodels_diagnostics
#' @aliases pithist qqrplot wormplot rootogram proscore procast
NULL

#' @export
pithist <- function(object, ...) {
  UseMethod("pithist")
}

#' @export
qqrplot <- function(object, ...) {
  UseMethod("qqrplot")
}

#' @export
wormplot <- function(object, ...) {
  UseMethod("wormplot")
}

#' @export
rootogram <- function(object, ...) {
  UseMethod("rootogram")
}

#' @export
proscore <- function(object, ...) {
  UseMethod("proscore")
}

#' @export
procast <- function(object, ...) {
  UseMethod("procast")
}

.gl_get_family_fun <- function(family_name, prefix) {
  fun_name <- paste0(prefix, family_name)
  if (exists(fun_name, envir = asNamespace("gamlss.dist"), mode = "function", inherits = FALSE)) {
    return(get(fun_name, envir = asNamespace("gamlss.dist"), mode = "function", inherits = FALSE))
  }
  if (exists(fun_name, mode = "function")) {
    return(get(fun_name, mode = "function"))
  }
  stop(
    "Distribution function '", fun_name,
    "' is not available in gamlss.dist or the current session."
  )
}

.gl_call_family_fun <- function(prefix, family_name, x, params, extra_args = list()) {
  fun <- .gl_get_family_fun(family_name, prefix)
  arg_name <- switch(prefix, p = "q", d = "x", q = "p", "x")
  args <- c(stats::setNames(list(x), arg_name), extra_args, params)
  args <- args[names(args) %in% formalArgs(fun)]
  do.call(fun, args)
}

.gl_diag_data <- function(object) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("Diagnostics are only available for objects of class 'gamlss.longitudinal'.")
  }

  copula_link <- get_copula_dist(object$copula_dist)$copula_link
  eta_out <- calc_eta(
    par_cov = object$par,
    mm = object$model_matrix,
    margin_dist = object$margin_dist,
    copula_link = copula_link,
    par_s = object$par_s
  )

  margin_params <- eta_out$eta_inv[names(object$margin_dist$parameters)]
  y <- object$response

  keep <- is.finite(y)
  for (par_name in names(margin_params)) {
    keep <- keep & is.finite(margin_params[[par_name]])
  }

  margin_params <- lapply(margin_params, function(x) x[keep])
  y <- y[keep]

  mu_hat <- if ("mu" %in% names(margin_params)) margin_params$mu else margin_params[[1]]
  sigma_hat <- if ("sigma" %in% names(margin_params)) margin_params$sigma else rep(stats::sd(y, na.rm = TRUE), length(y))

  sigma_hat <- pmax(as.numeric(sigma_hat), .Machine$double.eps)

  list(
    response = y,
    params = margin_params,
    mu_hat = as.numeric(mu_hat),
    sigma_hat = sigma_hat,
    family = object$margin_dist$family[1],
    subject = object$response_subject[keep],
    time = object$response_margin[keep]
  )
}

.gl_prepare_newdata_internal <- function(object, newdata, require_response = FALSE) {
  if (is.null(newdata)) return(NULL)

  if (is.null(object$formulas_int) || is.null(object$var_map)) {
    stop("newdata prediction requires a model fit created with stored formulas/variable map. Refit with the current package version.")
  }

  nd <- as.data.frame(newdata, stringsAsFactors = FALSE)

  # Translate user variable names to internal names used by model formulas.
  for (old_name in names(object$var_map)) {
    new_name <- object$var_map[[old_name]]
    if (old_name %in% names(nd) && !new_name %in% names(nd)) {
      names(nd)[names(nd) == old_name] <- new_name
    }
  }

  # Fitting keeps the original user time scale as `time_covariate` for
  # formulas, while `time` is used internally for ordering/pairing. Recreate
  # that column for prediction data supplied with either original or internal
  # names.
  if (!"time_covariate" %in% names(nd) && "time" %in% names(nd)) {
    nd$time_covariate <- nd$time
  }

  if (!"time" %in% names(nd) && "time" %in% names(object$model_matrix$x$mu)) {
    nd$time <- NA
  }
  if (!"subject" %in% names(nd)) {
    nd$subject <- seq_len(nrow(nd))
  }
  if (!"response" %in% names(nd)) {
    nd$response <- NA_real_
  }

  if (!is.null(object$dataset)) {
    factor_cols <- names(object$dataset)[vapply(object$dataset, is.factor, logical(1))]
    for (nm in intersect(factor_cols, names(nd))) {
      train_col <- object$dataset[[nm]]
      train_levels <- levels(train_col)
      nd_values <- as.character(nd[[nm]])
      unknown <- setdiff(unique(nd_values[!is.na(nd_values)]), train_levels)
      if (length(unknown) > 0L) {
        stop(
          "newdata column '", nm, "' contains level(s) not seen during fitting: ",
          paste(unknown, collapse = ", "),
          call. = FALSE
        )
      }
      nd[[nm]] <- factor(nd_values, levels = train_levels, ordered = is.ordered(train_col))
      if (is.ordered(train_col) && length(train_levels) > 1L) {
        contr <- contr.treatment(length(train_levels))
        colnames(contr) <- train_levels[-1]
        contrasts(nd[[nm]]) <- contr
      }
    }
  }

  if (require_response && all(is.na(nd$response))) {
    stop("newdata must include a response column (or mapped response variable) for this operation.")
  }

  nd
}

.gl_align_model_matrix_columns <- function(mm_use, mm_reference) {
  if (is.null(mm_use) || is.null(mm_reference)) {
    return(mm_use)
  }

  for (par_name in intersect(names(mm_reference$x), names(mm_use$x))) {
    ref_cols <- colnames(mm_reference$x[[par_name]])
    use_cols <- colnames(mm_use$x[[par_name]])
    missing_cols <- setdiff(ref_cols, use_cols)

    if (length(missing_cols) > 0L) {
      for (col_name in missing_cols) {
        mm_use$x[[par_name]][[col_name]] <- 0
      }
    }

    extra_cols <- setdiff(colnames(mm_use$x[[par_name]]), ref_cols)
    if (length(extra_cols) > 0L) {
      mm_use$x[[par_name]] <- mm_use$x[[par_name]][
        ,
        setdiff(colnames(mm_use$x[[par_name]]), extra_cols),
        drop = FALSE
      ]
    }

    mm_use$x[[par_name]] <- mm_use$x[[par_name]][, ref_cols, drop = FALSE]
  }

  mm_use
}

.gl_fitted_distribution <- function(object, newdata = NULL, require_response = TRUE) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("Diagnostics are only available for objects of class 'gamlss.longitudinal'.")
  }

  copula_link <- get_copula_dist(object$copula_dist)$copula_link

  if (is.null(newdata)) {
    mm_use <- object$model_matrix
    response <- object$response
    response_margin <- object$response_margin
    response_subject <- object$response_subject
  } else {
    nd <- .gl_prepare_newdata_internal(object, newdata, require_response = require_response)

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
        copula.link = copula_link,
        dataset = nd,
        quiet_gamlss2 = TRUE,
        preserve_factor_levels = TRUE
      )
    )
    mm_use <- .gl_align_model_matrix_columns(mm_use, object$model_matrix)

    response <- nd$response
    response_margin <- nd$time
    response_subject <- nd$subject
  }

  eta_out <- calc_eta(
    par_cov = object$par,
    mm = mm_use,
    margin_dist = object$margin_dist,
    copula_link = copula_link,
    par_s = object$par_s
  )

  margin_params <- eta_out$eta_inv[names(object$margin_dist$parameters)]

  keep <- rep(TRUE, length(response))
  if (require_response) {
    keep <- is.finite(response)
  }
  for (par_name in names(margin_params)) {
    keep <- keep & is.finite(margin_params[[par_name]])
  }

  common_n <- min(
    length(response),
    length(response_margin),
    length(response_subject),
    if (length(margin_params) > 0) min(vapply(margin_params, length, integer(1))) else length(response)
  )
  if (!is.finite(common_n) || common_n < 0) {
    common_n <- 0L
  }

  response <- response[seq_len(common_n)]
  response_margin <- response_margin[seq_len(common_n)]
  response_subject <- response_subject[seq_len(common_n)]
  margin_params <- lapply(margin_params, function(x) x[seq_len(common_n)])
  keep <- keep[seq_len(common_n)]

  margin_params <- lapply(margin_params, function(x) x[keep])
  response <- response[keep]
  response_margin <- response_margin[keep]
  response_subject <- response_subject[keep]

  mu_hat <- if ("mu" %in% names(margin_params)) margin_params$mu else margin_params[[1]]
  sigma_hat <- if ("sigma" %in% names(margin_params)) margin_params$sigma else rep(stats::sd(response, na.rm = TRUE), length(response))
  sigma_hat <- pmax(as.numeric(sigma_hat), .Machine$double.eps)

  list(
    response = response,
    params = margin_params,
    mu_hat = as.numeric(mu_hat),
    sigma_hat = sigma_hat,
    family = object$margin_dist$family[1],
    subject = response_subject,
    time = response_margin,
    keep_mask = keep,
    keep_index = which(keep)
  )
}

.gl_pit <- function(object, randomize = FALSE) {
  diag_data <- .gl_fitted_distribution(object, newdata = NULL, require_response = TRUE)
  y <- diag_data$response
  params <- diag_data$params

  pit_upper <- .gl_call_family_fun("p", diag_data$family, y, params)
  pit <- pit_upper

  if (randomize && .is_discrete_margin(object$margin_dist)) {
    pit_lower <- .gl_call_family_fun("p", diag_data$family, y - 1, params)
    pit_lower <- pmin(pmax(as.numeric(pit_lower), 0), 1)
    pit_upper <- pmin(pmax(as.numeric(pit_upper), 0), 1)
    interval_width <- pmax(pit_upper - pit_lower, 0)
    pit <- pit_lower + stats::runif(length(pit_upper)) * interval_width
  }

  pit <- pmin(pmax(as.numeric(pit), 0), 1)
  list(diag = diag_data, pit = pit)
}

#' @export
pithist.gamlss.longitudinal <- function(object, bins = 20, randomize = FALSE, plot = TRUE, by_time = FALSE, ...) {
  pit_out <- .gl_pit(object, randomize = randomize)
  pit <- pit_out$pit
  pit_df <- data.frame(pit = pit, time = as.factor(pit_out$diag$time))
  pit_df <- pit_df[is.finite(pit_df$pit), , drop = FALSE]

  if (!plot) {
    return(pit_df)
  }

  if (!by_time) {
    expected <- nrow(pit_df) / bins
    return(
      ggplot2::ggplot(pit_df, ggplot2::aes(x = pit)) +
        ggplot2::geom_histogram(breaks = seq(0, 1, length.out = bins + 1L), closed = "right", fill = "#2c7fb8", color = "white") +
        ggplot2::geom_hline(yintercept = expected, linetype = "dashed", linewidth = 0.4, color = "#444444") +
        ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +
        ggplot2::labs(x = "PIT", y = "Count", title = "PIT Histogram") +
        ggplot2::theme_minimal()
    )
  }

  expected_df <- aggregate(pit ~ time, data = pit_df, FUN = function(x) length(x) / bins)
  names(expected_df)[2] <- "expected"

  ggplot2::ggplot(pit_df, ggplot2::aes(x = pit)) +
    ggplot2::geom_histogram(breaks = seq(0, 1, length.out = bins + 1L), closed = "right", fill = "#2c7fb8", color = "white") +
    ggplot2::geom_hline(data = expected_df, ggplot2::aes(yintercept = expected), linetype = "dashed", linewidth = 0.4, color = "#444444") +
    ggplot2::facet_wrap(~time, scales = "free_y") +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +
    ggplot2::labs(x = "PIT", y = "Count", title = "PIT Histogram by Time") +
    ggplot2::theme_minimal()
}

#' @export
qqrplot.gamlss.longitudinal <- function(object, randomize = FALSE, plot = TRUE, by_time = FALSE, ...) {
  pit_out <- .gl_pit(object, randomize = randomize)
  z <- stats::qnorm(pmin(pmax(pit_out$pit, .Machine$double.eps), 1 - .Machine$double.eps))

  if (!by_time) {
    theo <- stats::qnorm(stats::ppoints(length(z)))
    qq_df <- data.frame(theoretical = theo, observed = sort(z))
  } else {
    split_z <- split(z, pit_out$diag$time)
    qq_list <- lapply(names(split_z), function(ti) {
      z_t <- split_z[[ti]]
      theo_t <- stats::qnorm(stats::ppoints(length(z_t)))
      data.frame(theoretical = theo_t, observed = sort(z_t), time = as.factor(ti))
    })
    qq_df <- do.call(rbind, qq_list)
  }

  if (!plot) {
    return(qq_df)
  }

  theoretical <- observed <- NULL
  p <- ggplot2::ggplot(qq_df, ggplot2::aes(x = theoretical, y = observed)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.7) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#666666") +
    ggplot2::labs(
      x = "Theoretical Normal Quantiles",
      y = "Randomized Quantiles",
      title = if (by_time) "QQ Plot by Time" else "QQ Plot"
    ) +
    ggplot2::theme_minimal()

  if (by_time) {
    p <- p + ggplot2::facet_wrap(~time, scales = "free")
  }

  p
}

#' @export
wormplot.gamlss.longitudinal <- function(object, randomize = FALSE, plot = TRUE, smooth = TRUE, by_time = FALSE, by = NULL, data = NULL, ...) {
  pit_out <- .gl_pit(object, randomize = randomize)

  resolve_split <- function(by, pit_out, object, data = NULL) {
    if (is.null(by) || (is.character(by) && length(by) == 1 && !nzchar(by))) {
      return(NULL)
    }
    if (!is.character(by) || length(by) != 1) {
      stop("'by' must be NULL or a single column name as a character string.")
    }

    if (by %in% c("time", "response_margin")) {
      return(as.factor(pit_out$diag$time))
    }
    if (by %in% c("subject", "response_subject")) {
      return(as.factor(pit_out$diag$subject))
    }

    if (is.null(data)) {
      stop("To split wormplot by '", by, "', provide data= containing that column.")
    }
    if (!is.data.frame(data)) {
      data <- as.data.frame(data, stringsAsFactors = FALSE)
    }
    if (!by %in% names(data)) {
      stop("Column '", by, "' not found in provided data.")
    }

    n_pit <- length(pit_out$pit)
    if (nrow(data) == n_pit) {
      return(as.factor(data[[by]]))
    }

    keep_mask <- pit_out$diag$keep_mask
    if (!is.null(keep_mask) && length(keep_mask) == nrow(data)) {
      vec <- data[[by]][keep_mask]
      if (length(vec) == n_pit) {
        return(as.factor(vec))
      }
    }

    keep_index <- pit_out$diag$keep_index
    if (!is.null(keep_index) && length(keep_index) == n_pit && max(keep_index) <= nrow(data)) {
      return(as.factor(data[[by]][keep_index]))
    }

    stop(
      "Could not align data rows with wormplot residual rows for by='", by, "'. ",
      "Provide data with row count equal to either length(pit) (", n_pit, ") or the original fit data rows."
    )
  }

  worm_band_frame <- function(theoretical, n, band_level = 0.95) {
    p <- stats::pnorm(theoretical)
    se <- sqrt(pmax(0, p * (1 - p) / n)) / stats::dnorm(theoretical)
    z <- stats::qnorm((1 + band_level) / 2)
    data.frame(
      theoretical = theoretical,
      band_lower = -z * se,
      band_upper = z * se
    )
  }

  if (isTRUE(by_time) && is.null(by)) {
    by <- "time"
  } else if (isTRUE(by_time) && !is.null(by)) {
    warning("Both by_time and by were provided; using by='", by, "'.", call. = FALSE)
  }

  split_group <- resolve_split(by, pit_out, object, data)
  split_by <- !is.null(split_group)

  if (!split_by) {
    z <- stats::qnorm(pmin(pmax(pit_out$pit, .Machine$double.eps), 1 - .Machine$double.eps))
    theo <- stats::qnorm(stats::ppoints(length(z)))
    worm_df <- data.frame(theoretical = theo, detrended = sort(z) - theo)
    worm_band <- worm_band_frame(theo, length(z))
  } else {
    split_pit <- split(pit_out$pit, split_group)
    worm_list <- lapply(names(split_pit), function(grp) {
      pit_t <- split_pit[[grp]]
      z_t <- stats::qnorm(pmin(pmax(pit_t, .Machine$double.eps), 1 - .Machine$double.eps))
      theo_t <- stats::qnorm(stats::ppoints(length(z_t)))
      band_t <- worm_band_frame(theo_t, length(z_t))
      data.frame(
        theoretical = theo_t,
        detrended = sort(z_t) - theo_t,
        band_lower = band_t$band_lower,
        band_upper = band_t$band_upper,
        split_group = as.factor(grp)
      )
    })
    worm_df <- do.call(rbind, worm_list)
  }

  if (!plot) {
    return(worm_df)
  }

  theoretical <- detrended <- band_lower <- band_upper <- NULL
  p <- ggplot2::ggplot(worm_df, ggplot2::aes(x = theoretical, y = detrended)) +
    ggplot2::geom_ribbon(
      data = if (split_by) worm_df else worm_band,
      ggplot2::aes(x = theoretical, ymin = band_lower, ymax = band_upper),
      inherit.aes = FALSE,
      fill = "#9ecae1",
      alpha = 0.25
    ) +
    ggplot2::geom_point(size = 1.2, alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#666666") +
    ggplot2::labs(
      x = "Theoretical Normal Quantiles",
      y = "Detrended Quantiles",
      title = if (split_by) paste0("Worm Plot by ", by) else "Worm Plot"
    ) +
    ggplot2::theme_minimal()

  if (split_by) {
    p <- p + ggplot2::facet_wrap(~split_group, scales = "free")
  }

  if (smooth && nrow(worm_df) > 5) {
    p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#d95f0e", linewidth = 0.7)
  }

  p
}

#' @export
rootogram.gamlss.longitudinal <- function(object, bins = 20, plot = TRUE, ...) {
  diag_data <- .gl_diag_data(object)
  y <- diag_data$response
  params <- diag_data$params

  if (length(y) == 0) {
    stop("No finite observations available for the rootogram.")
  }

  breaks <- pretty(range(y, na.rm = TRUE), n = bins)
  if (length(unique(breaks)) < 2) {
    breaks <- seq(min(y) - 0.5, max(y) + 0.5, length.out = bins + 1)
  }

  obs <- hist(y, breaks = breaks, plot = FALSE, include.lowest = TRUE, right = FALSE)$counts
  exp <- vapply(seq_len(length(breaks) - 1), function(i) {
    upper <- .gl_call_family_fun("p", diag_data$family, breaks[i + 1], params)
    lower <- .gl_call_family_fun("p", diag_data$family, breaks[i], params)
    sum(pmax(upper - lower, 0), na.rm = TRUE)
  }, numeric(1))

  root_df <- data.frame(
    lower = breaks[-length(breaks)],
    upper = breaks[-1],
    midpoint = (breaks[-length(breaks)] + breaks[-1]) / 2,
    observed = obs,
    expected = exp,
    root_diff = sqrt(obs) - sqrt(exp)
  )

  if (!plot) {
    return(root_df)
  }

  midpoint <- root_diff <- NULL
  ggplot2::ggplot(root_df, ggplot2::aes(x = midpoint, y = root_diff)) +
    ggplot2::geom_hline(yintercept = 0, color = "#666666") +
    ggplot2::geom_col(fill = "#2c7fb8", alpha = 0.8, width = diff(range(root_df$midpoint)) / max(length(root_df$midpoint), 1)) +
    ggplot2::labs(x = "Response", y = expression(sqrt(O) - sqrt(E)), title = "Rootogram") +
    ggplot2::theme_minimal()
}

.gl_crps_sample <- function(y, draws) {
  mean(abs(draws - y)) - 0.5 * mean(abs(outer(draws, draws, "-")))
}

#' @export
proscore.gamlss.longitudinal <- function(object, type = c("logs", "crps", "mae", "mse", "dss"), crps_grid = 25, ...) {
  type <- match.arg(type, several.ok = TRUE)
  diag_data <- .gl_diag_data(object)
  y <- diag_data$response
  params <- diag_data$params
  mu_hat <- diag_data$mu_hat
  sigma_hat <- diag_data$sigma_hat

  out <- setNames(numeric(length(type)), type)
  density_hat <- .gl_call_family_fun("d", diag_data$family, y, params)

  if ("logs" %in% type) {
    out["logs"] <- mean(-log(pmax(density_hat, .Machine$double.eps)), na.rm = TRUE)
  }
  if ("mae" %in% type) {
    out["mae"] <- mean(abs(y - mu_hat), na.rm = TRUE)
  }
  if ("mse" %in% type) {
    out["mse"] <- mean((y - mu_hat)^2, na.rm = TRUE)
  }
  if ("dss" %in% type) {
    out["dss"] <- mean(log(sigma_hat^2) + ((y - mu_hat)^2 / sigma_hat^2), na.rm = TRUE)
  }
  if ("crps" %in% type) {
    p_grid <- seq_len(crps_grid) / (crps_grid + 1)
    sample_mat <- vapply(p_grid, function(prob) {
      .gl_call_family_fun("q", diag_data$family, prob, params)
    }, numeric(length(y)))
    if (is.null(dim(sample_mat))) {
      sample_mat <- matrix(sample_mat, ncol = 1)
    }
    out["crps"] <- mean(vapply(seq_len(nrow(sample_mat)), function(i) {
      .gl_crps_sample(y[i], sample_mat[i, ])
    }, numeric(1)), na.rm = TRUE)
  }

  out
}

#' @export
procast.gamlss.longitudinal <- function(object, type = c("quantile", "cdf", "density"), at = c(0.025, 0.5, 0.975), newdata = NULL, ...) {
  type <- match.arg(type)
  require_response <- type %in% c("cdf", "density")
  diag_data <- .gl_fitted_distribution(object, newdata = newdata, require_response = require_response)
  y <- diag_data$response
  params <- diag_data$params

  if (type == "quantile") {
    quantile_df <- data.frame(response = y)
    for (prob in at) {
      quantile_df[[paste0("q", gsub("^0\\.", "", format(prob, trim = TRUE)) )]] <- .gl_call_family_fun("q", diag_data$family, prob, params)
    }
    return(quantile_df)
  }

  if (type == "cdf") {
    return(data.frame(response = y, cdf = .gl_call_family_fun("p", diag_data$family, y, params)))
  }

  data.frame(response = y, density = .gl_call_family_fun("d", diag_data$family, y, params))
}
