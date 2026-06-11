#' @export
pithist.gamlss.longitudinal <- function(object, bins = 20, randomize = FALSE, plot = TRUE, by_time = FALSE, by = NULL, data = NULL, ...) {

  pit_out <- .gl_pit(object, randomize = randomize)

  pit <- pit_out$pit

  split_info <- .gl_diag_split_info(by_time, by, pit_out$diag, data = data, plot_name = "PIT histogram")

  pit_df <- data.frame(pit = pit, time = as.factor(pit_out$diag$time))

  if (split_info$split_by) {

    pit_df$split_group <- split_info$group

  }

  pit_df <- pit_df[is.finite(pit_df$pit), , drop = FALSE]


  if (!plot) {

    return(pit_df)

  }


  if (!split_info$split_by) {

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


  expected_df <- aggregate(pit ~ split_group, data = pit_df, FUN = function(x) length(x) / bins)

  names(expected_df)[2] <- "expected"


  ggplot2::ggplot(pit_df, ggplot2::aes(x = pit)) +

    ggplot2::geom_histogram(breaks = seq(0, 1, length.out = bins + 1L), closed = "right", fill = "#2c7fb8", color = "white") +

    ggplot2::geom_hline(data = expected_df, ggplot2::aes(yintercept = expected), linetype = "dashed", linewidth = 0.4, color = "#444444") +

    ggplot2::facet_wrap(~split_group, scales = "free_y") +

    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +

    ggplot2::labs(x = "PIT", y = "Count", title = paste0("PIT Histogram by ", split_info$by)) +

    ggplot2::theme_minimal()

}


#' @export

qqrplot.gamlss.longitudinal <- function(object, randomize = FALSE, plot = TRUE, by_time = FALSE, by = NULL, data = NULL, ...) {

  pit_out <- .gl_pit(object, randomize = randomize)

  z <- stats::qnorm(pmin(pmax(pit_out$pit, .Machine$double.eps), 1 - .Machine$double.eps))

  split_info <- .gl_diag_split_info(by_time, by, pit_out$diag, data = data, plot_name = "QQ plot")


  if (!split_info$split_by) {

    theo <- stats::qnorm(stats::ppoints(length(z)))

    qq_df <- data.frame(theoretical = theo, observed = sort(z))

  } else {

    split_z <- split(z, split_info$group)

    qq_list <- lapply(names(split_z), function(grp) {

      z_t <- split_z[[grp]]

      theo_t <- stats::qnorm(stats::ppoints(length(z_t)))

      data.frame(theoretical = theo_t, observed = sort(z_t), split_group = as.factor(grp))

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

      title = if (split_info$split_by) paste0("QQ Plot by ", split_info$by) else "QQ Plot"

    ) +

    ggplot2::theme_minimal()


  if (split_info$split_by) {

    p <- p + ggplot2::facet_wrap(~split_group, scales = "free")

  }


  p

}


#' @export

wormplot.gamlss.longitudinal <- function(object, randomize = FALSE, plot = TRUE, smooth = TRUE, by_time = FALSE, by = NULL, data = NULL, ...) {

  pit_out <- .gl_pit(object, randomize = randomize)

  split_info <- .gl_diag_split_info(by_time, by, pit_out$diag, data = data, plot_name = "wormplot")


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


  if (!split_info$split_by) {

    z <- stats::qnorm(pmin(pmax(pit_out$pit, .Machine$double.eps), 1 - .Machine$double.eps))

    theo <- stats::qnorm(stats::ppoints(length(z)))

    worm_df <- data.frame(theoretical = theo, detrended = sort(z) - theo)

    worm_band <- worm_band_frame(theo, length(z))

  } else {

    split_pit <- split(pit_out$pit, split_info$group)

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

      data = if (split_info$split_by) worm_df else worm_band,

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

      title = if (split_info$split_by) paste0("Worm Plot by ", split_info$by) else "Worm Plot"

    ) +

    ggplot2::theme_minimal()


  if (split_info$split_by) {

    p <- p + ggplot2::facet_wrap(~split_group, scales = "free")

  }


  if (smooth && nrow(worm_df) > 5) {

    p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#d95f0e", linewidth = 0.7)

  }


  p

}


#' @export

rootogram.gamlss.longitudinal <- function(object, bins = 20, plot = TRUE, by_time = FALSE, by = NULL, data = NULL, ...) {

  diag_data <- .gl_fitted_distribution(object, newdata = NULL, require_response = TRUE)

  y <- diag_data$response

  params <- diag_data$params

  split_info <- .gl_diag_split_info(by_time, by, diag_data, data = data, plot_name = "rootogram")


  if (length(y) == 0) {

    stop("No finite observations available for the rootogram.")

  }


  breaks <- pretty(range(y, na.rm = TRUE), n = bins)

  if (length(unique(breaks)) < 2) {

    breaks <- seq(min(y) - 0.5, max(y) + 0.5, length.out = bins + 1)

  }


  root_frame <- function(y_i, params_i, split_group = NULL) {

    obs <- hist(y_i, breaks = breaks, plot = FALSE, include.lowest = TRUE, right = FALSE)$counts

    exp <- vapply(seq_len(length(breaks) - 1L), function(i) {

      upper <- .gl_call_family_fun("p", diag_data$family, breaks[i + 1L], params_i)

      lower <- .gl_call_family_fun("p", diag_data$family, breaks[i], params_i)

      sum(pmax(upper - lower, 0), na.rm = TRUE)

    }, numeric(1))


    out <- data.frame(

      lower = breaks[-length(breaks)],

      upper = breaks[-1],

      midpoint = (breaks[-length(breaks)] + breaks[-1]) / 2,

      observed = obs,

      expected = exp,

      root_diff = sqrt(obs) - sqrt(exp)

    )

    if (!is.null(split_group)) {

      out$split_group <- as.factor(split_group)

    }

    out

  }


  if (!split_info$split_by) {

    root_df <- root_frame(y, params)

  } else {

    group <- split_info$group

    root_groups <- unique(as.character(group[!is.na(group)]))

    root_list <- lapply(root_groups, function(grp) {

      idx <- group == grp

      params_i <- lapply(params, function(x) x[idx])

      root_frame(y[idx], params_i, split_group = grp)

    })

    root_df <- do.call(rbind, root_list)

  }


  if (!plot) {

    return(root_df)

  }


  midpoint <- root_diff <- NULL

  p <- ggplot2::ggplot(root_df, ggplot2::aes(x = midpoint, y = root_diff)) +

    ggplot2::geom_hline(yintercept = 0, color = "#666666") +

    ggplot2::geom_col(fill = "#2c7fb8", alpha = 0.8, width = diff(range(root_df$midpoint)) / max(length(root_df$midpoint), 1)) +

    ggplot2::labs(

      x = "Response",

      y = expression(sqrt(O) - sqrt(E)),

      title = if (split_info$split_by) paste0("Rootogram by ", split_info$by) else "Rootogram"

    ) +

    ggplot2::theme_minimal()


  if (split_info$split_by) {

    p <- p + ggplot2::facet_wrap(~split_group, scales = "free_y")

  }


  p

}


