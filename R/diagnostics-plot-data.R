#' Build QQ-plot data for randomized quantile residuals
#'
#' @noRd
.gl_qq_plot_frame <- function(z, split_info) {
  if (!split_info$split_by) {
    theo <- stats::qnorm(stats::ppoints(length(z)))
    return(data.frame(theoretical = theo, observed = sort(z)))
  }

  split_z <- split(z, split_info$group)
  qq_list <- lapply(names(split_z), function(grp) {
    z_t <- split_z[[grp]]
    theo_t <- stats::qnorm(stats::ppoints(length(z_t)))
    data.frame(theoretical = theo_t, observed = sort(z_t), split_group = as.factor(grp))
  })
  do.call(rbind, qq_list)
}

#' Build normal QQ worm-plot confidence bands
#'
#' @noRd
.gl_worm_band_frame <- function(theoretical, n, band_level = 0.95) {
  p <- stats::pnorm(theoretical)
  se <- sqrt(pmax(0, p * (1 - p) / n)) / stats::dnorm(theoretical)
  z <- stats::qnorm((1 + band_level) / 2)
  data.frame(
    theoretical = theoretical,
    band_lower = -z * se,
    band_upper = z * se
  )
}

#' Build worm-plot data for PIT values
#'
#' @noRd
.gl_worm_plot_frames <- function(pit, split_info) {
  if (!split_info$split_by) {
    z <- stats::qnorm(pmin(pmax(pit, .Machine$double.eps), 1 - .Machine$double.eps))
    theo <- stats::qnorm(stats::ppoints(length(z)))
    worm_df <- data.frame(theoretical = theo, detrended = sort(z) - theo)
    worm_band <- .gl_worm_band_frame(theo, length(z))
    return(list(worm_df = worm_df, worm_band = worm_band))
  }

  split_pit <- split(pit, split_info$group)
  worm_list <- lapply(names(split_pit), function(grp) {
    pit_t <- split_pit[[grp]]
    z_t <- stats::qnorm(pmin(pmax(pit_t, .Machine$double.eps), 1 - .Machine$double.eps))
    theo_t <- stats::qnorm(stats::ppoints(length(z_t)))
    band_t <- .gl_worm_band_frame(theo_t, length(z_t))
    data.frame(
      theoretical = theo_t,
      detrended = sort(z_t) - theo_t,
      band_lower = band_t$band_lower,
      band_upper = band_t$band_upper,
      split_group = as.factor(grp)
    )
  })
  list(worm_df = do.call(rbind, worm_list), worm_band = NULL)
}

#' Build one rootogram data frame for observed and expected bin counts
#'
#' @noRd
.gl_rootogram_frame <- function(y_i, params_i, breaks, family, split_group = NULL) {
  obs <- hist(y_i, breaks = breaks, plot = FALSE, include.lowest = TRUE, right = FALSE)$counts
  exp <- vapply(seq_len(length(breaks) - 1L), function(i) {
    upper <- .gl_call_family_fun("p", family, breaks[i + 1L], params_i)
    lower <- .gl_call_family_fun("p", family, breaks[i], params_i)
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

#' Build rootogram data, optionally split by diagnostic group
#'
#' @noRd
.gl_rootogram_plot_frame <- function(y, params, breaks, family, split_info) {
  if (!split_info$split_by) {
    return(.gl_rootogram_frame(y, params, breaks, family))
  }

  group <- split_info$group
  root_groups <- unique(as.character(group[!is.na(group)]))
  root_list <- lapply(root_groups, function(grp) {
    idx <- group == grp
    params_i <- lapply(params, function(x) x[idx])
    .gl_rootogram_frame(y[idx], params_i, breaks, family, split_group = grp)
  })
  do.call(rbind, root_list)
}
