#' Internal table builders for copula time summaries
#'
#' These helpers keep the public summary method focused on orchestration while
#' making the fitted-time and adjacent-pair aggregation logic easy to review.
#'
#' @noRd
.gl_copula_summary_agg_fun <- function(stat) {
  if (identical(stat, "median")) stats::median else mean
}

#' @noRd
.gl_copula_time_table <- function(fit_data, agg_fun, has_zeta) {
  time_summary <- do.call(rbind, lapply(split(fit_data, fit_data$time), function(x) {
    out <- data.frame(
      time = x$time[1],
      n_obs = nrow(x),
      theta_fit = agg_fun(x$theta_fit, na.rm = TRUE),
      tau_fit = agg_fun(x$tau_fit, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    if (has_zeta) {
      out$zeta_fit <- agg_fun(x$zeta_fit, na.rm = TRUE)
    }
    out
  }))

  time_summary[order(time_summary$time), , drop = FALSE]
}

#' @noRd
.gl_copula_pair_table <- function(pair_data, agg_fun, has_zeta) {
  do.call(rbind, lapply(split(pair_data, pair_data$time_pair), function(x) {
    out <- data.frame(
      time_pair = x$time_pair[1],
      n_pairs = nrow(x),
      theta_pair = agg_fun(x$theta_pair, na.rm = TRUE),
      tau_pair = agg_fun(x$tau_fit, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    if (has_zeta) {
      out$zeta_pair <- agg_fun(x$zeta_pair, na.rm = TRUE)
    }
    out
  }))
}

#' @noRd
.gl_copula_drop_zeta_if_absent <- function(fit_data, pair_data, has_zeta) {
  if (!has_zeta) {
    # Keep the returned data tidy for one-parameter copulas.
    if ("zeta_fit" %in% names(fit_data)) {
      fit_data$zeta_fit <- NULL
    }
    if ("zeta_pair" %in% names(pair_data)) {
      pair_data$zeta_pair <- NULL
    }
  }

  list(fit_data = fit_data, pair_data = pair_data)
}
