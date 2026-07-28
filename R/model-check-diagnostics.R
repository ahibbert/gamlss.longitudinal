.gl_residual_dependence_summary <- function(object, residual_lags = 1L) {
  out <- data.frame(lag = as.integer(residual_lags), normal_score_cor = NA_real_, n_pairs = 0L)

  if (length(residual_lags) == 0L) {
    return(data.frame(lag = integer(0), normal_score_cor = numeric(0), n_pairs = integer(0)))
  }

  copula_spec <- get_copula_dist(object$copula_dist)

  family_name <- .copula_family_code(copula_spec$copula_dist)

  family_num <- tryCatch(.copula_family_code(family_name), error = function(e) NA_character_)

  lag_summary <- tryCatch(
    {
      fit_data <- .copula_v2_fit_data(object)

      rosenblatt_df <- .copula_v2_rosenblatt_series(fit_data, family_num)

      .copula_v2_rosenblatt_lag_summary(rosenblatt_df, lag_values = residual_lags)
    },
    error = function(e) data.frame()
  )

  if (nrow(lag_summary) == 0L) {
    return(out)
  }

  data.frame(
    lag = as.integer(lag_summary$lag),
    normal_score_cor = as.numeric(lag_summary$cor_z),
    n_pairs = as.integer(lag_summary$n_pairs),
    stringsAsFactors = FALSE
  )
}

.gl_pit_tail_summary <- function(pit, thresholds = c(0.05, 0.10)) {
  n <- sum(is.finite(pit))

  rows <- lapply(thresholds, function(threshold) {
    lower <- mean(pit <= threshold, na.rm = TRUE)

    upper <- mean(pit >= 1 - threshold, na.rm = TRUE)

    data.frame(
      threshold = threshold,
      lower = lower,
      upper = upper,
      expected = threshold,
      lower_ratio = lower / threshold,
      upper_ratio = upper / threshold,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)

  if (n == 0L) {
    out[, c("lower", "upper", "lower_ratio", "upper_ratio")] <- NA_real_
  }

  tail_ratios <- c(out$lower_ratio, out$upper_ratio)

  attr(out, "tail_ratio_max") <- if (any(is.finite(tail_ratios))) {
    max(tail_ratios, na.rm = TRUE)
  } else {
    NA_real_
  }

  out
}
