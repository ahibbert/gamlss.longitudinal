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
