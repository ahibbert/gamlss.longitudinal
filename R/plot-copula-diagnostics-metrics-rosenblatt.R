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

.copula_v2_rosenblatt_qq_data <- function(rosenblatt_df) {
  z <- rosenblatt_df$z[is.finite(rosenblatt_df$z)]
  if (length(z) == 0) {
    return(data.frame())
  }

  data.frame(
    theoretical = stats::qnorm(stats::ppoints(length(z))),
    observed = sort(z),
    stringsAsFactors = FALSE
  )
}

.copula_v2_rosenblatt_lag_summary <- function(rosenblatt_df, lag_values = 1:3) {
  lag_values <- sort(unique(as.integer(lag_values)))

  lag_values <- lag_values[lag_values > 0]

  rows <- lapply(lag_values, function(lag_value) {
    pair_list <- lapply(split(rosenblatt_df, rosenblatt_df$subject), function(x) {
      x <- x[order(x$time), , drop = FALSE]

      if (nrow(x) <= lag_value) {
        return(NULL)
      }

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
