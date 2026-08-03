.select_copula_tail_thresholds <- function(tail_thresholds) {
  tail_thresholds <- sort(unique(as.numeric(tail_thresholds)))
  tail_thresholds <- tail_thresholds[is.finite(tail_thresholds)]
  if (length(tail_thresholds) < 1L || any(tail_thresholds <= 0 | tail_thresholds >= 0.5)) {
    stop("tail_thresholds must contain finite probabilities between 0 and 0.5.", call. = FALSE)
  }
  tail_thresholds
}

.select_copula_tail_diagnostics <- function(pairs, fits, thresholds, copula_time_intercepts = FALSE) {
  cooccurrence <- do.call(rbind, lapply(fits, function(fit) {
    .select_copula_tail_fit_rows(
      pairs = pairs,
      fit = fit,
      thresholds = thresholds,
      copula_time_intercepts = copula_time_intercepts
    )
  }))

  conditional <- cooccurrence
  conditional$empirical <- pmin(pmax(conditional$empirical / conditional$threshold, 0), 1)
  conditional$fitted <- pmin(pmax(conditional$fitted / conditional$threshold, 0), 1)
  conditional$diff <- conditional$fitted - conditional$empirical
  conditional$abs_diff <- abs(conditional$diff)
  conditional$ratio <- .select_copula_tail_ratio(conditional$fitted, conditional$empirical)

  list(cooccurrence = cooccurrence, conditional = conditional)
}

.select_copula_tail_fit_rows <- function(pairs, fit, thresholds, copula_time_intercepts = FALSE) {
  family <- fit$family[[1L]]
  par_values <- .select_copula_tail_pair_parameter(fit, pairs, "par", copula_time_intercepts)
  par2_values <- .select_copula_tail_pair_parameter(fit, pairs, "par2", copula_time_intercepts)

  rows <- lapply(thresholds, function(threshold) {
    upper_cut <- 1 - threshold
    lower_fitted <- .select_copula_mean_cdf(
      u1 = rep(threshold, nrow(pairs)),
      u2 = rep(threshold, nrow(pairs)),
      family = family,
      par = par_values,
      par2 = par2_values
    )
    upper_cdf <- .select_copula_mean_cdf(
      u1 = rep(upper_cut, nrow(pairs)),
      u2 = rep(upper_cut, nrow(pairs)),
      family = family,
      par = par_values,
      par2 = par2_values
    )
    upper_fitted <- 1 - 2 * upper_cut + upper_cdf

    empirical <- c(
      mean(pairs$u1 <= threshold & pairs$u2 <= threshold, na.rm = TRUE),
      mean(pairs$u1 >= upper_cut & pairs$u2 >= upper_cut, na.rm = TRUE)
    )
    fitted <- pmin(pmax(c(lower_fitted, upper_fitted), 0), 1)

    data.frame(
      family = family,
      threshold = threshold,
      tail = c("Lower", "Upper"),
      empirical = empirical,
      fitted = fitted,
      diff = fitted - empirical,
      abs_diff = abs(fitted - empirical),
      ratio = .select_copula_tail_ratio(fitted, empirical),
      n_pairs = nrow(pairs),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

.select_copula_tail_pair_parameter <- function(fit, pairs, parameter, copula_time_intercepts = FALSE) {
  if (!isTRUE(copula_time_intercepts)) {
    return(rep(fit[[parameter]][[1L]], nrow(pairs)))
  }

  time_fits <- attr(fit, "copula_time_fits")
  if (is.null(time_fits) || !parameter %in% names(time_fits)) {
    return(rep(NA_real_, nrow(pairs)))
  }

  time_fits[[parameter]][match(pairs$copula_time, time_fits$copula_time)]
}

.select_copula_mean_cdf <- function(u1, u2, family, par, par2) {
  par2 <- ifelse(is.finite(par2), par2, 0)
  cdf <- tryCatch(
    .copula_cdf(u1, u2, family = family, par = par, par2 = par2),
    error = function(e) {
      vapply(seq_along(u1), function(i) {
        if (!is.finite(par[i])) {
          return(NA_real_)
        }
        tryCatch(
          .copula_cdf(u1[i], u2[i], family = family, par = par[i], par2 = par2[i]),
          error = function(e) NA_real_
        )
      }, numeric(1), USE.NAMES = FALSE)
    }
  )
  if (!any(is.finite(cdf))) {
    return(NA_real_)
  }
  mean(cdf[is.finite(cdf)])
}

.select_copula_tail_ratio <- function(fitted, empirical) {
  ifelse(is.finite(empirical) & empirical > 0, fitted / empirical, NA_real_)
}

.select_copula_add_tail_summary <- function(out, tail_diagnostics) {
  cooccurrence <- tail_diagnostics$cooccurrence
  conditional <- tail_diagnostics$conditional
  threshold <- min(cooccurrence$threshold, na.rm = TRUE)

  out$tail_threshold <- threshold
  out$lower_tail_cooccurrence <- .select_copula_tail_value(out$family, cooccurrence, threshold, "Lower")
  out$upper_tail_cooccurrence <- .select_copula_tail_value(out$family, cooccurrence, threshold, "Upper")
  out$lower_tail_exceedance <- .select_copula_tail_value(out$family, conditional, threshold, "Lower")
  out$upper_tail_exceedance <- .select_copula_tail_value(out$family, conditional, threshold, "Upper")
  out
}

.select_copula_tail_value <- function(families, tail_df, threshold, tail) {
  vapply(families, function(family) {
    idx <- tail_df$family == family & tail_df$threshold == threshold & tail_df$tail == tail
    if (!any(idx)) {
      return(NA_real_)
    }
    tail_df$fitted[which(idx)[1L]]
  }, numeric(1), USE.NAMES = FALSE)
}
