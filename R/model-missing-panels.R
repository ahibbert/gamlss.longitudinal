#' Expand missing subject-time rows for fitting (we just add them and place NA in response and predictors)
#'
#' @noRd
.gl_expand_fit_panel <- function(
    dataset,
    time_covariate_is_factor,
    time_covariate_levels,
    time_covariate_ordered,
    verbose = 1) {
  # Mapping from internal time index back to preserved covariate values
  time_lookup <- .gl_fit_time_lookup(dataset)

  # Expand to full subject x time dataframe with NAs for missing rows, then merge back to original data to fill in observed values.
  observed_n <- nrow(dataset)
  full_grid <- .gl_fit_full_panel_grid(dataset)
  dataset <- merge(full_grid, dataset, by = c("subject", "time"), all.x = TRUE, sort = FALSE)
  dataset$time_covariate <- time_lookup$time_covariate[match(dataset$time, time_lookup$time)]

  if (time_covariate_is_factor) {
    dataset$time_covariate <- factor(as.character(dataset$time_covariate),
      levels = time_covariate_levels,
      ordered = time_covariate_ordered
    )

    if (time_covariate_ordered) {
      time_contr <- contr.treatment(length(time_covariate_levels))

      if (length(time_covariate_levels) > 1) {
        colnames(time_contr) <- time_covariate_levels[-1]
      }

      contrasts(dataset$time_covariate) <- time_contr
    }
  }

  dataset <- dataset[order(dataset$subject, dataset$time), , drop = FALSE]

  rownames(dataset) <- NULL

  inserted_n <- nrow(dataset) - observed_n

  if (verbose > 0 && inserted_n > 0) {
    cat("Inserted", inserted_n, "missing subject/time rows as NA entries.\n\n")
  }

  list(dataset = dataset, inserted_n = inserted_n)
}


#' Summarize response missingness by margin and margin pair for dependence
#'
#' @noRd
.gl_summarize_fit_missingness <- function(dataset, verbose = 1) {
  # Missingness summary by time and consecutive time pairs.

  time_levels <- sort(unique(dataset$time))
  n_time_levels <- length(time_levels)

  miss_by_time <- do.call(rbind, lapply(time_levels, function(ti) {
    idx <- dataset$time == ti

    n_total <- sum(idx)

    n_na <- sum(is.na(dataset$response[idx]))

    c(time = ti, n_total = n_total, n_na_response = n_na, n_observed_response = n_total - n_na)
  }))

  miss_by_time <- as.data.frame(miss_by_time)
  rownames(miss_by_time) <- NULL

  pair_summary <- data.frame(
    time1 = numeric(0),
    time2 = numeric(0),
    complete_pairs = integer(0),
    total_pairs = integer(0),
    total_observations = integer(0)
  )

  if (n_time_levels > 1) {
    for (i in seq_len(n_time_levels - 1)) {
      t1 <- time_levels[i]
      t2 <- time_levels[i + 1]
      d1 <- dataset[dataset$time == t1, c("subject", "response")]
      d2 <- dataset[dataset$time == t2, c("subject", "response")]
      names(d1) <- c("subject", "response_t1")
      names(d2) <- c("subject", "response_t2")
      merged_pair <- merge(d1, d2, by = "subject", all = FALSE)
      total_pairs <- nrow(merged_pair)
      complete_pairs <- sum(!is.na(merged_pair$response_t1) & !is.na(merged_pair$response_t2))

      pair_summary <- rbind(
        pair_summary,
        data.frame(
          time1 = t1,
          time2 = t2,
          complete_pairs = complete_pairs,
          total_pairs = total_pairs,
          total_observations = nrow(dataset)
        )
      )
    }
  }

  if (verbose > 0) {
    cat("Missingness Summary (by time):\n")
    print(miss_by_time)
    cat("\nConsecutive Pair Completeness:\n")

    if (nrow(pair_summary) > 0) {
      print(pair_summary)
    } else {
      cat("No consecutive time pairs available.\n")
    }

    cat("\n")
  }

  list(miss_by_time = miss_by_time, pair_summary = pair_summary)
}

#' Validate that each margin and each copula pair is not 100% missing any complete pairs
#'
#' @noRd
.gl_validate_fit_missingness_support <- function(miss_by_time, pair_summary) {
  # Hard stop if any margin is 100% missing.
  margin_all_missing <- miss_by_time$n_observed_response == 0
  if (any(margin_all_missing)) {
    bad_times <- miss_by_time$time[margin_all_missing]
    stop(
      "ERROR: 100% missing response values detected for margin time point(s): ",
      paste(bad_times, collapse = ", "),
      "\nModel fitting stopped because at least one margin has no observed outcomes."
    )
  }

  # Hard stop if any consecutive copula pair has 0 complete pairs while pairs exist.
  if (nrow(pair_summary) > 0) {
    pair_all_missing <- pair_summary$total_pairs > 0 & pair_summary$complete_pairs == 0
    if (any(pair_all_missing)) {
      bad_pairs <- apply(pair_summary[pair_all_missing, c("time1", "time2"), drop = FALSE], 1, function(x) {
        paste0("(", x[1], ",", x[2], ")")
      })
      stop(
        "ERROR: 100% missing complete copula pairs detected for consecutive time pair(s): ",
        paste(bad_pairs, collapse = ", "),
        "\nModel fitting stopped because at least one copula pair contributes no complete observations."
      )
    }
  }

  invisible(TRUE)
}
