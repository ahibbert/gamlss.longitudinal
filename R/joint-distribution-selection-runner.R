.joint_selection_fit_candidates <- function(
    combinations,
    data,
    response_var,
    time_var,
    subject_var,
    fit_args,
    time_intercepts,
    copula_time_intercepts,
    trace,
    progress,
    keep_fits) {
  fit_store <- vector("list", nrow(combinations))
  rows <- vector("list", nrow(combinations))
  screen_start <- Sys.time()
  n_combinations <- nrow(combinations)
  if (isTRUE(progress)) {
    message(
      "Joint distribution screen: fitting ",
      n_combinations,
      " model(s) across ",
      length(unique(combinations$margin_family)),
      " margin candidate(s) and ",
      length(unique(combinations$copula_family)),
      " copula candidate(s)."
    )
  }
  for (ii in seq_len(nrow(combinations))) {
    margin_family <- combinations$margin_family[[ii]]
    copula_family <- combinations$copula_family[[ii]]
    candidate_label <- paste(margin_family, copula_family, sep = "+")
    if (isTRUE(progress)) {
      message("[", ii, "/", n_combinations, "] Fitting ", candidate_label, "...")
    }
    margin_dist <- .margin_family_object(margin_family)
    fit_one <- .joint_selection_fit_one(
      data = data,
      response_var = response_var,
      time_var = time_var,
      subject_var = subject_var,
      margin_family = margin_family,
      margin_dist = margin_dist,
      copula_family = copula_family,
      fit_args = fit_args,
      time_intercepts = time_intercepts,
      copula_time_intercepts = copula_time_intercepts,
      trace = trace
    )
    rows[[ii]] <- fit_one$row
    if (isTRUE(keep_fits)) {
      fit_store[[ii]] <- fit_one$fit
    }
    if (isTRUE(progress)) {
      elapsed_total <- as.numeric(difftime(Sys.time(), screen_start, units = "secs"))
      elapsed_fit <- as.numeric(fit_one$row$elapsed_sec[[1L]])
      remaining <- if (ii < n_combinations && elapsed_total > 0) {
        elapsed_total / ii * (n_combinations - ii)
      } else {
        0
      }
      status <- if (!is.na(fit_one$row$error[[1L]])) {
        paste0("failed: ", fit_one$row$error[[1L]])
      } else if (isTRUE(fit_one$row$converged[[1L]])) {
        "completed"
      } else {
        "completed without confirmed convergence"
      }
      message(
        "[", ii, "/", n_combinations, "] ",
        candidate_label,
        " ", status,
        " in ", .joint_selection_format_seconds(elapsed_fit),
        ". Elapsed ", .joint_selection_format_seconds(elapsed_total),
        "; ETA ", .joint_selection_format_seconds(remaining),
        "."
      )
    }
  }

  list(rows = rows, fit_store = fit_store)
}

.joint_selection_finalize_result <- function(
    rows,
    fit_store,
    n_pairs,
    criterion,
    margin_selection,
    time_intercepts,
    time_var,
    copula_time_intercepts,
    keep_fits) {
  out <- do.call(rbind, rows)
  out$n_pairs <- n_pairs
  out$rank <- NA_integer_
  success <- is.na(out$error) & is.finite(out[[criterion]])
  ord_success <- order(out[success, criterion], decreasing = identical(criterion, "logLik"))
  success_idx <- which(success)[ord_success]
  failed_idx <- which(!success)
  out <- out[c(success_idx, failed_idx), , drop = FALSE]
  if (any(success)) {
    out$rank[seq_along(success_idx)] <- seq_along(success_idx)
  }
  out$delta_AIC <- if (any(is.finite(out$AIC))) out$AIC - min(out$AIC, na.rm = TRUE) else NA_real_
  rownames(out) <- NULL

  attr(out, "selected") <- if (length(success_idx) > 0L) {
    paste(out$margin_family[[1L]], out$copula_family[[1L]], sep = "+")
  } else {
    NA_character_
  }
  attr(out, "criterion") <- criterion
  attr(out, "margin_selection") <- margin_selection
  attr(out, "response_type") <- attr(margin_selection, "response_type")
  attr(out, "time_intercepts") <- isTRUE(time_intercepts)
  attr(out, "time_var") <- if (isTRUE(time_intercepts)) time_var else NULL
  attr(out, "copula_time_intercepts") <- isTRUE(copula_time_intercepts)
  attr(out, "copula_time_var") <- if (isTRUE(copula_time_intercepts)) time_var else NULL
  if (isTRUE(keep_fits)) {
    fit_store <- fit_store[c(success_idx, failed_idx)]
    attr(out, "fits") <- fit_store
  }
  class(out) <- c("joint_distribution_selection", "data.frame")
  out
}
