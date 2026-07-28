#' Finalize marginal-effects rows with contrasts and intervals
#'
#' @param rows List of one-row data frames from `.gl_effect_counterfactual_row()`.
#' @param reference Reference counterfactual value.
#' @param se.fit Logical; whether confidence intervals may be added.
#' @param level Confidence level.
#' @return Data frame with reference, contrast, and optional confidence limits.
#' @noRd
.gl_finalize_marginal_effects <- function(rows, reference, se.fit, level) {
  out <- do.call(rbind, rows)

  ref_idx <- match(as.character(reference), out$value)

  if (is.na(ref_idx)) {
    stop("'reference' must be one of the counterfactual values.", call. = FALSE)
  }

  out$reference <- out$value[[ref_idx]]

  out$contrast <- out$estimate - out$estimate[[ref_idx]]

  if (isTRUE(se.fit) && any(is.finite(out$std_error))) {
    alpha <- 1 - level

    z <- stats::qnorm(1 - alpha / 2)

    out$conf.low <- out$estimate - z * out$std_error

    out$conf.high <- out$estimate + z * out$std_error
  }

  rownames(out) <- NULL

  out
}
