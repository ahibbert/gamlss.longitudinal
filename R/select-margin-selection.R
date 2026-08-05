#' Select candidate marginal distributions

#'

#' `select_margin()` is a lightweight wrapper around [gamlss::fitDist()] for

#' the recommended longitudinal workflow: choose a plausible marginal family,

#' then fit dependence with [gamlss_longitudinal()].

#'

#' @param response Numeric response vector. For the common

#'   `select_margin(dat, response_var = "y")` call, a data frame supplied here

#'   is treated as `data`.

#' @param data Optional data frame containing the response.

#' @param response_var Optional response column name in `data`.

#' @param type Optional `gamlss::fitDist()` type. If `NULL`, a simple heuristic

#'   uses `"counts"` for non-negative integer responses, `"realplus"` for

#'   positive continuous responses, and `"realAll"` otherwise.

#' @param families Optional character vector used to filter the returned table.

#' @param time_intercepts Logical; if `TRUE`, rank retained families by a

#'   GAMLSS fit with factor time intercepts for each distribution parameter.

#'   This is useful for longitudinal screening where the final model is allowed

#'   to vary marginal location, scale, or shape over time. Finite-AIC fits are

#'   retained even if the temporary GAMLSS fit did not report convergence; check

#'   the returned `converged` column before selecting a final marginal family.

#' @param time_var Optional time column name in `data`, required when

#'   `time_intercepts = TRUE`.

#' @param try.gamlss,trace,... Passed to [gamlss::fitDist()].

#'

#' @return A data frame ordered by AIC and class `margin_selection`. With

#'   `time_intercepts = TRUE`, the table includes a `converged` column for the

#'   temporary time-intercept marginal fit. The selected family is also stored in

#'   the `"selected"` attribute for backward compatibility. Use [best_fit()] or

#'   [best_fit_family()] to extract the selected family in a fitting-friendly

#'   form.

#' @export

select_margin <- function(
    response = NULL,
    data = NULL,
    response_var = NULL,
    type = NULL,
    families = NULL,
    time_intercepts = FALSE,
    time_var = NULL,
    try.gamlss = FALSE,
    trace = FALSE,
    ...) {
  inputs <- .select_margin_inputs(
    response = response,
    data = data,
    response_var = response_var,
    type = type,
    time_intercepts = time_intercepts,
    time_var = time_var
  )
  response <- inputs$response
  response_all <- inputs$response_all
  data <- inputs$data
  type <- inputs$type
  time_intercepts <- inputs$time_intercepts
  time_var <- inputs$time_var

  fit <- NULL

  if (isTRUE(trace)) {
    fit <- gamlss::fitDist(

      response,
      type = type,
      try.gamlss = try.gamlss,
      trace = trace,
      ...
    )
  } else {
    invisible(utils::capture.output({
      invisible(utils::capture.output({
        fit <- suppressWarnings(suppressMessages(gamlss::fitDist(

          response,
          type = type,
          try.gamlss = try.gamlss,
          trace = trace,
          ...
        )))
      }, type = "message"))
    }))
  }

  if (is.null(fit$fits)) {
    stop("gamlss::fitDist() did not return a fits table.", call. = FALSE)
  }

  out <- data.frame(
    family = names(fit$fits),
    AIC = as.numeric(fit$fits),
    type = type,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  out <- out[is.finite(out$AIC), , drop = FALSE]

  if (!is.null(families)) {
    out <- out[out$family %in% families, , drop = FALSE]
  }

  out$screen_model <- "pooled_intercept"

  if (isTRUE(time_intercepts) && nrow(out) > 0L) {
    time_fits <- .select_margin_time_intercept_fits(
      response = response_all,
      time = data[[time_var]],
      families = out$family,
      trace = trace
    )

    out$pooled_AIC <- out$AIC

    out <- merge(

      out,
      time_fits,
      by = "family",
      all.x = TRUE,
      sort = FALSE
    )

    out$AIC <- out$time_intercept_AIC

    out$time_intercept_AIC <- NULL

    out$screen_model <- "time_intercepts"

    out <- out[is.finite(out$AIC), , drop = FALSE]

  }

  out <- out[order(out$AIC), , drop = FALSE]

  rownames(out) <- NULL

  out$rank <- seq_len(nrow(out))

  out$delta_AIC <- if (nrow(out) > 0L) out$AIC - min(out$AIC, na.rm = TRUE) else numeric(0)

  out$supported_by_longitudinal <- vapply(out$family, function(family) {
    all(vapply(

      paste0(c("d", "p", "q"), family),
      exists,
      logical(1),
      envir = asNamespace("gamlss.dist"),
      inherits = FALSE
    ))
  }, logical(1), USE.NAMES = FALSE)

  attr(out, "selected") <- if (nrow(out) > 0L) out$family[[1L]] else NA_character_

  attr(out, "response_type") <- type

  attr(out, "time_intercepts") <- isTRUE(time_intercepts)

  attr(out, "time_var") <- if (isTRUE(time_intercepts)) time_var else NULL

  class(out) <- c("margin_selection", "margin_screen", class(out))

  out
}
