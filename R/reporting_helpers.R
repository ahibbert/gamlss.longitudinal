#' Audit a fitted longitudinal GAMLSS-copula model specification
#'
#' `model_spec()` returns a compact, printable audit trail for a fitted model:
#' formulas, distributions, links, optimisation settings, missing-response
#' counts, likelihood components, and variance-covariance metadata.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#'
#' @return A list of class `gamlss_longitudinal_model_spec`.
#' @export
model_spec <- function(object) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }

  family <- object$margin_dist
  link_names <- names(family)[grepl("\\.link$", names(family))]
  margin_links <- data.frame(
    parameter = sub("\\.link$", "", link_names),
    link = as.character(unlist(family[link_names], use.names = FALSE)),
    stringsAsFactors = FALSE
  )
  if (nrow(margin_links) == 0L) {
    margin_links <- data.frame(parameter = character(0), link = character(0))
  }

  response <- object$response
  model_selection <- tryCatch(summary(object, include_vcov = FALSE)$fit$model_selection, error = function(e) NULL)
  dep <- tryCatch(copula_time_summary(object), error = function(e) NULL)

  out <- list(
    variables = list(
      response = object$response_var %||% "response",
      subject = object$subject_var %||% "subject",
      time = object$time_var %||% "time"
    ),
    distributions = list(
      margin = as.character(object$margin_dist$family[1]),
      margin_label = as.character(object$margin_dist$family[2] %||% object$margin_dist$family[1]),
      copula = as.character(object$copula_dist)
    ),
    formulas = object$formulas,
    margin_links = margin_links,
    optimisation = list(
      method = object$optim_method %||% object$convergence$method %||% NA_character_,
      converged = isTRUE(object$convergence$converged),
      stop_reason = object$convergence$stop_reason %||% NA_character_,
      outer_iterations = object$convergence$outer_iterations %||% NA_integer_,
      max_outer_iter = object$convergence$max_outer_iter %||% NA_integer_
    ),
    missingness = list(
      n_rows = length(response),
      n_missing_response = sum(is.na(response)),
      n_nonfinite_response = sum(!is.na(response) & !is.finite(response))
    ),
    likelihood = model_selection,
    vcov = c(object$vcov_meta %||% list(precomputed = FALSE), list(
      hessian_diagnostics = object$vcov$hessian_diagnostics %||% NULL
    )),
    dependence = if (!is.null(dep)) dep$time_summary else NULL
  )
  class(out) <- "gamlss_longitudinal_model_spec"
  out
}

#' @export
print.gamlss_longitudinal_model_spec <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nGAMLSS Longitudinal Model Specification\n")
  cat("--------------------------------------\n")
  cat("Response:", x$variables$response, " | Subject:", x$variables$subject, " | Time:", x$variables$time, "\n")
  cat("Margin:", x$distributions$margin, " | Copula:", x$distributions$copula, "\n")
  cat("Optimisation:", x$optimisation$method, " | Converged:", if (isTRUE(x$optimisation$converged)) "yes" else "no", "\n")
  cat("Missing responses:", x$missingness$n_missing_response, "of", x$missingness$n_rows, "\n")
  if (!is.null(x$vcov)) {
    cat("VCOV:", x$vcov$method_used %||% x$vcov$method %||% "not precomputed", "\n")
  }

  if (!is.null(x$formulas) && length(x$formulas) > 0L) {
    cat("\nFormulas\n")
    print(vapply(x$formulas, function(f) paste(deparse(f), collapse = " "), character(1)))
  }
  if (!is.null(x$margin_links) && nrow(x$margin_links) > 0L) {
    cat("\nMargin Links\n")
    print(x$margin_links, row.names = FALSE)
  }
  if (!is.null(x$likelihood)) {
    cat("\nLikelihood And Information Criteria\n")
    print(round(x$likelihood, digits))
  }
  if (!is.null(x$dependence) && nrow(x$dependence) > 0L) {
    cat("\nDependence By Time (includes Kendall's tau)\n")
    print(x$dependence, digits = digits, row.names = FALSE)
  }
  invisible(x)
}

#' Build an applied reporting table from a fitted model
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param newdata Data to summarize.
#' @param by Optional grouping columns in `newdata`, such as treatment and time.
#' @param probs Quantiles to include.
#' @param threshold Optional threshold for probabilities.
#' @param direction Probability direction when `threshold` is supplied.
#'
#' @return A data frame with grouped fitted means, medians, quantiles, and
#'   optional threshold probabilities.
#' @export
reporting_table <- function(
  object,
  newdata,
  by = NULL,
  probs = c(0.1, 0.5, 0.9),
  threshold = NULL,
  direction = c("above", "below")
) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }
  if (missing(newdata) || is.null(newdata)) {
    stop("'newdata' is required.", call. = FALSE)
  }
  newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)
  if (is.null(by)) {
    by <- character(0)
  }
  by <- as.character(by)
  missing_by <- setdiff(by, names(newdata))
  if (length(missing_by) > 0L) {
    stop("'by' column(s) not found in 'newdata': ", paste(missing_by, collapse = ", "), call. = FALSE)
  }

  direction <- match.arg(direction)
  pred <- data.frame(.row = seq_len(nrow(newdata)), stringsAsFactors = FALSE)
  pred$mean <- predict(object, newdata = newdata, type = "mean")
  pred$mu <- predict(object, newdata = newdata, type = "mu")
  pred$median <- predict(object, newdata = newdata, type = "median")
  q_pred <- predict(object, newdata = newdata, type = "quantile", probs = probs)
  q_cols <- setdiff(names(q_pred), c("subject", "time", "response"))
  pred <- cbind(pred, q_pred[q_cols])
  if (!is.null(threshold)) {
    p_pred <- predict(object, newdata = newdata, type = "probability", q = threshold, direction = direction)
    pred[[paste0("prob_", direction, "_", threshold)]] <- p_pred$probability
  }

  if (length(by) == 0L) {
    out <- as.data.frame(as.list(colMeans(pred[setdiff(names(pred), ".row")], na.rm = TRUE)), stringsAsFactors = FALSE)
    out$n <- nrow(newdata)
    return(out[c("n", setdiff(names(out), "n"))])
  }

  group_data <- newdata[by]
  agg <- stats::aggregate(pred[setdiff(names(pred), ".row")], group_data, mean, na.rm = TRUE)
  counts <- stats::aggregate(pred$.row, group_data, length)
  names(counts)[ncol(counts)] <- "n"
  merge(counts, agg, by = by, sort = FALSE)
}
