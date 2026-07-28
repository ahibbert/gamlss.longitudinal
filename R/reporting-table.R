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
    direction = c("above", "below")) {
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
