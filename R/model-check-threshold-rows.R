#' Map check status to reviewer-facing severity
#'
#' @param status Check status string.
#' @return Severity string used in `check.gamlss.longitudinal()` output.
#' @noRd
.gl_check_status_severity <- function(status) {
  switch(status,
    not_converged = "concern",
    flagged = "review",
    review = "review",
    unavailable = "review",
    "information"
  )
}

#' Build one reviewer-facing model-check threshold row
#'
#' @param area Diagnostic area.
#' @param quantity_checked Quantity or object slot being checked.
#' @param value Observed value formatted for display.
#' @param threshold_condition Threshold or rule being applied.
#' @param default Default threshold value.
#' @param status Check status.
#' @param message Reviewer-facing interpretation.
#' @param action Suggested next diagnostic or model action.
#' @return One-row data frame with the model-check table schema.
#' @noRd
.gl_check_row <- function(area, quantity_checked, value, threshold_condition, default,
                          status, message, action) {
  severity <- .gl_check_status_severity(status)

  data.frame(
    area = area,
    quantity_checked = quantity_checked,
    value = value,
    threshold_condition = threshold_condition,
    default = default,
    status = status,
    severity = severity,
    message = message,
    action = action,
    stringsAsFactors = FALSE
  )
}
