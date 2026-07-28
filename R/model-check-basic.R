.gl_basic_checks <- function(checks) {
  checks[, c("area", "status", "value", "threshold_condition", "message"), drop = FALSE]
}

.gl_basic_checks_result <- function(checks) {
  if (any(checks$status == "FAIL", na.rm = TRUE)) {
    return("failed")
  }

  if (any(checks$status == "REVIEW", na.rm = TRUE)) {
    return("review")
  }

  "passed"
}
