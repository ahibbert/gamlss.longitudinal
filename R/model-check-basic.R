.gl_basic_checks <- function(checks) {
  checks[, c("area", "status", "value", "threshold_condition", "message"), drop = FALSE]
}

.gl_basic_checks_result <- function(checks) {
  if (any(checks$status == "not_converged", na.rm = TRUE)) return("not_converged")
  if (any(checks$status %in% c("flagged", "review", "unavailable"), na.rm = TRUE)) return("review")
  "descriptive"
}
