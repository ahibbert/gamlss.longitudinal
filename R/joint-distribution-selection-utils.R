.joint_selection_format_seconds <- function(seconds) {
  seconds <- as.numeric(seconds)[1L]
  if (!is.finite(seconds) || seconds < 0) {
    return("unknown")
  }
  if (seconds < 60) {
    return(sprintf("%.1fs", seconds))
  }
  minutes <- floor(seconds / 60)
  remaining <- seconds - minutes * 60
  if (minutes < 60) {
    return(sprintf("%dm %.0fs", minutes, remaining))
  }
  hours <- floor(minutes / 60)
  minutes <- minutes - hours * 60
  sprintf("%dh %dm", hours, minutes)
}
