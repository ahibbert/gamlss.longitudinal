.sim_check_count <- function(x, label) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop(label, " must be a single positive integer.", call. = FALSE)
  }
  x <- as.integer(x)
  if (!is.finite(x) || x < 1L) {
    stop(label, " must be a single positive integer.", call. = FALSE)
  }
  x
}

.sim_base_long_data <- function(n, times, subject_var, time_var) {
  data.frame(
    .sim_subject_index = rep(seq_len(n), each = length(times)),
    .sim_time_index = rep(seq_along(times), times = n),
    stringsAsFactors = FALSE
  ) |>
    within({
      subject <- factor(.sim_subject_index)
      time <- times[.sim_time_index]
    }) |>
    .sim_rename_base_columns(subject_var, time_var)
}

.sim_rename_base_columns <- function(x, subject_var, time_var) {
  names(x)[names(x) == "subject"] <- subject_var
  names(x)[names(x) == "time"] <- time_var
  x
}
