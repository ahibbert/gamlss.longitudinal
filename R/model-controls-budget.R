#' Check whether a fit has exceeded its elapsed-time budget
#'
#' @noRd
.gl_check_elapsed_budget <- function(
    fit_start_time,
    max_elapsed_sec,
    stage = "optimisation",
    now = Sys.time()) {
  if (is.finite(max_elapsed_sec) && max_elapsed_sec > 0) {
    elapsed <- as.numeric(difftime(now, fit_start_time, units = "secs"))

    if (elapsed > max_elapsed_sec) {
      stop(
        sprintf(
          "Model exceeded max_elapsed_sec during %s (elapsed %.1f sec > %.1f sec).",
          stage,
          elapsed,
          max_elapsed_sec
        ),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

#' Build an elapsed-budget checker for one model fit
#'
#' @noRd
.gl_build_elapsed_budget_checker <- function(
    fit_start_time,
    max_elapsed_sec,
    check_fn = .gl_check_elapsed_budget) {
  function(stage = "optimisation") {
    check_fn(
      fit_start_time = fit_start_time,
      max_elapsed_sec = max_elapsed_sec,
      stage = stage
    )
  }
}
