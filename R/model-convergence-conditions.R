.gl_convergence_condition <- function(message, class, object = NULL, operation = NULL) {
  structure(
    list(
      message = message,
      call = NULL,
      stop_reason = object$convergence$stop_reason %||% NA_character_,
      operation = operation
    ),
    class = c(class, if (grepl("warning", class, fixed = TRUE)) "warning" else "error", "condition")
  )
}

.gl_warn_nonconverged_fit <- function(convergence) {
  reason <- convergence$stop_reason %||% "unknown"
  warning(.gl_convergence_condition(
    paste0(
      "Model returned without satisfying the optimizer convergence contract ",
      "(stop reason: ", reason, "). Point predictions remain available with a warning; ",
      "inference, uncertainty intervals, and model comparison are disabled."
    ),
    "gamlss.longitudinal_nonconvergence_warning",
    object = list(convergence = convergence),
    operation = "fit"
  ))
  invisible()
}

.gl_require_converged_fit <- function(object, operation = "inference") {
  if (!identical(object$convergence$converged, FALSE)) return(invisible(TRUE))
  reason <- object$convergence$stop_reason %||% "unknown"
  stop(.gl_convergence_condition(
    paste0(
      "Cannot perform ", operation, " because the model did not converge ",
      "(stop reason: ", reason, "). Refit successfully before using this result."
    ),
    "gamlss.longitudinal_nonconvergence_error",
    object = object,
    operation = operation
  ))
}

.gl_warn_nonconverged_prediction <- function(object) {
  if (!identical(object$convergence$converged, FALSE)) return(invisible(FALSE))
  reason <- object$convergence$stop_reason %||% "unknown"
  warning(.gl_convergence_condition(
    paste0(
      "Point prediction from a nonconverged model (stop reason: ", reason,
      "). Treat values as provisional; uncertainty is unavailable."
    ),
    "gamlss.longitudinal_nonconverged_prediction_warning",
    object = object,
    operation = "point prediction"
  ))
  invisible(TRUE)
}
