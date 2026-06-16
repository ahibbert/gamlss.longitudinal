#' Validate optimizer step controls input and make sure it's in the right format, or error.
#' 
#' Insert defaults if no value provided.
#'
#' @noRd
.gl_normalize_step_controls <- function(
    method,
    include_dlcopdpar,
    start_step_size,
    max_steps,
    step_adjustment,
    verbose) {
  if (length(start_step_size) != 1 || !is.numeric(start_step_size) ||
    !is.finite(start_step_size) || start_step_size <= 0) {
    stop("ERROR: start_step_size must be a single positive finite numeric value.")
  }

  if (length(max_steps) != 1 || !is.numeric(max_steps) || is.na(max_steps)) {
    stop("ERROR: max_steps must be a single non-negative integer.")
  }

  max_steps <- as.integer(max_steps)

  if (max_steps < 0) {
    stop("ERROR: max_steps must be a single non-negative integer.")
  }

  if (length(step_adjustment) != 1 || is.null(step_adjustment)) {
    stop("ERROR: step_adjustment must be a single positive numeric value, or NA for the method-specific default.")
  }

  step_adjustment <- as.numeric(step_adjustment)

  if (is.na(step_adjustment)) {
    rs_joint_step_adjustment_default <- 1
    rs_separate_step_adjustment_default <- 1

    step_adjustment <- if (method == "RS" && isTRUE(include_dlcopdpar)) {
      rs_joint_step_adjustment_default
    } else if (method == "RS") {
      rs_separate_step_adjustment_default
    } else {
      1
    }

    if (verbose > 0) {
      cat(
        "\nUsing automatic step_adjustment=",
        signif(step_adjustment, 4),
        " for ",
        if (method == "RS" && isTRUE(include_dlcopdpar)) "joint RS" else if (method == "RS") "separate RS" else method,
        ".\n",
        sep = ""
      )
    }
  } else if (!is.finite(step_adjustment) || step_adjustment <= 0) {
    stop("ERROR: step_adjustment must be a single positive numeric value, or NA for the method-specific default.")
  }

  list(
    start_step_size = start_step_size,
    max_steps = max_steps,
    step_adjustment = step_adjustment
  )
}
