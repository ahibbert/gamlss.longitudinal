#' Pass defaults to the warm start fit call
#'
#' @noRd
.gl_default_warm_start_info <- function() {
  list(
    used = FALSE,
    outer_iter = 0L,
    include_dlcopdpar = FALSE,
    log_lik = NULL
  )
}

#' Decide whether the separate RS warm-start phase should run based on inputs
#' and if so, run it and capture the output, warnings, and errors.
#' 
#' Don't run the separate RS warm-start phase if any of the following are true:
#' - The method is not RS
#' - The dlcopdpar is not included in the model
#' - The user has set warm_start_joint = FALSE
#' - The user has set warm_start_joint_iter = 0
#' - The user has supplied their own starting values
#'
#' @noRd
.gl_should_run_joint_warm_start <- function(
    method,
    include_dlcopdpar,
    warm_start_joint,
    warm_start_joint_iter,
    user_supplied_start) {
  method == "RS" &&
    isTRUE(include_dlcopdpar) &&
    isTRUE(warm_start_joint) &&
    warm_start_joint_iter > 0L &&
    !isTRUE(user_supplied_start)
}
