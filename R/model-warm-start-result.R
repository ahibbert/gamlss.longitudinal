#' Validate that the warm start has not failed and has actually returned starting values
#' If it does fail, tell user to try a cold start instead or CG.
#'
#' @noRd
.gl_validate_joint_warm_start_fit <- function(warm_fit, warm_err) {
  if (!is.null(warm_err)) {
    stop(
      "Separate RS warm-start phase failed: ",
      conditionMessage(warm_err),
      "\nSet warm_start_joint = FALSE to force a cold-start joint fit.",
      call. = FALSE
    )
  }

  if (is.null(warm_fit) || is.null(warm_fit$par)) {
    stop(
      "Separate RS warm-start phase did not return coefficient starting values.\n",
      "Set warm_start_joint = FALSE to force a cold-start joint fit or try method=\"CG\".",
      call. = FALSE
    )
  }

  invisible(warm_fit)
}

#' Output the result list from the warm start for use as the start argument in the main fit
#'
#' @noRd
.gl_joint_warm_start_result <- function(warm_fit, warm_output, warm_warnings, warm_start_joint_iter) {
  warm_start_par_s <- warm_fit$par_s
  list(
    start_from = warm_fit$par,
    warm_start_par_s = warm_start_par_s,
    warm_start_info = list(
      used = TRUE,
      outer_iter = warm_start_joint_iter,
      include_dlcopdpar = FALSE,
      log_lik = warm_fit$calc_lik_out_end$log_lik,
      carries_smooth = !is.null(warm_start_par_s) && any(vapply(warm_start_par_s, length, integer(1L)) > 0L),
      captured_output = warm_output,
      captured_warnings = warm_warnings
    )
  )
}
