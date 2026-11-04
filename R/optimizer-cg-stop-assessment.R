#' Assess CG stopping criteria after an outer iteration
#'
#' @noRd
.gl_assess_cg_stopping <- function(
    outer_log_lik_change,
    grad_inf,
    step_l2,
    stall_count,
    max_stall,
    raw_loglik_drop_from_best,
    raw_loglik_drop_tol,
    lambda_update_count,
    prevented_deterioration,
    outer_stop_crit,
    grad_tol,
    step_tol,
    stop_on_convergence = TRUE) {
  tolerance_met <- abs(outer_log_lik_change) <= outer_stop_crit &&
    is.finite(grad_inf) && grad_inf <= grad_tol &&
    is.finite(step_l2) && step_l2 <= step_tol
  max_stall_hit <- stall_count >= max_stall
  deterioration_hit <- is.finite(raw_loglik_drop_tol) &&
    lambda_update_count > 0L &&
    is.finite(raw_loglik_drop_from_best) &&
    raw_loglik_drop_from_best >= raw_loglik_drop_tol
  deterioration_hit <- isTRUE(deterioration_hit) || isTRUE(prevented_deterioration)
  stop_requested <- max_stall_hit || (isTRUE(stop_on_convergence) && tolerance_met) || deterioration_hit

  list(
    tolerance_met = isTRUE(tolerance_met),
    max_stall_hit = isTRUE(max_stall_hit),
    deterioration_hit = isTRUE(deterioration_hit),
    stop_requested = isTRUE(stop_requested)
  )
}

#' Select CG stop reason
#'
#' @noRd
.gl_cg_stop_reason <- function(tolerance_met, deterioration_hit) {
  if (isTRUE(tolerance_met)) {
    "tolerance"
  } else if (isTRUE(deterioration_hit)) {
    "raw_loglik_deterioration"
  } else {
    "max_stall"
  }
}

#' Determine whether CG convergence should wait for first lambda update
#'
#' @noRd
.gl_should_delay_cg_convergence_for_lambda_update <- function(update_lambda, has_smooths, lambda_update_count) {
  isTRUE(update_lambda) && isTRUE(has_smooths) && lambda_update_count == 0L
}

#' Handle a requested CG stop
#'
