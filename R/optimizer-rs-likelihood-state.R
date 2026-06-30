#' Assemble the RS likelihood state returned for one inner iteration
#'
#' @noRd
.gl_build_rs_iteration_likelihood_state <- function(
    eta_out,
    eta,
    eta_dr,
    eta_inv,
    likelihood_context,
    outer_start_state,
    timer,
    history_state) {
  c(
    list(
      eta_out = eta_out,
      eta = eta,
      eta_dr = eta_dr,
      eta_inv = eta_inv
    ),
    likelihood_context,
    list(
      first_outer_run = outer_start_state$first_outer_run,
      outer_start_log_lik = outer_start_state$outer_start_log_lik,
      timer = timer,
      log_lik_history = history_state$log_lik_history,
      par_history = history_state$par_history
    )
  )
}
