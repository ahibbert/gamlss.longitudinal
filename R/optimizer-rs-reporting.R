.gl_print_outer_iteration_summary <- function(
    outer_start_log_lik,
    outer_end_log_lik,
    outer_log_lik_change,
    verbose,
    cat_fn = cat,
    print_fn = print) {
  out_temp <- c(outer_start_log_lik, outer_end_log_lik, outer_log_lik_change)
  names(out_temp) <- c("Start LogLik", "End LogLik", "Change")

  if (verbose > 0) {
    cat_fn("\n")
    print_fn(out_temp)
  }

  invisible(out_temp)
}


#' Print the outer-convergence message when tolerance is satisfied
#'
#' @noRd
.gl_print_outer_convergence <- function(
    outer_log_lik_change,
    outer_stop_crit,
    verbose,
    cat_fn = cat,
    print_fn = print) {
  converged <- abs(outer_log_lik_change) <= outer_stop_crit

  if (isTRUE(converged) && verbose > 0) {
    print_fn(c(outer_log_lik_change))
    cat_fn("\nOUTER CONVERGED")
  }

  invisible(converged)
}
