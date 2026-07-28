#' Report CG optimizer startup controls
#'
#' @noRd
.gl_report_cg_optimizer_start <- function(
    verbose,
    cg_max_delta,
    cg_lambda_update_every,
    cg_update_lambda,
    cg_line_search,
    cg_gradient_method,
    cg_hessian_method) {
  if (verbose > 0) {
    cat("\nUsing optimization method: CG")

    cat(paste0(
      "\nCG controls: max_delta=", signif(cg_max_delta, 4),
      " | lambda_update_every=", cg_lambda_update_every,
      " | update_lambda=", isTRUE(cg_update_lambda),
      " | line_search=", cg_line_search,
      " | gradient=", cg_gradient_method,
      " | hessian=", cg_hessian_method,
      "\n"
    ))
  }

  invisible(NULL)
}
