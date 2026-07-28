#' Build CG lambda candidate trace rows
#'
#' @noRd
.gl_build_cg_lambda_trace_rows <- function(
    outer_iteration,
    parameter,
    smooth,
    lambda_before,
    candidates,
    lambda_scores,
    best) {
  data.frame(
    outer_iteration = as.integer(outer_iteration),
    parameter = parameter,
    smooth = smooth,
    lambda_before = as.numeric(lambda_before),
    lambda_candidate = candidates,
    raw_logLik_after_step = lambda_scores$raw_loglik,
    smooth_penalty_after_step = lambda_scores$penalty_value,
    penalized_logLik_after_step = lambda_scores$penalized_loglik,
    edf_after_step = lambda_scores$edf_values,
    gaic_score = lambda_scores$gaic_score,
    chosen = seq_along(candidates) == best,
    row.names = NULL
  )
}
