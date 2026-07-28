#' Detect whether CG smoothing parameters changed
#'
#' @noRd
.gl_cg_lambdas_changed <- function(
    lambda_before,
    lambda_after,
    tolerance = 1e-12) {
  !isTRUE(all.equal(
    unlist(lambda_before, use.names = TRUE),
    unlist(lambda_after, use.names = TRUE),
    tolerance = tolerance,
    check.attributes = FALSE
  ))
}
