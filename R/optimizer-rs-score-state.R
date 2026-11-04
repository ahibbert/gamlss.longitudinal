#' Assemble the RS parameter-score state returned for one inner iteration
#'
#' @noRd
.gl_build_rs_parameter_score_state <- function(
    copula_context,
    discrete_scores,
    score_assembly,
    timer) {
  c(
    copula_context,
    list(
      discrete_scores = discrete_scores,
      score_assembly = score_assembly,
      d1 = score_assembly$d1,
      d1_m = score_assembly$d1_m,
      d1_cop = score_assembly$d1_cop,
      timer = timer
    )
  )
}
