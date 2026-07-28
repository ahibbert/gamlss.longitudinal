#' Compute the CG gradient for one iteration
#'
#' @noRd
.gl_compute_cg_gradient <- function(
    beta_vec,
    mm_cg,
    eval_start,
    margin_dist,
    copula_dist,
    include_dlcopdpar,
    dataset,
    cg_gradient_method,
    finite_gradient_fn,
    analytical_gradient_fn = .cg_analytical_gradient) {
  if (identical(cg_gradient_method, "analytical")) {
    analytical_gradient_fn(
      beta_vec,
      mm_cg,
      eval_start$eta_out,
      eval_start$calc_lik,
      margin_dist,
      copula_dist,
      include_dlcopdpar,
      dataset$response,
      dataset$time,
      dataset$subject
    )
  } else {
    finite_gradient_fn(beta_vec, eval_start$loglik, mm_cg)
  }
}
