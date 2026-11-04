#' Build the model object passed to CG Hessian routines
#'
#' @noRd
.gl_build_cg_hessian_object <- function(
    dataset,
    margin_dist,
    copula_dist,
    mm_cg,
    beta_vec) {
  list(
    response = dataset$response,
    response_margin = dataset$time,
    response_subject = dataset$subject,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    model_matrix = mm_cg,
    par = beta_vec,
    par_s = setNames(lapply(names(mm_cg$x), function(x) list()), names(mm_cg$x))
  )
}
