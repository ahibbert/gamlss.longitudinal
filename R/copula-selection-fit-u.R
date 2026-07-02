.select_copula_u_from_fit <- function(object) {
  copula_link <- get_copula_dist(object$copula_dist)$copula_link
  eta_out <- calc_eta(
    object$par,
    object$model_matrix,
    object$margin_dist,
    copula_link,
    object$par_s
  )
  u <- calc_F_x(
    eta_out$eta_inv,
    object$model_matrix$x,
    object$margin_dist,
    object$response
  )
  .copula_clamp01(u)
}
