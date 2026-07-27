#' Prepare object fields and eta values for vcov calculation
#'
#' @noRd
.gl_prepare_vcov_evaluation <- function(
    object,
    par,
    numderiv,
    method,
    progress) {
  if (length(method) > 1L) {
    method <- method[[1L]]
  }
  method <- match.arg(method, c("analytical", "numderiv", "analytical_only", "sandwich"))
  if (isTRUE(numderiv)) method <- "numderiv"
  method_requested <- method
  method_used <- method

  progress <- isTRUE(progress)

  include_dlcopdpar <- TRUE

  response <- object$response

  response_margin <- object$response_margin

  response_subject <- object$response_subject

  margin_names <- unique(object$response_margin)

  num_margins <- length(margin_names)

  # se_out=object$par*0;

  margin_dist <- object$margin_dist
  copula_dist <- object$copula_dist
  copula_link <- get_copula_dist(copula_dist)$copula_link

  mm <- object$model_matrix

  if (all(is.na(par))) {
    par_cov <- object$par

    par_s <- object$par_s
  } else {
    par_cov <- par$par

    par_s <- par$par_s
  }

  eta_out <- calc_eta(par_cov, mm, margin_dist, copula_link, par_s = par_s)

  eta_inv <- eta_out$eta_inv
  eta_dr <- eta_out$eta_dr
  eta <- eta_out$eta
  eta_dr <- eta_out$eta_dr

  list(
    method = method,
    method_requested = method_requested,
    method_used = method_used,
    progress = progress,
    include_dlcopdpar = include_dlcopdpar,
    response = response,
    response_margin = response_margin,
    response_subject = response_subject,
    margin_names = margin_names,
    num_margins = num_margins,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    copula_link = copula_link,
    mm = mm,
    par_cov = par_cov,
    par_s = par_s,
    eta_out = eta_out,
    eta_inv = eta_inv,
    eta_dr = eta_dr,
    eta = eta
  )
}
