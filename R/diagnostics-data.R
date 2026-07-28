.gl_diag_data <- function(object) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("Diagnostics are only available for objects of class 'gamlss.longitudinal'.")
  }

  copula_link <- get_copula_dist(object$copula_dist)$copula_link

  eta_out <- calc_eta(
    par_cov = object$par,
    mm = object$model_matrix,
    margin_dist = object$margin_dist,
    copula_link = copula_link,
    par_s = object$par_s
  )

  margin_params <- eta_out$eta_inv[names(object$margin_dist$parameters)]

  y <- object$response

  keep <- is.finite(y)

  for (par_name in names(margin_params)) {
    keep <- keep & is.finite(margin_params[[par_name]])
  }

  margin_params <- lapply(margin_params, function(x) x[keep])

  y <- y[keep]

  mu_hat <- if ("mu" %in% names(margin_params)) margin_params$mu else margin_params[[1]]

  sigma_hat <- if ("sigma" %in% names(margin_params)) margin_params$sigma else rep(stats::sd(y, na.rm = TRUE), length(y))

  sigma_hat <- pmax(as.numeric(sigma_hat), .Machine$double.eps)

  list(
    response = y,
    params = margin_params,
    mu_hat = as.numeric(mu_hat),
    sigma_hat = sigma_hat,
    family = object$margin_dist$family[1],
    subject = object$response_subject[keep],
    time = object$response_margin[keep]
  )
}
