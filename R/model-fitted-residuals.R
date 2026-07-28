.gl_fitted_parameter_values <- function(object, parameter = "mu") {
  copula_link <- get_copula_dist(object$copula_dist)$copula_link
  eta_out <- calc_eta(
    par_cov = object$par,
    mm = object$model_matrix,
    margin_dist = object$margin_dist,
    copula_link = copula_link,
    par_s = object$par_s
  )
  params <- eta_out$eta_inv[names(object$margin_dist$parameters)]
  if (!parameter %in% names(params)) {
    stop("Parameter '", parameter, "' is not available in the fitted margin.", call. = FALSE)
  }
  as.numeric(params[[parameter]])
}

#' @rdname gamlss_longitudinal_accessors
#' @export
fitted.gamlss.longitudinal <- function(object, parameter = "mu", finite = FALSE, ...) {
  fit <- .gl_fitted_parameter_values(object, parameter = parameter)
  if (isTRUE(finite)) {
    keep <- is.finite(object$response) & is.finite(fit)
    fit <- fit[keep]
  }
  fit
}

#' @rdname gamlss_longitudinal_accessors
#' @export
residuals.gamlss.longitudinal <- function(object, type = c("response", "pearson", "quantile"), finite = TRUE, ...) {
  type <- match.arg(type)
  mu <- .gl_fitted_parameter_values(object, parameter = "mu")
  y <- as.numeric(object$response)
  keep <- is.finite(y) & is.finite(mu)

  if (identical(type, "response")) {
    out <- y - mu
  } else if (identical(type, "pearson")) {
    sigma <- if ("sigma" %in% names(object$margin_dist$parameters)) {
      .gl_fitted_parameter_values(object, parameter = "sigma")
    } else {
      rep(stats::sd(y, na.rm = TRUE), length(y))
    }
    sigma <- pmax(as.numeric(sigma), .Machine$double.eps)
    out <- (y - mu) / sigma
    keep <- keep & is.finite(sigma)
  } else {
    copula_link <- get_copula_dist(object$copula_dist)$copula_link
    eta_out <- calc_eta(
      par_cov = object$par,
      mm = object$model_matrix,
      margin_dist = object$margin_dist,
      copula_link = copula_link,
      par_s = object$par_s
    )
    params <- eta_out$eta_inv[names(object$margin_dist$parameters)]
    pit <- .gl_call_family_fun("p", object$margin_dist$family[1], y, params)
    pit <- pmin(pmax(pit, .Machine$double.eps), 1 - .Machine$double.eps)
    out <- stats::qnorm(pit)
    keep <- keep & is.finite(out)
  }

  if (isTRUE(finite)) {
    out <- out[keep]
  }
  out
}
