.gl_fitted_distribution <- function(object, newdata = NULL, require_response = TRUE) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("Diagnostics are only available for objects of class 'gamlss.longitudinal'.")
  }

  .gl_validate_capability_route(
    margin_dist = object$margin_dist,
    copula_dist = object$copula_dist,
    context = "diagnostics"
  )

  copula_link <- get_copula_dist(object$copula_dist)$copula_link

  if (is.null(newdata)) {
    mm_use <- object$model_matrix

    response <- object$response

    response_margin <- object$response_margin

    response_subject <- object$response_subject
  } else {
    nd <- .gl_prepare_newdata_internal(object, newdata, require_response = require_response)

    mm_use <- do.call(
      create_model_matrices,
      list(
        mu.formula = object$formulas_int$mu,
        sigma.formula = object$formulas_int$sigma,
        nu.formula = object$formulas_int$nu,
        tau.formula = object$formulas_int$tau,
        theta.formula = object$formulas_int$theta,
        zeta.formula = object$formulas_int$zeta,
        margin.family = object$margin_dist,
        copula.family = object$copula_dist,
        copula.link = copula_link,
        dataset = nd,
        quiet_gamlss2 = TRUE,
        preserve_factor_levels = TRUE
      )
    )

    mm_use <- .gl_align_model_matrix_columns(mm_use, object$model_matrix)

    response <- nd$response

    response_margin <- nd$time

    response_subject <- nd$subject
  }

  eta_out <- calc_eta(
    par_cov = object$par,
    mm = mm_use,
    margin_dist = object$margin_dist,
    copula_link = copula_link,
    par_s = object$par_s
  )

  margin_params <- eta_out$eta_inv[names(object$margin_dist$parameters)]

  keep <- rep(TRUE, length(response))

  if (require_response) {
    keep <- is.finite(response)
  }

  for (par_name in names(margin_params)) {
    keep <- keep & is.finite(margin_params[[par_name]])
  }

  common_n <- min(
    length(response),
    length(response_margin),
    length(response_subject),
    if (length(margin_params) > 0) min(vapply(margin_params, length, integer(1))) else length(response)
  )

  if (!is.finite(common_n) || common_n < 0) {
    common_n <- 0L
  }

  response <- response[seq_len(common_n)]

  response_margin <- response_margin[seq_len(common_n)]

  response_subject <- response_subject[seq_len(common_n)]

  margin_params <- lapply(margin_params, function(x) x[seq_len(common_n)])

  keep <- keep[seq_len(common_n)]

  margin_params <- lapply(margin_params, function(x) x[keep])

  response <- response[keep]

  response_margin <- response_margin[keep]

  response_subject <- response_subject[keep]

  mu_hat <- if ("mu" %in% names(margin_params)) margin_params$mu else margin_params[[1]]

  sigma_hat <- if ("sigma" %in% names(margin_params)) margin_params$sigma else rep(stats::sd(response, na.rm = TRUE), length(response))

  sigma_hat <- pmax(as.numeric(sigma_hat), .Machine$double.eps)

  list(
    response = response,
    params = margin_params,
    mu_hat = as.numeric(mu_hat),
    sigma_hat = sigma_hat,
    family = object$margin_dist$family[1],
    subject = response_subject,
    time = response_margin,
    keep_mask = keep,
    keep_index = which(keep)
  )
}
