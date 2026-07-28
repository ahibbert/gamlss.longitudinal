.gl_prediction_model_matrix <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    return(object$model_matrix)
  }

  copula_link <- get_copula_dist(object$copula_dist)$copula_link

  nd <- .gl_prepare_newdata_internal(object, newdata, require_response = FALSE)

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

  .gl_align_model_matrix_columns(mm_use, object$model_matrix)
}

.gl_predict_response_se <- function(object, newdata = NULL, method = "analytical", ...) {
  mm_use <- .gl_prediction_model_matrix(object, newdata = newdata)

  if (is.null(mm_use$x$mu) || ncol(mm_use$x$mu) == 0L) {
    return(rep(NA_real_, length(predict(object, newdata = newdata, type = "response"))))
  }

  X <- as.matrix(mm_use$x$mu)

  beta_names <- colnames(X)

  beta_names <- ifelse(startsWith(beta_names, "mu."), beta_names, paste0("mu.", beta_names))

  vc <- .resolve_vcov(object, extra_args = list(method = method, ...))

  V <- vc$vcov$overall

  V_names <- colnames(V) %||% rownames(V)

  idx <- match(beta_names, V_names)

  if (any(is.na(idx))) {
    return(rep(NA_real_, nrow(X)))
  }

  V_mu <- V[idx, idx, drop = FALSE]

  se_eta <- sqrt(pmax(0, rowSums((X %*% V_mu) * X)))

  copula_link <- get_copula_dist(object$copula_dist)$copula_link

  eta_out <- calc_eta(
    par_cov = object$par,
    mm = mm_use,
    margin_dist = object$margin_dist,
    copula_link = copula_link,
    par_s = object$par_s
  )

  mu_dr <- eta_out$eta_dr$mu %||% rep(1, length(se_eta))

  as.numeric(abs(mu_dr[seq_along(se_eta)]) * se_eta)
}
