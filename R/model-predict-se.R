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
    stop(.gl_inference_unavailable(list(
      status = "unavailable",
      failure_codes = "prediction_mu_design_unavailable",
      message = "Prediction inference is unavailable because the mu fixed-effect design is absent."
    )))
  }

  X <- as.matrix(mm_use$x$mu)

  beta_names <- colnames(X)

  beta_names <- ifelse(startsWith(beta_names, "mu."), beta_names, paste0("mu.", beta_names))

  vc <- .resolve_vcov(object, extra_args = list(method = method, ...))
  .gl_require_available_inference(vc)

  V <- vc$vcov$overall

  V_names <- colnames(V) %||% rownames(V)

  idx <- match(beta_names, V_names)

  if (any(is.na(idx))) {
    stop(.gl_inference_unavailable(list(
      status = "unavailable",
      failure_codes = "prediction_covariance_columns_unmatched",
      message = paste0(
        "Prediction inference is unavailable because mu design columns are missing ",
        "from the fixed-effect covariance matrix: ",
        paste(beta_names[is.na(idx)], collapse = ", "), "."
      ),
      unmatched_columns = beta_names[is.na(idx)]
    )))
  }

  V_mu <- V[idx, idx, drop = FALSE]

  se_eta <- .gl_sqrt_derived_variance(
    rowSums((X %*% V_mu) * X), "prediction covariance", allow_zero = TRUE
  )

  copula_link <- get_copula_dist(object$copula_dist)$copula_link

  eta_out <- calc_eta(
    par_cov = object$par,
    mm = mm_use,
    margin_dist = object$margin_dist,
    copula_link = copula_link,
    par_s = object$par_s
  )

  mu_dr <- eta_out$eta_dr$mu %||% rep(1, length(se_eta))

  out <- as.numeric(abs(mu_dr[seq_along(se_eta)]) * se_eta)
  if (!any(is.finite(out))) {
    stop(.gl_inference_unavailable(list(
      status = "unavailable",
      failure_codes = "prediction_standard_errors_nonfinite",
      message = "Prediction inference is unavailable because no finite standard errors were produced."
    )))
  }
  covariance_contract <- vc$inference_contract %||%
    .gl_fixed_inference_contract(vc, coefficient_names = beta_names)
  attr(out, "inference_contract") <- covariance_contract
  out
}
