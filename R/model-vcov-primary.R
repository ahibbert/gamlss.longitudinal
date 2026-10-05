#' Compute the primary fixed-effect vcov path
#'
#' @noRd
.gl_vcov_compute_primary <- function(object, par_cov, mm, margin_dist, response,
                                     response_margin, response_subject, method,
                                     progress, h,
                                     inference = inference_control("standard")) {
  vcov_final <- NULL
  se_final <- NULL
  hessian_diagnostics <- NULL
  hessian_nd <- NULL
  method_used <- method
  hessian_fallback <- NULL
  reference_hessian <- NULL

  if (method == "numderiv") {
    hessian_nd <- calc_true_SE_numderiv_only_covariates(
      object = object,
      par = par_cov,
      mm = mm$x,
      margin_dist = margin_dist,
      response = response,
      testing = FALSE,
      response_margin = response_margin,
      response_subject = response_subject,
      h = h,
      progress = progress
    )
  } else if (method %in% c("analytical", "analytical_only")) {
    .gl_source_analytical_hessian_helpers()

    analytical_hessian <- tryCatch(
      calc_analytical_hessian(object, progress = progress, h = h),
      error = function(e) structure(list(error = conditionMessage(e)), class = "vcov_hessian_error")
    )

    if (!inherits(analytical_hessian, "vcov_hessian_error") &&
        isTRUE(inference$check_agreement)) {
      reference_hessian <- calc_true_SE_numderiv_only_covariates(
        object = object, par = par_cov, mm = mm$x,
        margin_dist = margin_dist, response = response, testing = FALSE,
        response_margin = response_margin,
        response_subject = response_subject, h = h, progress = progress
      )
    }

    analytical_vcov <- if (inherits(analytical_hessian, "vcov_hessian_error")) {
      analytical_hessian
    } else {
      tryCatch(
        {
          gradient <- .gl_vcov_fitted_gradient(
            object, par_cov, analytical_hessian, inference$gradient_step
          )
          .gl_solve_hessian_vcov(
            analytical_hessian,
            parameter_names = names(par_cov),
            inference = inference,
            source = "analytical",
            gradient = gradient,
            reference_hessian = reference_hessian
          )
        },
        error = function(e) structure(
          list(error = conditionMessage(e), condition = e),
          class = "vcov_hessian_error"
        )
      )
    }

    if (inherits(analytical_vcov, "vcov_hessian_error")) {
      if (identical(method, "analytical_only")) {
        if (inherits(analytical_vcov$condition, "gamlss_longitudinal_inference_unavailable")) {
          stop(analytical_vcov$condition)
        }
        stop("Analytical Hessian vcov failed: ", analytical_vcov$error, call. = FALSE)
      }

      warning(
        "Analytical Hessian vcov failed; falling back to numerical Hessian. Reason: ",
        analytical_vcov$error,
        call. = FALSE
      )
      # For an analytical-Hessian validation failure, use a deliberately
      # less cancellation-prone finite-difference scale. The standard
      # numerical step (1e-4) can yield an indefinite observed-information
      # matrix for large, zero-heavy discrete models even when every diagonal
      # curvature term is positive. This is not an acceptance relaxation:
      # the full numerical Hessian is still passed through exactly the same
      # positive-definiteness, rank, covariance, and stationarity gates.
      numerical_hessian_step <- max(h, 1e-3)
      method_used <- "numderiv"
      hessian_fallback <- list(
        from = "analytical",
        to = "numderiv",
        reason = analytical_vcov$error,
        diagnostics = analytical_vcov$condition$diagnostics %||% NULL,
        numerical_hessian_step = numerical_hessian_step
      )
      hessian_nd <- calc_true_SE_numderiv_only_covariates(
        object = object,
        par = par_cov,
        mm = mm$x,
        margin_dist = margin_dist,
        response = response,
        testing = FALSE,
        response_margin = response_margin,
        response_subject = response_subject,
        h = numerical_hessian_step,
        progress = progress
      )
      numerical_gradient <- .gl_vcov_fitted_gradient(
        object, par_cov, hessian_nd, inference$gradient_step
      )
      numerical_vcov <- tryCatch(
        .gl_solve_hessian_vcov(
          hessian_nd,
          parameter_names = names(par_cov),
          inference = inference,
          source = "numderiv",
          fallback = hessian_fallback,
          gradient = numerical_gradient
        ),
        error = function(e) structure(
          list(error = conditionMessage(e), condition = e),
          class = "vcov_hessian_error"
        )
      )
      if (!inherits(numerical_vcov, "vcov_hessian_error")) {
        vcov_final <- numerical_vcov$vcov
        se_final <- numerical_vcov$se
        hessian_diagnostics <- numerical_vcov$hessian_diagnostics
      }
    } else {
      method_used <- "analytical"
      hessian_nd <- analytical_hessian
      vcov_final <- analytical_vcov$vcov
      se_final <- analytical_vcov$se
      hessian_diagnostics <- analytical_vcov$hessian_diagnostics
    }
  }

  list(
    method_used = method_used,
    hessian_fallback = hessian_fallback,
    reference_hessian = reference_hessian,
    hessian_nd = hessian_nd,
    vcov_final = vcov_final,
    se_final = se_final,
    hessian_diagnostics = hessian_diagnostics
  )
}
