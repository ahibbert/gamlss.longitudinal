#' Compute the primary fixed-effect vcov path
#'
#' @noRd
.gl_vcov_compute_primary <- function(object, par_cov, mm, margin_dist, response,
                                     response_margin, response_subject, method,
                                     progress, h) {
  vcov_final <- NULL
  se_final <- NULL
  hessian_diagnostics <- NULL
  hessian_nd <- NULL
  method_used <- method

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
      progress = progress
    )
  } else if (method %in% c("analytical", "analytical_only")) {
    .gl_source_analytical_hessian_helpers()

    analytical_hessian <- tryCatch(
      calc_analytical_hessian(object, progress = progress, h = h),
      error = function(e) structure(list(error = conditionMessage(e)), class = "vcov_hessian_error")
    )

    analytical_vcov <- if (inherits(analytical_hessian, "vcov_hessian_error")) {
      analytical_hessian
    } else {
      tryCatch(
        .gl_solve_hessian_vcov(analytical_hessian),
        error = function(e) structure(list(error = conditionMessage(e)), class = "vcov_hessian_error")
      )
    }

    if (inherits(analytical_vcov, "vcov_hessian_error")) {
      if (identical(method, "analytical_only")) {
        stop("Analytical Hessian vcov failed: ", analytical_vcov$error, call. = FALSE)
      }

      warning(
        "Analytical Hessian vcov failed; falling back to numerical Hessian. Reason: ",
        analytical_vcov$error,
        call. = FALSE
      )
      method_used <- "numderiv"
      hessian_nd <- calc_true_SE_numderiv_only_covariates(
        object = object,
        par = par_cov,
        mm = mm$x,
        margin_dist = margin_dist,
        response = response,
        testing = FALSE,
        response_margin = response_margin,
        response_subject = response_subject,
        progress = progress
      )
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
    hessian_nd = hessian_nd,
    vcov_final = vcov_final,
    se_final = se_final,
    hessian_diagnostics = hessian_diagnostics
  )
}
