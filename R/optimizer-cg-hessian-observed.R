.gl_cg_hessian_ok <- function(H, beta_names) {
  is.matrix(H) &&
    identical(dim(H), c(length(beta_names), length(beta_names))) &&
    all(is.finite(H))
}

#' Construct observed Hessian for CG optimization
#'
#' @noRd
.gl_cg_observed_hessian <- function(
    tmp_obj,
    beta_vec,
    mm_cg,
    analytical_fn,
    finite_fn,
    hessian_method,
    verbose,
    context = "CG iteration") {
  beta_names <- names(beta_vec)
  use_finite <- identical(hessian_method, "finite")
  H <- NULL
  analytical_error <- NULL

  if (!use_finite) {
    H <- tryCatch(
      analytical_fn(tmp_obj),
      error = function(e) {
        analytical_error <<- conditionMessage(e)
        NULL
      }
    )
    if (.gl_cg_hessian_ok(H, beta_names)) {
      return(0.5 * (H + t(H)))
    }
    if (verbose > 0) {
      msg <- if (!is.null(analytical_error)) analytical_error else "non-finite analytical Hessian"
      warning(
        "CG analytical Hessian failed during ", context,
        "; falling back to finite-difference Hessian. Reason: ", msg,
        call. = FALSE
      )
    }
  }

  H_fd <- finite_fn(beta_vec, beta_names, mm_cg)
  if (!.gl_cg_hessian_ok(H_fd, beta_names)) {
    msg <- if (!is.null(analytical_error)) analytical_error else "finite-difference Hessian was non-finite"
    stop("CG failed to construct a usable Hessian during ", context, ". Reason: ", msg, call. = FALSE)
  }
  0.5 * (H_fd + t(H_fd))
}

#' Compute CG smooth effective degrees of freedom
#'
