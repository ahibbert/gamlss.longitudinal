#' Check whether a fitted object can reuse its stored vcov output
#'
#' @noRd
.gl_vcov_cache_version <- function() 1L

#' @noRd
.can_use_cached_vcov <- function(object, numderiv = FALSE, method = NULL, extra_args = list()) {
  if (!inherits(object, "gamlss.longitudinal")) {
    return(FALSE)
  }

  extra_args_cache <- extra_args

  extra_args_cache$method <- NULL

  if (!is.null(extra_args_cache) && length(extra_args_cache) > 0) {
    return(FALSE)
  }

  if (is.null(object$vcov) || !is.list(object$vcov)) {
    return(FALSE)
  }

  if (!identical(object$vcov_meta$cache_version, .gl_vcov_cache_version())) {
    return(FALSE)
  }

  if (identical(object$vcov$hessian_diagnostics$status, "unavailable")) {
    return(FALSE)
  }

  if (is.null(object$vcov$vcov) || is.null(object$vcov$vcov$overall)) {
    return(FALSE)
  }

  if (!is.null(object$vcov_meta) && !is.null(object$vcov_meta$numderiv)) {
    numderiv_ok <- identical(isTRUE(object$vcov_meta$numderiv), isTRUE(numderiv))

    method_ok <- is.null(method) ||

      is.null(object$vcov_meta$method) ||

      identical(as.character(object$vcov_meta$method)[1], as.character(method)[1])

    return(numderiv_ok && method_ok)
  }

  TRUE
}

#' Resolve cached or recomputed vcov output for downstream summaries
#'
#' @noRd
.resolve_vcov <- function(object, numderiv = FALSE, extra_args = list()) {
  vcov_method <- extra_args$method %||% if (isTRUE(numderiv)) "numderiv" else "analytical"

  stored_unavailable <- identical(object$vcov_meta$inference_status, "unavailable")
  stored_method <- object$vcov_meta$method %||% "analytical"
  explicit_revalidation <- !is.null(extra_args$inference) ||
    !identical(as.character(vcov_method)[1], as.character(stored_method)[1])
  if (stored_unavailable && !explicit_revalidation) {
    stop(.gl_inference_unavailable(object$vcov_meta$hessian_diagnostics))
  }

  if (.can_use_cached_vcov(object, numderiv = numderiv, method = vcov_method, extra_args = extra_args)) {
    return(.gl_enrich_vcov_contract(object$vcov, object))
  }

  if (is.null(extra_args$method)) {
    extra_args$method <- vcov_method
  }

  do.call(
    vcov.gamlss.longitudinal,
    c(list(object = object, numderiv = numderiv, details = TRUE), extra_args)
  )
}
