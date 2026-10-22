#' Check whether a fitted object can reuse its stored vcov output
#'
#' @noRd
.gl_vcov_cache_version <- function() 2L

#' Continuous-curvature source version, separate from discrete caches
#' @noRd
.gl_vcov_continuous_curvature_version <- function() 2L

#' @noRd
.gl_vcov_continuous_curvature_applies <- function(object) {
  family <- object$copula_dist
  if (!is.character(family) || length(family) != 1L || is.na(family) ||
      !family %in% c("N", "C", "F", "G", "J", "t")) return(FALSE)
  discrete <- tryCatch(.is_discrete_margin(object$margin_dist), error = function(e) FALSE)
  !isTRUE(discrete)
}

#' @noRd
.gl_vcov_uses_analytical_curvature <- function(method, diagnostics = NULL) {
  key <- .gl_vcov_method_key(method)
  if (identical(key, "numderiv")) return(FALSE)
  if (identical(key, "sandwich") &&
      identical(.gl_vcov_method_key(diagnostics$bread_method), "numderiv")) return(FALSE)
  TRUE
}

#' @noRd
.gl_vcov_continuous_curvature_marker <- function(object) {
  marker <- list(version = .gl_vcov_continuous_curvature_version(),
    backend = .copula_backend(), copula = object$copula_dist,
    boundary_version = 2L)
  if (identical(object$copula_dist, "t") && identical(marker$backend, "native")) {
    marker$df_branch_version <- 1L
  }
  marker
}

#' @noRd
.gl_vcov_method_key <- function(method) {
  method <- as.character(method %||% "")[[1L]]
  if (grepl("sandwich", method, fixed = TRUE)) return("sandwich")
  if (grepl("numderiv|numerical", method)) return("numderiv")
  method
}

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

  stored_method <- object$vcov$method %||% object$vcov_meta$method_used %||%
    object$vcov_meta$method
  if (.gl_vcov_continuous_curvature_applies(object) &&
      .gl_vcov_uses_analytical_curvature(stored_method, object$vcov$hessian_diagnostics) &&
      !identical(object$vcov$continuous_curvature,
        .gl_vcov_continuous_curvature_marker(object))) {
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

    cached_methods <- c(
      object$vcov_meta$method,
      object$vcov_meta$method_used,
      object$vcov$method
    )
    cached_method_keys <- unique(vapply(
      cached_methods,
      .gl_vcov_method_key,
      character(1)
    ))
    method_ok <- is.null(method) ||
      is.null(object$vcov_meta$method) ||
      .gl_vcov_method_key(method) %in% cached_method_keys

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
  stored_diagnostics <- object$vcov_meta$hessian_diagnostics %||% list()
  stored_codes <- stored_diagnostics$failure_codes %||% character()
  stored_is_advisory <- length(stored_codes) > 0L &&
    length(.gl_hessian_blocking_failure_codes(stored_codes)) == 0L
  stored_method <- object$vcov_meta$method %||% "analytical"
  explicit_revalidation <- !is.null(extra_args$inference) ||
    !identical(as.character(vcov_method)[1], as.character(stored_method)[1])
  if (stored_unavailable && !stored_is_advisory && !explicit_revalidation) {
    stop(.gl_inference_unavailable(stored_diagnostics))
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
