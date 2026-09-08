#' Print a fitted longitudinal GAMLSS model
#'
#' @param x A fitted `gamlss.longitudinal` object.
#' @param digits Number of significant digits used for fit statistics.
#' @param ... Reserved for future methods.
#'
#' @return `x`, invisibly.
#' @export
print.gamlss.longitudinal <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!inherits(x, "gamlss.longitudinal")) {
    stop("`x` must be a fitted `gamlss.longitudinal` object.", call. = FALSE)
  }
  margin <- .gl_margin_display_label(x$margin_dist)
  copula <- .gl_copula_display_label(x$copula_dist)
  n_obs <- sum(is.finite(x$response))
  n_subjects <- length(unique(x$response_subject))
  n_time <- length(unique(x$response_margin))
  method <- x$optim_method %||% x$convergence$method
  if (is.null(method) && is.character(x$optimizer) && length(x$optimizer) == 1L) {
    method <- x$optimizer
  }
  method <- method %||% "unknown"
  optimization <- x$optimization %||%
    if (isTRUE(x$include_dlcopdpar)) "joint" else "separate"
  converged <- isTRUE(x$convergence$converged)
  iterations <- x$convergence$outer_iterations %||%
    x$convergence$outer_iter %||% max(0L, length(x$log_lik_history) - 1L)
  criteria <- x$model_selection
  metric <- function(name) {
    if (is.null(criteria) || !name %in% rownames(criteria) || !"joint" %in% colnames(criteria)) {
      return(NA_real_)
    }
    as.numeric(criteria[name, "joint"])
  }

  cat("\n<gamlss_longitudinal>\n")
  cat("Margin: ", margin, " | Copula: ", copula, "\n", sep = "")
  cat("Data: ", n_obs, " observations, ", n_subjects, " subjects, ", n_time,
      " time points\n", sep = "")
  cat("Optimizer: ", method, ", ", optimization, " | ",
      if (converged) "Converged" else "Not converged",
      if (is.finite(iterations)) paste0(" in ", iterations, " outer iterations") else "",
      "\n", sep = "")
  if (!converged) {
    cat("Stop reason: ", x$convergence$stop_reason %||% "unknown", "\n", sep = "")
  }
  cat(
    "logLik: ", .gl_format_console_number(metric("LogLik"), digits),
    " | AIC: ", .gl_format_console_number(metric("AIC"), digits),
    " | BIC: ", .gl_format_console_number(metric("BIC"), digits), "\n",
    sep = ""
  )
  inference <- if (isTRUE(x$vcov_meta$precomputed)) {
    x$vcov_meta$method_used %||% x$vcov_meta$method %||% "available"
  } else {
    x$vcov_meta$inference_status %||% "not computed"
  }
  cat("Inference: ", inference, "\n", sep = "")
  if (identical(x$likelihood_contract$objective, "segmented")) {
    cat("Missingness: segmented likelihood; observed segments are treated as independent.\n")
  }
  invisible(x)
}

#' @export
print.gamlss.longitudinal.control <- function(x, ...) {
  cat("\n<gamlss_longitudinal_control>\n")
  shared <- x$shared
  cat("Shared:\n")
  shared_display <- c(
    outer_tol = if (is.na(shared$outer_tol)) "auto" else shared$outer_tol,
    max_outer_iter = shared$max_outer_iter,
    max_elapsed_sec = if (is.infinite(shared$max_elapsed_sec)) "Inf" else shared$max_elapsed_sec,
    stop_on_convergence = shared$stop_on_convergence,
    lambda_start = if (is.null(shared$lambda_start) || is.na(shared$lambda_start)) "auto" else shared$lambda_start,
    lambda_penalty = shared$lambda_penalty %||% 2
  )
  cat(paste0("  ", names(shared_display), ": ", shared_display, collapse = "\n"), "\n", sep = "")
  specified <- attr(x, "specified")
  show_section <- function(label, values, keys) {
    if (!length(keys)) {
      cat(label, ": defaults\n", sep = "")
      return()
    }
    rendered <- vapply(keys, function(key) paste(values[[key]], collapse = ","), character(1))
    cat(label, " overrides:\n", sep = "")
    cat(paste0("  ", keys, ": ", rendered, collapse = "\n"), "\n", sep = "")
  }
  show_section("RS", x$rs, specified$rs)
  show_section("CG", x$cg, specified$cg)
  invisible(x)
}

#' @export
print.gamlss_longitudinal_inference_control <- function(x, ...) {
  cat("\n<gamlss_longitudinal_inference_control>\n")
  cat("Profile: ", x$profile, " | Defaults: ", x$defaults_version, "\n", sep = "")
  overrides <- x$expert_override %||% list()
  cat("Expert overrides: ", if (length(overrides)) paste(names(overrides), collapse = ", ") else "none", "\n", sep = "")
  invisible(x)
}
