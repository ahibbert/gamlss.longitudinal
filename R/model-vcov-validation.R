#' Construct a classed inference-unavailable error
#'
#' @noRd
.gl_inference_unavailable <- function(diagnostics) {
  codes <- diagnostics$failure_codes %||% "unknown_inference_failure"
  message <- diagnostics$message %||% paste0(
    "Hessian-based inference is unavailable (",
    paste(codes, collapse = ", "), "). The fitted model is unchanged."
  )
  structure(
    list(message = message, call = NULL, diagnostics = diagnostics),
    class = c(
      "gamlss_longitudinal_inference_unavailable",
      "error", "condition"
    )
  )
}

#' Curvature-scaled finite-difference score at the fitted coefficients
#'
#' @noRd
.gl_vcov_fitted_gradient <- function(object, par_cov, H, step) {
  if (!is.matrix(H) || nrow(H) != ncol(H) || any(!is.finite(H))) {
    return(list(
      gradient = stats::setNames(rep(NA_real_, length(par_cov)), names(par_cov)),
      steps = rep(NA_real_, length(par_cov)),
      scaled_max = Inf
    ))
  }
  info_diag <- diag(-(H + t(H)) / 2)
  if (length(info_diag) != length(par_cov) || any(!is.finite(info_diag)) ||
      any(info_diag <= 0)) {
    return(list(
      gradient = stats::setNames(rep(NA_real_, length(par_cov)), names(par_cov)),
      steps = rep(NA_real_, length(par_cov)),
      scaled_max = Inf
    ))
  }
  steps <- step / sqrt(info_diag)
  pair_cache <- build_copula_pair_cache(
    object$response, object$response_margin, object$response_subject
  )
  copula_link <- get_copula_dist(object$copula_dist)$copula_link
  eval_ll <- function(beta) {
    eta <- calc_eta(beta, object$model_matrix, object$margin_dist,
                    copula_link, par_s = object$par_s)
    out <- calc_likelihood_minimal(
      eta$eta_inv,
      mm = object$model_matrix$x,
      margin_dist = object$margin_dist,
      copula_dist = object$copula_dist,
      calc_d2 = FALSE,
      response = object$response,
      response_margin = object$response_margin,
      response_subject = object$response_subject,
      pair_cache = pair_cache,
      calc_margin_deriv = FALSE
    )
    as.numeric(out$log_lik["joint"])
  }
  gradient <- rep(NA_real_, length(par_cov))
  names(gradient) <- names(par_cov)
  for (j in seq_along(par_cov)) {
    plus <- minus <- par_cov
    plus[[j]] <- plus[[j]] + steps[[j]]
    minus[[j]] <- minus[[j]] - steps[[j]]
    gradient[[j]] <- (eval_ll(plus) - eval_ll(minus)) / (2 * steps[[j]])
  }
  list(
    gradient = gradient,
    steps = steps,
    scaled_max = max(abs(gradient) / sqrt(info_diag))
  )
}

#' Refuse downstream inference when validation has failed
#'
#' @noRd
.gl_require_available_inference <- function(x) {
  diagnostics <- x$hessian_diagnostics %||% NULL
  if (identical(diagnostics$status, "unavailable")) {
    stop(.gl_inference_unavailable(diagnostics))
  }
  invisible(TRUE)
}

#' Square root derived variances without silently masking negative values
#'
#' @noRd
.gl_sqrt_derived_variance <- function(x, context, allow_zero = TRUE) {
  bad <- !is.finite(x) | x < 0 | (!isTRUE(allow_zero) & x <= 0)
  if (any(bad)) {
    diagnostics <- list(
      status = "unavailable",
      failure_codes = "derived_variance_nonpositive",
      context = context,
      derived_variance = x
    )
    stop(.gl_inference_unavailable(diagnostics))
  }
  sqrt(x)
}

#' Validate log-likelihood curvature and construct covariance output
#'
#' @noRd
.gl_validate_hessian_inference <- function(
    H,
    parameter_names = NULL,
    control = inference_control("standard"),
    source = "unknown",
    fallback = NULL,
    gradient = NULL,
    reference_hessian = NULL) {
  control <- .gl_normalize_inference_control(control)
  failures <- character()
  add_failure <- function(code) failures <<- unique(c(failures, code))

  is_square <- is.matrix(H) && nrow(H) == ncol(H) && nrow(H) > 0L
  finite <- is_square && all(is.finite(H))
  conformable <- is_square && (is.null(parameter_names) || nrow(H) == length(parameter_names))
  matrix_names <- if (is_square) dimnames(H) else NULL
  named <- is_square && !is.null(matrix_names[[1]]) && !is.null(matrix_names[[2]]) &&
    !anyNA(matrix_names[[1]]) && !anyNA(matrix_names[[2]]) &&
    identical(matrix_names[[1]], matrix_names[[2]]) &&
    (is.null(parameter_names) || identical(matrix_names[[1]], parameter_names))

  if (!is_square) add_failure("hessian_not_square")
  if (!finite) add_failure("hessian_nonfinite")
  if (!conformable) add_failure("hessian_not_conformable")
  if (!named) add_failure("hessian_names_invalid")

  symmetry_relative <- NA_real_
  eigen_information <- numeric()
  numerical_rank <- NA_integer_
  scaled_condition <- NA_real_
  scale_factors <- numeric()
  H_sym <- NULL
  info <- NULL
  if (is_square && finite) {
    h_scale <- max(abs(H))
    symmetry_relative <- max(abs(H - t(H))) / max(h_scale, .Machine$double.eps)
    if (symmetry_relative > control$symmetry_tol) add_failure("hessian_asymmetric")
    H_sym <- (H + t(H)) / 2
    info <- -H_sym
    eigen_information <- tryCatch(
      eigen(info, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) rep(NA_real_, nrow(H))
    )
    if (any(!is.finite(eigen_information))) {
      add_failure("information_eigen_failure")
    } else if (any(eigen_information <= control$information_eigenvalue_min)) {
      add_failure("information_not_positive_definite")
    }
    info_diag <- diag(info)
    if (all(is.finite(info_diag)) && all(info_diag > 0)) {
      scale_factors <- 1 / sqrt(info_diag)
      scaled_info <- info * tcrossprod(scale_factors)
      numerical_rank <- qr(scaled_info, tol = control$rank_tol)$rank
      scaled_condition <- tryCatch(kappa(scaled_info, exact = TRUE), error = function(e) Inf)
      if (numerical_rank < nrow(H)) add_failure("information_rank_deficient")
      if (!is.finite(scaled_condition) || scaled_condition > control$condition_max) {
        add_failure("information_ill_conditioned")
      }
    } else {
      add_failure("information_nonpositive_diagonal")
    }
  }

  gradient_scaled_max <- if (is.null(gradient)) NA_real_ else gradient$scaled_max
  if (is.null(gradient)) {
    add_failure("fitted_gradient_not_checked")
  } else if (
      (!is.finite(gradient_scaled_max) || gradient_scaled_max > control$gradient_tol)) {
    add_failure("fitted_gradient_too_large")
  }

  agreement_relative <- NA_real_
  if (!is.null(reference_hessian) && is_square && finite &&
      is.matrix(reference_hessian) && identical(dim(reference_hessian), dim(H)) &&
      all(is.finite(reference_hessian)) && length(scale_factors) == nrow(H)) {
    delta_scaled <- (H_sym - (reference_hessian + t(reference_hessian)) / 2) *
      tcrossprod(scale_factors)
    base_scaled <- H_sym * tcrossprod(scale_factors)
    agreement_relative <- sqrt(sum(delta_scaled^2)) /
      max(sqrt(sum(base_scaled^2)), .Machine$double.eps)
    if (agreement_relative > control$agreement_tol) add_failure("hessian_method_disagreement")
  } else if (isTRUE(control$check_agreement) && identical(source, "analytical")) {
    add_failure("hessian_agreement_not_checked")
  }

  vc <- NULL
  covariance_diagonal <- numeric()
  if (length(failures) == 0L) {
    vc <- tryCatch(solve(info), error = function(e) NULL)
    if (is.null(vc) || any(!is.finite(vc))) {
      add_failure("covariance_nonfinite")
    } else {
      dimnames(vc) <- dimnames(H)
      covariance_diagonal <- diag(vc)
      if (any(!is.finite(covariance_diagonal)) || any(covariance_diagonal <= 0)) {
        add_failure("covariance_diagonal_nonpositive")
      }
    }
  }

  diagnostics <- list(
    status = if (length(failures)) "unavailable" else "available",
    failure_codes = failures,
    validation_profile = control$profile,
    validation_defaults_version = control$defaults_version,
    effective_thresholds = unclass(control[c(
      "information_eigenvalue_min", "symmetry_tol", "rank_tol", "condition_max", "gradient_tol",
      "agreement_tol", "gradient_step", "check_agreement"
    )]),
    expert_override = control$expert_override,
    hessian_source = source,
    hessian_fallback = fallback,
    dimension = if (is_square) nrow(H) else NA_integer_,
    finite = finite,
    conformable = conformable,
    named = named,
    relative_symmetry_error = symmetry_relative,
    min_information_eigenvalue = if (length(eigen_information)) min(eigen_information) else NA_real_,
    max_information_eigenvalue = if (length(eigen_information)) max(eigen_information) else NA_real_,
    numerical_rank = numerical_rank,
    scaled_condition_number = scaled_condition,
    fitted_gradient = gradient$gradient %||% NULL,
    fitted_gradient_steps = gradient$steps %||% NULL,
    scaled_gradient_max = gradient_scaled_max,
    hessian_agreement_relative = agreement_relative,
    covariance_diagonal = covariance_diagonal
  )
  if (length(failures)) stop(.gl_inference_unavailable(diagnostics))

  list(vcov = vc, se = sqrt(covariance_diagonal), hessian_diagnostics = diagnostics)
}
