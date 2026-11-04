# Explicit opt-in joint covariance. Fixed vcov and its advisory gates are unchanged.
.gl_joint_likelihood_evaluator <- function(object) {
  pairs <- build_copula_pair_cache(object$response, object$response_margin, object$response_subject)
  link <- get_copula_dist(object$copula_dist)$copula_link
  base_branch <- NULL
  function(q) {
    eta <- calc_eta(q, object$model_matrix, object$margin_dist, link, object$par_s)
    if (any(!is.finite(unlist(eta$eta_inv, use.names = FALSE)))) {
      .gl_joint_stop("joint_predictor_nonfinite", "Joint likelihood has nonfinite distribution parameters.")
    }
    likelihood <- calc_likelihood_minimal(
      eta$eta_inv, mm = object$model_matrix$x,
      margin_dist = object$margin_dist, copula_dist = object$copula_dist,
      calc_d2 = FALSE, response = object$response,
      response_margin = object$response_margin, response_subject = object$response_subject,
      pair_cache = pairs, calc_margin_deriv = FALSE)
    observed <- is.finite(object$response)
    active <- likelihood$pair_complete
    continuous <- !identical(likelihood$likelihood_type, "discrete_rectangle")
    floor <- if (continuous) .Machine$double.xmin * (1 + 1e-12) else 0
    if (identical(likelihood$valid, FALSE) ||
        any(!is.finite(likelihood$margin_d[observed]) | likelihood$margin_d[observed] <= floor) ||
        any(!is.finite(likelihood$copula_d[active]) | likelihood$copula_d[active] <= floor) ||
        !is.finite(unname(likelihood$log_lik["joint"]))) {
      .gl_joint_stop("joint_likelihood_invalid",
        "An included joint likelihood contribution is invalid, nonpositive or underflowed.")
    }
    if (identical(object$copula_dist, "t") && identical(.copula_backend(), "native")) {
      df <- likelihood$copula_par2[active]
      .copula_t_df_derivative_domain(df)  # Retains the exact-join refusal.
      branch <- df >= .copula_t_gaussian_limit_df
      if (is.null(base_branch)) base_branch <<- branch
      if (!identical(branch, base_branch)) {
        .gl_joint_stop("joint_df_branch_crossing", "A joint numerical stencil crosses the native Student-t approximation join.")
      }
    }
    list(loglik = unname(likelihood$log_lik["joint"]), eta = eta, likelihood = likelihood)
  }
}

.gl_vcov_joint <- function(object, par, method, progress, h, inference, details, cluster = NULL) {
  objective <- object$likelihood_contract$objective %||% object$missingness$objective %||% "ordinary"
  if (!identical(objective, "ordinary") || isTRUE(object$missingness$n_subjects_with_gaps > 0)) {
    .gl_joint_stop("joint_segmented_unsupported", "Joint penalized inference is unavailable for segmented or intermittent-gap likelihoods.")
  }
  if (identical(method, "sandwich") || !is.null(cluster)) {
    .gl_joint_stop("joint_sandwich_unsupported", "Joint cluster-sandwich inference is not implemented; fixed-only sandwich remains available.")
  }
  if (!(is.atomic(par) && length(par) == 1L && is.na(par))) {
    .gl_joint_stop("joint_par_override_unsupported", "Joint inference currently evaluates fitted coefficients only.")
  }
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0) {
    stop("'h' must be one positive finite number.", call. = FALSE)
  }
  control <- .gl_normalize_inference_control(inference)
  joint <- .gl_joint_design(object, control)
  evaluate <- .gl_joint_likelihood_evaluator(joint$object)
  base <- evaluate(joint$q)
  structured_evaluate <- .gl_joint_likelihood_evaluator(joint$structured)
  structured <- structured_evaluate(joint$structured$par)
  eta_equal <- isTRUE(all.equal(lapply(base$eta$eta, as.numeric),
                               lapply(structured$eta$eta, as.numeric), tolerance = 1e-10))
  if (!eta_equal || !isTRUE(all.equal(base$loglik, structured$loglik, tolerance = 1e-10))) {
    .gl_joint_stop("joint_representation_mismatch", "Augmented and structured fitted predictors/likelihood disagree.")
  }
  method_requested <- method
  preflight <- .gl_vcov_apply_margin_preflight(method, method, object$margin_dist, base$eta$eta_inv)
  method <- preflight$method
  .gl_vcov_apply_likelihood_preflight(method, method, base$likelihood, object$response)
  reference <- NULL
  if (identical(method, "numderiv")) {
    H <- .gl_joint_fd_hessian(joint$q, evaluate, h)
    used <- "numderiv"
  } else {
    .gl_source_analytical_hessian_helpers()
    # Do not replace a domain/underflow/curvature refusal with a numerical result.
    H <- tryCatch(calc_analytical_hessian(joint$object, progress = progress, h = h),
      error = function(e) .gl_joint_stop("joint_analytical_curvature_unavailable",
        paste("Joint analytical curvature is unavailable:", conditionMessage(e))))
    used <- "analytical"
    if (isTRUE(control$check_agreement)) reference <- .gl_joint_fd_hessian(joint$q, evaluate, h)
  }
  score <- .gl_joint_fd_gradient(joint$q, evaluate, control$gradient_step)
  solved <- .gl_joint_solve(H, joint$penalty, joint$q, score$gradient,
                           control, used, reference, score$steps)
  V <- solved$vcov
  smooth_vcov <- smooth_se <- list()
  for (pn in names(joint$smooth)) {
    if (!length(joint$smooth[[pn]])) next
    smooth_vcov[[pn]] <- smooth_se[[pn]] <- list()
    for (term in names(joint$smooth[[pn]])) {
      ids <- names(joint$smooth[[pn]][[term]])
      smooth_vcov[[pn]][[term]] <- V[ids, ids, drop = FALSE]
      smooth_se[[pn]][[term]] <- sqrt(diag(V[ids, ids, drop = FALSE]))
    }
  }
  smooth_names <- joint$map$name[joint$map$kind == "smooth"]
  contract <- .gl_inference_contract("joint_penalized_inverse", names(joint$q),
    method = used, validity_status = solved$hessian_diagnostics$status,
    failure_states = solved$hessian_diagnostics$failure_codes)
  contract$method_requested <- method_requested
  contract$method_used <- used
  contract$covariance_type <- "penalized_inverse"
  contract$diagnostics <- solved$hessian_diagnostics
  contract$fallback_used <- !identical(method_requested, method) # Includes margin preflight only.
  result <- list(
    vcov = list(overall = V, joint = V,
      fixed_joint = V[joint$fixed_names, joint$fixed_names, drop = FALSE],
      fixed_smooth = V[joint$fixed_names, smooth_names, drop = FALSE],
      smooth_vcov = smooth_vcov, smooth_se = smooth_se),
    se = list(overall = solved$se, smooth_se = smooth_se),
    coefficients = "joint", covariance_type = "penalized_inverse", joint_inference_version = 1L,
    coefficient_estimates = joint$q, coefficient_map = joint$map,
    conditioning = "fitted smoothing parameters and fitted constrained basis",
    method = used, method_requested = method_requested,
    hessian_diagnostics = solved$hessian_diagnostics,
    likelihood_hessian = solved$likelihood_hessian,
    penalty = solved$penalty, penalized_information = solved$penalized_information,
    smoothing_parameters = object$lambda_s,
    inference_contract = contract, smooth_inference_contract = contract,
    cache_policy = "joint covariance is recomputed; fixed and legacy caches are never reused")
  if (.gl_vcov_continuous_curvature_applies(object) && used == "analytical") {
    result$continuous_curvature <- .gl_vcov_continuous_curvature_marker(object)
  }
  result <- .gl_vcov_stamp_epoch(result)
  out <- .gl_vcov_format_result(result, return_details = details)
  if (!isTRUE(details)) {
    attr(out, "coefficients") <- "joint"
    attr(out, "covariance_type") <- "penalized_inverse"
    attr(out, "coefficient_estimates") <- joint$q
    attr(out, "coefficient_map") <- joint$map
  }
  out
}
