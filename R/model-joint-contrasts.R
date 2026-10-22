# Private next-slice candidate. No fitting, covariance recomputation or defaults changed.
.gl_joint_contrast_bundle <- function(joint) {
  fail <- function(message) .gl_joint_stop("joint_contrast_bundle_invalid", message)
  if (!is.list(joint) || !is.list(joint$inference_contract) || !is.list(joint$vcov) ||
      !identical(joint$coefficients, "joint") ||
      !identical(joint$joint_inference_version, 1L) ||
      !identical(joint$covariance_type, "penalized_inverse") ||
      !identical(joint$inference_contract$contract_id, "joint_penalized_inverse")) {
    fail("Contrasts require a detailed version-1 joint penalized-inverse bundle.")
  }
  q <- joint$coefficient_estimates
  map <- joint$coefficient_map
  V <- joint$vcov$overall
  if (!is.numeric(q) || !length(q) || any(!is.finite(q)) ||
      !.gl_joint_valid_names(names(q)) || !is.data.frame(map) ||
      !all(c("name", "parameter", "kind", "term", "basis_index") %in% names(map)) ||
      !identical(as.character(map$name), names(q)) || anyNA(map[c("name", "parameter", "kind", "term")]) ||
      any(!map$kind %in% c("fixed", "smooth"))) fail("Joint coefficient map/estimates are malformed.")
  if (!is.matrix(V) || !is.numeric(V) || any(!is.finite(V)) ||
      !identical(dim(V), c(length(q), length(q))) ||
      !identical(rownames(V), names(q)) || !identical(colnames(V), names(q))) {
    fail("Joint covariance must have the exact estimate coordinates.")
  }
  scale <- max(abs(V), .Machine$double.eps)
  if (max(abs(V - t(V))) > 1e-8 * scale ||
      is.null(tryCatch(chol((V + t(V)) / 2), error = function(e) NULL))) {
    fail("Joint contrast covariance must be symmetric positive definite.")
  }
  diagnostics <- joint$hessian_diagnostics
  if (!is.list(diagnostics) || length(diagnostics$status) != 1L ||
      !isTRUE(diagnostics$status %in% c("available", "available_with_warnings")) ||
      !identical(joint$inference_contract$validity_status, diagnostics$status) ||
      length(diagnostics$blocking_failure_codes)) fail("Joint inference is unavailable or has inconsistent validity metadata.")
  warnings <- unique(c(diagnostics$warning_codes, diagnostics$failure_codes,
                       joint$inference_contract$observed_failures))
  clean <- identical(diagnostics$status, "available") &&
    identical(diagnostics$stationarity_status, "passed") && !length(warnings)
  list(q = q, map = map, V = V, clean = clean, warnings = warnings,
       diagnostics = diagnostics, covariance_contract = joint$inference_contract)
}

.gl_joint_linear_contrast <- function(joint, L, level = 0.95,
    estimand = "named linear combination of joint coefficients") {
  bundle <- .gl_joint_contrast_bundle(joint)
  if (!is.matrix(L) || !is.numeric(L) || !nrow(L) || any(!is.finite(L)) ||
      !.gl_joint_valid_names(colnames(L)) || !.gl_joint_valid_names(rownames(L)) ||
      !setequal(colnames(L), names(bundle$q)) ||
      any(rowSums(abs(L)) == 0)) {
    .gl_joint_stop("joint_contrast_weights_invalid", "L must have unique named nonzero rows and exactly all joint coefficient columns; write explicit zeros.")
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1 ||
      !is.character(estimand) || length(estimand) != 1L || is.na(estimand) || !nzchar(estimand)) {
    .gl_joint_stop("joint_contrast_request_invalid", "A finite confidence level in (0,1) and one estimand description are required.")
  }
  L <- L[, names(bundle$q), drop = FALSE]
  estimate <- as.numeric(L %*% bundle$q)
  covariance <- L %*% bundle$V %*% t(L)
  variance <- diag(covariance)
  if (any(!is.finite(estimate)) || any(!is.finite(covariance)) || any(variance <= 0)) {
    .gl_joint_stop("joint_contrast_variance_invalid", "Contrast estimate/variance is nonfinite or nonpositive; no clipping or coefficient-SE substitution.")
  }
  se <- sqrt(variance)
  z <- stats::qnorm((1 + level) / 2)
  table <- data.frame(contrast = rownames(L), estimate = estimate, std_error = se,
    conf_low = estimate - z * se, conf_high = estimate + z * se,
    uncertainty_available = TRUE, clean_uncertainty_eligible = bundle$clean,
    eligibility_reason = if (bundle$clean) "joint_checks_passed" else "warning_qualified_joint_inference_requires_review",
    inference_status = bundle$diagnostics$status, stringsAsFactors = FALSE)
  contract <- list(contract_id = "joint_named_linear_contrast_v1", estimand = estimand,
    coefficient_scope = "fixed and constrained smooth coordinates",
    coefficient_names = names(bundle$q), covariance_type = "penalized_inverse",
    conditioning = joint$conditioning,
    omitted_uncertainty = bundle$covariance_contract$omitted_uncertainty,
    validity_status = bundle$diagnostics$status, observed_failures = bundle$warnings,
    covariance_contract = bundle$covariance_contract, level = level,
    interval = "pointwise normal approximation; not simultaneous or prediction interval",
    clean_uncertainty_eligible = bundle$clean,
    publication_acceptance = "requires the paper-specific validation/denominator policy; this flag is not scientific sign-off")
  list(table = table, L = L, covariance = covariance, inference_contract = contract)
}

.gl_joint_contrast_object <- function(object, joint) {
  .gl_require_converged_fit(object, "joint fixed-smooth paper contrasts")
  objective <- object$likelihood_contract$objective %||% object$missingness$objective %||% "ordinary"
  if (!identical(objective, "ordinary") || isTRUE(object$missingness$n_subjects_with_gaps > 0)) {
    .gl_joint_stop("joint_contrast_segmented_unsupported", "Segmented/intermittent joint model-Hessian contrasts remain unavailable.")
  }
  bundle <- .gl_joint_contrast_bundle(joint)
  design <- .gl_joint_design(object, inference_control())
  columns <- c("name", "parameter", "kind", "term", "basis_index")
  if (!isTRUE(all.equal(bundle$q, design$q, tolerance = 0)) ||
      !identical(bundle$map[columns], design$map[columns]) ||
      !isTRUE(all.equal(joint$penalty, design$penalty, tolerance = 0)) ||
      !isTRUE(all.equal(joint$smoothing_parameters, object$lambda_s, tolerance = 0))) {
    .gl_joint_stop("joint_contrast_object_mismatch", "Joint bundle coordinates, penalty or fitted lambda do not match the supplied object.")
  }
  # Callers must pass freshly computed details from this same object. q/P/map
  # matching is not a fingerprint of response data or proof of covariance origin.
  list(bundle = bundle, design = design)
}

.gl_joint_intercept_reference <- function(object, parameters,
    reference = c("bcpe_sample", "missingness_population")) {
  reference <- match.arg(reference)
  # A point estimand does not require available uncertainty. This mapping is
  # also valid for a segmented fit; only the later covariance request refuses it.
  design <- .gl_joint_design(object, inference_control())
  checked <- list(bundle = list(q = design$q, map = design$map))
  if (!.gl_joint_valid_names(parameters) || !length(parameters) ||
      any(!parameters %in% checked$bundle$map$parameter)) {
    .gl_joint_stop("joint_contrast_parameter_invalid", "Select unique existing parameter names.")
  }
  L <- matrix(0, length(parameters), length(checked$bundle$q),
    dimnames = list(paste0(parameters, ".", reference, "_intercept"), names(checked$bundle$q)))
  row_counts <- stats::setNames(integer(length(parameters)), parameters)
  terms_used <- stats::setNames(vector("list", length(parameters)), parameters)
  grid <- if (reference == "missingness_population") data.frame(s1 = seq(0, 1, length.out = 10001L)) else NULL
  for (parameter in parameters) {
    map <- checked$bundle$map
    intercept <- map$name[map$parameter == parameter & map$kind == "fixed" & map$term %in% c("intercept", "(Intercept)")]
    if (length(intercept) != 1L) .gl_joint_stop("joint_contrast_intercept_invalid", "Exactly one named fixed intercept is required per parameter.")
    row <- match(parameter, parameters)
    L[row, intercept] <- 1
    blocks <- object$model_matrix$s[[parameter]] %||% list()
    terms_used[[parameter]] <- names(blocks)
    row_counts[[parameter]] <- if (is.null(grid)) nrow(object$model_matrix$x[[parameter]]) else nrow(grid)
    if (!is.finite(row_counts[[parameter]]) || row_counts[[parameter]] < 1L) .gl_joint_stop("joint_contrast_rows_invalid", "The reference must contain at least one row.")
    for (term in names(blocks)) {
      if (is.null(grid)) {
        B <- blocks[[term]]
        colnames(B) <- paste(parameter, term, seq_len(ncol(B)), sep = ".")
      } else {
        spec <- .gl_joint_smooth_spec(blocks[[term]], paste(parameter, term, sep = ":"))
        if (!identical(as.character(spec$term), "s1") ||
            (!is.null(spec$by) && !identical(spec$by, "NA"))) {
          .gl_joint_stop("joint_population_reference_unsupported", "The registered missingness grid supports univariate s1 smooths without by variables only.")
        }
        B <- .gl_joint_smooth_prediction_matrix(object, parameter, term, grid)
      }
      expected <- map$name[map$parameter == parameter & map$kind == "smooth" & map$term == term]
      if (!setequal(colnames(B), expected) || nrow(B) != row_counts[[parameter]]) {
        .gl_joint_stop("joint_contrast_basis_mismatch", "Reference basis does not match all retained smooth coordinates/rows.")
      }
      L[row, colnames(B)] <- colMeans(B)
    }
  }
  list(L = L, coefficient_estimates = checked$bundle$q,
    parameters = parameters, terms = terms_used, row_counts = row_counts,
    reference = reference, grid = if (is.null(grid)) NULL else list(variable = "s1", lower = 0, upper = 1, points = 10001L, endpoints = TRUE, weights = "equal"),
    covariance_object_binding = "named q/map/penalty/lambda checked; caller must supply fresh details from the same response/design object")
}

.gl_joint_intercept_weights <- function(object, joint, parameters,
    reference = c("bcpe_sample", "missingness_population")) {
  .gl_joint_contrast_object(object, joint)
  .gl_joint_intercept_reference(object, parameters, match.arg(reference))
}

.gl_bcpe_sample_intercept_contrasts <- function(object, joint,
    parameters = c("mu", "sigma", "theta"), level = 0.95) {
  weights <- .gl_joint_intercept_weights(object, joint, parameters, "bcpe_sample")
  out <- .gl_joint_linear_contrast(joint, weights$L, level,
    "link-scale fixed intercept plus mean fitted smooth contribution over the fitted parameter-design rows")
  out$reference <- weights
  out$inference_contract$reference <- weights$reference
  out
}

.gl_missingness_population_intercept_contrasts <- function(object, joint,
    parameters = c("mu", "sigma", "theta"), level = 0.95) {
  weights <- .gl_joint_intercept_weights(object, joint, parameters, "missingness_population")
  out <- .gl_joint_linear_contrast(joint, weights$L, level,
    "link-scale fixed intercept plus fitted smooth mean on the fixed 10001-point Uniform(0,1) s1 reference grid")
  out$reference <- weights
  out$inference_contract$reference <- weights$grid
  out$inference_contract$point_estimate_epoch <- "retained PredictMat basis; historical training-value interpolation is not reused"
  out
}

.gl_joint_unvalidated_gamlss2_contrast <- function() {
  list(uncertainty_available = FALSE, clean_uncertainty_eligible = FALSE,
    std_error = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
    inference_status = "not_available_without_validated_comparator_joint_coordinates",
    reason = "Longitudinal joint covariance cannot certify gamlss2 fixed/smooth contrast covariance.")
}
