#' Inference scope and validity contracts
#'
#' Returns the versioned, machine-readable contracts used by covariance,
#' interval, test, bootstrap, plotting, and publication-table consumers. These
#' contracts describe the uncertainty that is actually computed; they are not
#' claims of unconditional or all-parameter inference.
#'
#' @return A data frame with one row per inference producer or consumer.
#' @export
inference_contracts <- function() {
  .gl_inference_contract_registry()
}

.gl_inference_contract_version <- function() "2026.1"

.gl_inference_contract_registry <- function() {
  row <- function(id, producer, method, estimand, coefficient_blocks,
                  conditioning, omitted_uncertainty, assumptions,
                  failure_states, validation_status, approximation) {
    data.frame(
      contract_version = .gl_inference_contract_version(),
      contract_id = id,
      producer = producer,
      method = method,
      estimand = estimand,
      coefficient_blocks = coefficient_blocks,
      conditioning = conditioning,
      omitted_uncertainty = omitted_uncertainty,
      assumptions = assumptions,
      failure_states = failure_states,
      validation_status = validation_status,
      approximation = approximation,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, list(
    row(
      "fixed_hessian_analytical", "vcov/summary/tidy/confint",
      "observed fixed-coefficient Hessian (semi-analytical)",
      "joint-likelihood fixed coefficients", "all fixed margin and copula coefficients",
      "conditional on fitted smooth coefficients, penalties, and smoothing parameters",
      "fixed-smooth cross-covariance; smooth-coefficient covariance; smoothing-parameter uncertainty",
      "local quadratic likelihood approximation and JSS-002 curvature validity",
      "invalid curvature, fitted score, rank, conditioning, derivative agreement, or nonpositive covariance",
      "JSS-002 validated per fit", "large-sample conditional Wald approximation"
    ),
    row(
      "fixed_hessian_numerical", "vcov/summary/tidy/confint",
      "observed fixed-coefficient Hessian (finite difference)",
      "joint-likelihood fixed coefficients", "all fixed margin and copula coefficients",
      "conditional on fitted smooth coefficients, penalties, and smoothing parameters",
      "fixed-smooth cross-covariance; smooth-coefficient covariance; smoothing-parameter uncertainty",
      "stable finite differences, local quadratic likelihood approximation, and JSS-002 curvature validity",
      "invalid curvature, fitted score, rank, conditioning, or nonpositive covariance",
      "JSS-002 validated per fit", "large-sample conditional numerical-Wald approximation"
    ),
    row(
      "fixed_sandwich_cluster", "vcov(method='sandwich')",
      "subject/cluster score sandwich with validated Hessian bread",
      "joint-likelihood fixed coefficients", "all fixed margin and copula coefficients",
      "conditional on fitted smooth coefficients, penalties, and smoothing parameters; robust over declared independent clusters",
      "fixed-smooth cross-covariance; smooth and smoothing-parameter uncertainty; robustness to dependence across clusters",
      "independent clusters, enough clusters, stable cluster scores, and valid JSS-002 bread",
      "cross-cluster copula pairs, invalid bread, nonfinite scores/covariance, or nonpositive covariance diagonal",
      "JSS-002 bread validation plus runtime cluster checks", "cluster-robust large-sample conditional approximation"
    ),
    row(
      "smooth_penalized_conditional", "vcov smooth_vcov/plot_smooth_terms",
      "penalized working-weight covariance",
      "fitted smooth contribution", "one smooth coefficient block at a time",
      "conditional on all fixed coefficients, other smooths, penalties, and selected smoothing parameters",
      "fixed-smooth and between-smooth covariance; smoothing-parameter uncertainty; full joint-likelihood curvature",
      "working-weight and penalty approximation is adequate",
      "singular/nonfinite penalized precision or nonpositive derived variance",
      "approximation only; not JSS-002 Hessian inference", "conditional pointwise smooth-band approximation"
    ),
    row(
      "bootstrap_parametric_fixed", "bootstrap_inference(type='parametric')",
      "parametric fitted-model refit bootstrap",
      "selected fixed coefficients under the fitted data-generating model", "selected fixed margin/copula coefficients",
      "conditions on the fitted model specification, design, missingness mask, and fitted smooth structure used by each refit",
      "model-selection uncertainty; uncertainty outside selected fixed coefficients",
      "fitted model is an adequate data-generating model and successful refits are representative",
      "refit error, explicit nonconvergence, missing/nonfinite target coefficient, or insufficient successful term replicates",
      "replicate failures and per-term successful counts reported; two successful replicates are only a computational minimum, not evidence of inferential adequacy", "percentile bootstrap"
    ),
    row(
      "bootstrap_cluster_fixed", "bootstrap_inference(type='cluster')",
      "subject cluster-resampling refit bootstrap",
      "selected fixed coefficients under empirical subject sampling", "selected fixed margin/copula coefficients",
      "conditions on model specification and within-subject records; resamples whole subjects",
      "model-selection uncertainty; uncertainty outside selected fixed coefficients",
      "subjects are independent sampling units and successful refits are representative",
      "refit error, explicit nonconvergence, missing/nonfinite target coefficient, or insufficient successful term replicates",
      "replicate failures and per-term successful counts reported; two successful replicates are only a computational minimum, not evidence of inferential adequacy", "cluster percentile bootstrap"
    ),
    row(
      "wald_fixed", "wald_test", "Wald z or chi-square test",
      "selected fixed-coefficient contrasts", "fixed coefficients touched by the contrast matrix",
      "inherits the selected fixed covariance contract",
      "inherits fixed-covariance omissions, including fixed-smooth and smoothing-parameter uncertainty",
      "contrast is identified and the selected covariance contract is valid",
      "unavailable covariance, unmatched coefficients, nonpositive contrast variance, or singular joint contrast covariance",
      "inherits JSS-002/sandwich validity", "large-sample conditional Wald test"
    ),
    row(
      "likelihood_ratio_nested", "likelihood_compare", "sequential likelihood-ratio chi-square reference test",
      "difference between nested joint-likelihood model specifications", "no coefficient covariance block",
      "conditions on the compared specifications, a common observed dataset, and the chosen row order",
      "model-selection uncertainty; boundary corrections; uncertainty from non-nested comparison or changed analysis samples",
      "models are correctly ordered and nested, use the same observations, and differ by identifiable interior parameters",
      "nonpositive degrees-of-freedom increment, negative/nonfinite likelihood ratio, changed observed sample or objective, unverified nesting, or boundary null",
      "conservative runtime sample/objective/design checks; boundary conditions still require user judgment", "asymptotic chi-square reference approximation"
    ),
    row(
      "confint_fixed", "confint", "Wald normal interval",
      "selected fixed coefficients", "selected fixed margin/copula coefficients",
      "inherits the selected fixed covariance contract",
      "inherits fixed-covariance omissions, including fixed-smooth and smoothing-parameter uncertainty",
      "selected covariance contract is valid and asymptotic normality is adequate",
      "unavailable covariance, unknown coefficient, or nonpositive variance",
      "inherits JSS-002/sandwich validity", "large-sample conditional Wald interval"
    ),
    row(
      "prediction_mu_delta", "predict(se.fit/interval)", "first-order delta method",
      "row-level fitted mu or exploratory response mean", "mu fixed coefficients only",
      "conditional on all non-mu coefficients, every smooth contribution, penalties, and smoothing parameters",
      "non-mu coefficient uncertainty; all smooth uncertainty; fixed-smooth covariance; future response variation; trajectory dependence",
      "mu-only linearization is adequate; response-mean use outside identity-mean families is exploratory",
      "unavailable fixed covariance, unmatched mu design columns, or nonpositive derived variance",
      "inherits fixed-covariance validity", "conditional confidence interval, not a prediction interval for a future observation"
    ),
    row(
      "fitted_distribution_plugin", "predict(quantile/cdf/density/probability)", "plug-in fitted marginal distribution",
      "row-level fitted marginal distribution summary", "no covariance block; all parameters fixed at point estimates",
      "conditional on fitted margin parameters and the supplied row; not conditioned on other observed responses",
      "parameter, smooth, model-selection, and smoothing-parameter uncertainty",
      "fitted marginal family is adequate; quantile pairs describe model-implied response variation",
      "family evaluation failure or nonfinite fitted parameters",
      "point-estimate distribution only", "plug-in distributional prediction, not a confidence interval"
    ),
    row(
      "marginal_effect_mu_delta", "marginal_effects(se.fit=TRUE)", "aggregated mu-only delta method",
      "average fitted mu counterfactual", "mu fixed coefficients only",
      "conditional on the supplied covariate distribution and all omitted model blocks",
      "cross-row covariance; non-mu coefficients; smooth uncertainty; fixed-smooth covariance; smoothing-parameter uncertainty",
      "rowwise standard-error aggregation is adequate for exploratory summaries",
      "unavailable prediction covariance or no finite rowwise standard errors",
      "inherits fixed-covariance validity", "exploratory conditional approximation"
    ),
    row(
      "fixed_term_pointwise", "plot_fixed_terms", "design-column Wald band",
      "displayed fixed design-column contribution(s)", "the fixed coefficient(s) represented by each displayed design column or grouped factor panel",
      "conditional on all other fixed and smooth coefficients",
      "covariance with other fixed terms; all smooth and smoothing-parameter uncertainty",
      "single-column term decomposition is the intended estimand",
      "unavailable fixed covariance, unmatched coefficient, or nonpositive variance",
      "inherits fixed-covariance validity", "conditional pointwise term band"
    ),
    row(
      "smooth_term_pointwise", "plot_smooth_terms", "penalized smooth covariance band",
      "one fitted smooth contribution", "one smooth coefficient block",
      "conditional on fixed coefficients, other smooths, penalties, and smoothing parameters",
      "fixed-smooth and between-smooth covariance; smoothing-parameter uncertainty",
      "penalized working-weight covariance is adequate",
      "missing/singular smooth covariance or nonpositive derived variance",
      "approximation only", "conditional pointwise smooth band"
    ),
    row(
      "publication_coefficients", "publication_table(table='coefficients')",
      "formatted inherited fixed-coefficient Wald inference",
      "fixed coefficient estimates and conditional Wald uncertainty", "all printed fixed coefficients",
      "inherits the selected fixed covariance contract",
      "inherits fixed covariance omissions; table formatting adds no uncertainty",
      "underlying summary covariance is valid",
      "any underlying inference failure is propagated",
      "inherits selected covariance validity", "formatted conditional inference"
    ),
    row(
      "publication_predictions_point", "reporting_table/publication prediction table",
      "grouped point prediction summaries", "grouped fitted means, parameters, quantiles, and probabilities",
      "no covariance block; point estimates only",
      "conditional on fitted model and supplied rows",
      "parameter uncertainty, smooth uncertainty, sampling uncertainty, and future-response variation",
      "point summaries answer the reporting estimand",
      "prediction evaluation failure or no finite group values",
      "point estimates only", "no inferential interval"
    )
  ))
}

.gl_inference_contract <- function(contract_id, coefficient_names = NULL,
                                   method = NULL, validity_status = NULL,
                                   failure_states = NULL) {
  registry <- .gl_inference_contract_registry()
  idx <- match(contract_id, registry$contract_id)
  if (is.na(idx)) stop("Unknown inference contract: ", contract_id, call. = FALSE)
  out <- as.list(registry[idx, , drop = FALSE])
  out <- lapply(out, `[[`, 1L)
  if (!is.null(coefficient_names)) out$coefficient_names <- as.character(coefficient_names)
  if (!is.null(method)) out$method_used <- as.character(method)[1L]
  if (!is.null(validity_status)) out$validity_status <- as.character(validity_status)[1L]
  if (!is.null(failure_states)) out$observed_failures <- as.character(failure_states)
  class(out) <- c("gamlss_longitudinal_inference_contract", "list")
  out
}

.gl_attach_inference_contract <- function(x, contract) {
  attr(x, "inference_contract") <- contract
  x
}

.gl_vcov_contract_id <- function(method) {
  method <- as.character(method %||% "analytical")[[1L]]
  if (grepl("sandwich", method, fixed = TRUE)) return("fixed_sandwich_cluster")
  if (grepl("numderiv|numerical", method)) return("fixed_hessian_numerical")
  "fixed_hessian_analytical"
}

.gl_fixed_inference_contract <- function(vcov_out, coefficient_names = NULL) {
  method <- vcov_out$method %||% vcov_out$method_requested %||% "analytical"
  diagnostics <- vcov_out$hessian_diagnostics %||% list()
  out <- .gl_inference_contract(
    .gl_vcov_contract_id(method),
    coefficient_names = coefficient_names,
    method = method,
    validity_status = diagnostics$status %||% "not_recorded",
    failure_states = diagnostics$failure_codes %||% character()
  )
  out$method_requested <- vcov_out$method_requested %||% method
  out$method_used <- vcov_out$method %||% method
  out$fallback_used <- !identical(
    as.character(out$method_requested)[1L], as.character(out$method_used)[1L]
  )
  out$diagnostics <- diagnostics
  out
}

.gl_smooth_coefficient_names <- function(object) {
  if (is.null(object$par_s) || length(object$par_s) == 0L) return(character())
  unlist(lapply(names(object$par_s), function(parameter) {
    blocks <- object$par_s[[parameter]]
    if (length(blocks) == 0L) return(character())
    unlist(lapply(names(blocks), function(term) {
      paste0(parameter, ":", term, "[", seq_along(blocks[[term]]), "]")
    }), use.names = FALSE)
  }), use.names = FALSE)
}

.gl_enrich_vcov_contract <- function(vcov_out, object) {
  if (is.null(vcov_out$inference_contract)) {
    vcov_out$inference_contract <- .gl_fixed_inference_contract(
      vcov_out, coefficient_names = names(object$par)
    )
  }

  expected_smooth <- list()
  for (parameter in names(object$par_s %||% list())) {
    for (term in names(object$par_s[[parameter]] %||% list())) {
      n_coef <- length(object$par_s[[parameter]][[term]])
      covariance <- vcov_out$vcov$smooth_vcov[[parameter]][[term]] %||% NULL
      standard_error <- vcov_out$vcov$smooth_se[[parameter]][[term]] %||% NULL
      covariance_ok <- is.matrix(covariance) &&
        identical(dim(covariance), c(n_coef, n_coef)) && all(is.finite(covariance))
      se_ok <- length(standard_error) == n_coef && all(is.finite(standard_error))
      expected_smooth[[length(expected_smooth) + 1L]] <- data.frame(
        parameter = parameter,
        term = term,
        status = if (covariance_ok) "complete" else if (se_ok) "diagonal_only" else "unavailable",
        reason = if (covariance_ok) "full_block_available" else if (se_ok) {
          "only_coefficient_standard_errors_available"
        } else {
          "smooth_covariance_block_missing_or_invalid"
        },
        stringsAsFactors = FALSE
      )
    }
  }
  block_status <- if (length(expected_smooth)) do.call(rbind, expected_smooth) else {
    data.frame(
      parameter = character(), term = character(), status = character(),
      reason = character(), stringsAsFactors = FALSE
    )
  }
  smooth_status <- if (!nrow(block_status)) {
    "not_applicable"
  } else if (all(block_status$status == "complete")) {
    "complete_approximate"
  } else if (all(block_status$status == "unavailable")) {
    "unavailable"
  } else {
    "partial_approximate"
  }
  if (is.null(vcov_out$smooth_inference_contract)) {
    vcov_out$smooth_inference_contract <- .gl_inference_contract(
      "smooth_penalized_conditional",
      coefficient_names = .gl_smooth_coefficient_names(object),
      validity_status = smooth_status
    )
  }
  vcov_out$smooth_inference_contract$validity_status <- smooth_status
  vcov_out$smooth_inference_contract$coefficient_names <-
    .gl_smooth_coefficient_names(object)
  vcov_out$smooth_inference_contract$observed_failures <- if (smooth_status == "unavailable") {
    "all_smooth_covariance_blocks_unavailable"
  } else if (smooth_status == "partial_approximate") {
    "one_or_more_smooth_covariance_blocks_partial_or_unavailable"
  } else character()
  vcov_out$smooth_inference_contract$block_status <- block_status
  if (is.matrix(vcov_out$vcov$overall)) {
    attr(vcov_out$vcov$overall, "inference_contract") <- vcov_out$inference_contract
  }
  for (parameter in names(vcov_out$vcov$smooth_vcov %||% list())) {
    for (term in names(vcov_out$vcov$smooth_vcov[[parameter]] %||% list())) {
      block_names <- paste0(
        parameter, ":", term, "[",
        seq_along(object$par_s[[parameter]][[term]] %||% numeric()), "]"
      )
      block_row <- block_status[
        block_status$parameter == parameter & block_status$term == term, , drop = FALSE
      ]
      block_contract <- .gl_inference_contract(
        "smooth_penalized_conditional", coefficient_names = block_names,
        validity_status = if (nrow(block_row)) block_row$status[[1L]] else "unavailable",
        failure_states = if (nrow(block_row) && block_row$status[[1L]] == "unavailable") {
          block_row$reason[[1L]]
        } else character()
      )
      attr(vcov_out$vcov$smooth_vcov[[parameter]][[term]], "inference_contract") <- block_contract
      if (!is.null(vcov_out$vcov$smooth_se[[parameter]][[term]])) {
        attr(vcov_out$vcov$smooth_se[[parameter]][[term]], "inference_contract") <- block_contract
      }
    }
  }
  if (!is.null(vcov_out$vcov$smooth_vcov)) {
    attr(vcov_out$vcov$smooth_vcov, "inference_contract") <- vcov_out$smooth_inference_contract
  }
  if (!is.null(vcov_out$vcov$smooth_se)) {
    attr(vcov_out$vcov$smooth_se, "inference_contract") <- vcov_out$smooth_inference_contract
  }
  vcov_out
}

.gl_inference_contract_note <- function(contract) {
  if (is.null(contract)) return(NULL)
  paste0(
    contract$estimand, "; ", contract$conditioning,
    ". Omitted: ", contract$omitted_uncertainty, "."
  )
}
