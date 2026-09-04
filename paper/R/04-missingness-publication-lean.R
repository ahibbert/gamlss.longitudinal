# Lean Phase 2 publication layer for the missingness study.
#
# The local checkpoint archive is producer infrastructure. Reviewers receive only
# portable CSV evidence and the two manuscript figures derived from those CSVs.

jss_missingness_lean_artifacts <- function() {
  c(
    "missingness_configuration.csv", "missingness_design_registry.csv",
    "missingness_estimand_registry.csv", "missingness_sensitivity_registry.csv",
    "fit_run_log.csv", "attempt_failure_summary.csv", "failure_reason_summary.csv",
    "missingness_warning_events.csv", "missingness_warning_audit.csv",
    "missingness_by_rep.csv", "missingness_pattern_by_subject_visit.csv",
    "fixed_effects_by_rep.csv", "smooth_estimates_by_rep.csv",
    "missingness_headline_summary.csv", "fixed_term_summary_by_missingness.csv",
    "smooth_irmse_summary.csv", "smooth_selected_plot_data.csv"
  )
}

jss_missingness_public_artifacts <- jss_missingness_lean_artifacts

jss_missingness_lean_bind <- function(x) {
  x <- x[vapply(x, is.data.frame, logical(1))]
  if (!length(x)) return(data.frame())
  fields <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(value) {
    for (field in setdiff(fields, names(value))) value[[field]] <- NA
    value[fields]
  })
  do.call(rbind, x)
}

jss_missingness_lean_order <- function(x) {
  keys <- intersect(c(
    "analysis_role", "missing_mechanism", "target_missing_rate", "scenario",
    "model", "parameter", "term", "rep", "id", "time", "s1", "warning_index"
  ), names(x))
  if (nrow(x) && length(keys)) {
    x <- x[do.call(order, c(lapply(x[keys], function(value) {
      if (is.numeric(value)) value else as.character(value)
    }), list(na.last = TRUE))), , drop = FALSE]
  }
  rownames(x) <- NULL
  x
}

jss_missingness_lean_metric <- function(values) {
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]
  n <- length(values)
  estimate <- if (n) mean(values) else NA_real_
  status <- if (!n) "no_success" else if (n == 1L) "one_success_no_mcse" else "estimated"
  mcse <- if (n > 1L) stats::sd(values) / sqrt(n) else NA_real_
  critical <- if (n > 1L) stats::qt(.975, n - 1L) else NA_real_
  list(n = as.integer(n), estimate = estimate, mcse = mcse,
    conf_low = if (n > 1L) estimate - critical * mcse else NA_real_,
    conf_high = if (n > 1L) estimate + critical * mcse else NA_real_,
    status = status)
}

jss_missingness_lean_wilson <- function(events, attempts) {
  if (!is.finite(events) || !is.finite(attempts) || attempts < 1L ||
      events < 0L || events > attempts) return(c(lower = NA_real_, upper = NA_real_))
  z <- stats::qnorm(.975); p <- events / attempts
  denominator <- 1 + z^2 / attempts
  centre <- (p + z^2 / (2 * attempts)) / denominator
  half <- z * sqrt(p * (1 - p) / attempts + z^2 / (4 * attempts^2)) / denominator
  c(lower = max(0, centre - half), upper = min(1, centre + half))
}

jss_missingness_lean_warning_classify <- function(message) {
  normalized <- trimws(gsub("[[:space:]]+", " ", enc2utf8(as.character(message))))
  inference <- c(
    "Analytical Hessian vcov failed; falling back to numerical Hessian. Reason: Hessian-based inference is unavailable (fitted_gradient_too_large). The fitted model is unchanged.",
    "Hessian-based inference is unavailable (fitted_gradient_too_large). The fitted model is unchanged.",
    "Hessian-based inference is unavailable (information_not_positive_definite, fitted_gradient_too_large). The fitted model is unchanged."
  )
  segmented <- grepl(
    paste0(
      "^Intermittent observation gaps or leading unobserved scheduled visits were detected for [0-9]+ subject\\(s\\)\\. ",
      "Using the segmented likelihood requested by missingness = \\\"segment\\\": observations in different contiguous segments are treated as independent\\. ",
      "AIC and BIC remain available under this explicit assumption\\. Ordinary model-Hessian inference is unavailable for this composite objective; ",
      "use cluster-sandwich or bootstrap sensitivity analyses\\. Numerical integration of dependence across gaps and a distinct delayed-entry likelihood ",
      "are not implemented; a future integration engine may add this support\\.$"
    ), normalized, perl = TRUE
  )
  if (normalized %in% inference) {
    return(c(classification = "inference_unavailable", expected = "TRUE",
      coverage_eligible = "FALSE"))
  }
  if (segmented) {
    return(c(classification = "segmented_likelihood_assumption", expected = "TRUE",
      coverage_eligible = "TRUE"))
  }
  c(classification = "unexpected", expected = "FALSE", coverage_eligible = "FALSE")
}

jss_missingness_lean_warning_events <- function(events) {
  required <- c("scenario", "rep", "missing_mechanism", "target_missing_rate",
    "model", "warning_index", "normalized_message")
  if (!is.data.frame(events) || !all(required %in% names(events))) {
    stop("Missingness warning events lack their lean publication schema.", call. = FALSE)
  }
  classified <- lapply(events$normalized_message, jss_missingness_lean_warning_classify)
  events$classification <- vapply(classified, `[[`, character(1), "classification")
  events$expected <- vapply(classified, function(x) identical(x[["expected"]], "TRUE"), logical(1))
  events$coverage_eligible <- vapply(classified,
    function(x) identical(x[["coverage_eligible"]], "TRUE"), logical(1))
  if (any(!events$expected)) {
    stop("Unexpected missingness warning text is fatal to publication.", call. = FALSE)
  }
  keep <- c(required, "classification", "expected", "coverage_eligible")
  jss_missingness_lean_order(events[keep])
}

jss_missingness_lean_subject_audit <- function(pattern) {
  by_subject <- split(pattern[order(pattern$id, pattern$time), , drop = FALSE], pattern$id)
  flags <- lapply(by_subject, function(x) {
    observed <- x$response_observed %in% TRUE
    observed_at <- which(observed)
    # A subject with no observed visit has neither an observed entry nor an
    # estimable within-subject segment; it is not labelled delayed entry.
    first <- if (length(observed_at)) min(observed_at) else NA_integer_
    last <- if (length(observed_at)) max(observed_at) else NA_integer_
    missing_at <- which(!observed)
    c(interior = length(observed_at) > 0L &&
        any(!observed & seq_along(observed) > first & seq_along(observed) < last),
      delayed = length(observed_at) > 0L && first > 1L,
      complete = all(observed),
      monotone = !length(missing_at) || all(!observed[min(missing_at):length(observed)]),
      observed = sum(observed), pairs = max(length(observed) - 1L, 0L),
      complete_pairs = if (length(observed) > 1L)
        sum(observed[-length(observed)] & observed[-1L]) else 0L)
  })
  flags <- do.call(rbind, flags)
  mechanism <- unique(as.character(pattern$missing_mechanism))
  c(
    observed_missing_rate = 1 - sum(flags[, "observed"]) / nrow(pattern),
    n_rows = nrow(pattern), n_observed_rows = sum(flags[, "observed"]),
    n_subjects = nrow(flags), n_complete_subjects = sum(flags[, "complete"]),
    n_dropout_subjects = if (identical(mechanism, "monotone_dropout"))
      sum(!flags[, "complete"] & flags[, "monotone"]) else 0L,
    n_monotone_incomplete_subjects = sum(!flags[, "complete"] & flags[, "monotone"]),
    n_subjects_with_interior_gaps = sum(flags[, "interior"]),
    n_subjects_with_delayed_entry = sum(flags[, "delayed"]),
    n_subjects_requiring_segmented_likelihood = sum(flags[, "interior"] | flags[, "delayed"]),
    no_observations_after_dropout = all(flags[, "monotone"]),
    total_adjacent_pairs = sum(flags[, "pairs"]),
    complete_adjacent_pairs = sum(flags[, "complete_pairs"]),
    complete_adjacent_pair_rate = sum(flags[, "complete_pairs"]) / sum(flags[, "pairs"])
  )
}

jss_missingness_lean_configuration <- function() {
  data.frame(
    field = c("design_version", "tasks", "fits", "subjects", "visits", "replicates",
      "missingness_rates", "missingness_mechanisms", "family", "copula", "compute_se",
      "max_outer_iter", "max_inner_iter", "max_elapsed_sec", "predictive_scores",
      "failure_penalty", "seed_rule"),
    value = c("missingness-v5", "240", "480", "500", "4", "20",
      "0,0.1,0.2,0.3,0.4,0.5",
      "monotone_dropout,time_dependent_intermittent_mar", "BCPE_log_mu", "t", "TRUE",
      "1000", "100", "180", "FALSE", "1",
      "simulation=100000+rep; missingness=simulation+mechanism_offset+1000*rate"),
    stringsAsFactors = FALSE
  )
}

jss_missingness_lean_validate_raw <- function(results) {
  if (!is.list(results) || length(results) != 240L) {
    stop("Lean missingness rebuild requires exactly 240 raw task payloads.", call. = FALSE)
  }
  if (any(!vapply(results, function(x) {
    is.list(x) && identical(x$checkpoint_schema_version, 5L) &&
      all(c("runs", "fixed", "smooth", "missingness", "missingness_pattern", "warning_events") %in% names(x))
  }, logical(1)))) stop("A raw missingness task lacks the registered v5 payload fields.", call. = FALSE)
  configurations <- lapply(results, function(x) x$checkpoint_spec$configuration)
  constant <- vapply(configurations[-1L], identical, logical(1), configurations[[1L]])
  if (!all(constant)) stop("Raw missingness tasks do not share one configuration.", call. = FALSE)
  config <- configurations[[1L]]
  checks <- list(
    design_version = identical(config$design_version, "missingness-v5"),
    family = identical(config$family, "BCPE_log_mu"), copula = identical(config$copula, "t"),
    smooth_k = identical(config$smooth_k, 10L), compute_se = isTRUE(config$compute_se),
    optimizer = identical(config$optim_method, "RS"),
    max_outer = identical(config$max_outer_iter, 1000L),
    max_inner = identical(config$max_inner_iter, 100L),
    max_elapsed = identical(config$max_elapsed_sec, 180),
    predictive = identical(config$predictive$compute, FALSE),
    seed_rule = identical(config$seed_formula_version,
      "simulation=100000+rep; mechanism offset + rate offset")
  )
  if (!all(unlist(checks))) stop("Raw missingness configuration is not the registered production configuration: ",
    paste(names(checks)[!unlist(checks)], collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

jss_missingness_lean_attempt_summary <- function(runs) {
  groups <- split(runs, interaction(runs$scenario, runs$model, drop = TRUE))
  out <- lapply(groups, function(x) {
    attempted <- nrow(x); fit_successful <- sum(x$success %in% TRUE)
    converged <- sum(x$converged %in% TRUE); retained <- sum(x$retained %in% TRUE)
    interval <- jss_missingness_lean_wilson(retained, attempted)
    data.frame(
      scenario = x$scenario[[1]], model = x$model[[1]],
      missing_mechanism = x$missing_mechanism[[1]], missingness_label = x$missingness_label[[1]],
      analysis_role = x$analysis_role[[1]], target_missing_rate = x$target_missing_rate[[1]],
      attempted = attempted, fit_successful = fit_successful, converged = converged,
      retained = retained, failed = attempted - retained,
      error_failures = sum(x$failure_type == "error", na.rm = TRUE),
      nonconvergence_failures = sum(!(x$retained %in% TRUE) & x$failure_type != "error", na.rm = TRUE),
      failure_inclusive_retention_rate = retained / attempted,
      retention_rate_mcse = sqrt((retained / attempted) * (1 - retained / attempted) / attempted),
      retention_rate_conf_low = interval[["lower"]], retention_rate_conf_high = interval[["upper"]],
      failure_rmse_penalty = 1, stringsAsFactors = FALSE)
  })
  jss_missingness_lean_order(do.call(rbind, out))
}

jss_missingness_lean_parent_retained <- function(table, runs) {
  key <- paste(runs$scenario, runs$rep, runs$model, sep = "\r")
  retained <- stats::setNames(runs$retained %in% TRUE, key)
  retained[paste(table$scenario, table$rep, table$model, sep = "\r")] %in% TRUE
}

jss_missingness_lean_fixed_summaries <- function(fixed, runs, attempts, warnings) {
  fixed$parent_retained <- jss_missingness_lean_parent_retained(fixed, runs)
  unavailable <- unique(paste(warnings$scenario[!warnings$coverage_eligible],
    warnings$rep[!warnings$coverage_eligible], warnings$model[!warnings$coverage_eligible], sep = "\r"))
  fixed$coverage_eligible <- fixed$parent_retained &
    !paste(fixed$scenario, fixed$rep, fixed$model, sep = "\r") %in% unavailable
  groups <- split(fixed, interaction(fixed$scenario, fixed$model, fixed$parameter, fixed$term, drop = TRUE))
  term <- lapply(groups, function(x) {
    ok <- x$parent_retained & is.finite(x$estimate) & is.finite(x$true_value)
    error <- x$estimate[ok] - x$true_value[ok]
    metric <- jss_missingness_lean_metric(error^2)
    n_metric <- sum(ok)
    n_attempted <- attempts$attempted[match(paste(x$scenario[[1]], x$model[[1]]),
      paste(attempts$scenario, attempts$model))]
    penalized <- c(error^2, rep(1, n_attempted - n_metric))
    inclusive <- jss_missingness_lean_metric(penalized)
    se_ok <- ok & x$coverage_eligible & is.finite(x$std_error) & x$std_error >= 0
    coverage <- if (any(se_ok)) abs(x$estimate[se_ok] - x$true_value[se_ok]) <= 1.96 * x$std_error[se_ok] else logical()
    coverage_metric <- jss_missingness_lean_metric(coverage)
    data.frame(
      scenario = x$scenario[[1]], model = x$model[[1]], parameter = x$parameter[[1]], term = x$term[[1]],
      true_value = x$true_value[[which.max(is.finite(x$true_value))]],
      mean_estimate = if (any(ok)) mean(x$estimate[ok]) else NA_real_,
      bias = if (length(error)) mean(error) else NA_real_, rmse = sqrt(as.numeric(metric[["estimate"]])),
      n_metric_successful = n_metric, metric_status = metric[["status"]],
      attempted_fit_count = n_attempted, metric_failures = n_attempted - n_metric,
      failure_inclusive_rmse = sqrt(as.numeric(inclusive[["estimate"]])),
      n_se = sum(se_ok), coverage_95 = as.numeric(coverage_metric[["estimate"]]),
      coverage_mcse = as.numeric(coverage_metric[["mcse"]]),
      coverage_conf_low = as.numeric(coverage_metric[["conf_low"]]),
      coverage_conf_high = as.numeric(coverage_metric[["conf_high"]]),
      coverage_status = coverage_metric[["status"]], stringsAsFactors = FALSE)
  })
  term <- jss_missingness_lean_order(do.call(rbind, term))
  groups <- split(term, interaction(term$scenario, term$model, term$parameter, drop = TRUE))
  summary <- lapply(groups, function(x) {
    a <- attempts[match(paste(x$scenario[[1]], x$model[[1]]), paste(attempts$scenario, attempts$model)), ]
    n_metric <- min(x$n_metric_successful)
    data.frame(
      scenario = x$scenario[[1]], model = x$model[[1]],
      missing_mechanism = a$missing_mechanism, target_missing_rate = a$target_missing_rate,
      parameter = x$parameter[[1]], mean_rmse = if (any(is.finite(x$rmse))) mean(x$rmse, na.rm = TRUE) else NA_real_,
      failure_inclusive_mean_rmse = if (any(is.finite(x$failure_inclusive_rmse)))
        mean(x$failure_inclusive_rmse, na.rm = TRUE) else NA_real_,
      mean_abs_bias = if (any(is.finite(x$bias))) mean(abs(x$bias), na.rm = TRUE) else NA_real_,
      mean_coverage_95 = if (any(is.finite(x$coverage_95))) mean(x$coverage_95, na.rm = TRUE) else NA_real_,
      metric_status = if (!n_metric) "no_success" else if (n_metric == 1L) "one_success_no_mcse" else "estimated",
      attempted_fit_count = a$attempted, fit_successful_count = a$fit_successful,
      converged_fit_count = a$converged, retained_fit_count = a$retained,
      fit_failure_count = a$failed, n_metric_successful = n_metric,
      metric_failures = a$attempted - n_metric,
      extraction_failure_count = a$retained - n_metric,
      target_missing_pct = 100 * a$target_missing_rate, stringsAsFactors = FALSE)
  })
  list(term = term, summary = jss_missingness_lean_order(do.call(rbind, summary)))
}

jss_missingness_lean_headline <- function(fixed, runs, attempts) {
  fixed$parent_retained <- jss_missingness_lean_parent_retained(fixed, runs)
  rep_groups <- split(fixed, interaction(fixed$scenario, fixed$model, fixed$parameter, fixed$rep, drop = TRUE))
  by_rep <- lapply(rep_groups, function(x) {
    ok <- x$parent_retained & is.finite(x$estimate) & is.finite(x$true_value)
    data.frame(scenario = x$scenario[[1]], model = x$model[[1]], parameter = x$parameter[[1]], rep = x$rep[[1]],
      rmse = if (any(ok)) sqrt(mean((x$estimate[ok] - x$true_value[ok])^2)) else NA_real_)
  })
  by_rep <- do.call(rbind, by_rep)
  groups <- split(by_rep, interaction(by_rep$scenario, by_rep$model, by_rep$parameter, drop = TRUE))
  out <- lapply(groups, function(x) {
    a <- attempts[match(paste(x$scenario[[1]], x$model[[1]]), paste(attempts$scenario, attempts$model)), ]
    conditional <- jss_missingness_lean_metric(x$rmse)
    n_metric <- sum(is.finite(x$rmse)); inclusive <- jss_missingness_lean_metric(c(x$rmse[is.finite(x$rmse)], rep(1, a$attempted - n_metric)))
    data.frame(
      scenario = x$scenario[[1]], model = x$model[[1]], parameter = x$parameter[[1]],
      missing_mechanism = a$missing_mechanism, missingness_label = a$missingness_label,
      analysis_role = a$analysis_role, target_missing_rate = a$target_missing_rate,
      conditional_mean_rmse = as.numeric(conditional[["estimate"]]),
      conditional_rmse_mcse = as.numeric(conditional[["mcse"]]),
      conditional_rmse_conf_low = as.numeric(conditional[["conf_low"]]),
      conditional_rmse_conf_high = as.numeric(conditional[["conf_high"]]),
      metric_status = conditional[["status"]], failure_penalty = 1,
      failure_inclusive_mean_rmse = as.numeric(inclusive[["estimate"]]),
      failure_inclusive_rmse_mcse = as.numeric(inclusive[["mcse"]]),
      failure_inclusive_rmse_conf_low = as.numeric(inclusive[["conf_low"]]),
      failure_inclusive_rmse_conf_high = as.numeric(inclusive[["conf_high"]]),
      attempted_fit_count = a$attempted, fit_successful_count = a$fit_successful,
      converged_fit_count = a$converged, retained_fit_count = a$retained,
      fit_failure_count = a$failed, n_metric_successful = n_metric,
      metric_failures = a$attempted - n_metric, extraction_failure_count = a$retained - n_metric,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  jss_missingness_lean_order(out[out$analysis_role == "headline", , drop = FALSE])
}

jss_missingness_lean_smooth <- function(smooth, runs, attempts) {
  smooth$parent_retained <- jss_missingness_lean_parent_retained(smooth, runs)
  rep_groups <- split(smooth, interaction(smooth$scenario, smooth$model, smooth$parameter, smooth$rep, drop = TRUE))
  by_rep <- lapply(rep_groups, function(x) {
    ok <- x$parent_retained & is.finite(x$smooth_hat) & is.finite(x$smooth_true)
    data.frame(scenario = x$scenario[[1]], model = x$model[[1]], parameter = x$parameter[[1]], rep = x$rep[[1]],
      irmse = if (any(ok)) sqrt(mean((x$smooth_hat[ok] - x$smooth_true[ok])^2)) else NA_real_)
  })
  by_rep <- do.call(rbind, by_rep)
  groups <- split(by_rep, interaction(by_rep$scenario, by_rep$model, by_rep$parameter, drop = TRUE))
  summary <- lapply(groups, function(x) {
    a <- attempts[match(paste(x$scenario[[1]], x$model[[1]]), paste(attempts$scenario, attempts$model)), ]
    metric <- jss_missingness_lean_metric(x$irmse); n_metric <- sum(is.finite(x$irmse))
    inclusive <- jss_missingness_lean_metric(c(x$irmse[is.finite(x$irmse)], rep(1, a$attempted - n_metric)))
    data.frame(
      scenario = x$scenario[[1]], model = x$model[[1]], parameter = x$parameter[[1]],
      missing_mechanism = a$missing_mechanism, target_missing_rate = a$target_missing_rate,
      mean_irmse = as.numeric(metric[["estimate"]]), irmse_mcse = as.numeric(metric[["mcse"]]),
      irmse_conf_low = as.numeric(metric[["conf_low"]]), irmse_conf_high = as.numeric(metric[["conf_high"]]),
      metric_status = metric[["status"]], n_metric_successful = n_metric,
      metric_failures = a$attempted - n_metric, failure_rmse_penalty = 1,
      failure_inclusive_mean_irmse = as.numeric(inclusive[["estimate"]]),
      failure_inclusive_irmse_mcse = as.numeric(inclusive[["mcse"]]),
      failure_inclusive_irmse_conf_low = as.numeric(inclusive[["conf_low"]]),
      failure_inclusive_irmse_conf_high = as.numeric(inclusive[["conf_high"]]),
      attempted_fit_count = a$attempted, fit_successful_count = a$fit_successful,
      converged_fit_count = a$converged, retained_fit_count = a$retained,
      fit_failure_count = a$failed, extraction_failure_count = a$retained - n_metric,
      stringsAsFactors = FALSE)
  })
  summary <- jss_missingness_lean_order(do.call(rbind, summary))
  selected <- smooth[smooth$target_missing_rate %in% c(0, .3, .5), , drop = FALSE]
  point_groups <- split(selected, interaction(selected$scenario, selected$model,
    selected$parameter, selected$s1, drop = TRUE))
  plot <- lapply(point_groups, function(x) {
    ok <- x$parent_retained & is.finite(x$smooth_hat)
    a <- attempts[match(paste(x$scenario[[1]], x$model[[1]]), paste(attempts$scenario, attempts$model)), ]
    data.frame(
      scenario = x$scenario[[1]], model = x$model[[1]], parameter = x$parameter[[1]], s1 = x$s1[[1]],
      smooth_true = x$smooth_true[[1]], smooth_median = if (any(ok)) stats::median(x$smooth_hat[ok]) else NA_real_,
      smooth_q05 = if (any(ok)) stats::quantile(x$smooth_hat[ok], .05, names = FALSE) else NA_real_,
      smooth_q95 = if (any(ok)) stats::quantile(x$smooth_hat[ok], .95, names = FALSE) else NA_real_,
      missing_mechanism = a$missing_mechanism, missingness_label = a$missingness_label,
      analysis_role = a$analysis_role, target_missing_rate = a$target_missing_rate,
      attempted_fit_count = a$attempted, retained_fit_count = a$retained,
      fit_failure_count = a$failed, missing_label = paste0(round(100 * a$target_missing_rate), "%"),
      stringsAsFactors = FALSE)
  })
  list(summary = summary, plot = jss_missingness_lean_order(do.call(rbind, plot)))
}

jss_missingness_materialize_lean_bundle <- function(raw_archive, output_dir) {
  results <- readRDS(raw_archive)
  jss_missingness_lean_validate_raw(results)
  runs <- jss_missingness_lean_order(jss_missingness_lean_bind(lapply(results, `[[`, "runs")))
  fixed <- jss_missingness_lean_order(jss_missingness_lean_bind(lapply(results, `[[`, "fixed")))
  smooth <- jss_missingness_lean_order(jss_missingness_lean_bind(lapply(results, `[[`, "smooth")))
  missingness <- jss_missingness_lean_order(jss_missingness_lean_bind(lapply(results, `[[`, "missingness")))
  pattern <- jss_missingness_lean_order(jss_missingness_lean_bind(lapply(results, `[[`, "missingness_pattern")))
  warning_events <- jss_missingness_lean_warning_events(
    jss_missingness_lean_bind(lapply(results, `[[`, "warning_events")))

  expected <- jss_missing_registered_task_grid()
  actual <- unique(runs[names(expected)])
  jss_missing_compare_frame(actual, expected, c("scenario", "rep"),
    "Lean missingness registered design", tolerance = 1e-15)
  if (nrow(runs) != 480L || anyDuplicated(paste(runs$scenario, runs$rep, runs$model)) ||
      !all(vapply(split(runs$model, interaction(runs$scenario, runs$rep, drop = TRUE)),
        function(x) identical(sort(as.character(x)), c("gamlss.longitudinal", "gamlss2")), logical(1)))) {
    stop("Lean missingness attempt rows do not form the exact 240-task/480-fit design.", call. = FALSE)
  }
  audits <- lapply(split(pattern, interaction(pattern$scenario, pattern$rep, drop = TRUE)),
    jss_missingness_lean_subject_audit)
  audit <- do.call(rbind, audits)
  audit_key <- names(audits)
  names(audit_key) <- audit_key
  key <- interaction(missingness$scenario, missingness$rep, drop = TRUE)
  audit <- audit[match(as.character(key), audit_key), , drop = FALSE]
  missingness$n_subjects_with_interior_gaps <- audit[, "n_subjects_with_interior_gaps"]
  missingness$n_subjects_with_delayed_entry <- audit[, "n_subjects_with_delayed_entry"]
  missingness$n_subjects_requiring_segmented_likelihood <- audit[, "n_subjects_requiring_segmented_likelihood"]
  missingness$analysis_missingness_mode <- ifelse(
    missingness$n_subjects_requiring_segmented_likelihood > 0L, "segment", "error")
  run_key <- paste(runs$scenario, runs$rep)
  missing_key <- paste(missingness$scenario, missingness$rep)
  for (field in c("n_subjects_with_interior_gaps", "n_subjects_with_delayed_entry",
      "n_subjects_requiring_segmented_likelihood", "analysis_missingness_mode")) {
    runs[[field]] <- missingness[[field]][match(run_key, missing_key)]
  }
  runs$n_subjects_with_interior_gaps_legacy <- NULL
  missingness$n_subjects_with_interior_gaps_legacy <- NULL
  runs$error <- NULL
  fixed$parent_retained <- NULL; smooth$parent_retained <- NULL

  attempts <- jss_missingness_lean_attempt_summary(runs)
  failures <- stats::aggregate(rep(1L, nrow(runs)),
    list(scenario = runs$scenario, model = runs$model,
      failure_reason = ifelse(runs$retained %in% TRUE, "none", runs$failure_type)), sum)
  names(failures)[[4L]] <- "attempts"
  warning_audit <- stats::aggregate(rep(1L, nrow(warning_events)),
    list(classification = warning_events$classification,
      expected = warning_events$expected, coverage_eligible = warning_events$coverage_eligible), sum)
  names(warning_audit)[[4L]] <- "events"
  affected_key <- paste(warning_events$scenario, warning_events$rep, warning_events$model)
  warning_audit$affected_fits <- vapply(warning_audit$classification, function(classification) {
    length(unique(affected_key[warning_events$classification == classification]))
  }, integer(1))
  warning_audit$attempted_fits <- nrow(runs)
  warning_audit$reconciliation_status <- "reconciled"

  fixed_outputs <- jss_missingness_lean_fixed_summaries(fixed, runs, attempts, warning_events)
  headline <- jss_missingness_lean_headline(fixed, runs, attempts)
  smooth_outputs <- jss_missingness_lean_smooth(smooth, runs, attempts)
  design <- unique(do.call(rbind, lapply(results, function(x)
    x$checkpoint_spec$configuration$mechanism_registry)))
  design <- design[design$missing_mechanism %in%
    c("monotone_dropout", "time_dependent_intermittent_mar"), , drop = FALSE]
  design <- design[match(c("monotone_dropout", "time_dependent_intermittent_mar"),
    design$missing_mechanism), , drop = FALSE]
  registries <- list(
    missingness_configuration.csv = jss_missingness_lean_configuration(),
    missingness_design_registry.csv = design,
    missingness_estimand_registry.csv = results[[1L]]$checkpoint_spec$configuration$estimand_registry,
    missingness_sensitivity_registry.csv = data.frame(
      penalty = 1, role = "primary", scale = "link-scale RMSE/IRMSE",
      application = "assigned to an attempted fit without a retained finite metric",
      justification = "Prespecified failure-inclusive sensitivity penalty; not an estimand.")
  )
  tables <- c(registries, list(
    fit_run_log.csv = runs, attempt_failure_summary.csv = attempts,
    failure_reason_summary.csv = jss_missingness_lean_order(failures),
    missingness_warning_events.csv = warning_events,
    missingness_warning_audit.csv = jss_missingness_lean_order(warning_audit),
    missingness_by_rep.csv = missingness,
    missingness_pattern_by_subject_visit.csv = pattern,
    fixed_effects_by_rep.csv = fixed, smooth_estimates_by_rep.csv = smooth,
    missingness_headline_summary.csv = headline,
    fixed_term_summary_by_missingness.csv = fixed_outputs$summary,
    smooth_irmse_summary.csv = smooth_outputs$summary,
    smooth_selected_plot_data.csv = smooth_outputs$plot
  ))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(tables)) {
    utils::write.csv(jss_missingness_lean_order(tables[[name]]),
      file.path(output_dir, name), row.names = FALSE, na = "")
  }
  jss_missingness_validate_public_bundle(output_dir)
}

jss_missingness_lean_no_private_location <- function(x) {
  text <- unlist(lapply(x[vapply(x, is.character, logical(1))], as.character), use.names = FALSE)
  if (!length(text)) return(TRUE)
  !any(grepl("(?:[A-Za-z]:[/\\\\]|[/](?:Users|home|mnt|tmp)[/]|\\\\\\\\|worker[_ -]?pid|hostname)",
    text, ignore.case = TRUE, perl = TRUE))
}

jss_missingness_lean_validate_metric_state <- function(x, n_field, estimate_field,
    interval_fields = character(), label = "Missingness metric") {
  required <- c(n_field, "metric_status", estimate_field, interval_fields)
  if (!all(required %in% names(x))) {
    stop(label, " lacks explicit metric-state fields.", call. = FALSE)
  }
  n <- suppressWarnings(as.integer(x[[n_field]]))
  expected <- ifelse(n == 0L, "no_success",
    ifelse(n == 1L, "one_success_no_mcse", "estimated"))
  if (anyNA(n) || any(n < 0L) || any(x$metric_status != expected)) {
    stop(label, " has inconsistent zero/one/estimated status.", call. = FALSE)
  }
  estimate <- suppressWarnings(as.numeric(x[[estimate_field]]))
  if (any(n == 0L & !is.na(estimate)) || any(n > 0L & !is.finite(estimate))) {
    stop(label, " has an estimate inconsistent with its successful-metric count.", call. = FALSE)
  }
  if (length(interval_fields)) {
    intervals <- lapply(x[interval_fields], function(value) suppressWarnings(as.numeric(value)))
    unavailable <- n < 2L
    if (any(vapply(intervals, function(value) any(unavailable & !is.na(value)), logical(1))) ||
        any(vapply(intervals, function(value) any(!unavailable & !is.finite(value)), logical(1)))) {
      stop(label, " has Monte Carlo uncertainty inconsistent with its metric state.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

jss_missingness_validate_public_bundle <- function(input, root = getwd(), ...) {
  paths <- file.path(input, jss_missingness_lean_artifacts())
  if (any(!file.exists(paths))) stop("Lean missingness public bundle is incomplete: ",
    paste(basename(paths[!file.exists(paths)]), collapse = ", "), call. = FALSE)
  read <- function(name) utils::read.csv(file.path(input, name), stringsAsFactors = FALSE,
    check.names = FALSE)
  tables <- stats::setNames(lapply(jss_missingness_lean_artifacts(), read),
    jss_missingness_lean_artifacts())
  if (!all(vapply(tables, jss_missingness_lean_no_private_location, logical(1)))) {
    stop("Lean missingness public bundle contains a host, path, or PID leak.", call. = FALSE)
  }
  runs <- tables[["fit_run_log.csv"]]; attempts <- tables[["attempt_failure_summary.csv"]]
  jss_missing_compare_frame(tables[["missingness_configuration.csv"]],
    jss_missingness_lean_configuration(), "field", "Lean missingness configuration")
  expected_registry <- jss_missing_mechanism_registry()
  expected_registry <- expected_registry[expected_registry$missing_mechanism %in%
    c("monotone_dropout", "time_dependent_intermittent_mar"), , drop = FALSE]
  expected_registry <- expected_registry[match(c("monotone_dropout",
    "time_dependent_intermittent_mar"), expected_registry$missing_mechanism), , drop = FALSE]
  jss_missing_compare_frame(tables[["missingness_design_registry.csv"]], expected_registry,
    "missing_mechanism", "Lean missingness mechanism registry")
  jss_missing_compare_frame(tables[["missingness_estimand_registry.csv"]],
    jss_missing_estimand_registry(), "component", "Lean missingness estimand registry")
  expected_sensitivity <- data.frame(
    penalty = 1, role = "primary", scale = "link-scale RMSE/IRMSE",
    application = "assigned to an attempted fit without a retained finite metric",
    justification = "Prespecified failure-inclusive sensitivity penalty; not an estimand.")
  jss_missing_compare_frame(tables[["missingness_sensitivity_registry.csv"]],
    expected_sensitivity, "role", "Lean missingness sensitivity registry")
  expected <- jss_missing_registered_task_grid()
  jss_missing_compare_frame(unique(runs[names(expected)]), expected,
    c("scenario", "rep"), "Lean missingness registered design", tolerance = 1e-15)
  if (nrow(runs) != 480L || nrow(attempts) != 24L ||
      nrow(tables[["missingness_by_rep.csv"]]) != 240L) {
    stop("Lean missingness public bundle does not contain 480 fits, 240 tasks, and 24 attempt cells.", call. = FALSE)
  }
  if (any(runs$retained %in% TRUE &
      (!(runs$success %in% TRUE) | !(runs$converged %in% TRUE))) ||
      any(runs$failure_type == "none" & !(runs$retained %in% TRUE))) {
    stop("Lean missingness fit-state semantics are inconsistent.", call. = FALSE)
  }
  expected_mechanisms <- c("monotone_dropout", "time_dependent_intermittent_mar")
  if (!setequal(runs$missing_mechanism, expected_mechanisms) ||
      any(runs$analysis_role[runs$missing_mechanism == "monotone_dropout"] != "headline") ||
      any(runs$analysis_role[runs$missing_mechanism == "time_dependent_intermittent_mar"] != "sensitivity")) {
    stop("Lean missingness headline/sensitivity labels are not registered.", call. = FALSE)
  }
  warning_events <- tables[["missingness_warning_events.csv"]]
  reclassified <- jss_missingness_lean_warning_events(warning_events)
  jss_missing_compare_frame(warning_events, reclassified,
    c("scenario", "rep", "model", "warning_index"), "Lean missingness warning classification")
  expected_warning_audit <- stats::aggregate(rep(1L, nrow(warning_events)),
    list(classification = warning_events$classification, expected = warning_events$expected,
      coverage_eligible = warning_events$coverage_eligible), sum)
  names(expected_warning_audit)[[4L]] <- "events"
  warning_fit_key <- paste(warning_events$scenario, warning_events$rep, warning_events$model)
  expected_warning_audit$affected_fits <- vapply(expected_warning_audit$classification,
    function(classification) length(unique(warning_fit_key[
      warning_events$classification == classification])), integer(1))
  expected_warning_audit$attempted_fits <- nrow(runs)
  expected_warning_audit$reconciliation_status <- "reconciled"
  expected_warning_audit <- jss_missingness_lean_order(expected_warning_audit)
  jss_missing_compare_frame(tables[["missingness_warning_audit.csv"]],
    expected_warning_audit, "classification", "Lean missingness warning audit")

  expected_attempts <- jss_missingness_lean_attempt_summary(runs)
  jss_missing_compare_frame(attempts, expected_attempts, c("scenario", "model"),
    "Lean missingness attempt summary")
  expected_failures <- stats::aggregate(rep(1L, nrow(runs)),
    list(scenario = runs$scenario, model = runs$model,
      failure_reason = ifelse(runs$retained %in% TRUE, "none", runs$failure_type)), sum)
  names(expected_failures)[[4L]] <- "attempts"
  expected_failures <- jss_missingness_lean_order(expected_failures)
  jss_missing_compare_frame(tables[["failure_reason_summary.csv"]], expected_failures,
    c("scenario", "model", "failure_reason"), "Lean missingness failure summary")

  pattern <- tables[["missingness_pattern_by_subject_visit.csv"]]
  missingness <- tables[["missingness_by_rep.csv"]]
  pattern_audits <- lapply(split(pattern,
    interaction(pattern$scenario, pattern$rep, drop = TRUE)), jss_missingness_lean_subject_audit)
  pattern_audit <- do.call(rbind, pattern_audits)
  pattern_key <- names(pattern_audits); names(pattern_key) <- pattern_key
  observed_key <- interaction(missingness$scenario, missingness$rep, drop = TRUE)
  pattern_audit <- pattern_audit[match(as.character(observed_key), pattern_key), , drop = FALSE]
  audit_fields <- colnames(pattern_audit)
  corrected <- data.frame(scenario = missingness$scenario, rep = missingness$rep,
    as.data.frame(pattern_audit, stringsAsFactors = FALSE, check.names = FALSE),
    analysis_missingness_mode = ifelse(
      pattern_audit[, "n_subjects_requiring_segmented_likelihood"] > 0L, "segment", "error"))
  reported_corrected <- missingness[c("scenario", "rep", audit_fields,
    "analysis_missingness_mode")]
  jss_missing_compare_frame(reported_corrected, corrected, c("scenario", "rep"),
    "Lean missingness interior-gap and delayed-entry audit")

  fixed <- tables[["fixed_effects_by_rep.csv"]]
  smooth <- tables[["smooth_estimates_by_rep.csv"]]
  reconstructed_fixed <- jss_missingness_lean_fixed_summaries(
    fixed, runs, attempts, warning_events)$summary
  reconstructed_headline <- jss_missingness_lean_headline(fixed, runs, attempts)
  reconstructed_smooth <- jss_missingness_lean_smooth(smooth, runs, attempts)
  jss_missing_compare_frame(tables[["fixed_term_summary_by_missingness.csv"]],
    reconstructed_fixed, c("scenario", "model", "parameter"),
    "Lean missingness fixed-effect summary")
  jss_missing_compare_frame(tables[["missingness_headline_summary.csv"]],
    reconstructed_headline, c("scenario", "model", "parameter"),
    "Lean missingness headline summary")
  jss_missing_compare_frame(tables[["smooth_irmse_summary.csv"]],
    reconstructed_smooth$summary, c("scenario", "model", "parameter"),
    "Lean missingness smooth summary")
  jss_missing_compare_frame(tables[["smooth_selected_plot_data.csv"]],
    reconstructed_smooth$plot, c("scenario", "model", "parameter", "s1"),
    "Lean missingness smooth plot data")
  status_tables <- tables[c("missingness_headline_summary.csv",
    "fixed_term_summary_by_missingness.csv", "smooth_irmse_summary.csv")]
  allowed_status <- c("no_success", "one_success_no_mcse", "estimated")
  if (any(!vapply(status_tables, function(x)
      "metric_status" %in% names(x) && all(x$metric_status %in% allowed_status), logical(1)))) {
    stop("Lean missingness summaries lack explicit zero/one/estimated metric states.", call. = FALSE)
  }
  headline <- tables[["missingness_headline_summary.csv"]]
  jss_missingness_lean_validate_metric_state(headline, "n_metric_successful",
    "conditional_mean_rmse", c("conditional_rmse_mcse", "conditional_rmse_conf_low",
      "conditional_rmse_conf_high"), "Lean missingness headline")
  jss_missingness_lean_validate_metric_state(
    tables[["fixed_term_summary_by_missingness.csv"]], "n_metric_successful",
    "mean_rmse", label = "Lean missingness fixed-effect summary")
  jss_missingness_lean_validate_metric_state(
    tables[["smooth_irmse_summary.csv"]], "n_metric_successful", "mean_irmse",
    c("irmse_mcse", "irmse_conf_low", "irmse_conf_high"),
    "Lean missingness smooth summary")
  if (any(headline$analysis_role != "headline") ||
      any(headline$missing_mechanism != "monotone_dropout")) {
    stop("Lean missingness headline table is not restricted to monotone dropout.", call. = FALSE)
  }
  invisible(list(
    registry = tables[["missingness_design_registry.csv"]], runs = runs,
    attempts = attempts, missingness = tables[["missingness_by_rep.csv"]],
    headline = tables[["missingness_headline_summary.csv"]], tables = tables
  ))
}

jss_missingness_validate_candidate_bundle <- jss_missingness_validate_public_bundle
jss_missingness_production_gate <- function(input, root = getwd())
  jss_missingness_validate_public_bundle(input, root = root)

jss_run_04_missingness_from_public_inputs <- function(settings) {
  input <- file.path(settings$public_data_dir, "missingness")
  validated <- jss_missingness_validate_public_bundle(input, settings$root)
  fixed <- validated$tables[["fixed_term_summary_by_missingness.csv"]]
  smooth <- validated$tables[["smooth_selected_plot_data.csv"]]
  labels <- c(
    monotone_dropout = "Monotone dropout (headline)",
    time_dependent_intermittent_mar = "Intermittent MAR (sensitivity)"
  )
  fixed$mechanism_panel <- unname(labels[fixed$missing_mechanism])
  smooth$mechanism_panel <- unname(labels[smooth$missing_mechanism])
  fixed$model <- factor(fixed$model, levels = c("gamlss2", "gamlss.longitudinal"))
  fixed_plot <- fixed[is.finite(fixed$mean_rmse) &
    fixed$parameter %in% c("mu", "sigma", "nu", "tau"), , drop = FALSE]
  p1 <- ggplot2::ggplot(fixed_plot,
    ggplot2::aes(target_missing_pct, mean_rmse, colour = model, shape = model, group = model)) +
    ggplot2::geom_line(linewidth = .7) + ggplot2::geom_point(size = 1.9) +
    ggplot2::facet_grid(parameter ~ mechanism_panel, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = seq(0, 50, 10)) +
    ggplot2::labs(x = "Target missingness (%)", y = "Conditional fixed-effect RMSE",
      colour = NULL, shape = NULL) +
    ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "top")
  dir.create(settings$figures_dir, recursive = TRUE, showWarnings = FALSE)
  f1 <- file.path(settings$figures_dir, "fixed_margin_rmse_by_missingness.png")
  ggplot2::ggsave(f1, p1, width = 8.5, height = 5.5, dpi = 180, bg = "white")

  smooth$missing_label <- factor(smooth$missing_label, levels = c("0%", "30%", "50%"))
  smooth$model <- factor(smooth$model, levels = c("gamlss2", "gamlss.longitudinal"))
  smooth$mechanism_panel <- unname(c(
    monotone_dropout = "Dropout headline",
    time_dependent_intermittent_mar = "MAR sensitivity"
  )[smooth$missing_mechanism])
  smooth_plot <- smooth[is.finite(smooth$smooth_true) & is.finite(smooth$smooth_median), , drop = FALSE]
  p2 <- ggplot2::ggplot(smooth_plot, ggplot2::aes(s1)) +
    ggplot2::geom_line(ggplot2::aes(y = smooth_true), colour = "black", linewidth = .75) +
    ggplot2::geom_line(ggplot2::aes(y = smooth_median, colour = model),
      linewidth = .75, linetype = "dashed") +
    ggplot2::facet_grid(parameter + mechanism_panel ~ missing_label, scales = "free_y") +
    ggplot2::labs(x = "s1", y = "Centered smooth", colour = NULL) +
    ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "top")
  f2 <- file.path(settings$figures_dir, "smooth_selected_recovery_curves.png")
  ggplot2::ggsave(f2, p2, width = 8.5, height = 7, dpi = 180, bg = "white")
  staged <- jss_stage_04_missingness_data(input, settings)
  list(module_id = "04-missingness-dropout-sensitivity", status = "regenerated",
    data = staged, tables = character(), figures = c(f1, f2),
    notes = "Monotone dropout is the headline; time-dependent intermittent MAR is a sensitivity analysis.")
}

jss_run_04_missingness_dropout <- function(settings) {
  if (identical(settings$profile, "paper")) {
    return(jss_run_04_missingness_from_public_inputs(settings))
  }
  out <- file.path(settings$out_dir, "missingness")
  smoke <- identical(settings$profile, "smoke")
  cores <- if (settings$workers > 1L) settings$workers else
    max(1L, parallel::detectCores() - 2L)
  jss_run_script(
    file.path(settings$root, "paper", "scripts", "final-simulations", "missingness",
      "run_missingness_study.R"),
    c(
      OUT_DIR = out, N_FITS = if (smoke) "1" else "20",
      N_CORES = as.character(cores), RESUME_CHECKPOINTS = "1",
      MISSING_MECHANISMS = "monotone_dropout,time_dependent_intermittent_mar",
      MISSING_LEVELS = if (smoke) "0,0.2" else "0,0.1,0.2,0.3,0.4,0.5",
      COMPUTE_SE = if (smoke) "0" else "1", COMPUTE_PREDICTIVE_SCORES = "0",
      MAX_OUTER_ITER = if (smoke) "3" else "1000",
      MAX_INNER_ITER = if (smoke) "3" else "100",
      MAX_ELAPSED_SEC = if (smoke) "20" else "180"
    ), settings$root)
  if (smoke) {
    figures <- file.path(out,
      c("fixed_margin_rmse_by_missingness.png", "smooth_selected_recovery_curves.png"))
    return(list(module_id = "04-missingness-dropout-sensitivity", status = "smoke_only",
      data = character(), tables = character(), figures = figures))
  }
  jss_missingness_materialize_lean_bundle(
    file.path(out, "missingness_checkpoint_payloads.rds"), out)
  paper_settings <- settings
  paper_settings$public_data_dir <- settings$out_dir
  jss_run_04_missingness_from_public_inputs(paper_settings)
}
