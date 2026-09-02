jss_run_04_missingness_dropout <- function(settings) {
  if (identical(settings$profile, "paper")) return(jss_run_04_missingness_from_public_inputs(settings))
  out <- file.path(settings$out_dir, "missingness")
  reps <- if (identical(settings$profile, "full")) "20" else if (identical(settings$profile, "paper")) "20" else "1"
  levels <- if (identical(settings$profile, "smoke")) "0,0.2" else "0,0.1,0.2,0.3,0.4,0.5"
  cores <- if (settings$workers > 1L) settings$workers else max(1L, parallel::detectCores() - 2L)
  jss_run_script(
    file.path(settings$root, "paper", "scripts", "final-simulations", "missingness", "run_missingness_study.R"),
    c(
      OUT_DIR = out, N_FITS = reps, N_CORES = as.character(cores),
      MISSING_MECHANISMS = "monotone_dropout,time_dependent_intermittent_mar",
      MISSING_LEVELS = levels, COMPUTE_SE = if (identical(settings$profile, "smoke")) "0" else "1",
      COMPUTE_PREDICTIVE_SCORES = "0", MAX_OUTER_ITER = if (identical(settings$profile, "smoke")) "3" else "1000",
      MAX_INNER_ITER = if (identical(settings$profile, "smoke")) "3" else "100",
      MAX_ELAPSED_SEC = if (identical(settings$profile, "smoke")) "20" else "180"
    ),
    settings$root
  )
  if (identical(settings$profile, "full")) {
    jss_missingness_production_gate(out, settings$root)
  }
  manuscript_figures <- file.path(
    settings$figures_dir,
    c("fixed_margin_rmse_by_missingness.png", "smooth_selected_recovery_curves.png")
  )
  source_figures <- file.path(out, basename(manuscript_figures))
  copied <- file.copy(source_figures, manuscript_figures, overwrite = TRUE)
  if (!all(copied)) {
    stop(
      "Missingness study completed but manuscript figures could not be staged: ",
      paste(source_figures[!copied], collapse = ", "),
      call. = FALSE
    )
  }
  staged_data <- jss_stage_04_missingness_data(out, settings)
  list(
    module_id = "04-missingness-dropout-sensitivity", status = "regenerated",
    data = staged_data, tables = character(),
    figures = manuscript_figures
  )
}

jss_missingness_public_artifacts <- function() {
  c(
    "missingness_design_registry.csv", "missingness_estimand_registry.csv",
    "missingness_sensitivity_registry.csv", "missingness_checkpoint_status.csv",
    "missingness_checkpoint_payloads.rds", "missingness_checkpoint_content_manifest.csv",
    "fit_run_log.csv", "attempt_failure_summary.csv", "failure_reason_summary.csv",
    "missingness_warning_events.csv", "missingness_warning_audit.csv",
    "missingness_by_rep.csv", "missingness_pattern_by_subject_visit.csv", "fixed_effects_by_rep.csv",
    "smooth_estimates_by_rep.csv", "missingness_headline_summary.csv",
    "fixed_term_summary_by_missingness.csv", "smooth_irmse_summary.csv",
    "smooth_selected_plot_data.csv"
  )
}

jss_missing_registered_task_grid <- function() {
  grid <- expand.grid(
    missing_mechanism = c("monotone_dropout", "time_dependent_intermittent_mar"),
    target_missing_rate = seq(0, .5, by = .1), rep = seq_len(20L),
    stringsAsFactors = FALSE)
  grid$n <- 500L; grid$d <- 4L
  grid$scenario <- sprintf("n500_d4_%s_miss%02d", grid$missing_mechanism,
    round(100 * grid$target_missing_rate))
  grid$simulation_seed <- 100000L + as.integer(grid$rep)
  grid$missingness_seed <- 100000L + as.integer(grid$rep) +
    vapply(grid$missing_mechanism, jss_missing_seed_offset, integer(1)) +
    as.integer(round(1000 * grid$target_missing_rate))
  grid[c("scenario", "n", "d", "rep", "missing_mechanism", "target_missing_rate",
    "simulation_seed", "missingness_seed")]
}

jss_missing_validate_registered_design <- function(runs, checkpoints) {
  expected <- jss_missing_registered_task_grid()
  task_fields <- names(expected)
  if (!is.data.frame(runs) || nrow(runs) != 480L ||
      !all(c(task_fields, "model") %in% names(runs)) ||
      !is.data.frame(checkpoints) || nrow(checkpoints) != 240L ||
      !all(c("scenario", "rep") %in% names(checkpoints))) {
    stop("Missingness bundle is not the registered 240-task/480-model design.", call. = FALSE)
  }
  actual_tasks <- unique(runs[task_fields])
  jss_missing_compare_frame(actual_tasks, expected, c("scenario", "rep"),
    "Missingness registered Cartesian task grid", tolerance = 1e-15)
  model_cells <- split(runs$model, interaction(runs$scenario, runs$rep, drop = TRUE))
  if (length(model_cells) != 240L || !all(vapply(model_cells, function(x)
      identical(sort(as.character(x)), c("gamlss.longitudinal", "gamlss2")), logical(1)))) {
    stop("Missingness task grid lacks exactly both registered model rows.", call. = FALSE)
  }
  checkpoint_tasks <- unique(checkpoints[c("scenario", "rep")])
  expected_checkpoint_tasks <- expected[c("scenario", "rep")]
  jss_missing_compare_frame(checkpoint_tasks, expected_checkpoint_tasks,
    c("scenario", "rep"), "Missingness registered checkpoint grid", tolerance = 0)
  invisible(TRUE)
}

jss_missing_compare_frame <- function(reported, expected, keys, label, tolerance = 1e-10) {
  if (!setequal(names(reported), names(expected))) {
    stop(label, " schema does not match its canonical reconstruction.", call. = FALSE)
  }
  order_frame <- function(x) {
    if (!nrow(x)) return(x[, names(expected), drop = FALSE])
    ord <- do.call(order, c(unname(lapply(x[keys], function(value) unname(as.character(value)))),
      list(na.last = TRUE)))
    x[ord, names(expected), drop = FALSE]
  }
  actual <- order_frame(reported); canonical <- order_frame(expected)
  rownames(actual) <- NULL; rownames(canonical) <- NULL
  if (nrow(actual) != nrow(canonical)) stop(label, " row count is not canonical.", call. = FALSE)
  for (nm in names(canonical)) {
    if (is.numeric(canonical[[nm]])) {
      a <- as.numeric(actual[[nm]]); e <- as.numeric(canonical[[nm]])
      a[is.nan(a)] <- NA_real_; e[is.nan(e)] <- NA_real_
      ok <- isTRUE(all.equal(a, e,
        tolerance = tolerance, check.attributes = FALSE))
    } else {
      ok <- identical(as.character(actual[[nm]]), as.character(canonical[[nm]]))
    }
    if (!ok) stop(label, " disagrees with raw rows in column ", nm, ".", call. = FALSE)
  }
  invisible(TRUE)
}

jss_missing_mc_interval <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]; n <- length(x)
  estimate <- if (n) mean(x) else NA_real_
  mcse <- if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_
  critical <- if (n > 1L) stats::qt(.975, n - 1L) else NA_real_
  c(n = n, estimate = estimate, mcse = mcse,
    conf_low = estimate - critical * mcse, conf_high = estimate + critical * mcse)
}

jss_missing_validate_intervals <- function(x, stems, label) {
  for (stem in stems) {
    fields <- paste0(stem, c("", "_mcse", "_conf_low", "_conf_high"))
    if (!all(fields %in% names(x))) stop(label, " lacks interval field ", stem, ".", call. = FALSE)
    est <- x[[fields[[1L]]]]; se <- x[[fields[[2L]]]]
    lo <- x[[fields[[3L]]]]; hi <- x[[fields[[4L]]]]
    complete <- is.finite(est) | is.finite(se) | is.finite(lo) | is.finite(hi)
    if (any(complete & (!is.finite(est) | !is.finite(se) | se < 0 |
        !is.finite(lo) | !is.finite(hi) | lo > est | est > hi))) {
      stop(label, " contains an impossible Monte Carlo interval for ", stem, ".", call. = FALSE)
    }
  }
  invisible(TRUE)
}

jss_missing_parse_rfc3339_utc <- function(x, label = "timestamp") {
  if (!is.character(x) || anyNA(x) ||
      any(!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$", x))) {
    stop(label, " is not an RFC3339 UTC instant.", call. = FALSE)
  }
  parsed <- as.POSIXct(strptime(x, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  if (anyNA(parsed) || !identical(format(parsed, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), x)) {
    stop(label, " is not a real RFC3339 UTC instant.", call. = FALSE)
  }
  parsed
}

jss_missing_require_external_attestation <- function(bundle_sha256, package_sha256,
    producer_sha256, checkpoints, root, attestation_path, signature_path) {
  if (!requireNamespace("sodium", quietly = TRUE)) stop("sodium is required for production promotion signatures.", call. = FALSE)
  if (!nzchar(attestation_path) || !nzchar(signature_path)) {
    stop("Missingness public bundle lacks a detached production promotion signature.", call. = FALSE)
  }
  root_prefix <- paste0(normalizePath(root, winslash = "/", mustWork = TRUE), "/")
  path <- normalizePath(attestation_path, winslash = "/", mustWork = TRUE)
  signature <- normalizePath(signature_path, winslash = "/", mustWork = TRUE)
  if (startsWith(tolower(paste0(path, "/")), tolower(root_prefix)) ||
      startsWith(tolower(paste0(signature, "/")), tolower(root_prefix))) {
    stop("Missingness promotion attestation/signature must be external to the checkout.", call. = FALSE)
  }
  public_key <- as.raw(c(66, 91, 71, 233, 21, 247, 172, 45, 215, 202, 170, 0, 64,
    43, 83, 206, 23, 50, 48, 154, 25, 217, 178, 37, 252, 59, 158, 195, 237, 0, 31, 216))
  message_raw <- readBin(path, "raw", n = file.info(path)$size)
  signature_raw <- readBin(signature, "raw", n = file.info(signature)$size)
  if (!isTRUE(tryCatch(sodium::sig_verify(message_raw, signature_raw, public_key),
      error = function(e) FALSE))) {
    stop("Missingness promotion attestation lacks a valid detached production signature.", call. = FALSE)
  }
  x <- tryCatch(unserialize(message_raw), error = function(e) NULL)
  expected_names <- c("schema_version", "study", "bundle_sha256", "package_source_sha256",
    "producer_sha256", "approved_at_utc", "approver", "checkpoint_manifest")
  if (!is.list(x) || !identical(names(x), expected_names) || !identical(x$schema_version, 1L) ||
      !identical(x$study, "missingness") || !identical(x$bundle_sha256, bundle_sha256) ||
      !identical(x$package_source_sha256, package_sha256) ||
      !identical(x$producer_sha256, producer_sha256) || !is.character(x$approver) ||
      length(x$approver) != 1L || !nzchar(x$approver) ||
      !identical(x$checkpoint_manifest, checkpoints)) {
    stop("Missingness external promotion attestation does not bind the canonical checkpoint manifest.", call. = FALSE)
  }
  approved <- jss_missing_parse_rfc3339_utc(x$approved_at_utc, "Missingness approval timestamp")
  executed <- jss_missing_parse_rfc3339_utc(checkpoints$timestamp_utc, "Missingness checkpoint timestamp")
  if (approved <= max(executed)) stop("Missingness approval must occur after every checkpoint execution.", call. = FALSE)
  invisible(TRUE)
}

jss_missingness_validate_candidate_bundle <- function(input, root = getwd(),
    attestation_path = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MISSINGNESS_ATTESTATION", unset = ""),
    signature_path = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MISSINGNESS_ATTESTATION_SIGNATURE", unset = ""),
    require_promotion = TRUE) {
  required_paths <- file.path(input, jss_missingness_public_artifacts())
  if (any(!file.exists(required_paths))) {
    stop("Missingness public bundle is incomplete: ",
      paste(basename(required_paths[!file.exists(required_paths)]), collapse = ", "), call. = FALSE)
  }
  bundle_manifest <- paste(sort(paste(basename(required_paths),
    vapply(required_paths, jss_missing_sha256_file, character(1)), sep = "\t")), collapse = "\n")
  bundle_sha256 <- digest::digest(bundle_manifest, algo = "sha256", serialize = FALSE)
  read <- function(name) utils::read.csv(file.path(input, name), stringsAsFactors = FALSE)
  registry <- read("missingness_design_registry.csv")
  estimands <- read("missingness_estimand_registry.csv")
  penalties <- read("missingness_sensitivity_registry.csv")
  checkpoints <- read("missingness_checkpoint_status.csv")
  checkpoint_payloads <- readRDS(file.path(input, "missingness_checkpoint_payloads.rds"))
  checkpoint_manifest <- read("missingness_checkpoint_content_manifest.csv")
  runs <- read("fit_run_log.csv")
  attempts <- read("attempt_failure_summary.csv")
  failures <- read("failure_reason_summary.csv")
  warning_events <- read("missingness_warning_events.csv")
  warning_audit <- read("missingness_warning_audit.csv")
  if (!nrow(warning_events)) {
    warning_events <- jss_missing_empty_warning_events()
    warning_events$public_payload_sha256 <- character()
  }
  missing_by_rep <- read("missingness_by_rep.csv")
  missingness_pattern <- read("missingness_pattern_by_subject_visit.csv")
  fixed_raw <- read("fixed_effects_by_rep.csv")
  smooth_raw <- read("smooth_estimates_by_rep.csv")
  headline <- read("missingness_headline_summary.csv")
  fixed_reported <- read("fixed_term_summary_by_missingness.csv")
  smooth_irmse_reported <- read("smooth_irmse_summary.csv")
  smooth_reported <- read("smooth_selected_plot_data.csv")
  exact_schema <- function(x, expected, label) {
    if (!setequal(names(x), expected)) stop(label, " has extra or missing columns.", call. = FALSE)
  }
  run_schema <- c("scenario", "n", "d", "rep", "simulation_seed", "missingness_seed",
    "missing_mechanism", "missingness_label", "missingness_pattern", "analysis_role",
    "target_missing_rate", "observed_missing_rate", "complete_adjacent_pair_rate",
    "complete_adjacent_pairs", "n_dropout_subjects", "n_subjects_with_interior_gaps",
    "no_observations_after_dropout", "model", "success", "retained", "logLik", "df",
    "converged", "hit_outer_limit", "hit_max_stall", "hit_raw_loglik_deterioration",
    "stop_reason", "grad_inf", "step_l2", "best_raw_loglik", "best_raw_loglik_iteration",
    "raw_loglik_drop_from_best", "raw_loglik_drop_tol", "outer_iterations",
    "outer_log_lik_change", "outer_stop_crit", "elapsed_sec", "error", "failure_type",
    "public_payload_sha256")
  fixed_schema <- c("model", "parameter", "term", "estimate", "std_error", "true_value",
    "intercept_includes_fitted_smooth_mean", "inference_status", "scenario", "n", "d", "rep",
    "missing_mechanism", "target_missing_rate", "observed_missing_rate",
    "complete_adjacent_pair_rate", "public_payload_sha256")
  smooth_schema <- c("scenario", "n", "d", "rep", "missing_mechanism", "target_missing_rate",
    "observed_missing_rate", "complete_adjacent_pair_rate", "model", "parameter", "s1",
    "smooth_hat", "smooth_true", "public_payload_sha256")
  missing_schema <- c("missing_mechanism", "missingness_label", "missingness_pattern", "analysis_role",
    "target_missing_rate", "observed_missing_rate", "n_rows", "n_observed_rows", "n_subjects",
    "n_complete_subjects", "n_dropout_subjects", "n_monotone_incomplete_subjects",
    "n_subjects_with_interior_gaps", "no_observations_after_dropout", "total_adjacent_pairs",
    "complete_adjacent_pairs", "complete_adjacent_pair_rate", "scenario", "n", "d", "rep",
    "public_payload_sha256")
  pattern_schema <- c("scenario", "n", "d", "rep", "missing_mechanism", "id", "time",
    "response_observed", "public_payload_sha256")
  warning_schema <- c(names(jss_missing_empty_warning_events()), "public_payload_sha256")
  warning_audit_schema <- c("audit_type", "classification", "policy_action", "expected",
    "events", "affected_fits", "attempted_fits", "omitted_rows", "reconciliation_status",
    "uncertainty")
  checkpoint_schema <- c("checkpoint_schema_version", "scenario", "rep", "checkpoint_content_sha256",
    "checkpoint", "public_payload_sha256", "package_source_sha256", "package_version", "package_fingerprint_scope",
    "package_source_file_count", "producer_sha256", "package_identity_verified", "timestamp_utc",
    "worker_pid", "host", "os", "platform", "r_version", "rng_kind", "blas", "lapack",
    "rlibs_user_sha256", "libpaths_sha256")
  exact_schema(runs, run_schema, "Missingness fit log")
  exact_schema(fixed_raw, fixed_schema, "Missingness fixed raw payload")
  exact_schema(smooth_raw, smooth_schema, "Missingness smooth raw payload")
  exact_schema(missing_by_rep, missing_schema, "Missingness replicate summary")
  exact_schema(missingness_pattern, pattern_schema, "Missingness subject-visit pattern")
  exact_schema(warning_events, warning_schema, "Missingness warning events")
  exact_schema(warning_audit, warning_audit_schema, "Missingness warning audit")
  exact_schema(checkpoints, checkpoint_schema, "Missingness checkpoint status")
  manifest_schema <- c("scenario", "rep", "checkpoint", "checkpoint_content_sha256",
    "public_payload_sha256")
  exact_schema(checkpoint_manifest, manifest_schema, "Missingness checkpoint content manifest")
  if (!is.list(checkpoint_payloads) || length(checkpoint_payloads) != nrow(checkpoints) ||
      !is.character(checkpoint_manifest$scenario) || !is.integer(checkpoint_manifest$rep) ||
      !all(vapply(checkpoint_manifest[c("checkpoint", "checkpoint_content_sha256",
        "public_payload_sha256")], is.character, logical(1))) ||
      anyDuplicated(checkpoints$checkpoint) || anyDuplicated(checkpoint_manifest$checkpoint) ||
      any(grepl("^([A-Za-z]:[/\\\\]|/)", checkpoints$checkpoint)) ||
      any(!grepl("^rep_results/.+[.]rds$", gsub("\\\\", "/", checkpoints$checkpoint)))) {
    stop("Missingness durable checkpoint archive or canonical paths are invalid.", call. = FALSE)
  }
  archive_manifest <- do.call(rbind, lapply(checkpoint_payloads, function(result) {
    if (!is.list(result) || !is.list(result$checkpoint_spec)) {
      return(data.frame(scenario = NA_character_, rep = NA_integer_, checkpoint = NA_character_,
        checkpoint_content_sha256 = NA_character_, public_payload_sha256 = NA_character_))
    }
    spec <- result$checkpoint_spec
    task <- data.frame(scenario_id = spec$scenario_id, n = spec$n, d = spec$d,
      rep = spec$replicate, missing_mechanism = spec$missing_mechanism,
      missing_rate = spec$missing_rate, stringsAsFactors = FALSE)
    if (!isTRUE(jss_missing_checkpoint_valid(result, task, spec$configuration))) {
      return(data.frame(scenario = spec$scenario, rep = spec$replicate, checkpoint = NA_character_,
        checkpoint_content_sha256 = NA_character_, public_payload_sha256 = NA_character_))
    }
    jss_missing_checkpoint_archive_record(result)
  }))
  jss_missing_compare_frame(checkpoint_manifest, archive_manifest,
    c("scenario", "rep"), "Missingness checkpoint content manifest")
  jss_missing_compare_frame(checkpoints[manifest_schema], archive_manifest,
    c("scenario", "rep"), "Missingness checkpoint status payload binding")
  archive_provenance <- do.call(rbind, lapply(checkpoint_payloads, function(result) {
    spec <- result$checkpoint_spec; provenance <- result$checkpoint_provenance
    checkpoint_name <- sprintf("%s_scenario%d_%s_miss%02d_n%d_d%d_rep%d.rds",
      result$checkpoint_configuration_key, spec$scenario_id, spec$missing_mechanism,
      round(100 * spec$missing_rate), spec$n, spec$d, spec$replicate)
    identity <- spec$configuration$package_identity
    data.frame(scenario = spec$scenario, rep = spec$replicate,
      checkpoint = gsub("\\\\", "/", file.path("rep_results", checkpoint_name)),
      checkpoint_content_sha256 = result$checkpoint_content_sha256,
      package_source_sha256 = identity$source_sha256, package_version = identity$version,
      package_fingerprint_scope = identity$fingerprint_scope,
      package_source_file_count = identity$source_file_count,
      producer_sha256 = spec$configuration$producer_sha256,
      package_identity_verified = isTRUE(provenance$package_identity_verified),
      timestamp_utc = provenance$timestamp_utc, worker_pid = provenance$pid,
      host = provenance$host, os = provenance$os, platform = provenance$platform,
      r_version = provenance$r_version, rng_kind = provenance$rng_kind,
      blas = provenance$blas, lapack = provenance$lapack,
      rlibs_user_sha256 = provenance$rlibs_user_sha256,
      libpaths_sha256 = provenance$libpaths_sha256, stringsAsFactors = FALSE)
  }))
  provenance_fields <- names(archive_provenance)
  jss_missing_compare_frame(checkpoints[provenance_fields], archive_provenance,
    c("scenario", "rep"), "Missingness durable checkpoint runtime provenance")
  normalize_declared_numeric <- function(x, fields) {
    for (field in fields) if (is.integer(x[[field]]) ||
        (is.logical(x[[field]]) && all(is.na(x[[field]])))) {
      x[[field]] <- as.numeric(x[[field]])
    }
    x
  }
  fixed_raw <- normalize_declared_numeric(fixed_raw,
    c("estimate", "std_error", "true_value", "target_missing_rate", "observed_missing_rate",
      "complete_adjacent_pair_rate"))
  smooth_raw <- normalize_declared_numeric(smooth_raw,
    c("target_missing_rate", "observed_missing_rate", "complete_adjacent_pair_rate", "s1",
      "smooth_hat", "smooth_true"))
  runs <- normalize_declared_numeric(runs,
    c("target_missing_rate", "observed_missing_rate", "complete_adjacent_pair_rate", "logLik", "df",
      "grad_inf", "step_l2", "best_raw_loglik", "raw_loglik_drop_from_best",
      "raw_loglik_drop_tol", "outer_log_lik_change", "outer_stop_crit", "elapsed_sec"))
  warning_events <- normalize_declared_numeric(warning_events, "target_missing_rate")
  if (!all(vapply(checkpoints[c("checkpoint_schema_version", "rep", "package_source_file_count",
      "worker_pid")], is.integer, logical(1))) || !is.logical(checkpoints$package_identity_verified)) {
    stop("Missingness checkpoint status violates exact integer/logical types.", call. = FALSE)
  }
  if (!all(vapply(warning_events[c("n", "d", "rep", "warning_index")], is.integer, logical(1))) ||
      !is.numeric(warning_events$target_missing_rate) || !is.logical(warning_events$expected) ||
      anyDuplicated(warning_events[c("scenario", "rep", "model", "warning_index")])) {
    stop("Missingness warning events violate exact declared types or unique event keys.", call. = FALSE)
  }
  if (is.logical(runs$error) && all(is.na(runs$error))) runs$error <- rep(NA_character_, nrow(runs))
  if (!is.character(runs$error) || !all(vapply(runs[c("success", "converged", "retained")],
      is.logical, logical(1))) || !all(vapply(runs[c("n", "d", "rep", "simulation_seed",
        "missingness_seed")], is.integer, logical(1))) || !is.numeric(runs$logLik) || !is.numeric(runs$df)) {
    stop("Missingness fit log violates canonical CSV boundary types.", call. = FALSE)
  }
  if (!all(vapply(fixed_raw[c("n", "d", "rep")], is.integer, logical(1))) ||
      !all(vapply(fixed_raw[c("estimate", "std_error", "true_value", "target_missing_rate",
        "observed_missing_rate", "complete_adjacent_pair_rate")], is.numeric, logical(1))) ||
      !is.logical(fixed_raw$intercept_includes_fitted_smooth_mean) ||
      !all(vapply(smooth_raw[c("n", "d", "rep")], is.integer, logical(1))) ||
      !all(vapply(smooth_raw[c("target_missing_rate", "observed_missing_rate",
        "complete_adjacent_pair_rate", "s1", "smooth_hat", "smooth_true")], is.numeric, logical(1))) ||
      !all(vapply(missing_by_rep[c("n", "d", "rep", "n_rows", "n_observed_rows", "n_subjects",
        "n_complete_subjects", "n_dropout_subjects", "n_monotone_incomplete_subjects",
        "n_subjects_with_interior_gaps", "total_adjacent_pairs", "complete_adjacent_pairs")],
        is.integer, logical(1))) || !is.logical(missing_by_rep$no_observations_after_dropout) ||
      !all(vapply(missingness_pattern[c("n", "d", "rep", "id")], is.integer, logical(1))) ||
      !is.numeric(missingness_pattern$time) || !is.logical(missingness_pattern$response_observed)) {
    stop("Missingness raw public payloads violate exact declared types.", call. = FALSE)
  }
  expected_mechanisms <- c("monotone_dropout", "time_dependent_intermittent_mar")
  expected_registry <- jss_missing_mechanism_registry()
  expected_registry <- expected_registry[expected_registry$missing_mechanism %in% expected_mechanisms, , drop = FALSE]
  rownames(expected_registry) <- NULL
  if (!setequal(registry$missing_mechanism, expected_mechanisms) ||
      !identical(as.character(registry$analysis_role[match(expected_mechanisms, registry$missing_mechanism)]),
        c("headline", "sensitivity")) || any(c("mar", "time_mar") %in% registry$missing_mechanism)) {
    stop("Missingness public registry has obsolete or unregistered mechanism roles.", call. = FALSE)
  }
  jss_missing_compare_frame(registry, expected_registry, "missing_mechanism", "Missingness mechanism registry")
  jss_missing_validate_registered_design(runs, checkpoints)
  expected_estimands <- jss_missing_estimand_registry()
  jss_missing_compare_frame(estimands, expected_estimands, "component", "Missingness estimand registry")
  if (!all(c("component", "estimand", "observed_data_projection", "reference_distribution") %in% names(estimands)) ||
      any(as.logical(estimands$observed_data_projection))) {
    stop("Missingness public bundle lacks the registered full-data population estimands.", call. = FALSE)
  }
  if (!all(c("penalty", "role", "scale", "application", "justification") %in% names(penalties)) ||
      !setequal(as.numeric(penalties$penalty), c(.5, 1, 2)) ||
      sum(penalties$role == "primary") != 1L || penalties$penalty[penalties$role == "primary"] != 1) {
    stop("Missingness failure-penalty sensitivity registry is not the registered 0.5/1/2 design.", call. = FALSE)
  }
  expected_penalties <- data.frame(
    penalty = c(.5, 1, 2), role = c("sensitivity", "primary", "sensitivity"),
    scale = "link-scale RMSE/IRMSE",
    application = "assigned only when an attempted fit has no retained estimate",
    justification = paste("Penalty values are not estimands; they show how conclusions change when failed or",
      "nonconverged attempts receive prespecified moderate, primary, and severe losses."),
    stringsAsFactors = FALSE)
  jss_missing_compare_frame(penalties, expected_penalties, "penalty",
    "Missingness sensitivity registry")
  required_run <- c("scenario", "model", "rep", "missing_mechanism", "analysis_role",
    "missingness_label", "missingness_pattern", "target_missing_rate", "observed_missing_rate",
    "complete_adjacent_pair_rate", "success", "converged", "retained", "stop_reason",
    "failure_type", "error", "no_observations_after_dropout")
  if (!all(required_run %in% names(runs)) || anyDuplicated(runs[c("scenario", "model", "rep")])) {
    stop("Missingness fit log has an invalid task schema or duplicate keys.", call. = FALSE)
  }
  eligible <- runs$success %in% TRUE & runs$converged %in% TRUE
  if (!identical(as.logical(runs$retained), eligible)) stop("Missingness fit log has inconsistent retention.", call. = FALSE)
  allowed_stops <- c("converged", "relative_deviance_tolerance", "max_iterations", "max_stall",
    "objective_deterioration", "invalid_likelihood", "numerical_failure", "time_limit",
    "outer_iteration_limit_or_invalid_loglik", "fit_error")
  expected_failure_type <- ifelse(eligible, "none", ifelse(runs$success %in% TRUE,
    paste0("optimizer_nonconvergence:", runs$stop_reason),
    ifelse(runs$stop_reason == "fit_error", "fit_error", paste0("fit_error:", runs$stop_reason))))
  if (anyNA(runs[c("success", "converged", "retained")]) ||
      any(!runs$stop_reason %in% allowed_stops) ||
      any(as.character(runs$failure_type) != expected_failure_type) ||
      any(!runs$success & runs$converged) ||
      any(!runs$success & !runs$stop_reason %in% c("fit_error", "time_limit")) ||
      any(!runs$success & (!is.na(runs$logLik) | !is.na(runs$df))) ||
      any(!runs$success & (is.na(runs$error) | !nzchar(runs$error))) ||
      any(runs$success & !is.na(runs$error) & nzchar(runs$error)) ||
      any(eligible & !runs$stop_reason %in% c("converged", "relative_deviance_tolerance"))) {
    stop("Missingness fit log violates the registered status/failure truth table.", call. = FALSE)
  }
  if (nrow(warning_events)) {
    classified <- lapply(warning_events$warning_message, jss_missing_classify_warning)
    classification_valid <- vapply(seq_len(nrow(warning_events)), function(i) {
      expected <- classified[[i]]
      identical(warning_events$normalized_message[[i]], expected$normalized_message) &&
        identical(warning_events$classification[[i]], expected$classification) &&
        identical(warning_events$policy_action[[i]], expected$policy_action) &&
        identical(warning_events$expected[[i]], expected$expected)
    }, logical(1))
    if (!all(classification_valid)) {
      stop("Missingness warning classifications disagree with the registered policy.", call. = FALSE)
    }
  }
  jss_missing_assert_warning_policy(warning_events)
  numeric_columns <- names(runs)[vapply(runs, is.numeric, logical(1))]
  if (any(vapply(runs[numeric_columns], function(x)
      any(!is.na(x) & (!is.finite(x) | abs(x) > 1e12)), logical(1)))) {
    stop("Missingness fit log contains nonfinite or out-of-range numeric values.", call. = FALSE)
  }
  registry_lookup <- registry[match(runs$missing_mechanism, registry$missing_mechanism), , drop = FALSE]
  if (anyNA(registry_lookup$missing_mechanism) ||
      any(runs$analysis_role != registry_lookup$analysis_role) ||
      any(runs$missingness_pattern != registry_lookup$pattern) ||
      any(runs$missingness_label != registry_lookup$display_name)) {
    stop("Missingness fit-log mechanism labels or roles disagree with the registry.", call. = FALSE)
  }
  per_task <- stats::aggregate(runs$model, runs[c("scenario", "rep")], length)
  if (any(per_task$x != 2L) || !all(vapply(split(runs$model, interaction(runs$scenario, runs$rep)),
      function(x) setequal(x, c("gamlss.longitudinal", "gamlss2")), logical(1)))) {
    stop("Missingness public bundle does not contain both models for every attempted task.", call. = FALSE)
  }
  reconstruct <- stats::aggregate(
    cbind(attempted = rep(1L, nrow(runs)), fit_successful = as.integer(runs$success),
      converged = as.integer(runs$converged), retained = as.integer(eligible), failed = as.integer(!eligible)),
    runs[c("scenario", "model")], sum
  )
  expected_attempts <- do.call(rbind, lapply(
    split(runs, interaction(runs$scenario, runs$model, drop = TRUE)),
    function(df) {
      retained <- df$success %in% TRUE & df$converged %in% TRUE
      rate <- mean(retained); se <- sqrt(rate * (1 - rate) / nrow(df))
      data.frame(
        scenario = df$scenario[[1L]], model = df$model[[1L]],
        missing_mechanism = df$missing_mechanism[[1L]],
        missingness_label = df$missingness_label[[1L]],
        missingness_pattern = df$missingness_pattern[[1L]],
        analysis_role = df$analysis_role[[1L]], target_missing_rate = df$target_missing_rate[[1L]],
        attempted = nrow(df), fit_successful = sum(df$success %in% TRUE),
        converged = sum(df$converged %in% TRUE), retained = sum(retained), failed = sum(!retained),
        error_failures = sum(!(df$success %in% TRUE)),
        nonconvergence_failures = sum(df$success %in% TRUE & df$converged %in% FALSE),
        failure_inclusive_retention_rate = rate, retention_rate_mcse = se,
        retention_rate_conf_low = max(0, rate - 1.96 * se),
        retention_rate_conf_high = min(1, rate + 1.96 * se),
        failure_rmse_penalty = penalties$penalty[penalties$role == "primary"],
        stringsAsFactors = FALSE
      )
    }
  ))
  jss_missing_compare_frame(attempts, expected_attempts, c("scenario", "model"),
    "Missingness attempt summary")
  checked <- merge(attempts, reconstruct, by = c("scenario", "model"), suffixes = c("_reported", "_raw"), all = TRUE)
  for (nm in c("attempted", "fit_successful", "converged", "retained", "failed")) {
    if (any(is.na(checked[[paste0(nm, "_reported")]]) |
        checked[[paste0(nm, "_reported")]] != checked[[paste0(nm, "_raw")]])) {
      stop("Missingness attempt summary disagrees with the raw fit log.", call. = FALSE)
    }
  }
  expected_reason <- ifelse(eligible, "retained", ifelse(!runs$success,
    runs$stop_reason,
    ifelse(runs$converged %in% FALSE, paste0("nonconverged: ", runs$stop_reason), "not_retained")))
  expected_failures <- stats::aggregate(rep(1L, nrow(runs)),
    by = data.frame(scenario = runs$scenario, model = runs$model,
      failure_reason = expected_reason, stringsAsFactors = FALSE), FUN = sum)
  names(expected_failures)[[4L]] <- "attempts"
  jss_missing_compare_frame(failures, expected_failures,
    c("scenario", "model", "failure_reason"), "Missingness failure-reason summary")
  failure_totals <- stats::aggregate(failures$attempts, failures[c("scenario", "model")], sum)
  names(failure_totals)[3L] <- "failure_total"
  failure_check <- merge(reconstruct, failure_totals, by = c("scenario", "model"), all = TRUE)
  if (any(failure_check$attempted != failure_check$failure_total)) {
    stop("Missingness failure reasons do not reconcile to attempts.", call. = FALSE)
  }
  monotone <- missing_by_rep$missing_mechanism == "monotone_dropout"
  pattern_required <- c("scenario", "n", "d", "rep", "missing_mechanism", "id", "time", "response_observed")
  if (!all(pattern_required %in% names(missingness_pattern)) || !any(monotone)) {
    stop("Missingness bundle lacks auditable subject-visit missingness patterns.", call. = FALSE)
  }
  monotone_rows <- missingness_pattern$missing_mechanism == "monotone_dropout"
  pattern_monotone <- vapply(
    split(missingness_pattern[monotone_rows, ], interaction(
      missingness_pattern$scenario[monotone_rows], missingness_pattern$rep[monotone_rows],
      missingness_pattern$id[monotone_rows], drop = TRUE
    )),
    function(x) {
      x <- x[order(x$time), , drop = FALSE]
      observed <- as.logical(x$response_observed)
      first_missing <- which(!observed)
      !length(first_missing) || all(!observed[min(first_missing):length(observed)])
    }, logical(1)
  )
  if (!length(pattern_monotone) || !all(pattern_monotone) ||
      any(!as.logical(missing_by_rep$no_observations_after_dropout[monotone]))) {
    stop("Missingness bundle lacks monotone no-post-dropout invariant evidence.", call. = FALSE)
  }
  pattern_tasks <- split(missingness_pattern,
    interaction(missingness_pattern$scenario, missingness_pattern$rep, drop = TRUE))
  reconstructed_patterns <- do.call(rbind, lapply(pattern_tasks, function(x) {
    key <- missing_by_rep$scenario == x$scenario[[1L]] & missing_by_rep$rep == x$rep[[1L]]
    if (sum(key) != 1L) stop("Subject-visit pattern has no unique replicate summary.", call. = FALSE)
    jss_missing_reconstruct_pattern_summary(x, missing_by_rep$target_missing_rate[key][[1L]])
  }))
  summary_fields <- setdiff(names(missing_by_rep), "public_payload_sha256")
  jss_missing_compare_frame(missing_by_rep[summary_fields], reconstructed_patterns[summary_fields],
    c("scenario", "rep"), "Missingness replicate pattern metrics")
  run_task_match <- match(paste(runs$scenario, runs$rep, sep = "\r"),
    paste(missing_by_rep$scenario, missing_by_rep$rep, sep = "\r"))
  run_pattern_fields <- intersect(c("n", "d", "missing_mechanism", "target_missing_rate",
    "observed_missing_rate", "complete_adjacent_pair_rate", "complete_adjacent_pairs",
    "n_dropout_subjects", "n_subjects_with_interior_gaps", "no_observations_after_dropout"),
    names(runs))
  if (anyNA(run_task_match) || any(vapply(run_pattern_fields, function(field) {
    actual <- runs[[field]]; expected <- missing_by_rep[[field]][run_task_match]
    if (is.numeric(expected)) any(is.na(actual) != is.na(expected) |
      (!is.na(actual) & abs(as.numeric(actual) - as.numeric(expected)) > 1e-12))
    else any(as.character(actual) != as.character(expected))
  }, logical(1)))) {
    stop("Missingness fit-log pattern metrics do not reconstruct from subject-visit rows.", call. = FALSE)
  }
  pattern_audit <- do.call(rbind, lapply(pattern_tasks, function(x) data.frame(
    scenario = x$scenario[[1L]], rep = x$rep[[1L]], n = x$n[[1L]], d = x$d[[1L]],
    missing_mechanism = x$missing_mechanism[[1L]], n_rows = nrow(x),
    observed_missing_rate = mean(!as.logical(x$response_observed)),
    duplicate_visits = anyDuplicated(x[c("id", "time")]), stringsAsFactors = FALSE)))
  pattern_check <- merge(missing_by_rep, pattern_audit,
    by = c("scenario", "rep", "n", "d", "missing_mechanism"), suffixes = c("_reported", "_raw"), all = TRUE)
  if (nrow(pattern_check) != nrow(missing_by_rep) || anyNA(pattern_check$duplicate_visits) ||
      any(pattern_check$duplicate_visits != 0L) || any(pattern_check$n_rows_raw != pattern_check$n * pattern_check$d) ||
      any(abs(pattern_check$observed_missing_rate_reported - pattern_check$observed_missing_rate_raw) > 1e-12)) {
    stop("Subject-visit missingness patterns do not reconcile with replicate summaries.", call. = FALSE)
  }
  eligible_keys <- do.call(paste, c(lapply(runs[eligible, c("scenario", "model", "rep")], as.character), sep = "\r"))
  for (payload in list(fixed_raw, smooth_raw)) {
    keys <- do.call(paste, c(lapply(payload[c("scenario", "model", "rep")], as.character), sep = "\r"))
    if (any(!keys %in% eligible_keys)) stop("Missingness metric payload leaks an ineligible task key.", call. = FALSE)
    run_match <- match(keys, do.call(paste, c(lapply(runs[c("scenario", "model", "rep")], as.character), sep = "\r")))
    for (field in intersect(c("n", "d", "missing_mechanism", "target_missing_rate"), names(payload))) {
      equal <- if (is.numeric(payload[[field]]))
        abs(payload[[field]] - runs[[field]][run_match]) <= 1e-12 else
        as.character(payload[[field]]) == as.character(runs[[field]][run_match])
      if (any(is.na(equal) | !equal)) stop("Missingness metric payload metadata disagree with the fit log.", call. = FALSE)
    }
    numeric_payload <- names(payload)[vapply(payload, is.numeric, logical(1))]
    if (any(vapply(payload[numeric_payload], function(x)
        any(!is.na(x) & (!is.finite(x) | abs(x) > 1e12)), logical(1)))) {
      stop("Missingness metric payload contains nonfinite or out-of-range values.", call. = FALSE)
    }
    for (field in intersect(c("s1", "target_missing_rate", "observed_missing_rate",
        "complete_adjacent_pair_rate"), names(payload))) {
      if (any(!is.na(payload[[field]]) & (payload[[field]] < 0 | payload[[field]] > 1)))
        stop("Missingness metric payload contains an impossible rate/grid value.", call. = FALSE)
    }
  }
  segmented_raw <- fixed_raw$model == "gamlss.longitudinal" &
    fixed_raw$missing_mechanism == "time_dependent_intermittent_mar"
  if (!all(c("std_error", "inference_status") %in% names(fixed_raw)) ||
      any(is.finite(fixed_raw$std_error[segmented_raw])) ||
      any(fixed_raw$inference_status[segmented_raw] != "not_applicable_segmented_model_hessian")) {
    stop("Segmented-model raw inference must be auditable as not applicable.", call. = FALSE)
  }
  population_intercept <- fixed_raw$term == "intercept" &
    as.logical(fixed_raw$intercept_includes_fitted_smooth_mean) & !segmented_raw
  if (any(is.finite(fixed_raw$std_error[population_intercept])) ||
      any(fixed_raw$inference_status[population_intercept] !=
        "not_available_without_joint_fixed_smooth_covariance")) {
    stop("Population-intercept inference improperly reuses a coefficient-only standard error.", call. = FALSE)
  }
  if (any(reconstruct$retained < 1L)) stop("Missingness public bundle contains an empty retained cell.", call. = FALSE)
  headline_check <- merge(headline, reconstruct[c("scenario", "model", "attempted")],
    by = c("scenario", "model"), suffixes = c("", "_raw"), all.x = TRUE)
  required_headline <- c("attempted", "retained", "failed", "conditional_rmse_mcse",
    "conditional_rmse_conf_low", "conditional_rmse_conf_high", "failure_inclusive_rmse_mcse",
    "failure_inclusive_rmse_conf_low", "failure_inclusive_rmse_conf_high", "inference_status")
  if (!all(required_headline %in% names(headline)) || any(headline_check$attempted != headline_check$attempted_raw) ||
      any(headline$retained <= 0 | headline$failed != headline$attempted - headline$retained) ||
      any(!headline$failure_penalty %in% penalties$penalty)) {
    stop("Missingness headline denominators or penalty sensitivity do not reconcile.", call. = FALSE)
  }
  fixed_nonintercept <- fixed_raw[fixed_raw$term != "intercept", , drop = FALSE]
  fixed_nonintercept$squared_error <- (fixed_nonintercept$estimate - fixed_nonintercept$true_value)^2
  rep_rmse <- stats::aggregate(fixed_nonintercept$squared_error,
    fixed_nonintercept[c("scenario", "model", "parameter", "rep")],
    function(x) sqrt(mean(x, na.rm = TRUE)))
  names(rep_rmse)[names(rep_rmse) == "x"] <- "replicate_rmse"
  rep_rmse <- rep_rmse[is.finite(rep_rmse$replicate_rmse), , drop = FALSE]
  expected_headline <- do.call(rbind, lapply(
    split(rep_rmse, interaction(rep_rmse$scenario, rep_rmse$model, rep_rmse$parameter, drop = TRUE)),
    function(df) {
      run <- runs[runs$scenario == df$scenario[[1L]] & runs$model == df$model[[1L]], , drop = FALSE]
      conditional <- jss_missing_mc_interval(df$replicate_rmse)
      do.call(rbind, lapply(penalties$penalty, function(penalty) {
        failures <- nrow(run) - conditional[["n"]]
        inclusive <- jss_missing_mc_interval(c(df$replicate_rmse, rep(penalty, failures)))
        data.frame(
          scenario = df$scenario[[1L]], model = df$model[[1L]], parameter = df$parameter[[1L]],
          missing_mechanism = run$missing_mechanism[[1L]], missingness_label = run$missingness_label[[1L]],
          missingness_pattern = run$missingness_pattern[[1L]], analysis_role = run$analysis_role[[1L]],
          target_missing_rate = run$target_missing_rate[[1L]], attempted = nrow(run),
          retained = unname(conditional[["n"]]), failed = failures,
          conditional_mean_rmse = unname(conditional[["estimate"]]),
          conditional_rmse_mcse = unname(conditional[["mcse"]]),
          conditional_rmse_conf_low = unname(conditional[["conf_low"]]),
          conditional_rmse_conf_high = unname(conditional[["conf_high"]]),
          failure_penalty = penalty,
          failure_inclusive_mean_rmse = unname(inclusive[["estimate"]]),
          failure_inclusive_rmse_mcse = unname(inclusive[["mcse"]]),
          failure_inclusive_rmse_conf_low = unname(inclusive[["conf_low"]]),
          failure_inclusive_rmse_conf_high = unname(inclusive[["conf_high"]]),
          inference_status = if (run$model[[1L]] == "gamlss.longitudinal" &&
            run$missingness_pattern[[1L]] == "intermittent")
            "not_applicable_segmented_model_hessian" else "conditional_on_retained_fits",
          stringsAsFactors = FALSE
        )
      }))
    }
  ))
  jss_missing_compare_frame(headline, expected_headline,
    c("scenario", "model", "parameter", "failure_penalty"), "Missingness headline summary")
  for (prefix in c("conditional", "failure_inclusive")) {
    estimate <- headline[[paste0(prefix, "_mean_rmse")]]
    mcse <- headline[[paste0(prefix, "_rmse_mcse")]]
    low <- headline[[paste0(prefix, "_rmse_conf_low")]]
    high <- headline[[paste0(prefix, "_rmse_conf_high")]]
    complete <- is.finite(estimate) | is.finite(mcse) | is.finite(low) | is.finite(high)
    if (any(complete & (!is.finite(estimate) | !is.finite(mcse) | mcse < 0 |
        !is.finite(low) | !is.finite(high) | low > estimate | estimate > high))) {
      stop("Missingness headline summary contains an impossible Monte Carlo interval.", call. = FALSE)
    }
  }
  segmented <- headline$model == "gamlss.longitudinal" & headline$missingness_pattern == "intermittent"
  if (any(headline$inference_status[segmented] != "not_applicable_segmented_model_hessian")) {
    stop("Segmented-model Hessian coverage must be reported as not applicable.", call. = FALSE)
  }
  fixed_term_raw <- fixed_raw[fixed_raw$term != "intercept", , drop = FALSE]
  fixed_term_level <- do.call(rbind, lapply(
    split(fixed_term_raw, interaction(fixed_term_raw$scenario, fixed_term_raw$model,
      fixed_term_raw$parameter, fixed_term_raw$term, drop = TRUE)),
    function(df) {
      err <- df$estimate - df$true_value
      covered <- df$true_value >= df$estimate - stats::qnorm(.975) * df$std_error &
        df$true_value <= df$estimate + stats::qnorm(.975) * df$std_error
      run <- runs[runs$scenario == df$scenario[[1L]] & runs$model == df$model[[1L]], , drop = FALSE]
      segmented_na <- run$model[[1L]] == "gamlss.longitudinal" && run$missingness_pattern[[1L]] == "intermittent"
      data.frame(scenario = df$scenario[[1L]], model = df$model[[1L]], parameter = df$parameter[[1L]],
        term = df$term[[1L]], bias = mean(err, na.rm = TRUE), rmse = sqrt(mean(err^2, na.rm = TRUE)),
        coverage_95 = if (segmented_na || !any(is.finite(df$std_error))) NA_real_ else mean(covered, na.rm = TRUE),
        n_successful_fits = sum(is.finite(df$estimate)), stringsAsFactors = FALSE)
    }
  ))
  expected_fixed <- do.call(rbind, lapply(
    split(fixed_term_level, interaction(fixed_term_level$scenario, fixed_term_level$model,
      fixed_term_level$parameter, drop = TRUE)),
    function(df) {
      run <- runs[runs$scenario == df$scenario[[1L]] & runs$model == df$model[[1L]], , drop = FALSE]
      data.frame(
        missing_mechanism = run$missing_mechanism[[1L]], target_missing_rate = run$target_missing_rate[[1L]],
        observed_missing_rate = mean(run$observed_missing_rate, na.rm = TRUE),
        complete_adjacent_pair_rate = mean(run$complete_adjacent_pair_rate, na.rm = TRUE),
        model = df$model[[1L]], parameter = df$parameter[[1L]], attempted = nrow(run),
        retained = min(df$n_successful_fits, na.rm = TRUE),
        failed = max(nrow(run) - df$n_successful_fits, na.rm = TRUE),
        mean_rmse = mean(df$rmse, na.rm = TRUE), mean_abs_bias = mean(abs(df$bias), na.rm = TRUE),
        mean_coverage_95 = mean(df$coverage_95, na.rm = TRUE),
        inference_status = if (run$model[[1L]] == "gamlss.longitudinal" &&
          run$missingness_pattern[[1L]] == "intermittent")
          "not_applicable_segmented_model_hessian" else "conditional_on_retained_fits",
        target_missing_pct = 100 * run$target_missing_rate[[1L]], stringsAsFactors = FALSE
      )
    }
  ))
  expected_fixed <- expected_fixed[is.finite(expected_fixed$mean_rmse), , drop = FALSE]
  jss_missing_compare_frame(fixed_reported, expected_fixed,
    c("missing_mechanism", "target_missing_rate", "model", "parameter"),
    "Missingness fixed-term public summary")

  smooth_integrated <- do.call(rbind, lapply(
    split(smooth_raw, interaction(smooth_raw$scenario, smooth_raw$model,
      smooth_raw$parameter, smooth_raw$rep, drop = TRUE)), function(df) {
      err <- df$smooth_hat - df$smooth_true
      data.frame(scenario = df$scenario[[1L]], model = df$model[[1L]], parameter = df$parameter[[1L]],
        rep = df$rep[[1L]], irmse = sqrt(mean(err^2, na.rm = TRUE)), stringsAsFactors = FALSE)
    }))
  expected_smooth_irmse <- do.call(rbind, lapply(
    split(smooth_integrated, interaction(smooth_integrated$scenario, smooth_integrated$model,
      smooth_integrated$parameter, drop = TRUE)), function(df) {
      run <- runs[runs$scenario == df$scenario[[1L]] & runs$model == df$model[[1L]], , drop = FALSE]
      conditional <- jss_missing_mc_interval(df$irmse)
      do.call(rbind, lapply(penalties$penalty, function(penalty) {
        failed <- nrow(run) - conditional[["n"]]
        inclusive <- jss_missing_mc_interval(c(df$irmse, rep(penalty, failed)))
        data.frame(
          missing_mechanism = run$missing_mechanism[[1L]], target_missing_rate = run$target_missing_rate[[1L]],
          observed_missing_rate = mean(run$observed_missing_rate, na.rm = TRUE),
          complete_adjacent_pair_rate = mean(run$complete_adjacent_pair_rate, na.rm = TRUE),
          model = df$model[[1L]], parameter = df$parameter[[1L]],
          mean_irmse = conditional[["estimate"]], irmse_mcse = conditional[["mcse"]],
          irmse_conf_low = conditional[["conf_low"]], irmse_conf_high = conditional[["conf_high"]],
          sd_irmse = stats::sd(df$irmse, na.rm = TRUE), attempted = nrow(run),
          n_metric_successful = conditional[["n"]], metric_failures = failed,
          failure_inclusive_mean_irmse = inclusive[["estimate"]],
          failure_inclusive_irmse_mcse = inclusive[["mcse"]],
          failure_inclusive_irmse_conf_low = inclusive[["conf_low"]],
          failure_inclusive_irmse_conf_high = inclusive[["conf_high"]],
          failure_rmse_penalty = penalty, stringsAsFactors = FALSE)
      }))
    }))
  jss_missing_compare_frame(smooth_irmse_reported, expected_smooth_irmse,
    c("missing_mechanism", "target_missing_rate", "model", "parameter", "failure_rmse_penalty"),
    "Missingness smooth IRMSE summary")
  for (prefix in c("", "failure_inclusive_")) {
    estimate_name <- if (nzchar(prefix)) paste0(prefix, "mean_irmse") else "mean_irmse"
    stem <- if (nzchar(prefix)) paste0(prefix, "irmse") else "irmse"
    estimate <- smooth_irmse_reported[[estimate_name]]; mcse <- smooth_irmse_reported[[paste0(stem, "_mcse")]]
    low <- smooth_irmse_reported[[paste0(stem, "_conf_low")]]; high <- smooth_irmse_reported[[paste0(stem, "_conf_high")]]
    if (any(!is.finite(estimate) | !is.finite(mcse) | mcse < 0 |
        !is.finite(low) | !is.finite(high) | low > estimate | estimate > high))
      stop("Missingness smooth IRMSE summary contains impossible Monte Carlo intervals.", call. = FALSE)
  }

  smooth_pointwise <- do.call(rbind, lapply(
    split(smooth_raw, interaction(smooth_raw$scenario, smooth_raw$model,
      smooth_raw$parameter, smooth_raw$s1, drop = TRUE)),
    function(df) {
      err <- df$smooth_hat - df$smooth_true
      data.frame(scenario = df$scenario[[1L]], model = df$model[[1L]], n = df$n[[1L]], d = df$d[[1L]],
        parameter = df$parameter[[1L]], s1 = df$s1[[1L]], smooth_true = df$smooth_true[[1L]],
        smooth_median = stats::median(df$smooth_hat, na.rm = TRUE),
        smooth_q05 = as.numeric(stats::quantile(df$smooth_hat, .05, na.rm = TRUE, names = FALSE)),
        smooth_q95 = as.numeric(stats::quantile(df$smooth_hat, .95, na.rm = TRUE, names = FALSE)),
        bias = mean(err, na.rm = TRUE), rmse = sqrt(mean(err^2, na.rm = TRUE)), stringsAsFactors = FALSE)
    }
  ))
  selected <- smooth_pointwise[smooth_pointwise$scenario %in%
    unique(runs$scenario[runs$target_missing_rate %in% c(0, .3, .5)]), , drop = FALSE]
  scenario_meta <- do.call(rbind, lapply(split(runs, runs$scenario), function(df) data.frame(
    scenario = df$scenario[[1L]], missing_mechanism = df$missing_mechanism[[1L]],
    missingness_label = df$missingness_label[[1L]], missingness_pattern = df$missingness_pattern[[1L]],
    analysis_role = df$analysis_role[[1L]], target_missing_rate = df$target_missing_rate[[1L]],
    observed_missing_rate = mean(df$observed_missing_rate, na.rm = TRUE),
    complete_adjacent_pair_rate = mean(df$complete_adjacent_pair_rate, na.rm = TRUE), stringsAsFactors = FALSE)))
  expected_smooth <- merge(selected, scenario_meta, by = "scenario", all = FALSE, sort = FALSE)
  expected_smooth <- merge(expected_smooth, expected_attempts[c("scenario", "model", "attempted", "retained", "failed",
    "failure_inclusive_retention_rate", "retention_rate_mcse", "retention_rate_conf_low", "retention_rate_conf_high")],
    by = c("scenario", "model"), all.x = TRUE, sort = FALSE)
  expected_smooth$missing_label <- paste0(round(100 * expected_smooth$target_missing_rate), "%")
  jss_missing_compare_frame(smooth_reported, expected_smooth,
    c("scenario", "model", "parameter", "s1"), "Missingness smooth public plot data")
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  expected_package <- jss_missing_checkout_identity(root)
  expected_producer <- jss_missing_producer_sha256(c(
    file.path(root, "paper", "R", "missingness-study-helpers.R"),
    file.path(root, "paper", "scripts", "final-simulations", "missingness", "run_missingness_study.R")
  ))
  observed_runtime <- jss_missing_runtime_identity()
  hex <- "^[0-9a-f]{64}$"
  if (!all(c("checkpoint_schema_version", "checkpoint_content_sha256", "package_source_sha256",
      "public_payload_sha256", "producer_sha256", "package_identity_verified", "timestamp_utc", "worker_pid", "host", "os", "platform", "r_version",
      "rng_kind", "blas", "lapack", "rlibs_user_sha256", "libpaths_sha256",
      "package_version", "package_fingerprint_scope",
      "package_source_file_count") %in% names(checkpoints)) || any(checkpoints$checkpoint_schema_version != 5L) ||
      any(!as.logical(checkpoints$package_identity_verified)) ||
      any(!grepl(hex, checkpoints$checkpoint_content_sha256)) || any(!grepl(hex, checkpoints$package_source_sha256)) ||
      any(!grepl(hex, checkpoints$producer_sha256)) || any(!grepl(hex, checkpoints$public_payload_sha256)) ||
      any(checkpoints$package_source_sha256 != expected_package$source_sha256) ||
      any(checkpoints$package_version != expected_package$version) ||
      any(checkpoints$package_fingerprint_scope != expected_package$fingerprint_scope) ||
      any(checkpoints$package_source_file_count != expected_package$source_file_count) ||
      any(checkpoints$producer_sha256 != expected_producer) ||
      any(checkpoints$rlibs_user_sha256 != observed_runtime$rlibs_user_sha256) ||
      any(checkpoints$libpaths_sha256 != observed_runtime$libpaths_sha256) ||
      any(checkpoints$host != observed_runtime$host) ||
      any(checkpoints$os != observed_runtime$os) ||
      any(checkpoints$platform != observed_runtime$platform) ||
      any(checkpoints$r_version != observed_runtime$r_version) ||
      any(checkpoints$rng_kind != observed_runtime$rng_kind) ||
      any(checkpoints$blas != observed_runtime$blas) ||
      any(checkpoints$lapack != observed_runtime$lapack) ||
      any(!is.finite(checkpoints$worker_pid) | checkpoints$worker_pid < 1 | checkpoints$worker_pid != floor(checkpoints$worker_pid)) ||
      any(grepl("^([A-Za-z]:[/\\\\]|/)", checkpoints$blas)) || any(grepl("^([A-Za-z]:[/\\\\]|/)", checkpoints$lapack)) ||
      any(vapply(checkpoints[c("timestamp_utc", "host", "os", "platform", "r_version", "rng_kind", "blas", "lapack",
        "rlibs_user_sha256", "libpaths_sha256")],
        function(x) any(is.na(x) | !nzchar(as.character(x))), logical(1)))) {
    stop("Missingness checkpoint provenance or content hashes are incomplete.", call. = FALSE)
  }
  jss_missing_parse_rfc3339_utc(checkpoints$timestamp_utc, "Missingness checkpoint timestamp")
  checkpoint_key <- paste(checkpoints$scenario, checkpoints$rep, sep = "\r")
  run_task_key <- unique(paste(runs$scenario, runs$rep, sep = "\r"))
  if (anyDuplicated(checkpoint_key) || !setequal(checkpoint_key, run_task_key)) {
    stop("Missingness checkpoint-status rows do not reconcile with attempted tasks.", call. = FALSE)
  }
  if (!exists("jss_missing_portable_task_sha256", mode = "function")) {
    stop("Missingness public payload hash validator is unavailable.", call. = FALSE)
  }
  for (i in seq_len(nrow(checkpoints))) {
    scenario <- checkpoints$scenario[[i]]; rep_id <- checkpoints$rep[[i]]
    select <- function(x) x[x$scenario == scenario & x$rep == rep_id, , drop = FALSE]
    actual <- jss_missing_portable_task_sha256(select(runs), select(fixed_raw), select(smooth_raw),
      select(missing_by_rep), select(missingness_pattern), select(warning_events))
    recorded <- checkpoints$public_payload_sha256[[i]]
    raw_recorded <- unique(c(select(runs)$public_payload_sha256, select(fixed_raw)$public_payload_sha256,
      select(smooth_raw)$public_payload_sha256, select(missing_by_rep)$public_payload_sha256,
      select(missingness_pattern)$public_payload_sha256,
      select(warning_events)$public_payload_sha256))
    if (!identical(actual, recorded) || length(raw_recorded) != 1L || !identical(as.character(raw_recorded), recorded)) {
      stop("Missingness raw public rows do not validate against checkpoint payload hashes.", call. = FALSE)
    }
  }
  coverage_payload <- jss_missing_coverage_plot_payload(fixed_reported,
    c("mu", "sigma", "nu", "tau"))
  expected_warning_audit <- jss_missing_warning_audit(warning_events, runs, coverage_payload)
  jss_missing_compare_frame(warning_audit, expected_warning_audit,
    c("audit_type", "classification"), "Missingness warning/plot omission audit", tolerance = 0)
  if (isTRUE(require_promotion)) {
    jss_missing_require_external_attestation(bundle_sha256, expected_package$source_sha256,
      expected_producer, checkpoints, root, attestation_path, signature_path)
  }
  final_manifest <- paste(sort(paste(basename(required_paths),
    vapply(required_paths, jss_missing_sha256_file, character(1)), sep = "\t")), collapse = "\n")
  final_package <- jss_missing_checkout_identity(root)
  final_producer <- jss_missing_producer_sha256(c(
    file.path(root, "paper", "R", "missingness-study-helpers.R"),
    file.path(root, "paper", "scripts", "final-simulations", "missingness", "run_missingness_study.R")
  ))
  if (!identical(digest::digest(final_manifest, "sha256", serialize = FALSE), bundle_sha256)) {
    stop("Missingness staged snapshot changed during validation.", call. = FALSE)
  }
  if (!identical(final_package$source_sha256, expected_package$source_sha256) ||
      !identical(final_producer, expected_producer)) {
    stop("Missingness checkout/producer identity changed during validation.", call. = FALSE)
  }
  list(registry = registry, estimands = estimands, runs = runs, attempts = attempts,
    missingness = missing_by_rep, missingness_pattern = missingness_pattern,
    fixed = fixed_raw, smooth = smooth_raw, headline = headline,
    warning_events = warning_events, warning_audit = warning_audit)
}

jss_missingness_validate_public_bundle <- function(input, root = getwd(),
    attestation_path = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MISSINGNESS_ATTESTATION", unset = ""),
    signature_path = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MISSINGNESS_ATTESTATION_SIGNATURE", unset = "")) {
  original <- file.path(input, jss_missingness_public_artifacts())
  if (any(!file.exists(original))) {
    stop("Missingness public bundle is incomplete.", call. = FALSE)
  }
  before <- vapply(original, jss_missing_sha256_file, character(1))
  stage <- tempfile("missingness-immutable-snapshot-")
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
  staged <- file.path(stage, basename(original))
  if (!all(file.copy(original, staged, overwrite = FALSE)) ||
      !identical(before, vapply(original, jss_missing_sha256_file, character(1))) ||
      !identical(unname(before), unname(vapply(staged, jss_missing_sha256_file, character(1))))) {
    stop("Missingness source bundle changed while creating its immutable validation snapshot.", call. = FALSE)
  }
  jss_missingness_validate_candidate_bundle(stage, root, attestation_path,
    signature_path, require_promotion = TRUE)
}

jss_missingness_production_gate <- function(input, root = getwd()) {
  validated <- jss_missingness_validate_public_bundle(input, root = root)
  if (nrow(validated$runs) != 480L || nrow(validated$missingness) != 240L ||
      nrow(validated$attempts) != 24L) {
    stop("Missingness production aggregate bundle gate did not reconcile registered dimensions.", call. = FALSE)
  }
  invisible(validated)
}

jss_stage_04_missingness_data <- function(input, settings) {
  names <- jss_missingness_public_artifacts()
  source <- file.path(input, names)
  missing <- source[!file.exists(source)]
  if (length(missing)) {
    stop("Missingness evidence bundle is incomplete: ", paste(basename(missing), collapse = ", "), call. = FALSE)
  }
  destination_dir <- file.path(settings$data_dir, "missingness")
  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
  destination <- file.path(destination_dir, names)
  copied <- file.copy(source, destination, overwrite = TRUE)
  if (any(!copied)) stop("Could not stage missingness evidence.", call. = FALSE)
  destination
}

jss_run_04_missingness_from_public_inputs <- function(settings) {
  input <- file.path(settings$public_data_dir, "missingness")
  required <- jss_missingness_public_artifacts()
  missing <- required[!file.exists(file.path(input, required))]
  if (length(missing)) {
    stop(
      "The tracked missingness public-input bundle is legacy and cannot be used for paper mode. ",
      "Run the current full missingness producer and promote all required outputs. Missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  validated <- jss_missingness_validate_public_bundle(input)
  fixed <- utils::read.csv(file.path(input, "fixed_term_summary_by_missingness.csv"), stringsAsFactors = FALSE)
  smooth <- utils::read.csv(file.path(input, "smooth_selected_plot_data.csv"), stringsAsFactors = FALSE)
  attempts <- validated$attempts
  headline <- validated$headline
  registry <- validated$registry
  expected_mechanisms <- c("monotone_dropout", "time_dependent_intermittent_mar")
  observed_mechanisms <- unique(c(
    as.character(fixed$missing_mechanism), as.character(smooth$missing_mechanism),
    as.character(attempts$missing_mechanism), as.character(headline$missing_mechanism)
  ))
  if (any(c("mar", "time_mar") %in% observed_mechanisms) ||
      !setequal(observed_mechanisms, expected_mechanisms) ||
      !identical(as.character(registry$missing_mechanism[1:2]), expected_mechanisms)) {
    stop(
      "The missingness public-input bundle uses obsolete mechanisms. Paper mode requires ",
      "monotone_dropout and time_dependent_intermittent_mar only.",
      call. = FALSE
    )
  }
  denominator_fields <- c("attempted", "retained", "failed")
  if (!all(denominator_fields %in% names(attempts)) ||
      !all(c(denominator_fields, "conditional_rmse_mcse", "failure_inclusive_rmse_mcse") %in% names(headline))) {
    stop("Missingness public inputs lack required denominators or Monte Carlo uncertainty.", call. = FALSE)
  }
  fixed$model <- factor(fixed$model, levels = c("gamlss2", "gamlss.longitudinal"))
  p1 <- ggplot2::ggplot(fixed[fixed$parameter %in% c("mu", "sigma", "nu", "tau"), ],
    ggplot2::aes(x = target_missing_pct, y = mean_rmse, colour = model, shape = model, group = model)) +
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 1.9) +
    ggplot2::facet_grid(parameter ~ missing_mechanism, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = seq(0, 50, 10)) +
    ggplot2::labs(x = "Target missingness (%)", y = "Fixed-effect RMSE", colour = NULL, shape = NULL) +
    ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "top")
  f1 <- file.path(settings$figures_dir, "fixed_margin_rmse_by_missingness.png")
  ggplot2::ggsave(f1, p1, width = 8.5, height = 5.5, dpi = 180, bg = "white")

  smooth$missing_label <- factor(smooth$missing_label, levels = c("0%", "30%", "50%"))
  smooth$model <- factor(smooth$model, levels = c("gamlss2", "gamlss.longitudinal"))
  p2 <- ggplot2::ggplot(smooth, ggplot2::aes(x = s1)) +
    ggplot2::geom_line(ggplot2::aes(y = smooth_true), colour = "black", linewidth = 0.75) +
    ggplot2::geom_line(ggplot2::aes(y = smooth_median, colour = model), linewidth = 0.75, linetype = "dashed") +
    ggplot2::facet_grid(parameter + missing_mechanism ~ missing_label, scales = "free_y") +
    ggplot2::labs(x = "s1", y = "Centered smooth", colour = NULL) +
    ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "top")
  f2 <- file.path(settings$figures_dir, "smooth_selected_recovery_curves.png")
  ggplot2::ggsave(f2, p2, width = 8.5, height = 7, dpi = 180, bg = "white")
  staged_data <- jss_stage_04_missingness_data(input, settings)
  list(module_id = "04-missingness-dropout-sensitivity", status = "regenerated",
    data = staged_data,
    tables = character(), figures = c(f1, f2))
}
