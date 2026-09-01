# Shared helpers for the JSS multivariate longitudinal simulation study.
#
# This module is intentionally script-local. It adds paper-workflow helpers
# without changing the package-facing API.

`%||%` <- function(a, b) if (!is.null(a)) a else b

mvt_find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "DESCRIPTION")) &&
        dir.exists(file.path(current, "R"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find repository root from ", start, call. = FALSE)
    }
    current <- parent
  }
}

mvt_repo_root <- mvt_find_repo_root()
mvt_script_dir <- file.path(
  mvt_repo_root,
  "paper",
  "R",
  "09-simulation-multivariate-longitudinal"
)
mvt_output_root <- file.path(
  mvt_repo_root,
  "results",
  "jss-exploratory",
  "09-simulation-multivariate-longitudinal"
)

mvt_env <- function(name, default = "") {
  value <- Sys.getenv(name, unset = default)
  if (length(value) == 0L || is.na(value) || !nzchar(value)) default else value
}

mvt_env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(mvt_env(name, as.character(default))))
  if (length(value) == 0L || is.na(value) || !is.finite(value)) default else value
}

mvt_env_num <- function(name, default) {
  value <- suppressWarnings(as.numeric(mvt_env(name, as.character(default))))
  if (length(value) == 0L || is.na(value) || !is.finite(value)) default else value
}

mvt_env_flag <- function(name, default = FALSE) {
  value <- tolower(mvt_env(name, if (isTRUE(default)) "true" else "false"))
  value %in% c("1", "true", "t", "yes", "y")
}

mvt_env_vector <- function(name, default = character()) {
  value <- mvt_env(name, "")
  if (!nzchar(value)) return(default)
  trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
}

mvt_default_comparators <- function() {
  c(
    "glm",
    "glmm",
    "gee_independence",
    "gee_exchangeable",
    "gee_ar1",
    "gee_unstructured",
    "gamlss.longitudinal",
    "gamCopula_markov",
    "gamCopula_vine_simplified"
  )
}

mvt_allowed_comparators <- function() {
  c(mvt_default_comparators(), "gamCopula", "gamCopula_vine", "glmm_slope")
}

mvt_active_comparators <- function() {
  requested <- mvt_env_vector("GAMLSS_LONGITUDINAL_MVT_COMPARATORS", mvt_default_comparators())
  if ("gee" %in% requested) {
    requested <- c(requested, grep("^gee_", mvt_default_comparators(), value = TRUE))
  }
  requested[requested == "gamCopula"] <- "gamCopula_markov"
  requested <- setdiff(unique(requested), "gee")
  unname(mvt_select_named(stats::setNames(mvt_allowed_comparators(), mvt_allowed_comparators()), requested))
}

mvt_standard_comparators <- function(active = mvt_active_comparators()) {
  intersect(c("glm", "glmm"), active)
}

mvt_gee_comparators <- function(active = mvt_active_comparators()) {
  sub("^gee_", "", grep("^gee_", active, value = TRUE))
}

mvt_env_int_vector <- function(name, default = integer()) {
  value <- mvt_env_vector(name, character())
  if (length(value) == 0L) return(default)
  out <- suppressWarnings(as.integer(value))
  out <- out[is.finite(out) & out > 0L]
  if (length(out) == 0L) default else unique(out)
}

mvt_rep_ids <- function(default_n) {
  explicit <- mvt_env_int_vector("GAMLSS_LONGITUDINAL_MVT_REP_IDS", integer())
  if (length(explicit) > 0L) return(sort(unique(explicit)))
  seq_len(as.integer(default_n))
}

mvt_timestamp <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

mvt_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

mvt_result_names <- function() {
  c(
    "fit_status",
    "benchmark_results",
    "coefficient_results",
    "dependence_recovery",
    "variogram_scores",
    "runtime"
  )
}

mvt_expected_output_files <- function() {
  c(
    "scenario_grid.csv",
    paste0(mvt_result_names(), "_by_rep.csv"),
    "benchmark_summary.csv",
    "coefficient_summary.csv",
    "dependence_recovery_summary.csv",
    "variogram_summary.csv",
    "case_method_completion_summary.csv",
    "pilot_feasibility_by_method.csv",
    "pilot_feasibility_by_scenario.csv",
    "pilot_feasibility_overall_method.csv",
    "pilot_feasibility.md",
    "preflight_checks.csv",
    "preflight_checks.md",
    "artifact_manifest.csv",
    "artifact_manifest.md",
    "run_metadata.csv",
    "package_versions.csv",
    "session_info.txt"
  )
}

mvt_study_protocol_tables <- function() {
  time_specs <- mvt_time_specs()
  time_table <- do.call(rbind, lapply(names(time_specs), function(name) {
    spec <- time_specs[[name]]
    data.frame(
      time_grid = name,
      n_time = spec$n_time,
      n_subject = spec$n_subject,
      total_rows = spec$total_rows,
      stringsAsFactors = FALSE
    )
  }))

  family_specs <- mvt_family_specs(include_special = TRUE)
  family_table <- do.call(rbind, lapply(names(family_specs), function(name) {
    spec <- family_specs[[name]]
    data.frame(
      family_name = name,
      label = spec$label %||% name,
      gamlss_family = spec$gamlss_family %||% "",
      standard_comparator_available = !is.null(spec$standard_family),
      special = isTRUE(spec$special),
      mean_formula = "response ~ time + x + z",
      stringsAsFactors = FALSE
    )
  }))

  dependence_specs <- mvt_dependence_specs(include_appendix = TRUE)
  dependence_table <- do.call(rbind, lapply(names(dependence_specs), function(name) {
    spec <- dependence_specs[[name]]
    data.frame(
      dependence_name = name,
      generator = spec$generator %||% "",
      dependence = spec$dependence %||% "",
      correlation_level = spec$correlation_level %||% "",
      appendix = isTRUE(spec$appendix),
      gee_correlations = paste(spec$gee_correlations %||% character(), collapse = ","),
      stringsAsFactors = FALSE
    )
  }))

  comparator_table <- data.frame(
    method = mvt_allowed_comparators(),
    role = c(
      rep("main", length(mvt_default_comparators())),
      "legacy_alias",
      "targeted_sensitivity",
      "appendix_sensitivity"
    ),
    stringsAsFactors = FALSE
  )
  comparator_table$description <- c(
    "GLM independence mean model",
    "Random-intercept GLMM primary comparator",
    "GEE independence working correlation",
    "GEE exchangeable working correlation",
    "GEE AR(1) working correlation",
    "GEE unstructured working correlation with timeout",
    "Joint GAMLSS-copula longitudinal model",
    "Two-stage gamCopula adjacent Markov comparator",
    "Two-stage gamCopula simplified-vine comparator",
    "Legacy shortcut for gamCopula_markov",
    "Two-stage gamCopula full-vine comparator for covariate-dependent dependence sensitivity",
    "Random-intercept plus random-time-slope GLMM sensitivity"
  )

  output_table <- data.frame(
    artifact = mvt_expected_output_files(),
    level = "run",
    stringsAsFactors = FALSE
  )
  suite_output <- data.frame(
    artifact = c(
      "README.md",
      "publication_suite_plan.csv",
      "publication_suite_plan.md",
      "publication_suite_artifacts.csv",
      "publication_suite_preflight.csv",
      "publication_suite_preflight.md",
      "publication_readiness_audit.csv",
      "publication_readiness_audit.md"
    ),
    level = "suite",
    stringsAsFactors = FALSE
  )
  output_table <- rbind(output_table, suite_output)

  metrics_table <- data.frame(
    metric_family = c(
      "fit_status",
      "mean_prediction",
      "distribution",
      "coefficients",
      "dependence",
      "variogram",
      "runtime"
    ),
    metrics = c(
      "success, warning, timeout, error classification",
      "MAE, RMSE, true-mean RMSE",
      "negative log score, PIT calibration, q90 error, interval coverage and width",
      "bias, reported SE calibration, 95% CI coverage",
      "adjacent and all-pair theta/tau/correlation recovery",
      "variogram scores at p=0.5 and p=2 where simulation is available",
      "elapsed seconds, timeout/failure rates, projected runtime from pilot"
    ),
    stringsAsFactors = FALSE
  )

  readiness <- rbind(
    transform(mvt_publication_readiness_spec(), evidence = "required_or_recommended_run"),
    transform(mvt_publication_readiness_optional_spec(), min_reps = NA_integer_, required_time = "", required_dependence = "", evidence = "optional_sensitivity_run")
  )

  review_checklist <- data.frame(
    requirement = c(
      "Total-row controlled multivariate time grids",
      "Four main marginal distributions",
      "GAMLSS-only distributional margin",
      "External Gaussian-copula/MVN generator",
      "Native gamlss.longitudinal simulator track",
      "Dependence scenarios cover exchangeable, AR(1), time-varying, and covariate-dependent structures",
      "Standard GLM/GEE/GLMM comparators",
      "gamCopula Markov, simplified-vine, and targeted full-vine comparators",
      "gamlss.longitudinal comparator",
      "GJRM excluded from high-dimensional comparator grid",
      "Primary coefficient, prediction, calibration, dependence, variogram, and runtime metrics",
      "Resume-safe per-replicate outputs",
      "Pilot feasibility before main run",
      "Publication readiness gate",
      "Review bundle and artifact manifest",
      "No missingness/dropout in this extension"
    ),
    implementation_evidence = c(
      "mvt_time_specs(); mvt_expand_grid()",
      "mvt_family_specs(include_special = TRUE)",
      "gg_continuous entry in mvt_family_specs() with sigma/nu covariate formulas",
      "mvt_simulate_external(); external_* dependence specs",
      "mvt_simulate_native(); native_* dependence specs",
      "mvt_dependence_specs(include_appendix = TRUE)",
      "mvt_fit_standard_models(); mvt_run_one_gee(); mvt_run_glmm_slope()",
      "mvt_fit_gamcopula(); mvt_fit_gamcopula_vine(); mvt_simulate_gamcopula_response(); gamCopula runtime dependency checks",
      "mvt_fit_gamlss_longitudinal(); mvt_variogram_score()",
      "Protocol scope and absence from mvt_default_comparators()",
      "mvt_distribution_metrics(); mvt_coef_table_one(); mvt_dependence_recovery_row(); mvt_variogram_score(); mvt_summarise_results()",
      "mvt_run_grid(); checkpoint files; mvt_completed_case_ids()",
      "01-run-pilot-grid.R; mvt_write_pilot_feasibility()",
      "08-publication-readiness-audit.R; mvt_publication_readiness_audit()",
      "07-review-audit.R; 09-make-review-bundle.R; mvt_write_artifact_manifest()",
      "No dropout generator or missingness mechanism is defined"
    ),
    proof_after_run = c(
      "scenario_grid.csv; publication_suite_plan.csv; time_design.csv",
      "scenario_grid.csv; family_design.csv; benchmark/coefficient summaries",
      "special_gamlss run directory; publication_readiness optional evidence",
      "scenario_grid.csv generator/dependence rows; truth columns in by-rep outputs",
      "scenario_grid.csv generator/dependence rows; truth columns in by-rep outputs",
      "scenario_grid.csv; dependence_recovery_by_rep.csv",
      "fit_status_by_rep.csv; benchmark_results_by_rep.csv; coefficient_results_by_rep.csv",
      "fit_status_by_rep.csv rows for gamCopula_markov and gamCopula_vine_simplified; dependence_recovery_by_rep.csv; variogram_scores_by_rep.csv; targeted full-vine sensitivity for covariate-dependent dependence",
      "fit_status_by_rep.csv rows for gamlss.longitudinal; dependence_recovery_by_rep.csv; variogram_scores_by_rep.csv",
      "comparators.csv; fit_status_by_rep.csv contains no GJRM method",
      "benchmark_results_by_rep.csv; coefficient_results_by_rep.csv; dependence_recovery_by_rep.csv; variogram_scores_by_rep.csv; runtime_by_rep.csv",
      "*_by_rep.csv files; *_checkpoint.csv files during interrupted runs",
      "pilot_feasibility.md/csv and publication_readiness pilot evidence",
      "publication_readiness_audit.md reports required evidence ready",
      "artifact_manifest.csv/md; review_bundle/README.md",
      "study_protocol.md; scenario_grid.csv has no missingness/dropout fields"
    ),
    status_before_full_run = c(
      rep("implemented_preflight_verifiable", 6),
      rep("implemented_requires_fit_evidence", 5),
      "implemented_preflight_verifiable",
      "implemented_requires_pilot_outputs",
      "implemented_requires_full_suite_outputs",
      "implemented_preflight_verifiable",
      "implemented_preflight_verifiable"
    ),
    stringsAsFactors = FALSE
  )

  list(
    time_design = time_table,
    family_design = family_table,
    dependence_design = dependence_table,
    comparators = comparator_table,
    metrics = metrics_table,
    expected_artifacts = output_table,
    publication_readiness = readiness,
    review_checklist = review_checklist
  )
}

mvt_write_study_protocol <- function(output_dir = mvt_output_root) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- mvt_study_protocol_tables()
  protocol_dir <- file.path(output_dir, "study_protocol")
  dir.create(protocol_dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(tables)) {
    mvt_write_csv(tables[[nm]], file.path(protocol_dir, paste0(nm, ".csv")))
  }

  table_lines <- function(title, dat, cols = names(dat), max_rows = Inf) {
    dat <- dat[seq_len(min(nrow(dat), max_rows)), cols, drop = FALSE]
    c(
      paste0("## ", title),
      "",
      paste0("- ", apply(dat, 1L, function(row) paste(paste(names(row), row, sep = "="), collapse = "; "))),
      ""
    )
  }

  lines <- c(
    "# Multivariate Longitudinal Simulation Study Protocol",
    "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "This protocol is generated from the workflow definitions in `00-multivariate-setup.R`.",
    "",
    "## Scope",
    "",
    "- Extends the JSS bivariate/trivariate workflow to multivariate longitudinal settings.",
    "- Uses total rows as the sample-size cap, with `T = 5, 20, 50` and no missingness/dropout.",
    "- Compares standard mean-regression approaches, `gamCopula`, and `gamlss.longitudinal`.",
    "- Treats GJRM as historical low-dimensional context rather than a high-dimensional comparator.",
    "",
    table_lines("Time Designs", tables$time_design),
    table_lines("Margins", tables$family_design, c("family_name", "label", "gamlss_family", "standard_comparator_available", "special")),
    table_lines("Dependence Scenarios", tables$dependence_design, c("dependence_name", "generator", "dependence", "appendix")),
    table_lines("Comparators", tables$comparators, c("method", "role", "description")),
    table_lines("Metric Families", tables$metrics),
    table_lines("Publication Readiness Evidence", tables$publication_readiness, c("role", "required", "min_reps", "required_time", "required_family", "required_dependence", "required_method"), max_rows = nrow(tables$publication_readiness)),
    table_lines("Requirement Checklist", tables$review_checklist, c("requirement", "implementation_evidence", "proof_after_run", "status_before_full_run"), max_rows = nrow(tables$review_checklist)),
    "## Generated CSV Tables",
    "",
    paste0("- `study_protocol/", names(tables), ".csv`"),
    "",
    "## Review Gate",
    "",
    "- Run `10-run-publication-suite.R` in dry-run mode to inspect planned case counts.",
    "- Run `10-run-publication-suite.R` with `GAMLSS_LONGITUDINAL_MVT_SUITE_PREFLIGHT_ONLY=true` before full compute.",
    "- Run the full suite with strict readiness enabled for paper-facing evidence.",
    "- Treat publication readiness as unproven until `publication_readiness_audit.md` reports required evidence ready."
  )
  writeLines(unlist(lines, use.names = FALSE), file.path(output_dir, "study_protocol.md"), useBytes = TRUE)
  invisible(list(md = file.path(output_dir, "study_protocol.md"), tables = tables))
}

mvt_status_row <- function(check, status, detail = "", evidence = "") {
  data.frame(
    check = check,
    status = status,
    detail = detail,
    evidence = evidence,
    stringsAsFactors = FALSE
  )
}

mvt_write_implementation_status <- function(
    output_dir = mvt_output_root,
    suite_dir = file.path(output_dir, "publication_suite_preflight_default"),
    smoke_dir = file.path(output_dir, "full_method_stack_smoke"),
    allfamily_smoke_dir = file.path(output_dir, "review_smoke_allfamilies"),
    t20_smoke_dir = file.path(output_dir, "t20_full_method_smoke"),
    pilot_shard_dir = file.path(output_dir, "pilot_real_size_rep1_shard"),
    main_shard_dir = file.path(output_dir, "main_core_rep1_shard"),
    appendix_review_dir = file.path(output_dir, "appendix_t5_t50_review"),
    special_review_dir = file.path(output_dir, "special_gamlss_review"),
    glmm_sensitivity_review_dir = file.path(output_dir, "glmm_sensitivity_review")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  script_files <- list.files(mvt_script_dir, pattern = "[.]R$", full.names = TRUE)
  parse_status <- tryCatch({
    for (script in script_files) invisible(parse(script))
    TRUE
  }, error = function(e) e)

  protocol_path <- file.path(output_dir, "study_protocol.md")
  checklist_path <- file.path(output_dir, "study_protocol", "review_checklist.csv")
  suite_plan_path <- file.path(suite_dir, "publication_suite_plan.csv")
  suite_preflight_path <- file.path(suite_dir, "publication_suite_preflight.csv")
  readiness_path <- file.path(output_dir, "publication_readiness_audit.csv")
  smoke_status <- mvt_read_optional_csv(file.path(smoke_dir, "fit_status_by_rep.csv"))
  smoke_audit <- mvt_read_optional_csv(file.path(smoke_dir, "review_audit.csv"))
  allfamily_grid <- mvt_read_optional_csv(file.path(allfamily_smoke_dir, "scenario_grid.csv"))
  allfamily_status <- mvt_read_optional_csv(file.path(allfamily_smoke_dir, "fit_status_by_rep.csv"))
  allfamily_audit <- mvt_read_optional_csv(file.path(allfamily_smoke_dir, "review_audit.csv"))
  t20_grid <- mvt_read_optional_csv(file.path(t20_smoke_dir, "scenario_grid.csv"))
  t20_status <- mvt_read_optional_csv(file.path(t20_smoke_dir, "fit_status_by_rep.csv"))
  t20_audit <- mvt_read_optional_csv(file.path(t20_smoke_dir, "review_audit.csv"))
  pilot_shard_grid <- mvt_read_optional_csv(file.path(pilot_shard_dir, "scenario_grid.csv"))
  pilot_shard_status <- mvt_read_optional_csv(file.path(pilot_shard_dir, "fit_status_by_rep.csv"))
  pilot_shard_audit <- mvt_read_optional_csv(file.path(pilot_shard_dir, "review_audit.csv"))
  main_shard_grid <- mvt_read_optional_csv(file.path(main_shard_dir, "scenario_grid.csv"))
  main_shard_status <- mvt_read_optional_csv(file.path(main_shard_dir, "fit_status_by_rep.csv"))
  main_shard_audit <- mvt_read_optional_csv(file.path(main_shard_dir, "review_audit.csv"))
  appendix_review_grid <- mvt_read_optional_csv(file.path(appendix_review_dir, "scenario_grid.csv"))
  appendix_review_status <- mvt_read_optional_csv(file.path(appendix_review_dir, "fit_status_by_rep.csv"))
  appendix_review_audit <- mvt_read_optional_csv(file.path(appendix_review_dir, "review_audit.csv"))
  special_review_grid <- mvt_read_optional_csv(file.path(special_review_dir, "scenario_grid.csv"))
  special_review_status <- mvt_read_optional_csv(file.path(special_review_dir, "fit_status_by_rep.csv"))
  special_review_audit <- mvt_read_optional_csv(file.path(special_review_dir, "review_audit.csv"))
  glmm_sensitivity_review_grid <- mvt_read_optional_csv(file.path(glmm_sensitivity_review_dir, "scenario_grid.csv"))
  glmm_sensitivity_review_status <- mvt_read_optional_csv(file.path(glmm_sensitivity_review_dir, "fit_status_by_rep.csv"))
  glmm_sensitivity_review_audit <- mvt_read_optional_csv(file.path(glmm_sensitivity_review_dir, "review_audit.csv"))

  suite_plan <- mvt_read_optional_csv(suite_plan_path)
  suite_preflight <- mvt_read_optional_csv(suite_preflight_path)
  readiness <- mvt_read_optional_csv(readiness_path)

  expected_roles <- c("pilot", "main_core", "appendix")
  expected_cases <- c(pilot = 90L, main_core = 1600L, appendix = 4800L)
  expected_rows <- c(pilot = 4000L, main_core = 4000L, appendix = 5000L)

  preflight_roles_ok <- nrow(suite_preflight) > 0L &&
    all(expected_roles %in% suite_preflight$role) &&
    all(suite_preflight$failures[suite_preflight$role %in% expected_roles] == 0L) &&
    all(suite_preflight$warnings[suite_preflight$role %in% expected_roles] == 0L)

  plan_cases_ok <- nrow(suite_plan) > 0L &&
    all(expected_roles %in% suite_plan$role) &&
    all(vapply(expected_roles, function(role) {
      row <- suite_plan[suite_plan$role == role, , drop = FALSE]
      nrow(row) == 1L &&
        isTRUE(row$estimated_cases[[1L]] == expected_cases[[role]]) &&
        isTRUE(row$max_total_rows[[1L]] == expected_rows[[role]])
    }, logical(1)))

  role_artifact_rows <- lapply(expected_roles, function(role) {
    run_dir <- file.path(suite_dir, role)
    required <- c(
      "scenario_grid.csv",
      "preflight_checks.csv",
      "preflight_checks.md",
      "run_metadata.csv",
      "package_versions.csv",
      "session_info.txt",
      "artifact_manifest.csv"
    )
    present <- file.exists(file.path(run_dir, required))
    mvt_status_row(
      paste0("preflight_artifacts_", role),
      if (all(present)) "pass" else "fail",
      paste(required, present, sep = "=", collapse = ", "),
      normalizePath(run_dir, winslash = "/", mustWork = FALSE)
    )
  })

  role_fit_rows <- lapply(expected_roles, function(role) {
    run_dir <- file.path(suite_dir, role)
    fit_outputs <- file.path(run_dir, paste0(mvt_result_names(), "_by_rep.csv"))
    mvt_status_row(
      paste0("preflight_has_no_fit_outputs_", role),
      if (!any(file.exists(fit_outputs))) "pass" else "warn",
      paste(basename(fit_outputs), file.exists(fit_outputs), sep = "=", collapse = ", "),
      normalizePath(run_dir, winslash = "/", mustWork = FALSE)
    )
  })

  readiness_ready <- nrow(readiness) > 0L &&
    "status" %in% names(readiness) &&
    all(readiness$status[readiness$role %in% c("pilot", "main_core")] == "pass")

  smoke_methods_ok <- nrow(smoke_status) > 0L &&
    all(mvt_default_comparators() %in% unique(smoke_status$method)) &&
    all(smoke_status$success[smoke_status$method %in% mvt_default_comparators()])
  smoke_audit_ok <- nrow(smoke_audit) > 0L &&
    "status" %in% names(smoke_audit) &&
    !any(smoke_audit$status == "fail")
  smoke_status_class <- if (isTRUE(smoke_methods_ok) && isTRUE(smoke_audit_ok)) {
    "pass"
  } else if (dir.exists(smoke_dir)) {
    "warn"
  } else {
    "incomplete"
  }
  allfamily_families_ok <- nrow(allfamily_grid) > 0L &&
    "family_name" %in% names(allfamily_grid) &&
    all(c("gaussian", "poisson", "gamma", "binomial") %in% unique(allfamily_grid$family_name))
  allfamily_dependence_ok <- nrow(allfamily_grid) > 0L &&
    "dependence_name" %in% names(allfamily_grid) &&
    all(c("external_exchangeable", "native_covariate_dependent_adjacent") %in% unique(allfamily_grid$dependence_name))
  allfamily_methods_ok <- nrow(allfamily_status) > 0L &&
    "method" %in% names(allfamily_status) &&
    all(c("glm", "gamCopula_markov", "gamlss.longitudinal") %in% unique(allfamily_status$method)) &&
    all(allfamily_status$success[allfamily_status$method %in% c("glm", "gamCopula_markov", "gamlss.longitudinal")])
  allfamily_audit_ok <- nrow(allfamily_audit) > 0L &&
    "status" %in% names(allfamily_audit) &&
    !any(allfamily_audit$status == "fail")
  allfamily_status_class <- if (
    isTRUE(allfamily_families_ok) &&
      isTRUE(allfamily_dependence_ok) &&
      isTRUE(allfamily_methods_ok) &&
      isTRUE(allfamily_audit_ok)
  ) {
    "pass"
  } else if (dir.exists(allfamily_smoke_dir)) {
    "warn"
  } else {
    "incomplete"
  }
  t20_grid_ok <- nrow(t20_grid) > 0L &&
    "n_time" %in% names(t20_grid) &&
    any(t20_grid$n_time == 20L)
  t20_methods_ok <- nrow(t20_status) > 0L &&
    "method" %in% names(t20_status) &&
    all(mvt_default_comparators() %in% unique(t20_status$method)) &&
    all(t20_status$success[t20_status$method %in% mvt_default_comparators()])
  t20_audit_ok <- nrow(t20_audit) > 0L &&
    "status" %in% names(t20_audit) &&
    !any(t20_audit$status == "fail")
  t20_status_class <- if (isTRUE(t20_grid_ok) && isTRUE(t20_methods_ok) && isTRUE(t20_audit_ok)) {
    "pass"
  } else if (dir.exists(t20_smoke_dir)) {
    "warn"
  } else {
    "incomplete"
  }
  pilot_shard_grid_ok <- nrow(pilot_shard_grid) > 0L &&
    "n_time" %in% names(pilot_shard_grid) &&
    "family_name" %in% names(pilot_shard_grid) &&
    "dependence_name" %in% names(pilot_shard_grid) &&
    nrow(pilot_shard_grid) == 18L &&
    all(c(5L, 20L) %in% unique(pilot_shard_grid$n_time)) &&
    all(c("gaussian", "gamma", "binomial") %in% unique(pilot_shard_grid$family_name)) &&
    all(c("external_exchangeable", "external_ar1", "native_covariate_dependent_adjacent") %in% unique(pilot_shard_grid$dependence_name))
  pilot_shard_methods_ok <- nrow(pilot_shard_status) > 0L &&
    "method" %in% names(pilot_shard_status) &&
    "status_class" %in% names(pilot_shard_status) &&
    all(mvt_default_comparators() %in% unique(pilot_shard_status$method)) &&
    all(nzchar(trimws(as.character(pilot_shard_status$status_class))))
  pilot_shard_audit_ok <- nrow(pilot_shard_audit) > 0L &&
    "status" %in% names(pilot_shard_audit) &&
    !any(pilot_shard_audit$status == "fail")
  pilot_shard_status_class <- if (
    isTRUE(pilot_shard_grid_ok) &&
      isTRUE(pilot_shard_methods_ok) &&
      isTRUE(pilot_shard_audit_ok)
  ) {
    "pass"
  } else if (dir.exists(pilot_shard_dir)) {
    "warn"
  } else {
    "incomplete"
  }
  main_shard_grid_ok <- nrow(main_shard_grid) > 0L &&
    "n_time" %in% names(main_shard_grid) &&
    "family_name" %in% names(main_shard_grid) &&
    "dependence_name" %in% names(main_shard_grid) &&
    nrow(main_shard_grid) == 16L &&
    all(main_shard_grid$n_time == 20L) &&
    all(c("gaussian", "poisson", "gamma", "binomial") %in% unique(main_shard_grid$family_name)) &&
    all(c("external_exchangeable", "external_ar1", "native_time_varying_adjacent", "native_covariate_dependent_adjacent") %in% unique(main_shard_grid$dependence_name))
  main_shard_methods_ok <- nrow(main_shard_status) > 0L &&
    "method" %in% names(main_shard_status) &&
    "status_class" %in% names(main_shard_status) &&
    all(mvt_default_comparators() %in% unique(main_shard_status$method)) &&
    all(nzchar(trimws(as.character(main_shard_status$status_class))))
  main_shard_audit_ok <- nrow(main_shard_audit) > 0L &&
    "status" %in% names(main_shard_audit) &&
    !any(main_shard_audit$status == "fail")
  main_shard_status_class <- if (
    isTRUE(main_shard_grid_ok) &&
      isTRUE(main_shard_methods_ok) &&
      isTRUE(main_shard_audit_ok)
  ) {
    "pass"
  } else if (dir.exists(main_shard_dir)) {
    "warn"
  } else {
    "incomplete"
  }
  appendix_review_ok <- nrow(appendix_review_grid) > 0L &&
    "n_time" %in% names(appendix_review_grid) &&
    "family_name" %in% names(appendix_review_grid) &&
    "dependence_name" %in% names(appendix_review_grid) &&
    all(c(5L, 50L) %in% unique(appendix_review_grid$n_time)) &&
    all(c("gaussian", "poisson", "gamma", "binomial") %in% unique(appendix_review_grid$family_name)) &&
    all(c("external_exchangeable", "external_ar1", "native_time_varying_adjacent", "native_covariate_dependent_adjacent") %in% unique(appendix_review_grid$dependence_name)) &&
    nrow(appendix_review_status) > 0L &&
    all(mvt_default_comparators() %in% unique(appendix_review_status$method)) &&
    nrow(appendix_review_audit) > 0L &&
    !any(appendix_review_audit$status == "fail")
  special_review_ok <- nrow(special_review_grid) > 0L &&
    "family_name" %in% names(special_review_grid) &&
    "method" %in% names(special_review_status) &&
    "gg_continuous" %in% unique(special_review_grid$family_name) &&
    "gamlss.longitudinal" %in% unique(special_review_status$method) &&
    nrow(special_review_audit) > 0L &&
    !any(special_review_audit$status == "fail")
  glmm_sensitivity_review_ok <- nrow(glmm_sensitivity_review_grid) > 0L &&
    "n_time" %in% names(glmm_sensitivity_review_grid) &&
    "method" %in% names(glmm_sensitivity_review_status) &&
    all(c(5L, 20L, 50L) %in% unique(glmm_sensitivity_review_grid$n_time)) &&
    "glmm_slope" %in% unique(glmm_sensitivity_review_status$method) &&
    nrow(glmm_sensitivity_review_audit) > 0L &&
    !any(glmm_sensitivity_review_audit$status == "fail")

  checks <- mvt_bind_rows_fill(
    list(
      mvt_status_row(
        "workflow_scripts_parse",
        if (isTRUE(parse_status)) "pass" else "fail",
        if (isTRUE(parse_status)) paste("scripts", length(script_files)) else conditionMessage(parse_status),
        normalizePath(mvt_script_dir, winslash = "/", mustWork = FALSE)
      ),
      mvt_status_row(
        "study_protocol_generated",
        if (file.exists(protocol_path) && file.exists(checklist_path)) "pass" else "fail",
        paste("study_protocol.md", file.exists(protocol_path), "| review_checklist.csv", file.exists(checklist_path)),
        normalizePath(output_dir, winslash = "/", mustWork = FALSE)
      ),
      mvt_status_row(
        "default_suite_preflight_clean",
        if (preflight_roles_ok) "pass" else "fail",
        if (nrow(suite_preflight) > 0L) {
          paste0(suite_preflight$role, ": failures=", suite_preflight$failures, ", warnings=", suite_preflight$warnings, collapse = "; ")
        } else {
          "publication_suite_preflight.csv missing"
        },
        normalizePath(suite_preflight_path, winslash = "/", mustWork = FALSE)
      ),
      mvt_status_row(
        "default_suite_plan_matches_protocol",
        if (plan_cases_ok) "pass" else "fail",
        if (nrow(suite_plan) > 0L) {
          paste0(suite_plan$role, ": cases=", suite_plan$estimated_cases, ", max_rows=", suite_plan$max_total_rows, collapse = "; ")
        } else {
          "publication_suite_plan.csv missing"
        },
        normalizePath(suite_plan_path, winslash = "/", mustWork = FALSE)
      )
    ),
    role_artifact_rows,
    role_fit_rows,
    list(
      mvt_status_row(
        "full_method_stack_smoke",
        smoke_status_class,
        if (nrow(smoke_status) > 0L) {
          paste0(smoke_status$method, "=", smoke_status$status_class, collapse = "; ")
        } else {
          "fit_status_by_rep.csv missing"
        },
        normalizePath(smoke_dir, winslash = "/", mustWork = FALSE)
      ),
      mvt_status_row(
        "all_family_review_smoke",
        allfamily_status_class,
        if (nrow(allfamily_grid) > 0L && nrow(allfamily_status) > 0L) {
          paste(
            "families", paste(sort(unique(allfamily_grid$family_name)), collapse = ","),
            "| dependence", paste(sort(unique(allfamily_grid$dependence_name)), collapse = ","),
            "| methods", paste(sort(unique(allfamily_status$method)), collapse = ",")
          )
        } else {
          "scenario_grid.csv or fit_status_by_rep.csv missing"
        },
        normalizePath(allfamily_smoke_dir, winslash = "/", mustWork = FALSE)
      ),
      mvt_status_row(
        "t20_full_method_smoke",
        t20_status_class,
        if (nrow(t20_grid) > 0L && nrow(t20_status) > 0L) {
          paste(
            "times", paste(sort(unique(t20_grid$n_time)), collapse = ","),
            "| methods", paste(sort(unique(t20_status$method)), collapse = ",")
          )
        } else {
          "scenario_grid.csv or fit_status_by_rep.csv missing"
        },
        normalizePath(t20_smoke_dir, winslash = "/", mustWork = FALSE)
      ),
      mvt_status_row(
        "pilot_real_size_rep1_shard",
        pilot_shard_status_class,
        if (nrow(pilot_shard_grid) > 0L && nrow(pilot_shard_status) > 0L) {
          paste(
            "cases", nrow(pilot_shard_grid),
            "| times", paste(sort(unique(pilot_shard_grid$n_time)), collapse = ","),
            "| statuses", paste(names(table(pilot_shard_status$status_class)), table(pilot_shard_status$status_class), sep = "=", collapse = ",")
          )
        } else {
          "scenario_grid.csv or fit_status_by_rep.csv missing"
        },
        normalizePath(pilot_shard_dir, winslash = "/", mustWork = FALSE)
      ),
      mvt_status_row(
        "main_core_rep1_shard",
        main_shard_status_class,
        if (nrow(main_shard_grid) > 0L && nrow(main_shard_status) > 0L) {
          paste(
            "cases", nrow(main_shard_grid),
            "| families", paste(sort(unique(main_shard_grid$family_name)), collapse = ","),
            "| statuses", paste(names(table(main_shard_status$status_class)), table(main_shard_status$status_class), sep = "=", collapse = ",")
          )
        } else {
          "scenario_grid.csv or fit_status_by_rep.csv missing"
        },
        normalizePath(main_shard_dir, winslash = "/", mustWork = FALSE)
      ),
      mvt_status_row(
        "appendix_t5_t50_review_evidence",
        if (isTRUE(appendix_review_ok)) "pass" else if (dir.exists(appendix_review_dir)) "warn" else "incomplete",
        if (nrow(appendix_review_grid) > 0L && nrow(appendix_review_status) > 0L) {
          paste(
            "cases", nrow(appendix_review_grid),
            "| times", paste(sort(unique(appendix_review_grid$n_time)), collapse = ","),
            "| statuses", paste(names(table(appendix_review_status$status_class)), table(appendix_review_status$status_class), sep = "=", collapse = ",")
          )
        } else {
          "appendix review scenario_grid.csv or fit_status_by_rep.csv missing"
        },
        normalizePath(appendix_review_dir, winslash = "/", mustWork = FALSE)
      ),
      mvt_status_row(
        "special_gamlss_review_evidence",
        if (isTRUE(special_review_ok)) "pass" else if (dir.exists(special_review_dir)) "warn" else "incomplete",
        if (nrow(special_review_grid) > 0L && nrow(special_review_status) > 0L) {
          paste(
            "cases", nrow(special_review_grid),
            "| families", paste(sort(unique(special_review_grid$family_name)), collapse = ","),
            "| statuses", paste(names(table(special_review_status$status_class)), table(special_review_status$status_class), sep = "=", collapse = ",")
          )
        } else {
          "special review scenario_grid.csv or fit_status_by_rep.csv missing"
        },
        normalizePath(special_review_dir, winslash = "/", mustWork = FALSE)
      ),
      mvt_status_row(
        "glmm_sensitivity_review_evidence",
        if (isTRUE(glmm_sensitivity_review_ok)) "pass" else if (dir.exists(glmm_sensitivity_review_dir)) "warn" else "incomplete",
        if (nrow(glmm_sensitivity_review_grid) > 0L && nrow(glmm_sensitivity_review_status) > 0L) {
          paste(
            "cases", nrow(glmm_sensitivity_review_grid),
            "| times", paste(sort(unique(glmm_sensitivity_review_grid$n_time)), collapse = ","),
            "| statuses", paste(names(table(glmm_sensitivity_review_status$status_class)), table(glmm_sensitivity_review_status$status_class), sep = "=", collapse = ",")
          )
        } else {
          "GLMM sensitivity scenario_grid.csv or fit_status_by_rep.csv missing"
        },
        normalizePath(glmm_sensitivity_review_dir, winslash = "/", mustWork = FALSE)
      ),
      mvt_status_row(
        "full_fit_evidence",
        if (isTRUE(readiness_ready)) "pass" else "incomplete",
        if (nrow(readiness) > 0L) {
          paste0(readiness$role, "/", readiness$check, "=", readiness$status, collapse = "; ")
        } else {
          "publication_readiness_audit.csv missing or no full run directories configured"
        },
        normalizePath(readiness_path, winslash = "/", mustWork = FALSE)
      ),
      mvt_status_row(
        "publication_review_ready",
        if (isTRUE(readiness_ready)) "pass" else "incomplete",
        if (isTRUE(readiness_ready)) "required pilot and main evidence passed readiness audit" else "requires full pilot and main model-fitting outputs",
        "publication_readiness_audit.md"
      )
    )
  )

  csv_path <- file.path(output_dir, "implementation_status.csv")
  md_path <- file.path(output_dir, "implementation_status.md")
  mvt_write_csv(checks, csv_path)
  non_pass <- checks[checks$status != "pass", , drop = FALSE]
  lines <- c(
    "# Multivariate Simulation Implementation Status",
    "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Source, preflight, and smoke evidence ready:", all(checks$status[!checks$check %in% c("full_fit_evidence", "publication_review_ready")] == "pass")),
    paste("Publication review ready:", isTRUE(readiness_ready)),
    paste("Checks passing:", sum(checks$status == "pass"), "of", nrow(checks)),
    paste("Incomplete:", sum(checks$status == "incomplete")),
    paste("Warnings:", sum(checks$status == "warn")),
    paste("Failures:", sum(checks$status == "fail")),
    "",
    "## Non-Passing Checks",
    ""
  )
  if (nrow(non_pass) == 0L) {
    lines <- c(lines, "None.")
  } else {
    lines <- c(lines, paste0("- ", non_pass$check, " [", non_pass$status, "]: ", non_pass$detail))
  }
  lines <- c(
    lines,
    "",
    "## Evidence Files",
    "",
    "- `study_protocol.md`",
    "- `study_protocol/review_checklist.csv`",
    "- `publication_suite_preflight_default/publication_suite_plan.csv`",
    "- `publication_suite_preflight_default/publication_suite_preflight.csv`",
    "- `full_method_stack_smoke/review_audit.md`",
    "- `review_smoke_allfamilies/review_audit.md`",
    "- `t20_full_method_smoke/review_audit.md`",
    "- `pilot_real_size_rep1_shard/review_audit.md`",
    "- `main_core_rep1_shard/review_audit.md`",
    "- `appendix_t5_t50_review/review_audit.md`",
    "- `special_gamlss_review/review_audit.md`",
    "- `glmm_sensitivity_review/review_audit.md`",
    "- `publication_readiness_audit.md`"
  )
  writeLines(lines, md_path, useBytes = TRUE)
  invisible(list(csv = csv_path, md = md_path, checks = checks, ready = isTRUE(readiness_ready)))
}

mvt_bind_rows_fill <- function(...) {
  flatten_piece <- function(x) {
    if (is.data.frame(x)) return(list(x))
    if (is.list(x)) return(unlist(lapply(x, flatten_piece), recursive = FALSE))
    list()
  }
  pieces <- unlist(lapply(list(...), flatten_piece), recursive = FALSE)
  pieces <- pieces[vapply(pieces, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1))]
  if (length(pieces) == 0L) return(data.frame())
  cols <- unique(unlist(lapply(pieces, names), use.names = FALSE))
  pieces <- lapply(pieces, function(x) {
    missing <- setdiff(cols, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, pieces)
}

mvt_select_named <- function(x, keep) {
  if (length(keep) == 0L) return(x)
  missing <- setdiff(keep, names(x))
  if (length(missing) > 0L) {
    stop("Unknown name(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  x[keep]
}

mvt_require_namespaces <- function(packages, strict = TRUE) {
  available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  if (isTRUE(strict) && any(!available)) {
    stop(
      "Missing required package(s): ",
      paste(packages[!available], collapse = ", "),
      call. = FALSE
    )
  }
  data.frame(package = packages, available = available, stringsAsFactors = FALSE)
}

mvt_require_runtime_packages <- function(include_gamcopula = TRUE) {
  packages <- c("gamlss", "gamlss.dist", "mvtnorm", "VineCopula", "callr")
  if (isTRUE(include_gamcopula)) packages <- c(packages, "gamCopula")
  mvt_require_namespaces(packages, strict = TRUE)
}

mvt_package_versions <- function(packages) {
  rows <- lapply(packages, function(pkg) {
    data.frame(
      package = pkg,
      version = if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(utils::packageVersion(pkg))
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mvt_load_package <- function() {
  source <- match.arg(
    mvt_env("GAMLSS_LONGITUDINAL_MVT_SOURCE", "installed"),
    c("installed", "local")
  )
  if (identical(source, "local")) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("Package 'pkgload' is required for local source loading.", call. = FALSE)
    }
    pkgload::load_all(mvt_repo_root, quiet = TRUE)
  } else {
    suppressPackageStartupMessages(library(gamlss.longitudinal))
  }
  invisible(source)
}

mvt_time_specs <- function() {
  list(
    t5 = list(n_time = 5L, n_subject = 400L, total_rows = 2000L),
    t20 = list(n_time = 20L, n_subject = 200L, total_rows = 4000L),
    t50 = list(n_time = 50L, n_subject = 100L, total_rows = 5000L)
  )
}

mvt_family_specs <- function(include_special = TRUE) {
  specs <- list(
    gaussian = list(
      family = "gaussian",
      gamlss_family = "NO",
      label = "Gaussian",
      margin_dist = function() gamlss.dist::NO(),
      standard_family = stats::gaussian(),
      linkinv = identity,
      eta = list(intercept = 0.1, time = 0.35, x = 0.55, z = -0.25),
      sigma = 1
    ),
    poisson = list(
      family = "poisson",
      gamlss_family = "PO",
      label = "Poisson",
      margin_dist = function() gamlss.dist::PO(mu.link = "log"),
      standard_family = stats::poisson(),
      linkinv = exp,
      eta = list(intercept = log(2.4), time = 0.25, x = 0.30, z = -0.18),
      sigma = NA_real_
    ),
    gamma = list(
      family = "gamma",
      gamlss_family = "GA",
      label = "Gamma",
      margin_dist = function() gamlss.dist::GA(mu.link = "log", sigma.link = "log"),
      standard_family = stats::Gamma(link = "log"),
      linkinv = exp,
      eta = list(intercept = log(2), time = 0.20, x = 0.38, z = -0.22),
      sigma = 0.45
    ),
    binomial = list(
      family = "binomial",
      gamlss_family = "BI",
      label = "Binomial",
      margin_dist = function() gamlss.dist::BI(mu.link = "logit"),
      standard_family = stats::binomial(),
      linkinv = stats::plogis,
      eta = list(intercept = stats::qlogis(0.35), time = 0.25, x = 0.65, z = -0.30),
      sigma = NA_real_
    )
  )

  if (isTRUE(include_special)) {
    specs$gg_continuous <- list(
      family = "gg_continuous",
      gamlss_family = "GG",
      label = "Generalized gamma",
      margin_dist = function() gamlss.dist::GG(mu.link = "log", sigma.link = "log", nu.link = "identity"),
      standard_family = NULL,
      linkinv = exp,
      eta = list(intercept = log(2), time = 0.25, x = 0.35, z = -0.20),
      sigma = function(data) exp(-0.65 + 0.25 * data$time + 0.20 * data$z),
      nu = function(data) 0.35 + 0.30 * data$x,
      special = TRUE
    )
  }

  specs
}

mvt_dependence_specs <- function(include_appendix = TRUE) {
  specs <- list(
    external_exchangeable = list(
      scenario = "external_exchangeable",
      generator = "external",
      dependence = "exchangeable",
      correlation_level = "moderate",
      rho = 0.45,
      theta_formula = ~1,
      gee_correlations = c("independence", "exchangeable", "ar1", "unstructured"),
      appendix = FALSE
    ),
    external_ar1 = list(
      scenario = "external_ar1",
      generator = "external",
      dependence = "ar1",
      correlation_level = "moderate_high",
      rho = 0.65,
      theta_formula = ~time,
      gee_correlations = c("independence", "ar1", "exchangeable", "unstructured"),
      appendix = FALSE
    ),
    native_time_varying_adjacent = list(
      scenario = "native_time_varying_adjacent",
      generator = "native",
      dependence = "time_varying_adjacent",
      correlation_level = "moderate_high",
      tau_edges = c(0.25, 0.65),
      theta_formula = ~time,
      gee_correlations = c("independence", "ar1", "exchangeable", "unstructured"),
      appendix = FALSE
    ),
    native_covariate_dependent_adjacent = list(
      scenario = "native_covariate_dependent_adjacent",
      generator = "native",
      dependence = "covariate_dependent_adjacent",
      correlation_level = "moderate_high",
      tau_base = 0.45,
      tau_effect = 0.90,
      theta_formula = ~x,
      gee_correlations = c("independence", "exchangeable", "ar1", "unstructured"),
      appendix = FALSE
    )
  )

  if (isTRUE(include_appendix)) {
    specs$external_block <- list(
      scenario = "external_block",
      generator = "external",
      dependence = "block",
      correlation_level = "mixed",
      rho_within = 0.65,
      rho_between = 0.20,
      theta_formula = ~time,
      gee_correlations = c("independence", "exchangeable", "ar1", "unstructured"),
      appendix = TRUE
    )
  }

  specs
}

mvt_main_scope <- function(value = mvt_env("GAMLSS_LONGITUDINAL_MVT_MAIN_SCOPE", "core")) {
  match.arg(value, c("core", "appendix", "all"))
}

mvt_default_main_time_names <- function(scope = mvt_main_scope()) {
  if (identical(scope, "core")) "t20" else c("t5", "t20", "t50")
}

mvt_resize_dependence <- function(dep, n_time) {
  dep$n_time <- as.integer(n_time)
  if (identical(dep$dependence, "time_varying_adjacent")) {
    rng <- dep$tau_edges
    dep$tau_edges <- seq(rng[[1L]], rng[[2L]], length.out = max(1L, n_time - 1L))
  }
  dep
}

mvt_expand_grid <- function(
    time_names = names(mvt_time_specs()),
    family_names = c("gaussian", "poisson", "gamma", "binomial"),
    dependence_names = names(mvt_dependence_specs(include_appendix = FALSE)),
    reps = 1L,
    include_special = FALSE,
    include_appendix = FALSE) {
  time_specs <- mvt_select_named(mvt_time_specs(), time_names)
  family_specs <- mvt_select_named(mvt_family_specs(include_special = include_special), family_names)
  dep_specs <- mvt_select_named(mvt_dependence_specs(include_appendix = include_appendix), dependence_names)
  reps <- if (length(reps) == 1L) seq_len(as.integer(reps)) else as.integer(reps)
  reps <- sort(unique(reps[is.finite(reps) & reps > 0L]))
  if (length(reps) == 0L) stop("At least one positive replicate id is required.", call. = FALSE)

  rows <- list()
  idx <- 1L
  for (time_name in names(time_specs)) {
    t_spec <- time_specs[[time_name]]
    for (dep_name in names(dep_specs)) {
      dep <- mvt_resize_dependence(dep_specs[[dep_name]], t_spec$n_time)
      for (family_name in names(family_specs)) {
        spec <- family_specs[[family_name]]
        for (rep_id in reps) {
          rows[[idx]] <- data.frame(
            case_id = paste(time_name, dep_name, family_name, sprintf("rep%03d", rep_id), sep = "__"),
            time_name = time_name,
            family_name = family_name,
            dependence_name = dep_name,
            rep = rep_id,
            n_time = t_spec$n_time,
            n_subject = t_spec$n_subject,
            total_rows = t_spec$total_rows,
            generator = dep$generator,
            scenario = paste(dep$scenario, time_name, sep = "_"),
            dependence = dep$dependence,
            correlation_level = dep$correlation_level,
            appendix = isTRUE(dep$appendix) || isTRUE(spec$special),
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
    }
  }
  if (length(rows) == 0L) return(data.frame())
  mvt_apply_size_overrides(do.call(rbind, rows))
}

mvt_apply_size_overrides <- function(grid) {
  if (nrow(grid) == 0L) return(grid)
  n_subject <- mvt_env_int("GAMLSS_LONGITUDINAL_MVT_N_SUBJECT", NA_integer_)
  if (is.finite(n_subject) && n_subject > 0L) {
    grid$n_subject <- as.integer(n_subject)
    grid$total_rows <- grid$n_subject * grid$n_time
  }
  grid
}

mvt_design_long <- function(n_subject, n_time) {
  subject <- rep(seq_len(n_subject), each = n_time)
  time_index <- rep(seq_len(n_time), times = n_subject)
  time <- if (n_time == 1L) rep(0, n_subject) else (time_index - 1) / (n_time - 1)
  subject_scaled <- gamlss.longitudinal::sim_rescale01(subject) - 0.5
  z <- as.numeric((subject %% 2L) == 0L)
  data.frame(
    subject = factor(subject),
    time_index = time_index,
    time = time,
    x = subject_scaled,
    z = z,
    stringsAsFactors = FALSE
  )
}

mvt_eta_mu <- function(data, spec) {
  spec$eta$intercept + spec$eta$time * data$time + spec$eta$x * data$x + spec$eta$z * data$z
}

mvt_margin_param_list <- function(spec) {
  params <- list(mu = function(data) spec$linkinv(mvt_eta_mu(data, spec)))
  if (identical(spec$gamlss_family, "NO") || identical(spec$gamlss_family, "GA")) {
    params$sigma <- spec$sigma
  } else if (identical(spec$gamlss_family, "GG")) {
    params$sigma <- spec$sigma
    params$nu <- spec$nu
  }
  params
}

mvt_eval_margin_param <- function(value, data) {
  if (is.function(value)) return(as.numeric(value(data)))
  rep(as.numeric(value)[[1L]], nrow(data))
}

mvt_add_truth_columns <- function(data, spec) {
  data$true_eta_mu <- mvt_eta_mu(data, spec)
  data$true_mu <- spec$linkinv(data$true_eta_mu)
  data$true_beta_intercept <- spec$eta$intercept
  data$true_beta_time <- spec$eta$time
  data$true_beta_x <- spec$eta$x
  data$true_beta_z <- spec$eta$z
  if (!is.null(spec$sigma) && (is.function(spec$sigma) || !is.na(spec$sigma)[[1L]])) {
    data$true_sigma <- mvt_eval_margin_param(spec$sigma, data)
  }
  if (!is.null(spec$nu)) {
    data$true_nu <- mvt_eval_margin_param(spec$nu, data)
  }
  data
}

mvt_repair_correlation <- function(R, label = "correlation matrix") {
  R <- as.matrix(R)
  R <- (R + t(R)) / 2
  diag(R) <- 1
  eig <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) > 1e-8) return(R)
  if (!mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_REPAIR_CORRELATION", TRUE)) {
    stop(label, " is not positive definite.", call. = FALSE)
  }
  jitter <- abs(min(eig)) + 1e-6
  R <- R + diag(jitter, nrow(R))
  D <- sqrt(diag(R))
  R <- R / outer(D, D)
  diag(R) <- 1
  R
}

mvt_external_correlation_matrix <- function(n_time, dep) {
  idx <- seq_len(n_time)
  if (identical(dep$dependence, "exchangeable")) {
    R <- matrix(dep$rho, n_time, n_time)
    diag(R) <- 1
    return(mvt_repair_correlation(R, dep$scenario))
  }
  if (identical(dep$dependence, "ar1")) {
    R <- dep$rho ^ abs(outer(idx, idx, "-"))
    return(mvt_repair_correlation(R, dep$scenario))
  }
  if (identical(dep$dependence, "block")) {
    block <- idx <= ceiling(n_time / 2)
    R <- matrix(dep$rho_between, n_time, n_time)
    R[outer(block, block, "==")] <- dep$rho_within
    diag(R) <- 1
    return(mvt_repair_correlation(R, dep$scenario))
  }
  stop("Unsupported external dependence: ", dep$dependence, call. = FALSE)
}

mvt_quantile_response <- function(u, data, spec) {
  qfun <- get(paste0("q", spec$gamlss_family), envir = asNamespace("gamlss.dist"))
  args <- list(p = u, mu = data$true_mu)
  if ("true_sigma" %in% names(data)) args$sigma <- data$true_sigma
  if ("true_nu" %in% names(data)) args$nu <- data$true_nu
  if (identical(spec$gamlss_family, "BI")) args$bd <- rep(1, nrow(data))
  args <- args[names(args) %in% formalArgs(qfun)]
  as.numeric(do.call(qfun, args))
}

mvt_simulate_external <- function(row, spec, dep, seed) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' is required for external simulation.", call. = FALSE)
  }
  set.seed(seed)
  dat <- mvt_design_long(row$n_subject, row$n_time)
  dat <- mvt_add_truth_columns(dat, spec)
  R <- mvt_external_correlation_matrix(row$n_time, dep)
  z <- mvtnorm::rmvnorm(row$n_subject, sigma = R)
  u <- pmin(pmax(as.numeric(t(stats::pnorm(z))), 1e-8), 1 - 1e-8)
  dat$response <- mvt_quantile_response(u, dat, spec)
  dat$u <- u
  dat$true_theta <- NA_real_
  dat$true_zeta <- 0
  edge <- dat$time_index > 1L
  dat$true_theta[edge] <- R[cbind(dat$time_index[edge] - 1L, dat$time_index[edge])]
  dat$true_tau <- (2 / pi) * asin(pmax(-1, pmin(1, dat$true_theta)))
  dat <- dat[c("subject", "time_index", "time", "response", setdiff(names(dat), c("subject", "time_index", "time", "response")))]
  attr(dat, "true_correlation_matrix") <- R
  dat
}

mvt_native_copula_params <- function(dep) {
  if (identical(dep$dependence, "time_varying_adjacent")) {
    return(list(tau = dep$tau_edges))
  }
  if (identical(dep$dependence, "covariate_dependent_adjacent")) {
    return(list(tau = function(edge_data) {
      stats::plogis(stats::qlogis(dep$tau_base) + dep$tau_effect * edge_data$x)
    }))
  }
  stop("Unsupported native dependence: ", dep$dependence, call. = FALSE)
}

mvt_simulate_native <- function(row, spec, dep, seed) {
  dat <- gamlss.longitudinal::simulate_longitudinal_dataset(
    n = row$n_subject,
    times = seq_len(row$n_time),
    margin_dist = spec$margin_dist(),
    copula_dist = "N",
    margin_params = mvt_margin_param_list(spec),
    copula_params = mvt_native_copula_params(dep),
    covariates = function(base) {
      out <- mvt_design_long(row$n_subject, row$n_time)
      out[c("time", "x", "z")]
    },
    seed = seed,
    subject_var = "subject",
    time_var = "time_index",
    response_var = "response",
    include_truth = TRUE,
    u_bounds = if (spec$gamlss_family %in% c("PO", "BI")) c(1e-8, 1 - 1e-8) else NULL
  )
  dat$time <- if (row$n_time == 1L) 0 else (dat$time_index - 1) / (row$n_time - 1)
  dat$subject <- factor(dat$subject)
  dat <- mvt_add_truth_columns(dat, spec)
  dat$true_tau <- (2 / pi) * asin(pmax(-1, pmin(1, dat$true_theta)))
  dat[c("subject", "time_index", "time", "response", setdiff(names(dat), c("subject", "time_index", "time", "response")))]
}

mvt_simulate_case <- function(row, seed = NULL) {
  families <- mvt_family_specs(include_special = TRUE)
  deps <- mvt_dependence_specs(include_appendix = TRUE)
  spec <- families[[row$family_name]]
  dep <- mvt_resize_dependence(deps[[row$dependence_name]], row$n_time)
  seed <- seed %||% (20260818L + as.integer(row$rep))
  if (identical(dep$generator, "external")) {
    dat <- mvt_simulate_external(row, spec, dep, seed)
  } else {
    dat <- mvt_simulate_native(row, spec, dep, seed)
  }
  attr(dat, "mvt_spec") <- spec
  attr(dat, "mvt_dependence") <- dep
  dat
}

mvt_capture <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

mvt_elapsed_capture <- function(expr) {
  start <- Sys.time()
  out <- mvt_capture(expr)
  out$elapsed_sec <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  out
}

mvt_elapsed_do_call <- function(fun, args = list()) {
  start <- Sys.time()
  out <- mvt_capture(do.call(fun, args))
  out$elapsed_sec <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  out
}

mvt_run_fit_with_timeout <- function(fun_name, args = list(), timeout = Inf) {
  start <- Sys.time()
  if (is.finite(timeout) && timeout > 0 && requireNamespace("callr", quietly = TRUE)) {
    value <- tryCatch(
      callr::r(
        function(fun_name, args, repo_root) {
          setwd(repo_root)
          source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))
          mvt_load_package()
          mvt_elapsed_do_call(get(fun_name, envir = globalenv()), args)
        },
        args = list(fun_name = fun_name, args = args, repo_root = mvt_repo_root),
        timeout = timeout
      ),
      error = function(e) e
    )
    if (inherits(value, "error")) {
      return(list(
        value = value,
        warnings = character(),
        elapsed_sec = as.numeric(difftime(Sys.time(), start, units = "secs"))
      ))
    }
    return(value)
  }
  mvt_elapsed_do_call(get(fun_name, envir = parent.frame()), args)
}

mvt_score_predictions <- function(y, fitted) {
  ok <- is.finite(y) & is.finite(fitted)
  if (!any(ok)) return(c(mae = NA_real_, rmse = NA_real_, true_mean_rmse = NA_real_))
  err <- fitted[ok] - y[ok]
  c(mae = mean(abs(err)), rmse = sqrt(mean(err^2)))
}

mvt_empty_distribution_metrics <- function() {
  c(
    benchmark_mean_bias = NA_real_,
    benchmark_mean_mae = NA_real_,
    benchmark_mean_rmse = NA_real_,
    benchmark_neg_log_score = NA_real_,
    benchmark_pit_ks_p_value = NA_real_,
    benchmark_pit_mean_abs_error = NA_real_,
    benchmark_q90_mae = NA_real_,
    benchmark_upper_tail_error_90 = NA_real_,
    benchmark_interval_coverage_95 = NA_real_,
    benchmark_interval_width_95 = NA_real_,
    benchmark_tail_error_lower_05 = NA_real_,
    benchmark_tail_error_upper_05 = NA_real_
  )
}

mvt_fit_status_class <- function(success, error = NA_character_, warnings = character()) {
  if (isTRUE(success)) {
    warnings <- as.character(warnings)
    warnings <- warnings[!is.na(warnings) & nzchar(trimws(warnings))]
    warn <- paste(warnings, collapse = " | ")
    if (nzchar(warn)) return("warning")
    return("ok")
  }
  error <- paste(error %||% "", collapse = " | ")
  if (grepl("timeout|timed out|elapsed|time limit|reached maximum", error, ignore.case = TRUE)) {
    return("timeout")
  }
  "error"
}

mvt_failure_reason_short <- function(success, status_class, error = NA_character_, warning = NA_character_) {
  if (mvt_status_success(success) && identical(as.character(status_class), "ok")) return("none")
  text <- if (!is.na(error) && nzchar(trimws(as.character(error)))) {
    as.character(error)
  } else if (!is.na(warning) && nzchar(trimws(as.character(warning)))) {
    as.character(warning)
  } else {
    as.character(status_class %||% "unknown")
  }
  text <- gsub("[\r\n\t]+", " ", text)
  text <- gsub("\\s+", " ", trimws(text))
  if (grepl("timeout|timed out|elapsed|time limit", text, ignore.case = TRUE)) return("timeout")
  if (grepl("Hessian|singular|non-positive|positive definite", text, ignore.case = TRUE)) return("vcov_or_hessian_warning")
  if (grepl("converg|iteration|maximum", text, ignore.case = TRUE)) return("convergence_warning")
  if (grepl("package|namespace|there is no package", text, ignore.case = TRUE)) return("missing_dependency")
  if (!nzchar(text)) return("unknown")
  substr(text, 1L, 120L)
}

mvt_scalar_extra <- function(value) {
  if (length(value) == 0L) return(NA)
  if (length(value) == 1L) return(value)
  value[[1L]]
}

mvt_loglik_value <- function(fit) {
  value <- tryCatch(suppressWarnings(stats::logLik(fit)), error = function(e) NA_real_)
  value <- as.numeric(value)
  value_names <- names(value)
  if (!is.null(value_names) && "joint" %in% value_names) return(value[["joint"]])
  if (length(value) > 1L) return(value[[length(value)]])
  if (length(value) == 0L) NA_real_ else value[[1L]]
}

mvt_loglik_df <- function(fit) {
  value <- tryCatch(suppressWarnings(stats::logLik(fit)), error = function(e) NULL)
  df <- attr(value, "df")
  if (length(df) > 0L && is.finite(df[[1L]])) return(as.numeric(df[[1L]]))
  coefs <- tryCatch(stats::coef(fit), error = function(e) numeric())
  sum(is.finite(as.numeric(coefs)))
}

mvt_loglik_nobs <- function(fit) {
  value <- tryCatch(suppressWarnings(stats::logLik(fit)), error = function(e) NULL)
  n <- attr(value, "nobs")
  if (length(n) > 0L && is.finite(n[[1L]])) return(as.numeric(n[[1L]]))
  n <- tryCatch(stats::nobs(fit), error = function(e) NA_real_)
  n <- suppressWarnings(as.numeric(n))
  if (length(n) > 0L && is.finite(n[[1L]])) n[[1L]] else NA_real_
}

mvt_aic_value <- function(fit) {
  value <- tryCatch(suppressWarnings(stats::AIC(fit)), error = function(e) NA_real_)
  value <- as.numeric(value)
  if (length(value) > 0L && is.finite(value[[1L]])) return(value[[1L]])
  ll <- mvt_loglik_value(fit)
  df <- mvt_loglik_df(fit)
  if (is.finite(ll) && is.finite(df)) -2 * ll + 2 * df else NA_real_
}

mvt_bic_value <- function(fit) {
  value <- tryCatch(suppressWarnings(stats::BIC(fit)), error = function(e) NA_real_)
  value <- as.numeric(value)
  if (length(value) > 0L && is.finite(value[[1L]])) return(value[[1L]])
  ll <- mvt_loglik_value(fit)
  df <- mvt_loglik_df(fit)
  n <- mvt_loglik_nobs(fit)
  if (is.finite(ll) && is.finite(df) && is.finite(n) && n > 1) -2 * ll + log(n) * df else NA_real_
}

mvt_distribution_args <- function(y, params, spec, data = NULL, include_q = TRUE) {
  args <- list()
  if (!is.null(y)) args$x <- y
  args$mu <- params$mu
  if (!is.null(params$sigma)) args$sigma <- params$sigma
  if (!is.null(params$nu)) args$nu <- params$nu
  if (!is.null(params$tau)) args$tau <- params$tau
  if (identical(spec$gamlss_family, "BI")) args$bd <- rep(1, length(params$mu))
  args
}

mvt_call_gamlss_dist <- function(prefix, spec, args) {
  fn <- get(paste0(prefix, spec$gamlss_family), envir = asNamespace("gamlss.dist"))
  args <- args[names(args) %in% formalArgs(fn)]
  as.numeric(do.call(fn, args))
}

mvt_param_vector <- function(params, name, n) {
  value <- params[[name]]
  if (is.null(value)) return(rep(NA_real_, n))
  value <- as.numeric(value)
  if (length(value) == 1L) return(rep(value, n))
  if (length(value) == n) return(value)
  rep(NA_real_, n)
}

mvt_truth_params <- function(data) {
  params <- list(mu = data$true_mu)
  if ("true_sigma" %in% names(data)) params$sigma <- data$true_sigma
  if ("true_nu" %in% names(data)) params$nu <- data$true_nu
  params
}

mvt_distribution_metrics <- function(data, params, spec, p = 0.9, interval_level = 0.95) {
  out <- mvt_empty_distribution_metrics()
  y <- as.numeric(data$response)
  n <- length(y)
  params$mu <- mvt_param_vector(params, "mu", n)
  params$sigma <- mvt_param_vector(params, "sigma", n)
  params$nu <- mvt_param_vector(params, "nu", n)
  params$tau <- mvt_param_vector(params, "tau", n)
  params <- params[vapply(params, function(x) any(is.finite(x)), logical(1))]

  if ("true_mu" %in% names(data)) {
    ok_truth <- is.finite(data$true_mu) & is.finite(params$mu)
    if (any(ok_truth)) {
      err_mu <- params$mu[ok_truth] - data$true_mu[ok_truth]
      out["benchmark_mean_bias"] <- mean(err_mu, na.rm = TRUE)
      out["benchmark_mean_mae"] <- mean(abs(err_mu), na.rm = TRUE)
      out["benchmark_mean_rmse"] <- sqrt(mean(err_mu^2, na.rm = TRUE))
    }
  }

  dens_args <- mvt_distribution_args(y, params, spec)
  dens <- tryCatch(mvt_call_gamlss_dist("d", spec, dens_args), error = function(e) rep(NA_real_, n))
  ok_density <- is.finite(dens) & dens > 0
  if (any(ok_density)) {
    out["benchmark_neg_log_score"] <- mean(-log(pmax(dens[ok_density], .Machine$double.xmin)), na.rm = TRUE)
  }

  cdf_args <- dens_args
  names(cdf_args)[names(cdf_args) == "x"] <- "q"
  cdf <- tryCatch(mvt_call_gamlss_dist("p", spec, cdf_args), error = function(e) rep(NA_real_, n))
  pit <- cdf
  if (spec$gamlss_family %in% c("PO", "BI")) {
    pmass <- dens
    pit <- cdf - 0.5 * pmass
  }
  pit <- pmin(pmax(pit, 0), 1)
  ok_pit <- is.finite(pit)
  if (sum(ok_pit) >= 3L) {
    pit_ok <- pit[ok_pit]
    out["benchmark_pit_ks_p_value"] <- tryCatch(
      suppressWarnings(stats::ks.test(pit_ok, "punif")$p.value),
      error = function(e) NA_real_
    )
    out["benchmark_pit_mean_abs_error"] <- abs(mean(pit_ok, na.rm = TRUE) - 0.5)
    out["benchmark_tail_error_lower_05"] <- mean(pit_ok <= 0.05, na.rm = TRUE) - 0.05
    out["benchmark_tail_error_upper_05"] <- mean(pit_ok >= 0.95, na.rm = TRUE) - 0.05
  }

  alpha <- (1 - interval_level) / 2
  q_args <- mvt_distribution_args(NULL, params, spec)
  q_args$p <- rep(c(alpha, 1 - alpha, p), each = n)
  q_args$mu <- rep(params$mu, times = 3L)
  if (!is.null(params$sigma)) q_args$sigma <- rep(params$sigma, times = 3L)
  if (!is.null(params$nu)) q_args$nu <- rep(params$nu, times = 3L)
  if (!is.null(params$tau)) q_args$tau <- rep(params$tau, times = 3L)
  if (identical(spec$gamlss_family, "BI")) q_args$bd <- rep(1, 3L * n)
  q_pred <- tryCatch(matrix(mvt_call_gamlss_dist("q", spec, q_args), nrow = n, ncol = 3L), error = function(e) matrix(NA_real_, n, 3L))
  lower <- q_pred[, 1L]
  upper <- q_pred[, 2L]
  q_p <- q_pred[, 3L]

  ok_interval <- is.finite(y) & is.finite(lower) & is.finite(upper)
  if (any(ok_interval)) {
    out["benchmark_interval_coverage_95"] <- mean(y[ok_interval] >= lower[ok_interval] & y[ok_interval] <= upper[ok_interval])
    out["benchmark_interval_width_95"] <- mean(upper[ok_interval] - lower[ok_interval], na.rm = TRUE)
  }

  truth <- mvt_truth_params(data)
  truth_q_args <- mvt_distribution_args(NULL, truth, spec)
  truth_q_args$p <- rep(p, n)
  truth_q <- tryCatch(mvt_call_gamlss_dist("q", spec, truth_q_args), error = function(e) rep(NA_real_, n))
  ok_q <- is.finite(q_p) & is.finite(truth_q)
  if (any(ok_q)) {
    out["benchmark_q90_mae"] <- mean(abs(q_p[ok_q] - truth_q[ok_q]), na.rm = TRUE)
  }

  cdf_truth_args <- mvt_distribution_args(truth_q, params, spec)
  names(cdf_truth_args)[names(cdf_truth_args) == "x"] <- "q"
  cdf_at_truth <- tryCatch(mvt_call_gamlss_dist("p", spec, cdf_truth_args), error = function(e) rep(NA_real_, n))
  ok_tail <- is.finite(cdf_at_truth) & is.finite(truth_q)
  if (any(ok_tail)) {
    out["benchmark_upper_tail_error_90"] <- mean((1 - cdf_at_truth[ok_tail]) - (1 - p), na.rm = TRUE)
  }
  out
}

mvt_method_result_row <- function(row, spec, method, fit_result, pred = NULL, extra = list()) {
  success <- !inherits(fit_result$value, "error") && !is.null(fit_result$value)
  error_msg <- if (success) {
    NA_character_
  } else if (inherits(fit_result$value, "error")) {
    conditionMessage(fit_result$value)
  } else {
    "Fit returned NULL."
  }
  y <- as.numeric(extra$y %||% numeric())
  scores <- if (!is.null(pred) && length(y) == length(pred)) {
    mvt_score_predictions(y, pred)
  } else {
    c(mae = NA_real_, rmse = NA_real_)
  }
  out <- data.frame(
    case_id = row$case_id,
    scenario = row$scenario,
    generator = row$generator,
    dependence = row$dependence,
    correlation_level = row$correlation_level,
    n_time = row$n_time,
    n_subject = row$n_subject,
    total_rows = row$total_rows,
    family = spec$family,
    gamlss_family = spec$gamlss_family,
    rep = row$rep,
    method = method,
    available = TRUE,
    success = success,
    status_class = mvt_fit_status_class(
      success,
      error_msg,
      fit_result$warnings %||% character()
    ),
    elapsed_sec = fit_result$elapsed_sec %||% NA_real_,
    nobs = row$total_rows,
    mae = unname(scores[["mae"]]),
    rmse = unname(scores[["rmse"]]),
    as.list(mvt_empty_distribution_metrics()),
    warning = if (length(fit_result$warnings)) paste(fit_result$warnings, collapse = " | ") else NA_character_,
    error = error_msg,
    stringsAsFactors = FALSE
  )
  out$failure_reason_short <- mvt_failure_reason_short(
    out$success,
    out$status_class,
    out$error,
    out$warning
  )
  for (nm in names(extra)) {
    if (!nm %in% c("y")) out[[nm]] <- mvt_scalar_extra(extra[[nm]])
  }
  out
}

mvt_fit_gamlss_longitudinal <- function(dat, spec, dep) {
  fit_dat <- dat[c("subject", "time", "response", "x", "z")]
  sigma_formula <- if (identical(spec$gamlss_family, "GG")) ~ time + z else ~1
  nu_formula <- if (identical(spec$gamlss_family, "GG")) ~ x else ~1
  gamlss.longitudinal::gamlss_longitudinal(
    dataset = fit_dat,
    margin_dist = spec$margin_dist(),
    copula_dist = "N",
    time_var = "time",
    subject_var = "subject",
    mu.formula = response ~ time + x + z,
    sigma.formula = sigma_formula,
    nu.formula = nu_formula,
    tau.formula = ~1,
    theta.formula = dep$theta_formula,
    zeta.formula = ~1,
    include_dlcopdpar = TRUE,
    method = "RS",
    start_from = NA,
    warm_start_joint = TRUE,
    compute_vcov = mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_COMPUTE_VCOV", TRUE),
    max_elapsed_sec = mvt_env_num("GAMLSS_LONGITUDINAL_MVT_PRIMARY_MAX_ELAPSED_SEC", 180),
    max_outer_iter = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_MAX_OUTER_ITER", 60L),
    max_inner_iter = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_MAX_INNER_ITER", 60L),
    verbose = 0
  )
}

mvt_case_spec_dep <- function(row) {
  families <- mvt_family_specs(include_special = TRUE)
  deps <- mvt_dependence_specs(include_appendix = TRUE)
  list(
    spec = families[[row$family_name]],
    dep = mvt_resize_dependence(deps[[row$dependence_name]], row$n_time)
  )
}

mvt_fit_gamlss_longitudinal_for_case <- function(dat, row) {
  parts <- mvt_case_spec_dep(row)
  mvt_fit_gamlss_longitudinal(dat, parts$spec, parts$dep)
}

mvt_fit_standard_models <- function(dat, spec, row, dep) {
  if (is.null(spec$standard_family)) {
    return(list(results = data.frame(), coefficients = data.frame(), fits = list()))
  }
  comparators <- mvt_standard_comparators()
  if (length(comparators) == 0L) {
    return(list(results = data.frame(), coefficients = data.frame(), fits = list()))
  }
  fit <- gamlss.longitudinal::benchmark_standard_models(
    data = dat,
    formula = response ~ time + x + z,
    subject_var = "subject",
    family = spec$standard_family,
    comparators = comparators,
    correlation = "independence",
    truth_family = spec$gamlss_family,
    distributional_metrics = TRUE
  )
  fit$results$case_id <- row$case_id
  fit$results$scenario <- row$scenario
  fit$results$generator <- row$generator
  fit$results$dependence <- row$dependence
  fit$results$correlation_level <- row$correlation_level
  fit$results$n_time <- row$n_time
  fit$results$n_subject <- row$n_subject
  fit$results$total_rows <- row$total_rows
  fit$results$family <- spec$family
  fit$results$gamlss_family <- spec$gamlss_family
  fit$results$rep <- row$rep
  coefs <- fit$coefficients$long
  if (nrow(coefs) > 0L) {
    coefs <- cbind(
      data.frame(
        case_id = row$case_id,
        scenario = row$scenario,
        generator = row$generator,
        dependence = row$dependence,
        correlation_level = row$correlation_level,
        n_time = row$n_time,
        n_subject = row$n_subject,
        total_rows = row$total_rows,
        family = spec$family,
        gamlss_family = spec$gamlss_family,
        rep = row$rep,
        stringsAsFactors = FALSE
      ),
      coefs,
      row.names = NULL
    )
  }
  list(results = fit$results, coefficients = coefs, fits = fit$fits)
}

mvt_run_one_gee <- function(dat, spec, row, corstr) {
  if (is.null(spec$standard_family)) return(list(results = data.frame(), coefficients = data.frame(), fit = NULL))
  timeout <- if (identical(corstr, "unstructured")) {
    mvt_env_num("GAMLSS_LONGITUDINAL_MVT_GEE_UNSTRUCTURED_TIMEOUT_SEC", if (row$n_time >= 50L) 30 else Inf)
  } else {
    mvt_env_num("GAMLSS_LONGITUDINAL_MVT_GEE_TIMEOUT_SEC", Inf)
  }
  run <- function() {
    gamlss.longitudinal::benchmark_standard_models(
      data = dat,
      formula = response ~ time + x + z,
      subject_var = "subject",
      family = spec$standard_family,
      comparators = "gee",
      correlation = corstr,
      truth_family = spec$gamlss_family,
      distributional_metrics = TRUE,
      waves = dat$time_index
    )
  }
  start <- Sys.time()
  value <- if (is.finite(timeout) && timeout > 0 && requireNamespace("callr", quietly = TRUE)) {
    tryCatch(callr::r(function(dat, spec, corstr) {
      gamlss.longitudinal::benchmark_standard_models(
        data = dat,
        formula = response ~ time + x + z,
        subject_var = "subject",
        family = spec$standard_family,
        comparators = "gee",
        correlation = corstr,
        truth_family = spec$gamlss_family,
        distributional_metrics = TRUE,
        waves = dat$time_index
      )
    }, args = list(dat = dat, spec = spec, corstr = corstr), timeout = timeout), error = function(e) e)
  } else {
    tryCatch(run(), error = function(e) e)
  }
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  method <- paste0("gee_", corstr)
  if (inherits(value, "error")) {
    results <- mvt_method_result_row(
      row,
      spec,
      method,
      list(value = value, warnings = character(), elapsed_sec = elapsed),
      extra = list(y = dat$response)
    )
    return(list(results = results, coefficients = data.frame(), fit = NULL))
  }
  value$results$method <- method
  value$results$comparator <- method
  value$results$elapsed_sec <- elapsed
  prefix <- data.frame(
    case_id = row$case_id,
    scenario = row$scenario,
    generator = row$generator,
    dependence = row$dependence,
    correlation_level = row$correlation_level,
    n_time = row$n_time,
    n_subject = row$n_subject,
    total_rows = row$total_rows,
    family = spec$family,
    gamlss_family = spec$gamlss_family,
    rep = row$rep,
    stringsAsFactors = FALSE
  )
  results <- cbind(prefix[rep(1L, nrow(value$results)), , drop = FALSE], value$results, row.names = NULL)
  coefs <- value$coefficients$long
  if (nrow(coefs) > 0L) {
    coefs$method <- method
    coefs <- cbind(prefix[rep(1L, nrow(coefs)), , drop = FALSE], coefs, row.names = NULL)
  }
  list(results = results, coefficients = coefs, fit = value$fits$gee)
}

mvt_fit_glmm_slope <- function(dat, spec) {
  if (is.null(spec$standard_family)) {
    stop("Random-slope GLMM sensitivity is unavailable for this non-standard margin.", call. = FALSE)
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for the glmm_slope comparator.", call. = FALSE)
  }
  form <- response ~ time + x + z + (1 + time | subject)
  control <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000))
  if (identical(spec$standard_family$family, "gaussian")) {
    lme4::lmer(form, data = dat, control = lme4::lmerControl(optimizer = "bobyqa"))
  } else {
    lme4::glmer(form, data = dat, family = spec$standard_family, control = control)
  }
}

mvt_run_glmm_slope <- function(dat, spec, row) {
  fit <- mvt_elapsed_capture(mvt_fit_glmm_slope(dat, spec))
  pred <- if (!inherits(fit$value, "error") && !is.null(fit$value)) {
    tryCatch(
      as.numeric(stats::predict(fit$value, newdata = dat, type = "response", re.form = NA)),
      error = function(e) rep(NA_real_, nrow(dat))
    )
  } else {
    NULL
  }
  result <- mvt_method_result_row(row, spec, "glmm_slope", fit, pred = pred, extra = list(y = dat$response))
  coefs <- if (!inherits(fit$value, "error") && !is.null(fit$value)) {
    mvt_coef_table_one(fit$value, "glmm_slope", row, spec)
  } else {
    data.frame()
  }
  list(results = result, coefficients = coefs, fit = if (inherits(fit$value, "error")) NULL else fit$value)
}

mvt_predict_gamlss_params <- function(fit, data, spec) {
  params <- list()
  for (what in c("mu", "sigma", "nu", "tau")) {
    pred <- tryCatch(stats::fitted(fit, parameter = what), error = function(e) NULL)
    if (!is.null(pred) && length(pred) == nrow(data)) params[[what]] <- as.numeric(pred)
  }
  if (length(params) == 0L) {
    params$mu <- as.numeric(stats::fitted(fit))
  }
  if (!"sigma" %in% names(params) && !is.null(spec$sigma) && !is.function(spec$sigma) && !is.na(spec$sigma)[[1L]]) {
    params$sigma <- rep(as.numeric(spec$sigma)[[1L]], nrow(data))
  }
  params
}

mvt_margin_pit <- function(fit, data, spec, randomized_discrete = FALSE, seed = NULL) {
  pfun <- get(paste0("p", spec$gamlss_family), envir = asNamespace("gamlss.dist"))
  params <- mvt_predict_gamlss_params(fit, data, spec)
  args <- c(list(q = data$response), params)
  if (identical(spec$gamlss_family, "BI")) args$bd <- rep(1, nrow(data))
  args <- args[names(args) %in% formalArgs(pfun)]
  upper <- pmin(pmax(as.numeric(do.call(pfun, args)), 1e-6), 1 - 1e-6)
  if (!isTRUE(randomized_discrete) || !identical(spec$family, "poisson") && !identical(spec$family, "binomial")) {
    return(upper)
  }
  lower_args <- args
  lower_args$q <- data$response - 1
  lower <- pmin(pmax(as.numeric(do.call(pfun, lower_args)), 0), 1 - 1e-6)
  lower <- pmin(lower, upper)
  if (!is.null(seed)) set.seed(seed)
  pmin(pmax(lower + stats::runif(length(upper)) * pmax(upper - lower, 0), 1e-6), 1 - 1e-6)
}

mvt_gamcopula_independence_pair <- function(error = NA_character_) {
  structure(
    list(error = error),
    class = "mvt_gamcopula_independence_pair"
  )
}

mvt_fit_gamcopula <- function(dat, spec, dep) {
  if (!requireNamespace("gamCopula", quietly = TRUE)) {
    stop("Package 'gamCopula' is required for the multivariate simulation module.", call. = FALSE)
  }
  if (!requireNamespace("gamlss", quietly = TRUE)) {
    stop("Package 'gamlss' is required for the gamCopula two-stage margin.", call. = FALSE)
  }
  fit_dat <- dat[c("subject", "time_index", "time", "response", "x", "z")]
  margin_fit <- gamlss::gamlss(
    response ~ time + x + z,
    sigma.formula = if (identical(spec$gamlss_family, "GG")) ~ time + z else ~1,
    nu.formula = if (identical(spec$gamlss_family, "GG")) ~ x else ~1,
    family = spec$margin_dist(),
    data = fit_dat,
    trace = FALSE
  )
  dat$u_hat <- mvt_margin_pit(
    margin_fit,
    fit_dat,
    spec,
    randomized_discrete = TRUE,
    seed = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_GAMCOPULA_PIT_SEED", 7349L)
  )
  families <- c(1, 3, 4, 5)
  fits <- list()
  pair_covariates <- list()
  time_pairs <- sort(unique(dat$time_index))
  for (j in seq_len(length(time_pairs) - 1L)) {
    left <- dat[dat$time_index == time_pairs[[j]], c("subject", "u_hat", "time", "x", "z"), drop = FALSE]
    right <- dat[dat$time_index == time_pairs[[j + 1L]], c("subject", "u_hat"), drop = FALSE]
    names(left)[names(left) == "u_hat"] <- "u1"
    names(right)[names(right) == "u_hat"] <- "u2"
    pair <- merge(left, right, by = "subject", sort = FALSE)
    pair$z_binary <- pair$z
    covariates <- if (identical(dep$dependence, "covariate_dependent_adjacent")) {
      c("x", "z_binary")
    } else if (identical(dep$dependence, "time_varying_adjacent") || identical(dep$dependence, "ar1")) {
      "time"
    } else {
      character()
    }
    form <- if (length(covariates) > 0L) {
      stats::as.formula(paste("~", paste(covariates, collapse = " + ")))
    } else {
      ~1
    }
    udata <- pair[c("u1", "u2")]
    lin_covs <- if (length(covariates) > 0L) pair[covariates] else NULL
    fit_data <- if (is.null(lin_covs)) udata else cbind(udata, lin_covs)
    selected <- tryCatch(
      gamCopula::gamBiCopSelect(
        udata = udata,
        lin.covs = lin_covs,
        familyset = families,
        rotations = FALSE
      ),
      error = function(e) NULL
    )
    fit <- if (!is.null(selected) && !is.null(selected$res)) {
      selected$res
    } else {
      fallback <- tryCatch(
        gamCopula::gamBiCopFit(fit_data, form, family = 1)$res,
        error = function(e) e
      )
      if (inherits(fallback, "error")) {
        warning(
          "gamCopula pair ",
          j,
          " failed; using independence fallback for this pair: ",
          conditionMessage(fallback),
          call. = FALSE
        )
        mvt_gamcopula_independence_pair(conditionMessage(fallback))
      } else {
        fallback
      }
    }
    if (length(covariates) > 0L && !inherits(fit, "mvt_gamcopula_independence_pair")) {
      selected_family <- tryCatch(as.integer(methods::slot(fit, "family")), error = function(e) NA_integer_)
      if (length(selected_family) > 0L && is.finite(selected_family[[1L]]) && selected_family[[1L]] > 0L) {
        formula_fit <- tryCatch(
          gamCopula::gamBiCopFit(fit_data, form, family = selected_family[[1L]])$res,
          error = function(e) e
        )
        if (!inherits(formula_fit, "error")) {
          fit <- formula_fit
        } else {
          warning(
            "gamCopula pair ",
            j,
            " selected family refit failed; retaining selected fit: ",
            conditionMessage(formula_fit),
            call. = FALSE
          )
        }
      }
    }
    fits[[j]] <- fit
    pair_covariates[[j]] <- pair[c("subject", "time", "x", "z_binary")]
  }
  structure(
    list(
      margin_fit = margin_fit,
      pair_fits = fits,
      pair_covariates = pair_covariates,
      spec = spec,
      dep = dep
    ),
    class = "mvt_gamcopula_fit"
  )
}

mvt_fit_gamcopula_for_case <- function(dat, row) {
  parts <- mvt_case_spec_dep(row)
  mvt_fit_gamcopula(dat, parts$spec, parts$dep)
}

mvt_fit_gamcopula_vine <- function(dat, spec, dep, simplified = FALSE) {
  if (!requireNamespace("gamCopula", quietly = TRUE)) {
    stop("Package 'gamCopula' is required for the multivariate simulation module.", call. = FALSE)
  }
  if (!requireNamespace("gamlss", quietly = TRUE)) {
    stop("Package 'gamlss' is required for the gamCopula two-stage margin.", call. = FALSE)
  }
  fit_dat <- dat[c("subject", "time_index", "time", "response", "x", "z")]
  margin_fit <- gamlss::gamlss(
    response ~ time + x + z,
    sigma.formula = if (identical(spec$gamlss_family, "GG")) ~ time + z else ~1,
    nu.formula = if (identical(spec$gamlss_family, "GG")) ~ x else ~1,
    family = spec$margin_dist(),
    data = fit_dat,
    trace = FALSE
  )
  dat$u_hat <- mvt_margin_pit(
    margin_fit,
    fit_dat,
    spec,
    randomized_discrete = TRUE,
    seed = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_GAMCOPULA_PIT_SEED", 7349L)
  )
  dat <- dat[order(dat$subject, dat$time_index), , drop = FALSE]
  time_values <- sort(unique(dat$time_index))
  udata <- do.call(cbind, lapply(time_values, function(t) {
    dat$u_hat[dat$time_index == t]
  }))
  colnames(udata) <- paste0("t", time_values)
  if (any(!is.finite(udata))) {
    stop("gamCopula vine pseudo-observations contain non-finite values.", call. = FALSE)
  }
  vine_fit <- tryCatch(
    gamCopula::gamVineStructureSelect(
      udata = udata,
      simplified = isTRUE(simplified),
      familyset = c(1, 2, 5),
      rotations = TRUE,
      verbose = FALSE
    ),
    error = function(e) e
  )
  vine_engine <- "gamCopula"
  if (inherits(vine_fit, "error")) {
    if (!requireNamespace("VineCopula", quietly = TRUE)) {
      stop("gamCopula vine fit failed and VineCopula fallback is unavailable: ", conditionMessage(vine_fit), call. = FALSE)
    }
    vine_fit <- VineCopula::RVineStructureSelect(
      udata,
      familyset = c(1, 2, 3, 4, 5),
      rotations = TRUE,
      progress = FALSE
    )
    vine_engine <- "VineCopula"
  }
  vine_simplified <- isTRUE(simplified)
  if (identical(vine_engine, "VineCopula")) {
    vine_simplified <- TRUE
  }
  structure(
    list(
      margin_fit = margin_fit,
      vine_fit = vine_fit,
      vine_engine = vine_engine,
      vine_simplified = vine_simplified,
      vine_simplified_requested = isTRUE(simplified),
      udata = udata,
      udata_names = colnames(udata),
      spec = spec,
      dep = dep
    ),
    class = c("mvt_gamcopula_vine_fit", "mvt_gamcopula_fit")
  )
}

mvt_fit_gamcopula_vine_for_case <- function(dat, row) {
  parts <- mvt_case_spec_dep(row)
  mvt_fit_gamcopula_vine(dat, parts$spec, parts$dep, simplified = FALSE)
}

mvt_fit_gamcopula_vine_simplified_for_case <- function(dat, row) {
  parts <- mvt_case_spec_dep(row)
  mvt_fit_gamcopula_vine(dat, parts$spec, parts$dep, simplified = TRUE)
}

mvt_predict_gamcopula_mean <- function(fit, data) {
  as.numeric(stats::fitted(fit$margin_fit, what = "mu"))
}

mvt_predict_gamcopula_params <- function(fit, data, spec) {
  params <- list()
  for (what in c("mu", "sigma", "nu", "tau")) {
    pred <- tryCatch(stats::fitted(fit$margin_fit, what = what), error = function(e) NULL)
    if (!is.null(pred) && length(pred) == nrow(data)) params[[what]] <- as.numeric(pred)
  }
  if (!"mu" %in% names(params)) {
    params$mu <- mvt_predict_gamcopula_mean(fit, data)
  }
  if (!"sigma" %in% names(params) && !is.null(spec$sigma) && !is.function(spec$sigma) && !is.na(spec$sigma)[[1L]]) {
    params$sigma <- rep(as.numeric(spec$sigma)[[1L]], nrow(data))
  }
  params
}

mvt_vinecopula_family_df <- function(family, par = NA_real_, par2 = NA_real_) {
  family <- abs(suppressWarnings(as.integer(family)))
  if (length(family) == 0L || is.na(family[[1L]]) || family[[1L]] == 0L) return(0)
  if (family[[1L]] %in% c(2L, 7L, 8L, 9L, 10L, 17L, 18L, 19L, 20L, 27L, 28L, 29L, 30L, 37L, 38L, 39L, 40L)) {
    return(2)
  }
  if (is.finite(suppressWarnings(as.numeric(par2))) && abs(suppressWarnings(as.numeric(par2))) > 0) return(2)
  if (is.finite(suppressWarnings(as.numeric(par))) || family[[1L]] != 0L) 1 else 0
}

mvt_gamcopula_bicop_loglik <- function(object) {
  if (inherits(object, "gamBiCop")) {
    value <- tryCatch(suppressWarnings(stats::logLik(object)), error = function(e) NULL)
    if (is.null(value) && exists("logLik.gamBiCop", envir = asNamespace("gamCopula"), mode = "function")) {
      value <- tryCatch(suppressWarnings(get("logLik.gamBiCop", envir = asNamespace("gamCopula"))(object)), error = function(e) NULL)
    }
    value <- as.numeric(value)
    if (length(value) > 0L && is.finite(value[[1L]])) return(value[[1L]])
    return(NA_real_)
  }
  if (is.list(object) && !is.null(object$family)) return(0)
  NA_real_
}

mvt_gamcopula_bicop_df <- function(object) {
  if (inherits(object, "gamBiCop")) {
    value <- tryCatch(suppressWarnings(stats::logLik(object)), error = function(e) NULL)
    if (is.null(value) && exists("logLik.gamBiCop", envir = asNamespace("gamCopula"), mode = "function")) {
      value <- tryCatch(suppressWarnings(get("logLik.gamBiCop", envir = asNamespace("gamCopula"))(object)), error = function(e) NULL)
    }
    df <- attr(value, "df")
    if (length(df) > 0L && is.finite(df[[1L]])) return(as.numeric(df[[1L]]))
    edf <- tryCatch(sum(object@model$edf), error = function(e) NA_real_)
    if (is.finite(edf)) return(edf + if (identical(as.integer(object@family), 2L)) 1 else 0)
    return(NA_real_)
  }
  if (is.list(object) && !is.null(object$family)) {
    return(mvt_vinecopula_family_df(object$family, object$par %||% NA_real_, object$par2 %||% NA_real_))
  }
  NA_real_
}

mvt_gamcopula_bicop_nobs <- function(object) {
  if (inherits(object, "gamBiCop")) {
    n <- tryCatch(stats::nobs(object), error = function(e) NA_real_)
    n <- suppressWarnings(as.numeric(n))
    if (length(n) > 0L && is.finite(n[[1L]])) return(n[[1L]])
    n <- tryCatch(nrow(stats::na.omit(object@model$data)), error = function(e) NA_real_)
    if (is.finite(n)) return(n)
  }
  NA_real_
}

mvt_gamcopula_bicop_bic <- function(object) {
  ll <- mvt_gamcopula_bicop_loglik(object)
  df <- mvt_gamcopula_bicop_df(object)
  n <- mvt_gamcopula_bicop_nobs(object)
  if (is.list(object) && !inherits(object, "gamBiCop")) n <- NA_real_
  if (is.finite(ll) && is.finite(df) && is.finite(n) && n > 1) -2 * ll + log(n) * df else NA_real_
}

mvt_vinecopula_df <- function(vine_fit) {
  family <- suppressWarnings(as.integer(vine_fit$family))
  par <- suppressWarnings(as.numeric(vine_fit$par))
  par2 <- suppressWarnings(as.numeric(vine_fit$par2))
  sum(mapply(mvt_vinecopula_family_df, family, par, par2), na.rm = TRUE)
}

mvt_copula_loglik_gamcopula <- function(fit) {
  if (inherits(fit, "mvt_gamcopula_vine_fit")) {
    if (identical(fit$vine_engine, "VineCopula")) {
      value <- suppressWarnings(as.numeric(fit$vine_fit$logLik))
      if (length(value) > 0L && is.finite(value[[1L]])) return(value[[1L]])
      vine_pdf <- tryCatch(VineCopula::RVinePDF(fit$udata, fit$vine_fit), error = function(e) NA_real_)
      return(sum(log(pmax(as.numeric(vine_pdf), .Machine$double.xmin)), na.rm = TRUE))
    }
    return(sum(vapply(fit$vine_fit@model, mvt_gamcopula_bicop_loglik, numeric(1)), na.rm = TRUE))
  }
  sum(vapply(fit$pair_fits, mvt_gamcopula_bicop_loglik, numeric(1)), na.rm = TRUE)
}

mvt_gamcopula_copula_df <- function(fit) {
  if (inherits(fit, "mvt_gamcopula_vine_fit")) {
    if (identical(fit$vine_engine, "VineCopula")) return(mvt_vinecopula_df(fit$vine_fit))
    return(sum(vapply(fit$vine_fit@model, mvt_gamcopula_bicop_df, numeric(1)), na.rm = TRUE))
  }
  sum(vapply(fit$pair_fits, mvt_gamcopula_bicop_df, numeric(1)), na.rm = TRUE)
}

mvt_gamcopula_copula_bic <- function(fit) {
  if (inherits(fit, "mvt_gamcopula_vine_fit")) {
    if (identical(fit$vine_engine, "VineCopula")) {
      value <- suppressWarnings(as.numeric(fit$vine_fit$BIC))
      if (length(value) > 0L && is.finite(value[[1L]])) return(value[[1L]])
      ll <- mvt_copula_loglik_gamcopula(fit)
      df <- mvt_gamcopula_copula_df(fit)
      n <- nrow(fit$udata)
      return(if (is.finite(ll) && is.finite(df) && is.finite(n) && n > 1) -2 * ll + log(n) * df else NA_real_)
    }
    pair_bic <- vapply(fit$vine_fit@model, mvt_gamcopula_bicop_bic, numeric(1))
    missing <- !is.finite(pair_bic)
    if (any(missing)) {
      pair_ll <- vapply(fit$vine_fit@model[missing], mvt_gamcopula_bicop_loglik, numeric(1))
      pair_df <- vapply(fit$vine_fit@model[missing], mvt_gamcopula_bicop_df, numeric(1))
      pair_bic[missing] <- -2 * pair_ll + log(nrow(fit$udata)) * pair_df
    }
    return(sum(pair_bic, na.rm = TRUE))
  }
  pair_bic <- vapply(fit$pair_fits, mvt_gamcopula_bicop_bic, numeric(1))
  sum(pair_bic, na.rm = TRUE)
}

mvt_loglik_gamcopula <- function(fit) {
  margin_ll <- suppressWarnings(as.numeric(stats::logLik(fit$margin_fit)))
  margin_ll + mvt_copula_loglik_gamcopula(fit)
}

mvt_loglik_df_gamcopula <- function(fit) {
  mvt_loglik_df(fit$margin_fit) + mvt_gamcopula_copula_df(fit)
}

mvt_aic_gamcopula <- function(fit) {
  ll <- mvt_loglik_gamcopula(fit)
  df <- mvt_loglik_df_gamcopula(fit)
  if (is.finite(ll) && is.finite(df)) -2 * ll + 2 * df else NA_real_
}

mvt_bic_gamcopula <- function(fit) {
  margin_bic <- mvt_bic_value(fit$margin_fit)
  copula_bic <- mvt_gamcopula_copula_bic(fit)
  if (is.finite(margin_bic) && is.finite(copula_bic)) margin_bic + copula_bic else NA_real_
}

mvt_coef_table_gamlss_margin <- function(fit, method) {
  coef_mu <- tryCatch(stats::coef(fit, what = "mu"), error = function(e) numeric())
  if (length(coef_mu) == 0L) return(data.frame())
  se_mu <- rep(NA_real_, length(coef_mu))
  names(se_mu) <- names(coef_mu)
  coef_table <- tryCatch({
    tmp <- NULL
    invisible(utils::capture.output(tmp <- summary(fit, type = "qr")))
    as.matrix(tmp)
  }, error = function(e) NULL)
  if (!is.null(coef_table) && nrow(coef_table) >= length(coef_mu)) {
    se_col <- match("Std. Error", colnames(coef_table))
    if (is.na(se_col)) se_col <- 2L
    se_mu <- as.numeric(coef_table[seq_along(coef_mu), se_col])
  }
  data.frame(
    method = method,
    parameter = "mu",
    term = names(coef_mu),
    estimate = as.numeric(coef_mu),
    std_error = as.numeric(se_mu),
    stringsAsFactors = FALSE
  )
}

mvt_coef_table_one <- function(fit, method, row, spec) {
  out <- tryCatch({
    if (inherits(fit, "gamlss.longitudinal")) {
      s <- summary(
        fit,
        include_vcov = mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_SUMMARY_VCOV", TRUE)
      )
      tbl <- as.data.frame(s$coefficients, stringsAsFactors = FALSE)
      tbl <- tbl[tbl$parameter == "mu", , drop = FALSE]
      data.frame(
        method = method,
        parameter = "mu",
        term = sub("^mu[.]", "", tbl$term),
        estimate = tbl$estimate,
        std_error = tbl$std_error,
        stringsAsFactors = FALSE
      )
    } else if (inherits(fit, "mvt_gamcopula_fit")) {
      mvt_coef_table_gamlss_margin(fit$margin_fit, method)
    } else {
      mat <- coef(summary(fit))
      se_col <- match("Std. Error", colnames(mat))
      if (is.na(se_col)) se_col <- 2L
      data.frame(
        method = method,
        parameter = "mu",
        term = rownames(mat),
        estimate = mat[, 1],
        std_error = mat[, se_col],
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) data.frame())
  if (nrow(out) == 0L) return(out)
  out$term <- mvt_normalize_term(out$term)
  out <- mvt_annotate_coefficients(out, spec)
  cbind(
    data.frame(
      case_id = row$case_id,
      scenario = row$scenario,
      generator = row$generator,
      dependence = row$dependence,
      correlation_level = row$correlation_level,
      n_time = row$n_time,
      n_subject = row$n_subject,
      total_rows = row$total_rows,
      family = spec$family,
      gamlss_family = spec$gamlss_family,
      rep = row$rep,
      stringsAsFactors = FALSE
    ),
    out,
    row.names = NULL
  )
}

mvt_normalize_term <- function(term) {
  term <- as.character(term)
  term[term %in% c("(Intercept)", "Intercept", "mu.(Intercept)", "mu.Intercept")] <- "intercept"
  term <- sub("^mu[.]", "", term)
  term[term %in% c("time_covariate", "time.covariate", "time_index", "timeindex")] <- "time"
  term
}

mvt_truth_for_term <- function(term, spec) {
  out <- rep(NA_real_, length(term))
  out[term == "intercept"] <- spec$eta$intercept
  out[term == "time"] <- spec$eta$time
  out[term == "x"] <- spec$eta$x
  out[term == "z"] <- spec$eta$z
  out
}

mvt_annotate_coefficients <- function(coefs, spec, level = 0.95) {
  if (nrow(coefs) == 0L) return(coefs)
  coefs$term <- mvt_normalize_term(coefs$term)
  coefs$truth <- mvt_truth_for_term(coefs$term, spec)
  z <- stats::qnorm((1 + level) / 2)
  coefs$conf.low <- coefs$estimate - z * coefs$std_error
  coefs$conf.high <- coefs$estimate + z * coefs$std_error
  coefs$bias <- coefs$estimate - coefs$truth
  coefs$ci_width <- coefs$conf.high - coefs$conf.low
  coefs$ci_covers_truth <- is.finite(coefs$truth) & coefs$conf.low <= coefs$truth & coefs$conf.high >= coefs$truth
  coefs$p_value <- 2 * stats::pnorm(abs(coefs$estimate / coefs$std_error), lower.tail = FALSE)
  coefs$false_positive <- is.finite(coefs$truth) & abs(coefs$truth) < 1e-12 & coefs$p_value < 0.05
  coefs
}

mvt_true_pair_dependence <- function(dat) {
  time_values <- sort(unique(dat$time_index))
  R <- attr(dat, "true_correlation_matrix")
  rows <- list()
  k <- 1L
  for (left in seq_len(length(time_values) - 1L)) {
    for (right in seq.int(left + 1L, length(time_values))) {
      lag <- right - left
      if (!is.null(R) && all(c(left, right) <= nrow(R))) {
        theta <- R[left, right]
      } else if ("true_theta" %in% names(dat)) {
        subject_theta <- vapply(split(dat, dat$subject), function(subject_data) {
          subject_data <- subject_data[order(subject_data$time_index), , drop = FALSE]
          edge <- subject_data$time_index > left & subject_data$time_index <= right
          vals <- subject_data$true_theta[edge]
          vals <- vals[is.finite(vals)]
          if (length(vals) == lag) prod(vals) else NA_real_
        }, numeric(1))
        theta <- mean(subject_theta, na.rm = TRUE)
      } else {
        theta <- NA_real_
      }
      rows[[k]] <- data.frame(
        time_left_idx = left,
        time_right_idx = right,
        lag = lag,
        true_theta = theta,
        true_tau = (2 / pi) * asin(pmax(-1, pmin(1, theta))),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
}

mvt_gamcopula_pair_tau_values <- function(pair_fit, newdata = NULL, n = NULL) {
  n_out <- if (!is.null(n)) {
    as.integer(n)
  } else if (is.data.frame(newdata)) {
    nrow(newdata)
  } else {
    1L
  }
  if (length(n_out) == 0L || is.na(n_out) || n_out < 1L) n_out <- 1L
  finish <- function(x) {
    x <- as.numeric(x)
    if (length(x) == 0L) return(rep(NA_real_, n_out))
    if (length(x) == 1L && n_out > 1L) x <- rep(x, n_out)
    if (length(x) != n_out) x <- rep(x, length.out = n_out)
    x
  }
  if (inherits(pair_fit, "mvt_gamcopula_independence_pair")) return(rep(0, n_out))
  if (is.list(pair_fit) && all(c("family", "par", "par2") %in% names(pair_fit))) {
    tau <- tryCatch(
      VineCopula::BiCopPar2Tau(pair_fit$family, pair_fit$par, pair_fit$par2),
      error = function(e) NA_real_
    )
    return(finish(tau))
  }
  if (inherits(pair_fit, "gamBiCop") && !is.null(newdata)) {
    if ("subject" %in% names(newdata)) newdata <- newdata[setdiff(names(newdata), "subject")]
    pred <- tryCatch(
      gamCopula::gamBiCopPredict(pair_fit, newdata = newdata, target = "tau", type = "response"),
      error = function(e) NULL
    )
    tau <- if (is.list(pred) && "tau" %in% names(pred)) pred$tau else pred
    tau <- as.numeric(tau)
    if (length(tau) > 0L && any(is.finite(tau))) return(finish(tau))

    pred_par <- tryCatch(
      gamCopula::gamBiCopPredict(pair_fit, newdata = newdata, target = "par", type = "response"),
      error = function(e) NULL
    )
    par <- if (is.list(pred_par) && "par" %in% names(pred_par)) pred_par$par else pred_par
    par <- as.numeric(par)
    if (length(par) > 0L && any(is.finite(par)) && requireNamespace("VineCopula", quietly = TRUE)) {
      par2 <- as.numeric(methods::slot(pair_fit, "par2"))
      if (length(par2) == 0L || !is.finite(par2[[1L]])) par2 <- 0
      tau <- VineCopula::BiCopPar2Tau(
        family = as.integer(methods::slot(pair_fit, "family")),
        par = par,
        par2 = rep(par2[[1L]], length(par))
      )
      return(finish(tau))
    }
  }
  pred <- tryCatch(
    gamCopula::gamBiCopPredict(pair_fit, target = "tau", type = "response"),
    error = function(e) NULL
  )
  if (is.list(pred) && "tau" %in% names(pred)) {
    return(finish(pred$tau))
  }
  if (is.numeric(pred)) {
    return(finish(pred))
  }
  rep(NA_real_, n_out)
}

mvt_gamcopula_pair_tau <- function(pair_fit, newdata = NULL) {
  tau <- mvt_gamcopula_pair_tau_values(pair_fit, newdata)
  mean(tau, na.rm = TRUE)
}

mvt_gamcopula_pair_par <- function(pair_fit, newdata = NULL) {
  if (inherits(pair_fit, "mvt_gamcopula_independence_pair")) {
    return(list(independent = TRUE, family = 0L, par = 0, par2 = 0))
  }
  if (is.list(pair_fit) && all(c("family", "par", "par2") %in% names(pair_fit))) {
    return(list(
      independent = FALSE,
      family = as.integer(pair_fit$family[[1L]]),
      par = as.numeric(pair_fit$par[[1L]]),
      par2 = as.numeric(pair_fit$par2[[1L]])
    ))
  }
  if (!inherits(pair_fit, "gamBiCop")) {
    stop("Unsupported gamCopula pair fit class: ", paste(class(pair_fit), collapse = "/"), call. = FALSE)
  }
  pred <- tryCatch(
    gamCopula::gamBiCopPredict(pair_fit, newdata = newdata, target = "par", type = "response"),
    error = function(e) e
  )
  if (inherits(pred, "error")) {
    stop("Could not predict gamCopula pair parameter: ", conditionMessage(pred), call. = FALSE)
  }
  par <- if (is.list(pred) && "par" %in% names(pred)) pred$par else pred
  par <- as.numeric(par)
  par <- par[is.finite(par)]
  if (length(par) == 0L) {
    stop("gamCopula pair parameter prediction returned no finite values.", call. = FALSE)
  }
  par2 <- as.numeric(methods::slot(pair_fit, "par2"))
  if (length(par2) == 0L || !is.finite(par2[[1L]])) par2 <- 0
  list(
    independent = FALSE,
    family = as.integer(methods::slot(pair_fit, "family")),
    par = par[[1L]],
    par2 = par2[[1L]]
  )
}

mvt_simulate_gamcopula_uniforms <- function(fit, data, nsim = 50L, seed = NULL) {
  if (!requireNamespace("VineCopula", quietly = TRUE)) {
    stop("Package 'VineCopula' is required to simulate gamCopula joint paths.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  data <- data[order(data$subject, data$time_index), , drop = FALSE]
  subjects <- unique(data$subject)
  n_time <- length(unique(data$time_index))
  if (inherits(fit, "mvt_gamcopula_vine_fit")) {
    out <- matrix(NA_real_, nrow = nrow(data), ncol = nsim)
    row_map <- split(seq_len(nrow(data)), data$subject)
    for (sim_id in seq_len(nsim)) {
      sim <- if (identical(fit$vine_engine, "VineCopula")) {
        VineCopula::RVineSim(length(subjects), fit$vine_fit)
      } else {
        gamCopula::gamVineSimulate(length(subjects), fit$vine_fit)
      }
      sim <- as.matrix(sim)
      if (ncol(sim) != n_time) {
        stop("gamCopula vine simulation dimension did not match the number of time points.", call. = FALSE)
      }
      for (subject_pos in seq_along(subjects)) {
        idx <- row_map[[as.character(subjects[[subject_pos]])]]
        idx <- idx[order(data$time_index[idx])]
        out[idx, sim_id] <- sim[subject_pos, seq_len(n_time)]
      }
    }
    return(pmin(pmax(out, 1e-6), 1 - 1e-6))
  }
  if (length(fit$pair_fits) != n_time - 1L) {
    stop("gamCopula pair-fit count does not match the number of adjacent time intervals.", call. = FALSE)
  }
  out <- matrix(NA_real_, nrow = nrow(data), ncol = nsim)
  row_map <- split(seq_len(nrow(data)), data$subject)
  eps <- 1e-6
  for (sim_id in seq_len(nsim)) {
    for (subject_id in subjects) {
      idx <- row_map[[as.character(subject_id)]]
      if (length(idx) != n_time) {
        stop("gamCopula simulation expects a balanced longitudinal panel.", call. = FALSE)
      }
      idx <- idx[order(data$time_index[idx])]
      u_prev <- stats::runif(1L, eps, 1 - eps)
      out[idx[[1L]], sim_id] <- u_prev
      for (j in seq_along(fit$pair_fits)) {
        covs <- fit$pair_covariates[[j]]
        newdata <- covs[covs$subject == subject_id, , drop = FALSE]
        if (nrow(newdata) == 0L) newdata <- data.frame(time = data$time[idx[[j]]], x = data$x[idx[[j]]], z_binary = data$z[idx[[j]]])
        newdata <- newdata[1L, setdiff(names(newdata), "subject"), drop = FALSE]
        par <- mvt_gamcopula_pair_par(fit$pair_fits[[j]], newdata = newdata)
        u_next <- if (isTRUE(par$independent)) {
          stats::runif(1L, eps, 1 - eps)
        } else {
          VineCopula::BiCopCondSim(
            N = 1L,
            cond.val = pmin(pmax(u_prev, eps), 1 - eps),
            cond.var = 1L,
            family = par$family,
            par = par$par,
            par2 = par$par2
          )
        }
        u_prev <- pmin(pmax(as.numeric(u_next)[[1L]], eps), 1 - eps)
        out[idx[[j + 1L]], sim_id] <- u_prev
      }
    }
  }
  out
}

mvt_quantile_matrix <- function(u, params, spec) {
  n <- nrow(u)
  nsim <- ncol(u)
  out <- matrix(NA_real_, nrow = n, ncol = nsim)
  for (sim_id in seq_len(nsim)) {
    q_args <- mvt_distribution_args(NULL, params, spec)
    q_args$p <- as.numeric(u[, sim_id])
    if (identical(spec$gamlss_family, "BI")) q_args$bd <- rep(1, n)
    out[, sim_id] <- tryCatch(
      mvt_call_gamlss_dist("q", spec, q_args),
      error = function(e) rep(NA_real_, n)
    )
  }
  out
}

mvt_simulate_gamcopula_response <- function(fit, data, spec, nsim = 50L, seed = NULL) {
  u <- mvt_simulate_gamcopula_uniforms(fit, data, nsim = nsim, seed = seed)
  params <- mvt_predict_gamcopula_params(fit, data, spec)
  params$mu <- mvt_param_vector(params, "mu", nrow(data))
  params$sigma <- mvt_param_vector(params, "sigma", nrow(data))
  params$nu <- mvt_param_vector(params, "nu", nrow(data))
  params$tau <- mvt_param_vector(params, "tau", nrow(data))
  params <- params[vapply(params, function(x) any(is.finite(x)), logical(1))]
  mvt_quantile_matrix(u, params, spec)
}

mvt_dependence_recovery_row <- function(dat, fit, method, row, spec) {
  truth <- mvt_true_pair_dependence(dat)
  if (nrow(truth) == 0L) return(data.frame())
  fitted_theta <- rep(NA_real_, nrow(truth))
  fitted_tau <- rep(NA_real_, nrow(truth))
  if (inherits(fit, "gamlss.longitudinal")) {
    dep <- tryCatch(gamlss.longitudinal::copula_time_summary(fit), error = function(e) NULL)
    if (!is.null(dep) && is.data.frame(dep$fit_data) && "theta_fit" %in% names(dep$fit_data)) {
      theta_by_time <- stats::aggregate(theta_fit ~ time, data = dep$fit_data, FUN = mean, na.rm = TRUE)
      theta_by_time <- theta_by_time[order(theta_by_time$time), , drop = FALSE]
      adjacent_theta <- theta_by_time$theta_fit
      for (i in seq_len(nrow(truth))) {
        vals <- adjacent_theta[seq.int(truth$time_left_idx[i], truth$time_right_idx[i] - 1L)]
        if (length(vals) > 0L && all(is.finite(vals))) fitted_theta[i] <- prod(vals)
      }
    }
  } else if (inherits(fit, "mvt_gamcopula_vine_fit")) {
    nsim <- mvt_env_int("GAMLSS_LONGITUDINAL_MVT_VINE_DEPENDENCE_NSIM", 500L)
    sim <- tryCatch({
      if (identical(fit$vine_engine, "VineCopula")) {
        VineCopula::RVineSim(nsim, fit$vine_fit)
      } else {
        gamCopula::gamVineSimulate(nsim, fit$vine_fit)
      }
    }, error = function(e) NULL)
    if (!is.null(sim)) {
      sim <- as.matrix(sim)
      for (i in seq_len(nrow(truth))) {
        left <- truth$time_left_idx[[i]]
        right <- truth$time_right_idx[[i]]
        if (all(c(left, right) <= ncol(sim))) {
          fitted_tau[i] <- suppressWarnings(stats::cor(sim[, left], sim[, right], method = "kendall", use = "complete.obs"))
          fitted_theta[i] <- sin(pi * pmax(-1, pmin(1, fitted_tau[i])) / 2)
        }
      }
    }
  } else if (inherits(fit, "mvt_gamcopula_fit")) {
    adjacent_theta <- vector("list", length(fit$pair_fits))
    for (j in seq_along(fit$pair_fits)) {
      covs <- fit$pair_covariates[[j]]
      tau <- mvt_gamcopula_pair_tau_values(fit$pair_fits[[j]], covs, n = nrow(covs))
      adjacent_theta[[j]] <- data.frame(
        subject = as.character(covs$subject),
        theta = sin(pi * tau / 2),
        stringsAsFactors = FALSE
      )
    }
    for (i in seq_len(nrow(truth))) {
      edge_ids <- seq.int(truth$time_left_idx[i], truth$time_right_idx[i] - 1L)
      edge_ids <- edge_ids[edge_ids >= 1L & edge_ids <= length(adjacent_theta)]
      if (length(edge_ids) > 0L) {
        subjects <- Reduce(intersect, lapply(adjacent_theta[edge_ids], function(x) x$subject))
        if (length(subjects) > 0L) {
          subject_theta <- vapply(subjects, function(subject_id) {
            vals <- vapply(adjacent_theta[edge_ids], function(edge_data) {
              edge_data$theta[match(subject_id, edge_data$subject)]
            }, numeric(1))
            if (all(is.finite(vals))) prod(vals) else NA_real_
          }, numeric(1))
          if (any(is.finite(subject_theta))) fitted_theta[i] <- mean(subject_theta, na.rm = TRUE)
        }
      }
    }
  }
  missing_tau <- !is.finite(fitted_tau) & is.finite(fitted_theta)
  fitted_tau[missing_tau] <- (2 / pi) * asin(pmax(-1, pmin(1, fitted_theta[missing_tau])))
  theta_err <- fitted_theta - truth$true_theta
  tau_err <- fitted_tau - truth$true_tau
  data.frame(
    case_id = row$case_id,
    scenario = row$scenario,
    generator = row$generator,
    dependence = row$dependence,
    correlation_level = row$correlation_level,
    n_time = row$n_time,
    n_subject = row$n_subject,
    total_rows = row$total_rows,
    family = spec$family,
    gamlss_family = spec$gamlss_family,
    rep = row$rep,
    method = method,
    dependence_scope = "all_pairs",
    dependence_n = sum(is.finite(theta_err) | is.finite(tau_err)),
    theta_mean_error = mean(theta_err, na.rm = TRUE),
    theta_mae = mean(abs(theta_err), na.rm = TRUE),
    theta_rmse = sqrt(mean(theta_err^2, na.rm = TRUE)),
    tau_mean_error = mean(tau_err, na.rm = TRUE),
    tau_mae = mean(abs(tau_err), na.rm = TRUE),
    tau_rmse = sqrt(mean(tau_err^2, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

mvt_variogram_empty <- function(row, spec, method, error = NA_character_) {
  data.frame(
    case_id = row$case_id,
    scenario = row$scenario,
    generator = row$generator,
    dependence = row$dependence,
    correlation_level = row$correlation_level,
    n_time = row$n_time,
    n_subject = row$n_subject,
    total_rows = row$total_rows,
    family = spec$family,
    gamlss_family = spec$gamlss_family,
    rep = row$rep,
    method = method,
    variogram_nsim = NA_integer_,
    variogram_score_p05 = NA_real_,
    variogram_score_p2 = NA_real_,
    error = error,
    stringsAsFactors = FALSE
  )
}

mvt_standard_dispersion <- function(fit, dat, mu, spec) {
  if (identical(spec$family, "poisson") || identical(spec$family, "binomial")) return(1)
  value <- tryCatch({
    if (requireNamespace("lme4", quietly = TRUE) && inherits(fit, "merMod")) {
      lme4::sigma(fit)^2
    } else {
      summary(fit)$dispersion
    }
  }, error = function(e) NA_real_)
  value <- as.numeric(value)[[1L]]
  if (is.finite(value) && value > 0) return(value)
  residual_var <- stats::var(as.numeric(dat$response) - as.numeric(mu), na.rm = TRUE)
  if (is.finite(residual_var) && residual_var > 0) residual_var else 1
}

mvt_simulate_standard_independent_response <- function(fit, data, spec, nsim = 50L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mu <- tryCatch({
    if (requireNamespace("lme4", quietly = TRUE) && inherits(fit, "merMod")) {
      stats::predict(fit, newdata = data, type = "response", re.form = NA, allow.new.levels = TRUE)
    } else {
      stats::predict(fit, newdata = data, type = "response")
    }
  }, error = function(e) e)
  if (inherits(mu, "error")) {
    stop("Could not predict fitted means for independent variogram simulation: ", conditionMessage(mu), call. = FALSE)
  }
  mu <- as.numeric(mu)
  if (length(mu) != nrow(data) || !any(is.finite(mu))) {
    stop("Independent variogram simulation received invalid fitted means.", call. = FALSE)
  }
  phi <- mvt_standard_dispersion(fit, data, mu, spec)
  out <- matrix(NA_real_, nrow = nrow(data), ncol = nsim)
  for (sim_id in seq_len(nsim)) {
    out[, sim_id] <- switch(
      spec$family,
      gaussian = stats::rnorm(nrow(data), mean = mu, sd = sqrt(phi)),
      poisson = stats::rpois(nrow(data), lambda = pmax(mu, 1e-8)),
      gamma = {
        shape <- 1 / pmax(phi, 1e-8)
        stats::rgamma(nrow(data), shape = shape, scale = pmax(mu, 1e-8) / shape)
      },
      binomial = stats::rbinom(nrow(data), size = 1L, prob = pmin(pmax(mu, 1e-8), 1 - 1e-8)),
      stop("Independent variogram simulation is not implemented for family: ", spec$family, call. = FALSE)
    )
  }
  out
}

mvt_is_standard_variogram_fit <- function(fit) {
  inherits(fit, c("lm", "glm", "geeglm")) ||
    (requireNamespace("lme4", quietly = TRUE) && inherits(fit, "merMod"))
}

mvt_variogram_score <- function(dat, fit, method, row, spec, nsim = 50L, seed = NULL) {
  if (!requireNamespace("scoringRules", quietly = TRUE)) {
    return(mvt_variogram_empty(row, spec, method, "Package 'scoringRules' is not installed."))
  }
  if (inherits(fit, "gamlss.longitudinal")) {
    sim <- tryCatch(stats::simulate(fit, nsim = nsim, seed = seed), error = function(e) e)
  } else if (inherits(fit, "mvt_gamcopula_fit")) {
    sim <- tryCatch(mvt_simulate_gamcopula_response(fit, dat, spec, nsim = nsim, seed = seed), error = function(e) e)
  } else if (mvt_is_standard_variogram_fit(fit)) {
    sim <- tryCatch(mvt_simulate_standard_independent_response(fit, dat, spec, nsim = nsim, seed = seed), error = function(e) e)
  } else {
    return(mvt_variogram_empty(row, spec, method, "Joint simulation variogram scoring is not implemented for this method."))
  }
  if (inherits(sim, "error")) return(mvt_variogram_empty(row, spec, method, conditionMessage(sim)))
  if (nrow(as.data.frame(sim)) != nrow(dat)) {
    return(mvt_variogram_empty(row, spec, method, "Simulation output row count did not match the scored data."))
  }
  subject_scores <- list(p05 = numeric(), p2 = numeric())
  for (id in unique(dat$subject)) {
    idx <- which(dat$subject == id)
    idx <- idx[order(dat$time_index[idx])]
    y <- as.numeric(dat$response[idx])
    draws <- as.matrix(sim[idx, , drop = FALSE])
    ok <- is.finite(y) & rowSums(is.finite(draws)) == ncol(draws)
    if (sum(ok) < 2L) next
    subject_scores$p05 <- c(subject_scores$p05, scoringRules::vs_sample(y[ok], draws[ok, , drop = FALSE], p = 0.5))
    subject_scores$p2 <- c(subject_scores$p2, scoringRules::vs_sample(y[ok], draws[ok, , drop = FALSE], p = 2))
  }
  out <- mvt_variogram_empty(row, spec, method)
  out$variogram_nsim <- nsim
  out$variogram_score_p05 <- mean(subject_scores$p05, na.rm = TRUE)
  out$variogram_score_p2 <- mean(subject_scores$p2, na.rm = TRUE)
  out$error <- NA_character_
  out
}

mvt_status_from_results <- function(results) {
  if (!is.data.frame(results) || nrow(results) == 0L) return(data.frame())
  out <- results
  if (!"available" %in% names(out)) out$available <- TRUE
  if (!"success" %in% names(out)) out$success <- NA
  if (!"warning" %in% names(out)) out$warning <- NA_character_
  if (!"error" %in% names(out)) out$error <- NA_character_
  if (!"status_class" %in% names(out)) out$status_class <- NA_character_
  missing_status <- is.na(out$status_class) | !nzchar(trimws(as.character(out$status_class)))
  if (any(missing_status)) {
    out$status_class[missing_status] <- vapply(which(missing_status), function(i) {
      success <- isTRUE(out$success[[i]]) || identical(out$success[[i]], "TRUE")
      mvt_fit_status_class(success, out$error[[i]], out$warning[[i]])
    }, character(1))
  }
  if (!"failure_reason_short" %in% names(out)) {
    out$failure_reason_short <- mapply(
      mvt_failure_reason_short,
      out$success,
      out$status_class,
      out$error,
      out$warning,
      USE.NAMES = FALSE
    )
  }
  keep <- intersect(
    c(
      "case_id", "scenario", "generator", "dependence", "correlation_level",
      "n_time", "n_subject", "total_rows", "family", "gamlss_family", "rep",
      "method", "available", "success", "status_class", "elapsed_sec", "nobs",
      "failure_reason_short", "warning", "error"
    ),
    names(out)
  )
  out[, keep, drop = FALSE]
}

mvt_run_case <- function(row, seed_base = 20260818L, require_gamcopula = TRUE) {
  if (isTRUE(require_gamcopula)) {
    mvt_require_namespaces("gamCopula", strict = TRUE)
  }
  families <- mvt_family_specs(include_special = TRUE)
  deps <- mvt_dependence_specs(include_appendix = TRUE)
  spec <- families[[row$family_name]]
  dep <- mvt_resize_dependence(deps[[row$dependence_name]], row$n_time)
  seed <- seed_base + match(row$family_name, names(families)) * 100000L +
    match(row$dependence_name, names(deps)) * 1000L + row$n_time * 10L + row$rep

  sim_capture <- mvt_elapsed_capture(mvt_simulate_case(row, seed = seed))
  if (inherits(sim_capture$value, "error")) {
    status <- mvt_method_result_row(row, spec, "simulation", sim_capture)
    return(list(
      fit_status = status,
      benchmark_results = status,
      coefficient_results = data.frame(),
      dependence_recovery = data.frame(),
      variogram_scores = data.frame(),
      runtime = status[c("case_id", "method", "elapsed_sec", "success", "error")]
    ))
  }
  dat <- sim_capture$value
  model_dat <- dat[c("subject", "time_index", "time", "response", "x", "z")]
  y <- dat$response

  result_rows <- list()
  coef_rows <- list()
  dep_rows <- list()
  vario_rows <- list()
  status_rows <- list()
  active_comparators <- mvt_active_comparators()

  if (length(mvt_standard_comparators(active_comparators)) > 0L) {
    standard <- mvt_capture(mvt_fit_standard_models(dat, spec, row, dep))
    if (!inherits(standard$value, "error")) {
      result_rows[[length(result_rows) + 1L]] <- standard$value$results
      coef_rows[[length(coef_rows) + 1L]] <- mvt_annotate_coefficients(standard$value$coefficients, spec)
      standard_fits <- standard$value$fits
      standard_fits <- standard_fits[names(standard_fits) %in% mvt_standard_comparators(active_comparators)]
      for (fit_name in names(standard_fits)) {
        if (!is.null(standard_fits[[fit_name]])) {
          vario_rows[[length(vario_rows) + 1L]] <- mvt_variogram_score(
            model_dat, standard_fits[[fit_name]], fit_name, row, spec,
            nsim = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM", 20L),
            seed = seed + 7000L + match(fit_name, names(standard_fits))
          )
        }
      }
    } else {
      result_rows[[length(result_rows) + 1L]] <- mvt_method_result_row(
        row,
        spec,
        "standard_models",
        list(value = standard$value, warnings = standard$warnings, elapsed_sec = NA_real_),
        extra = list(y = y)
      )
    }
  }

  for (corstr in intersect(dep$gee_correlations, mvt_gee_comparators(active_comparators))) {
    gee <- mvt_run_one_gee(dat, spec, row, corstr)
    result_rows[[length(result_rows) + 1L]] <- gee$results
    coef_rows[[length(coef_rows) + 1L]] <- mvt_annotate_coefficients(gee$coefficients, spec)
    if (!is.null(gee$fit)) {
      vario_rows[[length(vario_rows) + 1L]] <- mvt_variogram_score(
        model_dat, gee$fit, paste0("gee_", corstr), row, spec,
        nsim = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM", 20L),
        seed = seed + 7200L + match(corstr, dep$gee_correlations)
      )
    }
  }

  if ("glmm_slope" %in% active_comparators) {
    glmm_slope <- mvt_run_glmm_slope(dat, spec, row)
    result_rows[[length(result_rows) + 1L]] <- glmm_slope$results
    coef_rows[[length(coef_rows) + 1L]] <- mvt_annotate_coefficients(glmm_slope$coefficients, spec)
    if (!is.null(glmm_slope$fit)) {
      vario_rows[[length(vario_rows) + 1L]] <- mvt_variogram_score(
        model_dat, glmm_slope$fit, "glmm_slope", row, spec,
        nsim = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM", 20L),
        seed = seed + 7400L
      )
    }
  }

  if ("gamlss.longitudinal" %in% active_comparators) {
    gl <- mvt_run_fit_with_timeout(
      "mvt_fit_gamlss_longitudinal_for_case",
      args = list(dat = dat, row = row),
      timeout = mvt_env_num("GAMLSS_LONGITUDINAL_MVT_PRIMARY_MAX_ELAPSED_SEC", 180)
    )
    extra_gl <- list(y = y)
    pred_gl <- NULL
    if (inherits(gl$value, "gamlss.longitudinal")) {
      gl_params <- mvt_predict_gamlss_params(gl$value, model_dat, spec)
      pred_gl <- gl_params$mu
      extra_gl <- c(extra_gl, as.list(mvt_distribution_metrics(dat, gl_params, spec)))
      extra_gl$logLik <- mvt_loglik_value(gl$value)
      extra_gl$logLik_df <- mvt_loglik_df(gl$value)
      extra_gl$AIC <- mvt_aic_value(gl$value)
      extra_gl$BIC <- mvt_bic_value(gl$value)
    } else {
      pred_gl <- NULL
    }
    gl_row <- mvt_method_result_row(row, spec, "gamlss.longitudinal", gl, pred = pred_gl, extra = extra_gl)
    result_rows[[length(result_rows) + 1L]] <- gl_row
    status_rows[[length(status_rows) + 1L]] <- gl_row
    if (inherits(gl$value, "gamlss.longitudinal")) {
      coef_rows[[length(coef_rows) + 1L]] <- mvt_coef_table_one(gl$value, "gamlss.longitudinal", row, spec)
      dep_rows[[length(dep_rows) + 1L]] <- mvt_dependence_recovery_row(dat, gl$value, "gamlss.longitudinal", row, spec)
      vario_rows[[length(vario_rows) + 1L]] <- mvt_variogram_score(
        model_dat, gl$value, "gamlss.longitudinal", row, spec,
        nsim = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM", 20L),
        seed = seed + 9000L
      )
    }
  }

  gamcopula_specs <- list(
    gamCopula_markov = "mvt_fit_gamcopula_for_case",
    gamCopula_vine_simplified = "mvt_fit_gamcopula_vine_simplified_for_case",
    gamCopula_vine = "mvt_fit_gamcopula_vine_for_case"
  )
  for (gc_method in intersect(names(gamcopula_specs), active_comparators)) {
    timeout_name <- if (gc_method %in% c("gamCopula_vine", "gamCopula_vine_simplified")) {
      "GAMLSS_LONGITUDINAL_MVT_GAMCOPULA_VINE_TIMEOUT_SEC"
    } else {
      "GAMLSS_LONGITUDINAL_MVT_GAMCOPULA_MARKOV_TIMEOUT_SEC"
    }
    gc_fit <- mvt_run_fit_with_timeout(
      gamcopula_specs[[gc_method]],
      args = list(dat = dat, row = row),
      timeout = mvt_env_num(
        timeout_name,
        mvt_env_num(
          "GAMLSS_LONGITUDINAL_MVT_GAMCOPULA_TIMEOUT_SEC",
          mvt_env_num("GAMLSS_LONGITUDINAL_MVT_PRIMARY_MAX_ELAPSED_SEC", 180)
        )
      )
    )
    extra_gc <- list(y = y)
    pred_gc <- NULL
    if (inherits(gc_fit$value, "mvt_gamcopula_fit")) {
      gc_params <- mvt_predict_gamcopula_params(gc_fit$value, dat, spec)
      pred_gc <- gc_params$mu
      extra_gc <- c(extra_gc, as.list(mvt_distribution_metrics(dat, gc_params, spec)))
      extra_gc$logLik <- mvt_loglik_gamcopula(gc_fit$value)
      extra_gc$logLik_df <- mvt_loglik_df_gamcopula(gc_fit$value)
      extra_gc$margin_df <- mvt_loglik_df(gc_fit$value$margin_fit)
      extra_gc$copula_df <- mvt_gamcopula_copula_df(gc_fit$value)
      extra_gc$AIC <- mvt_aic_gamcopula(gc_fit$value)
      extra_gc$BIC <- mvt_bic_gamcopula(gc_fit$value)
      extra_gc$copula_engine <- gc_fit$value$vine_engine %||% "gamCopula"
      extra_gc$copula_simplified <- gc_fit$value$vine_simplified %||% NA
    }
    gc_row <- mvt_method_result_row(row, spec, gc_method, gc_fit, pred = pred_gc, extra = extra_gc)
    result_rows[[length(result_rows) + 1L]] <- gc_row
    status_rows[[length(status_rows) + 1L]] <- gc_row
    if (inherits(gc_fit$value, "mvt_gamcopula_fit")) {
      coef_rows[[length(coef_rows) + 1L]] <- mvt_coef_table_one(gc_fit$value, gc_method, row, spec)
      dep_rows[[length(dep_rows) + 1L]] <- mvt_dependence_recovery_row(dat, gc_fit$value, gc_method, row, spec)
      vario_rows[[length(vario_rows) + 1L]] <- mvt_variogram_score(
        model_dat, gc_fit$value, gc_method, row, spec,
        nsim = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM", 20L),
        seed = seed + if (identical(gc_method, "gamCopula_vine")) 9700L else if (identical(gc_method, "gamCopula_vine_simplified")) 9650L else 9500L
      )
    }
  }

  benchmark_results <- mvt_bind_rows_fill(result_rows)
  fit_status <- mvt_status_from_results(benchmark_results)
  coefficient_results <- mvt_bind_rows_fill(coef_rows)
  dependence_recovery <- mvt_bind_rows_fill(dep_rows)
  variogram_scores <- mvt_bind_rows_fill(vario_rows)
  if (nrow(benchmark_results) > 0L && "method" %in% names(benchmark_results)) {
    scored_methods <- if (nrow(variogram_scores) > 0L && "method" %in% names(variogram_scores)) {
      unique(variogram_scores$method)
    } else {
      character()
    }
    missing_vario <- setdiff(unique(benchmark_results$method), scored_methods)
    if (length(missing_vario) > 0L) {
      empty_vario <- lapply(
        missing_vario,
        function(method) mvt_variogram_empty(
          row,
          spec,
          method,
          "Joint simulation variogram scoring is not implemented for this method in the paper workflow."
        )
      )
      variogram_scores <- mvt_bind_rows_fill(variogram_scores, empty_vario)
    }
  }
  runtime <- benchmark_results[, intersect(c("case_id", "scenario", "family", "method", "elapsed_sec", "success", "error"), names(benchmark_results)), drop = FALSE]

  list(
    fit_status = fit_status,
    benchmark_results = benchmark_results,
    coefficient_results = coefficient_results,
    dependence_recovery = dependence_recovery,
    variogram_scores = variogram_scores,
    runtime = runtime
  )
}

mvt_default_run_dir <- function(stage) {
  file.path(mvt_output_root, paste0(stage, "_", mvt_timestamp()))
}

mvt_read_result_table <- function(run_dir, name) {
  candidates <- file.path(run_dir, paste0(name, c("_by_rep.csv", "_checkpoint.csv")))
  path <- candidates[file.exists(candidates)]
  if (length(path) == 0L) return(data.frame())
  path <- path[[1L]]
  if (is.na(file.info(path)$size) || file.info(path)$size == 0L) return(data.frame())
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

mvt_read_existing_results <- function(run_dir) {
  stats::setNames(
    lapply(mvt_result_names(), function(name) mvt_read_result_table(run_dir, name)),
    mvt_result_names()
  )
}

mvt_completed_case_ids <- function(existing) {
  sources <- existing[c("fit_status", "benchmark_results", "runtime")]
  ids <- unique(unlist(lapply(sources, function(x) {
    if (is.data.frame(x) && "case_id" %in% names(x)) unique(x$case_id) else character()
  }), use.names = FALSE))
  ids[nzchar(ids)]
}

mvt_preflight_check <- function(check, status, detail = "", n = NA_integer_) {
  data.frame(
    check = check,
    status = status,
    n = n,
    detail = detail,
    stringsAsFactors = FALSE
  )
}

mvt_write_preflight <- function(grid, run_dir, require_gamcopula = TRUE, resume = TRUE) {
  active <- mvt_active_comparators()
  packages <- c("gamlss", "gamlss.dist", "mvtnorm", "VineCopula", "callr")
  if (isTRUE(require_gamcopula) || any(grepl("^gamCopula", active))) packages <- c(packages, "gamCopula")
  if (any(grepl("^gee_", active))) packages <- c(packages, "geepack")
  if (any(c("glmm", "glmm_slope") %in% active)) packages <- c(packages, "lme4")
  packages <- unique(packages)
  package_status <- mvt_require_namespaces(packages, strict = FALSE)
  existing_outputs <- file.exists(file.path(run_dir, paste0(mvt_result_names(), "_by_rep.csv"))) |
    file.exists(file.path(run_dir, paste0(mvt_result_names(), "_checkpoint.csv")))
  row_cap_ok <- nrow(grid) > 0L && "total_rows" %in% names(grid) && all(grid$total_rows <= 5000L)
  unstructured_t50 <- "gee_unstructured" %in% active && "n_time" %in% names(grid) && any(grid$n_time >= 50L)
  unstructured_timeout <- mvt_env_num("GAMLSS_LONGITUDINAL_MVT_GEE_UNSTRUCTURED_TIMEOUT_SEC", if (unstructured_t50) 30 else Inf)
  checks <- list(
    mvt_preflight_check(
      "grid_nonempty",
      if (nrow(grid) > 0L) "pass" else "fail",
      paste("cases", nrow(grid)),
      nrow(grid)
    ),
    mvt_preflight_check(
      "total_rows_within_cap",
      if (isTRUE(row_cap_ok)) "pass" else "fail",
      paste("max total_rows", if ("total_rows" %in% names(grid) && nrow(grid) > 0L) max(grid$total_rows, na.rm = TRUE) else NA),
      if ("total_rows" %in% names(grid) && nrow(grid) > 0L) max(grid$total_rows, na.rm = TRUE) else NA_integer_
    ),
    mvt_preflight_check(
      "active_comparators_valid",
      if (length(active) > 0L) "pass" else "fail",
      paste(active, collapse = ", "),
      length(active)
    ),
    mvt_preflight_check(
      "required_packages_available",
      if (all(package_status$available)) "pass" else "fail",
      paste(package_status$package, package_status$available, sep = "=", collapse = ", "),
      sum(package_status$available)
    ),
    mvt_preflight_check(
      "resume_configuration",
      if (isTRUE(resume) || !any(existing_outputs)) "pass" else "warn",
      paste("resume", resume, "| existing result/checkpoint files", sum(existing_outputs)),
      sum(existing_outputs)
    ),
    mvt_preflight_check(
      "gee_unstructured_timeout",
      if (!unstructured_t50 || is.finite(unstructured_timeout)) "pass" else "warn",
      paste("T>=50 unstructured GEE", unstructured_t50, "| timeout", unstructured_timeout),
      if (is.finite(unstructured_timeout)) unstructured_timeout else NA_integer_
    )
  )
  out <- mvt_bind_rows_fill(checks)
  mvt_write_csv(out, file.path(run_dir, "preflight_checks.csv"))
  lines <- c(
    "# Multivariate Simulation Preflight",
    "",
    paste("Run directory:", normalizePath(run_dir, winslash = "/", mustWork = FALSE)),
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Cases:", nrow(grid)),
    paste("Active comparators:", paste(active, collapse = ", ")),
    "",
    "## Checks",
    "",
    paste0("- ", out$check, " [", out$status, "]: ", out$detail)
  )
  writeLines(lines, file.path(run_dir, "preflight_checks.md"), useBytes = TRUE)
  invisible(out)
}

mvt_write_run_metadata <- function(run_dir, grid, seed_base, require_gamcopula, resume) {
  metadata <- data.frame(
    name = c(
      "created_at",
      "repo_root",
      "seed_base",
      "require_gamcopula",
      "resume",
      "active_comparators",
      "n_cases",
      "r_version",
      "platform"
    ),
    value = c(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      mvt_repo_root,
      as.character(seed_base),
      as.character(require_gamcopula),
      as.character(resume),
      paste(mvt_active_comparators(), collapse = ","),
      as.character(nrow(grid)),
      R.version.string,
      R.version$platform
    ),
    stringsAsFactors = FALSE
  )
  env_names <- sort(grep("^GAMLSS_LONGITUDINAL_MVT_", names(Sys.getenv()), value = TRUE))
  if (length(env_names) > 0L) {
    metadata <- rbind(
      metadata,
      data.frame(
        name = paste0("env.", env_names),
        value = Sys.getenv(env_names),
        stringsAsFactors = FALSE
      )
    )
  }
  mvt_write_csv(metadata, file.path(run_dir, "run_metadata.csv"))
  mvt_write_csv(
    mvt_package_versions(c(
      "gamlss.longitudinal", "gamlss", "gamlss.dist", "gamCopula",
      "VineCopula", "mvtnorm", "geepack", "lme4", "mgcv",
      "callr", "scoringRules", "ggplot2"
    )),
    file.path(run_dir, "package_versions.csv")
  )
  writeLines(capture.output(utils::sessionInfo()), file.path(run_dir, "session_info.txt"), useBytes = TRUE)
  invisible(run_dir)
}

mvt_run_grid <- function(grid, run_dir, seed_base = 20260818L, checkpoint_every = 5L, require_gamcopula = TRUE) {
  dir.create(mvt_output_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  resume <- mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_RESUME", TRUE)
  existing <- if (isTRUE(resume)) mvt_read_existing_results(run_dir) else stats::setNames(vector("list", length(mvt_result_names())), mvt_result_names())
  completed <- mvt_completed_case_ids(existing)
  mvt_write_csv(grid, file.path(run_dir, "scenario_grid.csv"))
  mvt_write_preflight(grid, run_dir, require_gamcopula = require_gamcopula, resume = resume)
  mvt_write_run_metadata(run_dir, grid, seed_base, require_gamcopula, resume)
  rows <- lapply(mvt_result_names(), function(name) {
    current <- existing[[name]]
    if (is.data.frame(current) && nrow(current) > 0L) list(current) else list()
  })
  names(rows) <- mvt_result_names()
  for (i in seq_len(nrow(grid))) {
    if (grid$case_id[[i]] %in% completed) {
      message("[", i, "/", nrow(grid), "] ", grid$case_id[[i]], " (resume skip)")
      next
    }
    message("[", i, "/", nrow(grid), "] ", grid$case_id[[i]])
    out <- mvt_run_case(grid[i, , drop = FALSE], seed_base = seed_base, require_gamcopula = require_gamcopula)
    for (nm in names(rows)) rows[[nm]][[length(rows[[nm]]) + 1L]] <- out[[nm]]
    if (checkpoint_every > 0L && i %% checkpoint_every == 0L) {
      for (nm in names(rows)) {
        mvt_write_csv(mvt_bind_rows_fill(rows[[nm]]), file.path(run_dir, paste0(nm, "_checkpoint.csv")))
      }
    }
  }
  for (nm in names(rows)) {
    mvt_write_csv(mvt_bind_rows_fill(rows[[nm]]), file.path(run_dir, paste0(nm, "_by_rep.csv")))
  }
  run_dir_norm <- normalizePath(run_dir, winslash = "/", mustWork = FALSE)
  output_root_norm <- normalizePath(mvt_output_root, winslash = "/", mustWork = FALSE)
  if (startsWith(run_dir_norm, output_root_norm)) {
    writeLines(run_dir, file.path(mvt_output_root, "latest_run_dir.txt"), useBytes = TRUE)
  }
  invisible(run_dir)
}

mvt_read_run_dir <- function() {
  run_dir <- mvt_env("GAMLSS_LONGITUDINAL_MVT_RUN_DIR", "")
  if (!nzchar(run_dir)) {
    latest <- file.path(mvt_output_root, "latest_run_dir.txt")
    if (!file.exists(latest)) {
      stop("No run directory supplied and latest_run_dir.txt does not exist.", call. = FALSE)
    }
    run_dir <- readLines(latest, warn = FALSE)[[1L]]
  }
  normalizePath(run_dir, winslash = "/", mustWork = TRUE)
}

mvt_merge_run_shards <- function(source_dirs, target_dir) {
  source_dirs <- normalizePath(source_dirs[nzchar(source_dirs)], winslash = "/", mustWork = TRUE)
  if (length(source_dirs) == 0L) stop("No source shard directories supplied.", call. = FALSE)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

  read_table <- function(dir, file) {
    path <- file.path(dir, file)
    out <- mvt_read_optional_csv(path)
    if (nrow(out) == 0L) return(data.frame())
    out$source_shard <- basename(dir)
    out
  }
  read_result_table <- function(dir, name) {
    files <- paste0(name, c("_by_rep.csv", "_checkpoint.csv"))
    paths <- file.path(dir, files)
    path <- paths[file.exists(paths) & !is.na(file.info(paths)$size) & file.info(paths)$size > 0L]
    if (length(path) == 0L) return(data.frame())
    out <- mvt_read_optional_csv(path[[1L]])
    if (nrow(out) == 0L) return(data.frame())
    out$source_shard <- basename(dir)
    out
  }

  grid <- unique(mvt_bind_rows_fill(lapply(source_dirs, read_table, file = "scenario_grid.csv")))
  if ("source_shard" %in% names(grid)) grid$source_shard <- NULL
  if (nrow(grid) == 0L) stop("No scenario_grid.csv rows found in shard directories.", call. = FALSE)
  grid <- grid[order(grid$n_time, grid$dependence_name, grid$family_name, grid$rep), , drop = FALSE]
  mvt_write_csv(grid, file.path(target_dir, "scenario_grid.csv"))

  for (name in mvt_result_names()) {
    combined <- mvt_bind_rows_fill(lapply(source_dirs, read_result_table, name = name))
    if ("source_shard" %in% names(combined)) combined$source_shard <- NULL
    if (nrow(combined) > 0L) {
      key_cols <- intersect(c("case_id", "method", "parameter", "term", "dependence_scope"), names(combined))
      if (length(key_cols) > 0L) combined <- unique(combined)
      mvt_write_csv(combined, file.path(target_dir, paste0(name, "_by_rep.csv")))
    }
  }

  preflight <- mvt_bind_rows_fill(lapply(source_dirs, read_table, file = "preflight_checks.csv"))
  if (nrow(preflight) > 0L) mvt_write_csv(preflight, file.path(target_dir, "preflight_checks.csv"))
  metadata <- mvt_bind_rows_fill(lapply(source_dirs, read_table, file = "run_metadata.csv"))
  if (nrow(metadata) > 0L) mvt_write_csv(metadata, file.path(target_dir, "run_metadata.csv"))
  versions <- unique(mvt_bind_rows_fill(lapply(source_dirs, read_table, file = "package_versions.csv")))
  if ("source_shard" %in% names(versions)) versions$source_shard <- NULL
  if (nrow(versions) > 0L) mvt_write_csv(versions, file.path(target_dir, "package_versions.csv"))
  writeLines(
    c(
      "Merged multivariate longitudinal simulation shards.",
      paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      "Source shards:",
      paste0("- ", source_dirs)
    ),
    file.path(target_dir, "session_info.txt"),
    useBytes = TRUE
  )
  writeLines(
    c(
      "# Merged Preflight",
      "",
      paste("Source shards:", length(source_dirs)),
      "",
      paste0("- ", basename(source_dirs))
    ),
    file.path(target_dir, "preflight_checks.md"),
    useBytes = TRUE
  )
  mvt_summarise_results(target_dir)
  mvt_write_artifact_manifest(target_dir)
  invisible(target_dir)
}

mvt_summarise_results <- function(run_dir = mvt_read_run_dir()) {
  results <- mvt_read_optional_csv(file.path(run_dir, "benchmark_results_by_rep.csv"))
  coefs <- mvt_read_optional_csv(file.path(run_dir, "coefficient_results_by_rep.csv"))
  if (nrow(results) == 0L) stop("benchmark_results_by_rep.csv is missing or empty.", call. = FALSE)
  if (nrow(coefs) == 0L) stop("coefficient_results_by_rep.csv is missing or empty.", call. = FALSE)
  dep <- mvt_read_optional_csv(file.path(run_dir, "dependence_recovery_by_rep.csv"))
  vario <- mvt_read_optional_csv(file.path(run_dir, "variogram_scores_by_rep.csv"))

  result_metrics <- intersect(
    c("mae", "rmse", "benchmark_mean_rmse", "benchmark_neg_log_score", "benchmark_pit_mean_abs_error", "logLik", "logLik_df", "AIC", "BIC", "elapsed_sec"),
    names(results)
  )
  result_summary <- mvt_group_metric_summary(
    results,
    group_cols = intersect(c("scenario", "n_time", "family", "method"), names(results)),
    metrics = result_metrics
  )
  coef_summary <- if (nrow(coefs) > 0L) {
    mvt_group_metric_summary(
      coefs,
      group_cols = intersect(c("scenario", "n_time", "family", "method", "term"), names(coefs)),
      metrics = intersect(c("bias", "ci_width"), names(coefs)),
      extra_bool = intersect(c("ci_covers_truth", "false_positive"), names(coefs))
    )
  } else {
    data.frame()
  }
  dep_summary <- if (nrow(dep) > 0L) {
    mvt_group_metric_summary(
      dep,
      group_cols = intersect(c("scenario", "n_time", "family", "method"), names(dep)),
      metrics = intersect(c("theta_mae", "theta_rmse", "tau_mae", "tau_rmse"), names(dep))
    )
  } else {
    data.frame()
  }
  vario_summary <- if (nrow(vario) > 0L) {
    mvt_group_metric_summary(
      vario,
      group_cols = intersect(c("scenario", "n_time", "family", "method"), names(vario)),
      metrics = intersect(c("variogram_score_p05", "variogram_score_p2"), names(vario))
    )
  } else {
    data.frame()
  }
  mvt_write_csv(result_summary, file.path(run_dir, "benchmark_summary.csv"))
  mvt_write_csv(coef_summary, file.path(run_dir, "coefficient_summary.csv"))
  mvt_write_csv(dep_summary, file.path(run_dir, "dependence_recovery_summary.csv"))
  mvt_write_csv(vario_summary, file.path(run_dir, "variogram_summary.csv"))
  status_path <- file.path(run_dir, "fit_status_by_rep.csv")
  if (file.exists(status_path)) {
    status <- utils::read.csv(status_path, stringsAsFactors = FALSE, check.names = FALSE)
    status <- mvt_status_from_results(status)
    mvt_write_csv(status, status_path)
    mvt_write_csv(mvt_case_method_completion_summary(status), file.path(run_dir, "case_method_completion_summary.csv"))
  }
  feasibility <- mvt_write_pilot_feasibility(run_dir)
  invisible(list(
    results = result_summary,
    coefficients = coef_summary,
    dependence = dep_summary,
    variogram = vario_summary,
    feasibility = feasibility
  ))
}

mvt_status_success <- function(x) {
  x %in% c(TRUE, "TRUE", "true", "1", 1L)
}

mvt_collapse_nonempty <- function(x, empty = "") {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) empty else paste(x, collapse = "; ")
}

mvt_case_method_completion_summary <- function(status) {
  if (!is.data.frame(status) || nrow(status) == 0L || !"case_id" %in% names(status)) return(data.frame())
  if (!"success" %in% names(status)) status$success <- NA
  if (!"status_class" %in% names(status)) status$status_class <- NA_character_
  if (!"method" %in% names(status)) status$method <- NA_character_
  if (!"failure_reason_short" %in% names(status)) {
    if (!"error" %in% names(status)) status$error <- NA_character_
    if (!"warning" %in% names(status)) status$warning <- NA_character_
    status$failure_reason_short <- mapply(
      mvt_failure_reason_short,
      status$success,
      status$status_class,
      status$error,
      status$warning,
      USE.NAMES = FALSE
    )
  }
  group_cols <- intersect(
    c(
      "case_id", "scenario", "generator", "dependence", "correlation_level",
      "n_time", "n_subject", "total_rows", "family", "gamlss_family", "rep"
    ),
    names(status)
  )
  groups <- unique(status[group_cols])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    idx <- rep(TRUE, nrow(status))
    for (col in group_cols) idx <- idx & status[[col]] == groups[[col]][i]
    sub <- status[idx, , drop = FALSE]
    success <- mvt_status_success(sub$success)
    status_class <- as.character(sub$status_class)
    failed <- !success
    warned <- success & status_class == "warning"
    attention <- failed | warned
    attention[is.na(attention)] <- FALSE
    failure_pairs <- if (any(attention)) {
      paste0(sub$method[attention], "=", sub$failure_reason_short[attention])
    } else {
      character()
    }
    rows[[i]] <- cbind(
      groups[i, , drop = FALSE],
      data.frame(
        n_methods_attempted = length(unique(sub$method)),
        n_methods_success = sum(success, na.rm = TRUE),
        n_methods_ok = sum(status_class == "ok", na.rm = TRUE),
        n_methods_warning = sum(warned, na.rm = TRUE),
        n_methods_timeout = sum(status_class == "timeout", na.rm = TRUE),
        n_methods_error = sum(status_class == "error", na.rm = TRUE),
        effective_methods_completed = mvt_collapse_nonempty(sub$method[success], empty = "none"),
        methods_with_warnings = mvt_collapse_nonempty(sub$method[warned], empty = "none"),
        methods_failed = mvt_collapse_nonempty(sub$method[failed], empty = "none"),
        failure_reasons_short = mvt_collapse_nonempty(failure_pairs, empty = "none"),
        stringsAsFactors = FALSE
      )
    )
  }
  mvt_bind_rows_fill(rows)
}

mvt_status_counts <- function(status, group_cols) {
  if (nrow(status) == 0L || length(group_cols) == 0L) return(data.frame())
  groups <- unique(status[group_cols])
  rows <- list()
  for (i in seq_len(nrow(groups))) {
    idx <- rep(TRUE, nrow(status))
    for (col in group_cols) idx <- idx & status[[col]] == groups[[col]][i]
    sub <- status[idx, , drop = FALSE]
    elapsed <- suppressWarnings(as.numeric(sub$elapsed_sec))
    success <- mvt_status_success(sub$success)
    status_class <- as.character(sub$status_class %||% rep(NA_character_, nrow(sub)))
    rows[[i]] <- cbind(
      groups[i, , drop = FALSE],
      data.frame(
        n_attempts = nrow(sub),
        n_success = sum(success, na.rm = TRUE),
        n_warning = sum(status_class == "warning", na.rm = TRUE),
        n_timeout = sum(status_class == "timeout", na.rm = TRUE),
        n_error = sum(status_class == "error", na.rm = TRUE),
        success_rate = mean(success, na.rm = TRUE),
        warning_rate = mean(status_class == "warning", na.rm = TRUE),
        timeout_rate = mean(status_class == "timeout", na.rm = TRUE),
        error_rate = mean(status_class == "error", na.rm = TRUE),
        mean_elapsed_sec = if (any(is.finite(elapsed))) mean(elapsed[is.finite(elapsed)]) else NA_real_,
        median_elapsed_sec = if (any(is.finite(elapsed))) stats::median(elapsed[is.finite(elapsed)]) else NA_real_,
        p90_elapsed_sec = if (sum(is.finite(elapsed)) > 0L) as.numeric(stats::quantile(elapsed[is.finite(elapsed)], 0.9, names = FALSE)) else NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }
  mvt_bind_rows_fill(rows)
}

mvt_write_pilot_feasibility <- function(run_dir = mvt_read_run_dir(), target_reps = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_TARGET_REPS", 100L)) {
  status_path <- file.path(run_dir, "fit_status_by_rep.csv")
  grid_path <- file.path(run_dir, "scenario_grid.csv")
  if (!file.exists(status_path) || !file.exists(grid_path)) return(data.frame())
  status <- utils::read.csv(status_path, stringsAsFactors = FALSE, check.names = FALSE)
  grid <- utils::read.csv(grid_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(status) == 0L || nrow(grid) == 0L) return(data.frame())

  method_summary <- mvt_status_counts(
    status,
    group_cols = intersect(c("method", "family", "n_time", "dependence"), names(status))
  )
  scenario_summary <- mvt_status_counts(
    status,
    group_cols = intersect(c("case_id", "scenario", "family", "n_time", "dependence"), names(status))
  )
  overall_method <- mvt_status_counts(status, group_cols = intersect("method", names(status)))

  n_pilot_reps <- length(unique(grid$rep))
  projection_multiplier <- if (n_pilot_reps > 0L) target_reps / n_pilot_reps else NA_real_
  if (nrow(method_summary) > 0L) {
    method_summary$target_reps <- target_reps
    method_summary$projected_elapsed_sec <- method_summary$mean_elapsed_sec * projection_multiplier
  }
  if (nrow(scenario_summary) > 0L) {
    scenario_summary$target_reps <- target_reps
    scenario_summary$projected_elapsed_sec <- scenario_summary$mean_elapsed_sec * projection_multiplier
  }
  if (nrow(overall_method) > 0L) {
    overall_method$target_reps <- target_reps
    overall_method$projected_elapsed_sec <- overall_method$mean_elapsed_sec * projection_multiplier
    overall_method$recommendation <- ifelse(
      overall_method$success_rate >= 0.9 & overall_method$error_rate <= 0.1,
      "main_or_appendix_candidate",
      ifelse(overall_method$timeout_rate > 0, "report_feasibility_or_appendix", "investigate_before_main")
    )
  }

  mvt_write_csv(method_summary, file.path(run_dir, "pilot_feasibility_by_method.csv"))
  mvt_write_csv(scenario_summary, file.path(run_dir, "pilot_feasibility_by_scenario.csv"))
  mvt_write_csv(overall_method, file.path(run_dir, "pilot_feasibility_overall_method.csv"))

  total_mean_sec <- sum(overall_method$mean_elapsed_sec, na.rm = TRUE)
  projected_core_sec <- total_mean_sec * projection_multiplier
  non_candidate <- overall_method[
    overall_method$recommendation != "main_or_appendix_candidate",
    ,
    drop = FALSE
  ]
  lines <- c(
    "# Pilot Feasibility Report",
    "",
    paste("Run directory:", normalizePath(run_dir, winslash = "/", mustWork = FALSE)),
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Pilot replicates:", n_pilot_reps),
    paste("Target replicates:", target_reps),
    paste("Projected elapsed seconds for target replicates over this scenario set:", round(projected_core_sec, 1)),
    "",
    "## Method Recommendations",
    ""
  )
  if (nrow(overall_method) == 0L) {
    lines <- c(lines, "No method status rows were available.")
  } else {
    lines <- c(
      lines,
      paste0(
        "- ", overall_method$method,
        ": success=", round(overall_method$success_rate, 3),
        ", warning=", round(overall_method$warning_rate, 3),
        ", timeout=", round(overall_method$timeout_rate, 3),
        ", error=", round(overall_method$error_rate, 3),
        ", recommendation=", overall_method$recommendation
      )
    )
  }
  lines <- c(lines, "", "## Methods Needing Attention", "")
  if (nrow(non_candidate) == 0L) {
    lines <- c(lines, "None.")
  } else {
    lines <- c(lines, paste0("- ", non_candidate$method, ": ", non_candidate$recommendation))
  }
  writeLines(lines, file.path(run_dir, "pilot_feasibility.md"), useBytes = TRUE)
  invisible(list(by_method = method_summary, by_scenario = scenario_summary, overall_method = overall_method))
}

mvt_audit_check <- function(check, status, detail = "", n = NA_integer_) {
  data.frame(
    check = check,
    status = status,
    n = n,
    detail = detail,
    stringsAsFactors = FALSE
  )
}

mvt_artifact_type <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv", "tsv")) return("data")
  if (ext %in% c("tex")) return("latex")
  if (ext %in% c("png", "jpg", "jpeg", "pdf", "svg")) return("figure")
  if (ext %in% c("md", "txt")) return("report")
  if (ext %in% c("r", "rmd")) return("source")
  "other"
}

mvt_write_artifact_manifest <- function(run_dir = mvt_read_run_dir()) {
  files <- list.files(run_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  if (length(files) == 0L) {
    manifest <- data.frame()
  } else {
    rel <- gsub("\\\\", "/", substring(normalizePath(files, winslash = "/", mustWork = FALSE), nchar(normalizePath(run_dir, winslash = "/", mustWork = FALSE)) + 2L))
    info <- file.info(files)
    manifest <- data.frame(
      path = rel,
      artifact_type = vapply(rel, mvt_artifact_type, character(1)),
      bytes = as.numeric(info$size),
      modified = format(info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
      md5 = unname(tools::md5sum(files)),
      stringsAsFactors = FALSE
    )
    manifest <- manifest[order(manifest$artifact_type, manifest$path), , drop = FALSE]
  }
  mvt_write_csv(manifest, file.path(run_dir, "artifact_manifest.csv"))
  summary <- if (nrow(manifest) > 0L) {
    aggregate(path ~ artifact_type, data = manifest, FUN = length)
  } else {
    data.frame(artifact_type = character(), path = integer())
  }
  lines <- c(
    "# Artifact Manifest",
    "",
    paste("Run directory:", normalizePath(run_dir, winslash = "/", mustWork = FALSE)),
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Files:", nrow(manifest)),
    paste("Total bytes:", sum(manifest$bytes, na.rm = TRUE)),
    "",
    "## Counts By Type",
    ""
  )
  if (nrow(summary) == 0L) {
    lines <- c(lines, "No files found.")
  } else {
    lines <- c(lines, paste0("- ", summary$artifact_type, ": ", summary$path))
  }
  lines <- c(lines, "", "See `artifact_manifest.csv` for file paths, sizes, timestamps, and MD5 hashes.")
  writeLines(lines, file.path(run_dir, "artifact_manifest.md"), useBytes = TRUE)
  invisible(manifest)
}

mvt_write_review_bundle <- function(run_dir = mvt_read_run_dir()) {
  bundle_dir <- file.path(run_dir, "review_bundle")
  dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- mvt_write_artifact_manifest(run_dir)
  index_path <- file.path(bundle_dir, "README.md")

  build_index <- function(manifest) {
    key_files <- c(
      "preflight_checks.md",
      "review_audit.md",
      "publication_readiness_audit.md",
      "pilot_feasibility.md",
      "artifact_manifest.md",
      "scenario_grid.csv",
      "fit_status_by_rep.csv",
      "benchmark_summary.csv",
      "coefficient_summary.csv",
      "dependence_recovery_summary.csv",
      "case_method_completion_summary.csv",
      "variogram_summary.csv"
    )
    present_key <- key_files[file.exists(file.path(run_dir, key_files))]
    table_files <- if (dir.exists(file.path(run_dir, "paper_tables"))) {
      gsub("\\\\", "/", file.path("paper_tables", list.files(file.path(run_dir, "paper_tables"))))
    } else {
      character()
    }
    figure_files <- if (dir.exists(file.path(run_dir, "figures"))) {
      gsub("\\\\", "/", file.path("figures", list.files(file.path(run_dir, "figures"))))
    } else {
      character()
    }
    lines <- c(
      "# Multivariate Simulation Review Bundle",
      "",
      paste("Run directory:", normalizePath(run_dir, winslash = "/", mustWork = FALSE)),
      paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      paste("Manifest files:", nrow(manifest)),
      "",
      "## Start Here",
      "",
      paste0("- ", present_key),
      "",
      "## Paper Tables",
      ""
    )
    lines <- c(lines, if (length(table_files) == 0L) "No table files found." else paste0("- ", table_files))
    lines <- c(lines, "", "## Diagnostic Figures", "")
    c(lines, if (length(figure_files) == 0L) "No figure files found." else paste0("- ", figure_files))
  }

  writeLines(build_index(manifest), index_path, useBytes = TRUE)
  manifest <- mvt_write_artifact_manifest(run_dir)
  writeLines(build_index(manifest), index_path, useBytes = TRUE)
  mvt_write_artifact_manifest(run_dir)
  invisible(index_path)
}

mvt_success_value <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  value <- tolower(trimws(as.character(x)))
  value %in% c("true", "t", "1", "yes", "y")
}

mvt_copula_variogram_coverage <- function(status, vario, copula_methods) {
  if (length(copula_methods) == 0L) {
    return(list(status = "pass", detail = "no active copula methods", n = 0L))
  }
  required_cols <- c("case_id", "method", "success")
  if (nrow(status) == 0L || !all(required_cols %in% names(status))) {
    return(list(status = "fail", detail = "fit_status_by_rep.csv missing case_id/method/success columns", n = 0L))
  }
  successful <- status[
    status$method %in% copula_methods & mvt_success_value(status$success),
    ,
    drop = FALSE
  ]
  if (nrow(successful) == 0L) {
    return(list(status = "warn", detail = "no successful copula fits to score", n = 0L))
  }
  if (nrow(vario) == 0L || !all(c("case_id", "method", "variogram_score_p05") %in% names(vario))) {
    return(list(status = "fail", detail = paste("successful copula fits", nrow(successful), "| variogram table missing"), n = 0L))
  }
  required_keys <- unique(paste(successful$case_id, successful$method, sep = "\r"))
  finite_vario <- vario[
    vario$method %in% copula_methods & is.finite(vario$variogram_score_p05),
    ,
    drop = FALSE
  ]
  finite_keys <- unique(paste(finite_vario$case_id, finite_vario$method, sep = "\r"))
  missing_keys <- setdiff(required_keys, finite_keys)
  finite_n <- length(intersect(required_keys, finite_keys))
  list(
    status = if (length(missing_keys) == 0L) "pass" else "fail",
    detail = paste(
      "successful copula fits", length(required_keys),
      "| finite variograms", finite_n,
      "| missing after successful fit", length(missing_keys)
    ),
    n = finite_n
  )
}

mvt_audit_run_dir <- function(run_dir = mvt_read_run_dir()) {
  mvt_write_artifact_manifest(run_dir)
  files <- mvt_expected_output_files()
  file_rows <- lapply(files, function(path) {
    full <- file.path(run_dir, path)
    mvt_audit_check(
      paste0("file:", path),
      if (file.exists(full)) "pass" else "fail",
      if (file.exists(full)) "present" else "missing"
    )
  })

  read_optional <- function(path) {
    full <- file.path(run_dir, path)
    mvt_read_optional_csv(full)
  }

  grid <- read_optional("scenario_grid.csv")
  status <- read_optional("fit_status_by_rep.csv")
  bench <- read_optional("benchmark_results_by_rep.csv")
  coefs <- read_optional("coefficient_results_by_rep.csv")
  dep <- read_optional("dependence_recovery_by_rep.csv")
  vario <- read_optional("variogram_scores_by_rep.csv")
  metadata <- read_optional("run_metadata.csv")
  preflight <- read_optional("preflight_checks.csv")

  metadata_value <- function(name, default = NA_character_) {
    if (nrow(metadata) == 0L || !all(c("name", "value") %in% names(metadata))) return(default)
    value <- metadata$value[metadata$name == name]
    if (length(value) == 0L || !nzchar(value[[1L]])) default else value[[1L]]
  }
  active_methods <- strsplit(metadata_value("active_comparators", paste(mvt_default_comparators(), collapse = ",")), ",", fixed = TRUE)[[1L]]
  active_methods <- trimws(active_methods[nzchar(active_methods)])
  primary_methods <- intersect(c("glm", "gamCopula_markov", "gamCopula_vine_simplified", "gamCopula_vine", "gamlss.longitudinal"), active_methods)
  copula_methods <- intersect(c("gamCopula_markov", "gamCopula_vine_simplified", "gamCopula_vine", "gamlss.longitudinal"), active_methods)
  successful_methods <- if (nrow(status) > 0L && all(c("method", "success") %in% names(status))) {
    unique(status$method[mvt_status_success(status$success)])
  } else {
    character()
  }
  primary_methods_successful <- intersect(primary_methods, successful_methods)
  copula_methods_successful <- intersect(copula_methods, successful_methods)

  checks <- c(file_rows, list(
    mvt_audit_check(
      "scenario_grid_nonempty",
      if (nrow(grid) > 0L) "pass" else "fail",
      paste("rows", nrow(grid)),
      nrow(grid)
    ),
    mvt_audit_check(
      "fit_status_all_cases_present",
      if (nrow(grid) > 0L && all(unique(grid$case_id) %in% unique(status$case_id))) "pass" else "fail",
      paste("cases with status", length(unique(status$case_id))),
      length(unique(status$case_id))
    ),
    mvt_audit_check(
      "fit_status_classified",
      if (nrow(status) > 0L && all(nzchar(trimws(as.character(status$status_class))))) "pass" else "fail",
      paste("classified rows", sum(nzchar(trimws(as.character(status$status_class))))),
      nrow(status)
    ),
    mvt_audit_check(
      "preflight_no_failures",
      if (nrow(preflight) > 0L && !any(preflight$status == "fail")) "pass" else "fail",
      if (nrow(preflight) > 0L) {
        paste("failures", sum(preflight$status == "fail"), "| warnings", sum(preflight$status == "warn"))
      } else {
        "preflight rows missing"
      },
      if (nrow(preflight) > 0L) sum(preflight$status == "fail") else NA_integer_
    ),
    mvt_audit_check(
      "benchmark_main_methods_present",
      if (nrow(bench) > 0L && all(primary_methods %in% unique(bench$method))) "pass" else "fail",
      paste("expected", paste(primary_methods, collapse = ", "), "| methods", paste(sort(unique(bench$method)), collapse = ", ")),
      length(unique(bench$method))
    ),
    mvt_audit_check(
      "coefficient_truth_for_mean_terms",
      if (nrow(coefs) > 0L) {
        focal <- coefs[coefs$term %in% c("intercept", "time", "x", "z"), , drop = FALSE]
        if (nrow(focal) > 0L && all(is.finite(focal$truth))) "pass" else "fail"
      } else {
        "fail"
      },
      paste("coefficient rows", nrow(coefs)),
      nrow(coefs)
    ),
    mvt_audit_check(
      "coefficient_rows_present_for_primary_methods",
      if (length(primary_methods_successful) == 0L) {
        "pass"
      } else if (nrow(coefs) > 0L) {
        focal <- coefs[coefs$parameter == "mu" & coefs$term %in% c("time", "x", "z"), , drop = FALSE]
        ok <- all(primary_methods_successful %in% unique(focal$method))
        if (ok) "pass" else "fail"
      } else {
        "fail"
      },
      paste("expected successful", paste(primary_methods_successful, collapse = ", "), "| methods", paste(sort(unique(coefs$method)), collapse = ", ")),
      length(unique(coefs$method))
    ),
    mvt_audit_check(
      "finite_standard_errors_for_primary_mean_terms",
      if (length(primary_methods_successful) == 0L) {
        "pass"
      } else if (nrow(coefs) > 0L && "std_error" %in% names(coefs)) {
        focal <- coefs[
          coefs$method %in% primary_methods_successful &
            coefs$parameter == "mu" &
            coefs$term %in% c("time", "x", "z"),
          ,
          drop = FALSE
        ]
        expected_n <- length(primary_methods_successful) * 3L
        if (nrow(focal) >= expected_n && all(is.finite(focal$std_error))) "pass" else "fail"
      } else {
        "fail"
      },
      "requires SE calibration and CI coverage for primary mean terms",
      if (nrow(coefs) > 0L && "std_error" %in% names(coefs)) sum(is.finite(coefs$std_error)) else 0L
    ),
    mvt_audit_check(
      "dependence_recovery_present_for_copula_methods",
      if (length(copula_methods_successful) == 0L) {
        "pass"
      } else if (nrow(dep) > 0L && all(copula_methods_successful %in% unique(dep$method))) {
        "pass"
      } else {
        "fail"
      },
      paste("expected successful", paste(copula_methods_successful, collapse = ", "), "| methods", paste(sort(unique(dep$method)), collapse = ", ")),
      length(unique(dep$method))
    ),
    mvt_audit_check(
      "finite_dependence_recovery_for_copula_methods",
      if (length(copula_methods_successful) == 0L) {
        "pass"
      } else if (nrow(dep) > 0L && all(c("theta_mae", "tau_mae") %in% names(dep))) {
        focal <- dep[dep$method %in% copula_methods_successful, , drop = FALSE]
        if (all(copula_methods_successful %in% unique(focal$method)) &&
            all(is.finite(focal$theta_mae)) &&
            all(is.finite(focal$tau_mae))) "pass" else "fail"
      } else {
        "fail"
      },
      "requires finite theta/tau recovery metrics for active copula comparators",
      if (nrow(dep) > 0L && "theta_mae" %in% names(dep)) sum(is.finite(dep$theta_mae)) else 0L
    ),
    mvt_audit_check(
      "variogram_rows_cover_benchmark_methods",
      if (nrow(bench) > 0L && nrow(vario) > 0L && all(unique(bench$method) %in% unique(vario$method))) "pass" else "fail",
      paste("variogram methods", paste(sort(unique(vario$method)), collapse = ", ")),
      length(unique(vario$method))
    ),
    mvt_audit_check(
      "finite_gamlss_longitudinal_variogram",
      if (!"gamlss.longitudinal" %in% active_methods) {
        "pass"
      } else if (nrow(vario) > 0L) {
        focal <- vario[vario$method == "gamlss.longitudinal", , drop = FALSE]
        if (nrow(focal) > 0L && any(is.finite(focal$variogram_score_p05))) "pass" else "warn"
      } else {
        "fail"
      },
      "requires successful simulation-capable fitted joint model",
      if (nrow(vario) > 0L) sum(is.finite(vario$variogram_score_p05)) else 0L
    ),
    mvt_audit_check(
      "finite_copula_variograms",
      mvt_copula_variogram_coverage(status, vario, copula_methods)$status,
      mvt_copula_variogram_coverage(status, vario, copula_methods)$detail,
      mvt_copula_variogram_coverage(status, vario, copula_methods)$n
    )
  ))

  audit <- mvt_bind_rows_fill(checks)
  mvt_write_csv(audit, file.path(run_dir, "review_audit.csv"))

  summary_lines <- c(
    "# Multivariate Simulation Review Audit",
    "",
    paste("Run directory:", run_dir),
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    paste("Checks passing:", sum(audit$status == "pass"), "of", nrow(audit)),
    paste("Warnings:", sum(audit$status == "warn")),
    paste("Failures:", sum(audit$status == "fail")),
    "",
    "## Non-Passing Checks",
    ""
  )
  non_pass <- audit[audit$status != "pass", , drop = FALSE]
  if (nrow(non_pass) == 0L) {
    summary_lines <- c(summary_lines, "None.")
  } else {
    summary_lines <- c(
      summary_lines,
      paste0("- ", non_pass$check, " [", non_pass$status, "]: ", non_pass$detail)
    )
  }
  writeLines(summary_lines, file.path(run_dir, "review_audit.md"), useBytes = TRUE)
  mvt_write_artifact_manifest(run_dir)
  invisible(audit)
}

mvt_read_optional_csv <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0L) {
    return(data.frame())
  }
  tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) data.frame()
  )
}

mvt_publication_readiness_spec <- function() {
  data.frame(
    role = c("pilot", "main_core", "appendix"),
    env_var = c(
      "GAMLSS_LONGITUDINAL_MVT_PILOT_RUN_DIR",
      "GAMLSS_LONGITUDINAL_MVT_MAIN_RUN_DIR",
      "GAMLSS_LONGITUDINAL_MVT_APPENDIX_RUN_DIR"
    ),
    required = c(TRUE, TRUE, FALSE),
    min_reps = c(5L, 100L, 100L),
    required_time = c("5,20", "20", "5,20,50"),
    required_family = c("gaussian,gamma,binomial", "gaussian,poisson,gamma,binomial", "gaussian,poisson,gamma,binomial"),
    required_dependence = c(
      "exchangeable,ar1,covariate_dependent_adjacent",
      "exchangeable,ar1,time_varying_adjacent,covariate_dependent_adjacent",
      "exchangeable,ar1,time_varying_adjacent,covariate_dependent_adjacent"
    ),
    required_method = c(
      "glm,glmm,gee_independence,gee_exchangeable,gee_ar1,gee_unstructured,gamCopula_markov,gamCopula_vine_simplified,gamlss.longitudinal",
      "glm,glmm,gee_independence,gee_exchangeable,gee_ar1,gee_unstructured,gamCopula_markov,gamCopula_vine_simplified,gamlss.longitudinal",
      "glm,glmm,gee_independence,gee_exchangeable,gee_ar1,gee_unstructured,gamCopula_markov,gamCopula_vine_simplified,gamlss.longitudinal"
    ),
    stringsAsFactors = FALSE
  )
}

mvt_publication_readiness_optional_spec <- function() {
  data.frame(
    role = c("special_gamlss", "glmm_sensitivity"),
    env_var = c(
      "GAMLSS_LONGITUDINAL_MVT_SPECIAL_RUN_DIR",
      "GAMLSS_LONGITUDINAL_MVT_GLMM_SENSITIVITY_RUN_DIR"
    ),
    required = c(FALSE, FALSE),
    required_family = c("gg_continuous", ""),
    required_method = c("gamlss.longitudinal", "glmm_slope"),
    stringsAsFactors = FALSE
  )
}

mvt_split_csv_value <- function(value) {
  value <- as.character(value %||% "")
  if (!nzchar(value)) return(character())
  trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
}

mvt_readiness_check <- function(role, check, status, detail = "", n = NA_integer_, run_dir = NA_character_) {
  data.frame(
    role = role,
    check = check,
    status = status,
    n = n,
    detail = detail,
    run_dir = run_dir,
    stringsAsFactors = FALSE
  )
}

mvt_publication_readiness_one <- function(spec_row) {
  role <- spec_row$role[[1L]]
  run_dir <- mvt_env(spec_row$env_var[[1L]], "")
  required <- isTRUE(spec_row$required[[1L]])
  if (!nzchar(run_dir)) {
    return(mvt_readiness_check(
      role,
      "run_dir_configured",
      if (required) "fail" else "warn",
      paste("set", spec_row$env_var[[1L]], "to include this evidence"),
      run_dir = NA_character_
    ))
  }
  if (!dir.exists(run_dir)) {
    return(mvt_readiness_check(role, "run_dir_exists", "fail", "directory not found", run_dir = run_dir))
  }

  grid <- mvt_read_optional_csv(file.path(run_dir, "scenario_grid.csv"))
  status <- mvt_read_optional_csv(file.path(run_dir, "fit_status_by_rep.csv"))
  bench <- mvt_read_optional_csv(file.path(run_dir, "benchmark_results_by_rep.csv"))
  coef <- mvt_read_optional_csv(file.path(run_dir, "coefficient_results_by_rep.csv"))
  dep <- mvt_read_optional_csv(file.path(run_dir, "dependence_recovery_by_rep.csv"))
  vario <- mvt_read_optional_csv(file.path(run_dir, "variogram_scores_by_rep.csv"))
  audit <- mvt_read_optional_csv(file.path(run_dir, "review_audit.csv"))

  required_time <- suppressWarnings(as.integer(mvt_split_csv_value(spec_row$required_time[[1L]])))
  required_family <- mvt_split_csv_value(spec_row$required_family[[1L]])
  required_dependence <- mvt_split_csv_value(spec_row$required_dependence[[1L]])
  required_method <- mvt_split_csv_value(spec_row$required_method[[1L]])

  rows <- list(
    mvt_readiness_check(role, "run_dir_configured", "pass", spec_row$env_var[[1L]], run_dir = run_dir),
    mvt_readiness_check(role, "scenario_grid_present", if (nrow(grid) > 0L) "pass" else "fail", paste("rows", nrow(grid)), nrow(grid), run_dir),
    mvt_readiness_check(
      role,
      "review_audit_clean",
      if (nrow(audit) > 0L && !any(audit$status == "fail")) "pass" else "fail",
      if (nrow(audit) > 0L) paste("failures", sum(audit$status == "fail"), "| warnings", sum(audit$status == "warn")) else "review_audit.csv missing",
      if (nrow(audit) > 0L) sum(audit$status == "fail") else NA_integer_,
      run_dir
    ),
    mvt_readiness_check(
      role,
      "minimum_replicates",
      if (nrow(grid) > 0L && "rep" %in% names(grid) && length(unique(grid$rep)) >= spec_row$min_reps[[1L]]) "pass" else "fail",
      paste("found", if (nrow(grid) > 0L && "rep" %in% names(grid)) length(unique(grid$rep)) else 0L, "| required", spec_row$min_reps[[1L]]),
      if (nrow(grid) > 0L && "rep" %in% names(grid)) length(unique(grid$rep)) else 0L,
      run_dir
    ),
    mvt_readiness_check(
      role,
      "time_grid_complete",
      if (nrow(grid) > 0L && "n_time" %in% names(grid) && all(required_time %in% unique(grid$n_time))) "pass" else "fail",
      paste("required", paste(required_time, collapse = ","), "| found", if ("n_time" %in% names(grid)) paste(sort(unique(grid$n_time)), collapse = ",") else ""),
      if ("n_time" %in% names(grid)) length(unique(grid$n_time)) else 0L,
      run_dir
    ),
    mvt_readiness_check(
      role,
      "family_grid_complete",
      if (nrow(grid) > 0L && "family_name" %in% names(grid) && all(required_family %in% unique(grid$family_name))) "pass" else "fail",
      paste("required", paste(required_family, collapse = ","), "| found", if ("family_name" %in% names(grid)) paste(sort(unique(grid$family_name)), collapse = ",") else ""),
      if ("family_name" %in% names(grid)) length(unique(grid$family_name)) else 0L,
      run_dir
    ),
    mvt_readiness_check(
      role,
      "dependence_grid_complete",
      if (nrow(grid) > 0L && "dependence" %in% names(grid) && all(required_dependence %in% unique(grid$dependence))) "pass" else "fail",
      paste("required", paste(required_dependence, collapse = ","), "| found", if ("dependence" %in% names(grid)) paste(sort(unique(grid$dependence)), collapse = ",") else ""),
      if ("dependence" %in% names(grid)) length(unique(grid$dependence)) else 0L,
      run_dir
    ),
    mvt_readiness_check(
      role,
      "status_all_cases_present",
      if (nrow(grid) > 0L && nrow(status) > 0L && all(unique(grid$case_id) %in% unique(status$case_id))) "pass" else "fail",
      paste("cases with status", if ("case_id" %in% names(status)) length(unique(status$case_id)) else 0L),
      if ("case_id" %in% names(status)) length(unique(status$case_id)) else 0L,
      run_dir
    ),
    mvt_readiness_check(
      role,
      "methods_complete",
      if (nrow(bench) > 0L && "method" %in% names(bench) && all(required_method %in% unique(bench$method))) "pass" else "fail",
      paste("required", paste(required_method, collapse = ","), "| found", if ("method" %in% names(bench)) paste(sort(unique(bench$method)), collapse = ",") else ""),
      if ("method" %in% names(bench)) length(unique(bench$method)) else 0L,
      run_dir
    ),
    mvt_readiness_check(
      role,
      "coefficient_truth_present",
      if (nrow(coef) > 0L && all(c("term", "truth") %in% names(coef)) && any(coef$term %in% c("intercept", "time", "x", "z")) && all(is.finite(coef$truth[coef$term %in% c("intercept", "time", "x", "z")]))) "pass" else "fail",
      paste("coefficient rows", nrow(coef)),
      nrow(coef),
      run_dir
    ),
    mvt_readiness_check(
      role,
      "dependence_recovery_for_copulas",
      if (nrow(dep) > 0L && "method" %in% names(dep) && all(intersect(c("gamCopula_markov", "gamCopula_vine_simplified", "gamCopula_vine", "gamlss.longitudinal"), required_method) %in% unique(dep$method))) "pass" else "fail",
      paste("methods", if ("method" %in% names(dep)) paste(sort(unique(dep$method)), collapse = ",") else ""),
      if ("method" %in% names(dep)) length(unique(dep$method)) else 0L,
      run_dir
    ),
    mvt_readiness_check(
      role,
      "finite_joint_variogram",
      if (nrow(vario) > 0L && "method" %in% names(vario) && "variogram_score_p05" %in% names(vario)) {
        focal <- vario[vario$method == "gamlss.longitudinal", , drop = FALSE]
        if (nrow(focal) > 0L && any(is.finite(focal$variogram_score_p05))) "pass" else "fail"
      } else {
        "fail"
      },
      paste("finite scores", if ("variogram_score_p05" %in% names(vario)) sum(is.finite(vario$variogram_score_p05)) else 0L),
      if ("variogram_score_p05" %in% names(vario)) sum(is.finite(vario$variogram_score_p05)) else 0L,
      run_dir
    ),
    mvt_readiness_check(
      role,
      "finite_copula_variograms",
      mvt_copula_variogram_coverage(
        status,
        vario,
        intersect(c("gamCopula_markov", "gamCopula_vine_simplified", "gamCopula_vine", "gamlss.longitudinal"), required_method)
      )$status,
      mvt_copula_variogram_coverage(
        status,
        vario,
        intersect(c("gamCopula_markov", "gamCopula_vine_simplified", "gamCopula_vine", "gamlss.longitudinal"), required_method)
      )$detail,
      mvt_copula_variogram_coverage(
        status,
        vario,
        intersect(c("gamCopula_markov", "gamCopula_vine_simplified", "gamCopula_vine", "gamlss.longitudinal"), required_method)
      )$n,
      run_dir
    )
  )
  mvt_bind_rows_fill(rows)
}

mvt_publication_readiness_optional_one <- function(spec_row) {
  role <- spec_row$role[[1L]]
  run_dir <- mvt_env(spec_row$env_var[[1L]], "")
  if (!nzchar(run_dir)) {
    return(mvt_readiness_check(role, "run_dir_configured", "warn", paste("optional evidence; set", spec_row$env_var[[1L]]), run_dir = NA_character_))
  }
  if (!dir.exists(run_dir)) {
    return(mvt_readiness_check(role, "run_dir_exists", "fail", "directory not found", run_dir = run_dir))
  }
  grid <- mvt_read_optional_csv(file.path(run_dir, "scenario_grid.csv"))
  bench <- mvt_read_optional_csv(file.path(run_dir, "benchmark_results_by_rep.csv"))
  audit <- mvt_read_optional_csv(file.path(run_dir, "review_audit.csv"))
  required_family <- mvt_split_csv_value(spec_row$required_family[[1L]])
  required_method <- mvt_split_csv_value(spec_row$required_method[[1L]])
  rows <- list(
    mvt_readiness_check(role, "run_dir_configured", "pass", spec_row$env_var[[1L]], run_dir = run_dir),
    mvt_readiness_check(role, "review_audit_clean", if (nrow(audit) > 0L && !any(audit$status == "fail")) "pass" else "fail", if (nrow(audit) > 0L) paste("failures", sum(audit$status == "fail")) else "review_audit.csv missing", run_dir = run_dir),
    mvt_readiness_check(role, "required_family_present", if (length(required_family) == 0L || (nrow(grid) > 0L && "family_name" %in% names(grid) && all(required_family %in% unique(grid$family_name)))) "pass" else "fail", paste(required_family, collapse = ","), run_dir = run_dir),
    mvt_readiness_check(role, "required_method_present", if (length(required_method) == 0L || (nrow(bench) > 0L && "method" %in% names(bench) && all(required_method %in% unique(bench$method)))) "pass" else "fail", paste(required_method, collapse = ","), run_dir = run_dir)
  )
  mvt_bind_rows_fill(rows)
}

mvt_publication_readiness_audit <- function(output_dir = mvt_output_root) {
  main_spec <- mvt_publication_readiness_spec()
  optional_spec <- mvt_publication_readiness_optional_spec()
  audit <- mvt_bind_rows_fill(
    lapply(seq_len(nrow(main_spec)), function(i) mvt_publication_readiness_one(main_spec[i, , drop = FALSE])),
    lapply(seq_len(nrow(optional_spec)), function(i) mvt_publication_readiness_optional_one(optional_spec[i, , drop = FALSE]))
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  csv_path <- file.path(output_dir, "publication_readiness_audit.csv")
  md_path <- file.path(output_dir, "publication_readiness_audit.md")
  mvt_write_csv(audit, csv_path)
  required_roles <- main_spec$role[main_spec$required]
  required_audit <- audit[audit$role %in% required_roles, , drop = FALSE]
  ready <- nrow(required_audit) > 0L && all(required_audit$status == "pass")
  non_pass <- audit[audit$status != "pass", , drop = FALSE]
  lines <- c(
    "# Publication Readiness Audit",
    "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Required evidence ready:", ready),
    paste("Required checks passing:", sum(required_audit$status == "pass"), "of", nrow(required_audit)),
    paste("Warnings:", sum(audit$status == "warn")),
    paste("Failures:", sum(audit$status == "fail")),
    "",
    "## Required Run Directory Environment Variables",
    "",
    paste0("- ", main_spec$role, ": `", main_spec$env_var, "`", ifelse(main_spec$required, " (required)", " (appendix/recommended)")),
    "",
    "## Non-Passing Checks",
    ""
  )
  if (nrow(non_pass) == 0L) {
    lines <- c(lines, "None.")
  } else {
    lines <- c(lines, paste0("- ", non_pass$role, "/", non_pass$check, " [", non_pass$status, "]: ", non_pass$detail))
  }
  writeLines(lines, md_path, useBytes = TRUE)
  invisible(list(audit = audit, ready = ready, csv = csv_path, md = md_path))
}

mvt_group_metric_summary <- function(data, group_cols, metrics, extra_bool = character()) {
  if (nrow(data) == 0L || length(metrics) == 0L) return(data.frame())
  groups <- unique(data[group_cols])
  rows <- list()
  k <- 1L
  for (i in seq_len(nrow(groups))) {
    idx <- rep(TRUE, nrow(data))
    for (col in group_cols) idx <- idx & data[[col]] == groups[[col]][i]
    sub <- data[idx, , drop = FALSE]
    for (metric in metrics) {
      vals <- suppressWarnings(as.numeric(sub[[metric]]))
      finite <- is.finite(vals)
      rows[[k]] <- cbind(
        groups[i, , drop = FALSE],
        data.frame(
          metric = metric,
          n = length(vals),
          n_finite = sum(finite),
          mean = if (any(finite)) mean(vals[finite]) else NA_real_,
          median = if (any(finite)) stats::median(vals[finite]) else NA_real_,
          sd = if (sum(finite) > 1L) stats::sd(vals[finite]) else NA_real_,
          stringsAsFactors = FALSE
        )
      )
      k <- k + 1L
    }
    for (metric in extra_bool) {
      vals <- sub[[metric]] %in% c(TRUE, "TRUE", "true", "1")
      known <- !is.na(sub[[metric]])
      rows[[k]] <- cbind(
        groups[i, , drop = FALSE],
        data.frame(
          metric = metric,
          n = length(vals),
          n_finite = sum(known),
          mean = if (any(known)) mean(vals[known]) else NA_real_,
          median = if (any(known)) stats::median(as.numeric(vals[known])) else NA_real_,
          sd = if (sum(known) > 1L) stats::sd(as.numeric(vals[known])) else NA_real_,
          stringsAsFactors = FALSE
        )
      )
      k <- k + 1L
    }
  }
  if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
}
