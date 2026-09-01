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

mvt_registered_timeout <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(as.numeric(default))
  value <- suppressWarnings(as.numeric(raw))
  if (length(value) != 1L || is.na(value) || value <= 0) return(Inf)
  value
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

mvt_nearest_neighbor_comparators <- function() {
  # Prespecified bound for JSS-014: exactly two tractable, task-matched
  # copula comparators. Full-vine and low-dimensional alternatives remain
  # documented sensitivities/exclusions, not additional headline benchmarks.
  c("gamCopula_markov", "gamCopula_vine_simplified")
}

mvt_capability_snapshot <- function(as_of = as.Date("2026-09-01")) {
  software <- c("gamlss.longitudinal", "gamCopula", "GJRM", "geepack", "glmmTMB", "gamlss", "VineCopula")
  installed_versions <- vapply(software, function(package) {
    package_path <- suppressWarnings(find.package(package, quiet = TRUE))
    if (!length(package_path) || !nzchar(package_path)) return("not_installed")
    description <- tryCatch(read.dcf(file.path(package_path, "DESCRIPTION"), fields = "Version"), error = function(e) NULL)
    if (is.null(description) || !length(description) || !nzchar(description[[1L]])) {
      stop("Installed version unavailable for ", package, call. = FALSE)
    }
    as.character(description[[1L]])
  }, character(1L))
  source_versions <- c(
    gamlss.longitudinal = as.character(read.dcf(file.path(mvt_repo_root, "DESCRIPTION"), fields = "Version")[[1L]]),
    gamCopula = "0.0-8", GJRM = "0.2-6.9", geepack = "1.3.13",
    glmmTMB = "1.1.14", gamlss = "5.5-0", VineCopula = "2.6.1"
  )
  claim_software <- c(
    "gamlss.longitudinal", "gamCopula", "gamCopula", "gamCopula", "GJRM",
    "geepack", "geepack", "glmmTMB", "gamlss", "VineCopula"
  )
  version_sources <- c(
    gamlss.longitudinal = "DESCRIPTION",
    gamCopula = "https://cran.r-project.org/web/packages/gamCopula/DESCRIPTION",
    GJRM = "https://cran.r-project.org/web/packages/GJRM/DESCRIPTION",
    geepack = "https://cran.r-project.org/web/packages/geepack/DESCRIPTION",
    glmmTMB = "https://cran.r-project.org/web/packages/glmmTMB/DESCRIPTION",
    gamlss = "https://cran.r-project.org/web/packages/gamlss/DESCRIPTION",
    VineCopula = "https://cran.r-project.org/web/packages/VineCopula/DESCRIPTION"
  )
  data.frame(
    as_of = rep(format(as_of, "%Y-%m-%d"), length(claim_software)),
    retrieved_on = rep(format(as_of, "%Y-%m-%d"), length(claim_software)),
    claim_id = sprintf("CAP-%03d", seq_along(claim_software)),
    software = claim_software,
    source_version = unname(source_versions[claim_software]),
    version_source_url = unname(version_sources[claim_software]),
    installed_version = unname(installed_versions[claim_software]),
    claim_type = c(
      "evaluation_role", "headline_workflow", "headline_workflow", "tractability",
      "task_equivalence", "standard_benchmark", "stress_test", "context_only",
      "context_only", "backend_only"
    ),
    claim_value = c(
      "focal", "retain_markov", "retain_simplified_vine", "full_vine_sensitivity_only",
      "exclude_non_equivalent", "family_specific_GEE", "unstructured_T20plus",
      "mixed_model_context", "independence_margin_context", "copula_backend"
    ),
    claim = c(
      "The focal workflow fits longitudinal distributional regression with adjacent copulas.",
      "The two-stage adjacent Markov gamCopula workflow is a retained headline comparator.",
      "The two-stage simplified-vine gamCopula workflow is a retained headline comparator.",
      "The full-vine gamCopula workflow is excluded from the registered headline grid by an a priori bounded-comparator scope decision.",
      "GJRM is excluded from the empirical grid because visit-as-equation encoding changes the model dimension and estimand.",
      "geepack GEE results are reported separately for every response family in the scenario grid.",
      "Unstructured geepack GEE at T at least 20 is treated as a feasibility stress test.",
      "glmmTMB is capability context rather than an additional empirical method because the retained lme4 fit already represents the mixed-model class.",
      "gamlss is independence-margin context and cannot estimate longitudinal dependence.",
      "VineCopula is a copula backend and does not itself fit matched longitudinal marginal regressions."
    ),
    documentation_topic = c(
      "method-traceability registry model fitting wrapper",
      "gamCopula::gamBiCopFit maximum penalized likelihood",
      "gamCopula::gamVineSeqFit sequential GAM-vine estimation",
      "gamCopula::gamVineSeqFit sequential GAM-vine estimation",
      "GJRM::gjrm formula-list interface",
      "geepack::geeglm family and corstr arguments",
      "geepack::geeglm unstructured corstr",
      "glmmTMB vignette section The unstructured covariance",
      "gamlss::gamlss location scale and shape interface",
      "VineCopula::RVineCopSelect pair-copula selection"
    ),
    documentation_url = c(
      "inst/standards/method-traceability.csv",
      "https://search.r-project.org/CRAN/refmans/gamCopula/html/gamBiCopFit.html",
      "https://search.r-project.org/CRAN/refmans/gamCopula/html/gamVineSeqFit.html",
      "https://search.r-project.org/CRAN/refmans/gamCopula/html/gamVineSeqFit.html",
      "https://search.r-project.org/CRAN/refmans/GJRM/html/gjrm.html",
      "https://search.r-project.org/CRAN/refmans/geepack/html/geeglm.html",
      "https://search.r-project.org/CRAN/refmans/geepack/html/geeglm.html",
      "https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html#the-unstructured-covariance",
      "https://search.r-project.org/CRAN/refmans/gamlss/html/gamlss.html",
      "https://search.r-project.org/CRAN/refmans/VineCopula/html/RVineCopSelect.html"
    ),
    empirical_evidence_link = c(
      "paper/R/09-simulation-multivariate-longitudinal/00-multivariate-setup.R#mvt_fit_gamlss_longitudinal",
      "nearest_neighbor_results.csv#gamCopula_markov",
      "nearest_neighbor_results.csv#gamCopula_vine_simplified",
      "comparator_scope_registry.csv#gamCopula_vine",
      "comparator_scope_registry.csv#GJRM",
      "gee_family_results.csv",
      "gee_unstructured_stress_test.csv",
      "comparator_scope_registry.csv#glmmTMB",
      "comparator_scope_registry.csv#gamlss_independence",
      "comparator_scope_registry.csv#VineCopula"
    ),
    empirical_evidence = c(
      "Implementation hook for the fitted focal method.",
      "Attempt and uncertainty rows for the retained Markov workflow.",
      "Attempt and uncertainty rows for the retained simplified-vine workflow.",
      "A priori qualitative scope decision; no unfrozen legacy runtime is used as evidence.",
      "A priori non-equivalence decision rather than an empirical infeasibility claim.",
      "Family-specific GEE attempt and uncertainty rows.",
      "T=20 unstructured GEE attempt and failure rows.",
      "Capability-only decision recorded in the scope registry.",
      "Independence-only decision recorded in the scope registry.",
      "Backend-only decision recorded in the scope registry."
    ),
    stringsAsFactors = FALSE
  )
}

mvt_capability_registered_urls <- function() {
  list(
    version_source_url = c(
      "DESCRIPTION",
      "https://cran.r-project.org/web/packages/gamCopula/DESCRIPTION",
      "https://cran.r-project.org/web/packages/GJRM/DESCRIPTION",
      "https://cran.r-project.org/web/packages/geepack/DESCRIPTION",
      "https://cran.r-project.org/web/packages/glmmTMB/DESCRIPTION",
      "https://cran.r-project.org/web/packages/gamlss/DESCRIPTION",
      "https://cran.r-project.org/web/packages/VineCopula/DESCRIPTION"
    ),
    documentation_url = c(
      "inst/standards/method-traceability.csv",
      "https://search.r-project.org/CRAN/refmans/gamCopula/html/gamBiCopFit.html",
      "https://search.r-project.org/CRAN/refmans/gamCopula/html/gamVineSeqFit.html",
      "https://search.r-project.org/CRAN/refmans/GJRM/html/gjrm.html",
      "https://search.r-project.org/CRAN/refmans/geepack/html/geeglm.html",
      "https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html#the-unstructured-covariance",
      "https://search.r-project.org/CRAN/refmans/gamlss/html/gamlss.html",
      "https://search.r-project.org/CRAN/refmans/VineCopula/html/RVineCopSelect.html"
    )
  )
}

mvt_capability_provenance_validation <- function(capabilities) {
  required <- c(
    "claim_id", "as_of", "retrieved_on", "software", "source_version", "version_source_url",
    "installed_version", "claim_type", "claim_value", "claim", "documentation_topic",
    "documentation_url", "empirical_evidence_link", "empirical_evidence"
  )
  problems <- character()
  if (!all(required %in% names(capabilities))) {
    problems <- c(problems, paste("missing columns", paste(setdiff(required, names(capabilities)), collapse = ",")))
    return(list(valid = FALSE, problems = problems, required = required))
  }
  if (!nrow(capabilities)) problems <- c(problems, "no claims")
  if (anyDuplicated(capabilities$claim_id)) problems <- c(problems, "duplicate claim_id")
  if (any(!grepl("^CAP-[0-9]{3}$", capabilities$claim_id))) problems <- c(problems, "invalid claim_id")
  atomic_key <- paste(capabilities$software, capabilities$claim_type, capabilities$claim_value, capabilities$claim, sep = "\r")
  if (anyDuplicated(atomic_key)) problems <- c(problems, "duplicate atomic claim")
  if (any(!nzchar(trimws(capabilities$claim)))) problems <- c(problems, "blank claim")
  valid_date <- function(x) {
    parsed <- suppressWarnings(as.Date(as.character(x), format = "%Y-%m-%d"))
    grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", x) & !is.na(parsed)
  }
  if (any(!valid_date(capabilities$as_of))) problems <- c(problems, "invalid as_of")
  if (any(!valid_date(capabilities$retrieved_on))) problems <- c(problems, "invalid retrieved_on")
  pinned <- "^[0-9]+([.-][0-9]+)+$"
  if (any(!grepl(pinned, capabilities$source_version))) problems <- c(problems, "unknown or unpinned source_version")
  if (any(!(grepl(pinned, capabilities$installed_version) | capabilities$installed_version == "not_installed"))) {
    problems <- c(problems, "unknown installed_version")
  }
  registered_urls <- mvt_capability_registered_urls()
  if (any(!capabilities$version_source_url %in% registered_urls$version_source_url)) {
    problems <- c(problems, "unregistered or unresolved version_source_url")
  }
  if (any(!capabilities$documentation_url %in% registered_urls$documentation_url)) {
    problems <- c(problems, "unregistered or unresolved documentation_url")
  }
  if (any(!nzchar(trimws(capabilities$claim_type))) || any(!nzchar(trimws(capabilities$claim_value)))) {
    problems <- c(problems, "non-atomic claim encoding")
  }
  if (any(!nzchar(trimws(capabilities$empirical_evidence_link)))) problems <- c(problems, "blank empirical evidence link")
  list(valid = length(problems) == 0L, problems = unique(problems), required = required)
}

mvt_comparator_scope_registry <- function() {
  data.frame(
    method = c(
      "gamCopula_markov", "gamCopula_vine_simplified", "gamCopula_vine",
      "GJRM", "geepack_GEE", "glmmTMB", "gamlss_independence", "VineCopula"
    ),
    decision = c(
      "retain", "retain", "targeted_sensitivity_only", "exclude_empirical",
      "standard_model_track", "context_only", "context_only", "backend_only"
    ),
    headline_empirical_comparator = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    matched_task = c(
      "same data, marginal mean formulas, adjacent dependence target, and replicate seeds",
      "same data, marginal mean formulas, multivariate copula target, and replicate seeds",
      "same task but not tractable across the prespecified production grid",
      "equation-oriented low-dimensional joint models are not equivalent to an arbitrary-T longitudinal first-order fit",
      "same marginal-mean estimand but no full distributional joint-likelihood estimand",
      "mixed-model capability context; the empirical GLMM benchmark is implemented with lme4",
      "same marginal distribution class without a repeated-outcome dependence model",
      "copula backend only; no fitted longitudinal marginal regression"
    ),
    limitation_or_infeasibility = c(
      "two-stage fit; margin and dependence uncertainty is not equivalent to joint estimation",
      "simplified-vine restriction; two-stage uncertainty is not equivalent to joint estimation",
      "excluded a priori to keep the headline comparison to two bounded nearest-neighbor workflows; reserve for a separately registered sensitivity",
      "forcing repeated visits into separate equations changes dimension, formulas, and estimand",
      "working-correlation and robust mean inference answer a narrower question; unstructured correlation scales quadratically in T",
      "does not provide a task-equivalent copula estimand; avoid duplicating the retained random-intercept GLMM class benchmark",
      "independence cannot assess recovery of longitudinal dependence",
      "dependency used by gamCopula workflows; adding it as a separate method would not define matched fitted margins"
    ),
    stringsAsFactors = FALSE
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

mvt_replace_file_atomic <- function(temporary, path, backup = NULL) {
  if (!file.exists(temporary)) stop("Atomic replacement source does not exist: ", temporary, call. = FALSE)
  if (!file.exists(path)) {
    if (!file.rename(temporary, path)) stop("Could not atomically install file: ", path, call. = FALSE)
    return(invisible(path))
  }
  if (is.null(backup)) backup <- tempfile(paste0(".", basename(path), "-previous-"), tmpdir = dirname(path))
  if (file.exists(backup)) stop("Atomic replacement backup path unexpectedly exists: ", backup, call. = FALSE)
  if (.Platform$OS.type == "windows") {
    ps_literal <- function(value) paste0("'", gsub("'", "''", normalizePath(value, winslash = "/", mustWork = FALSE), fixed = TRUE), "'")
    command <- paste0(
      "[System.IO.File]::Replace(", ps_literal(temporary), ",",
      ps_literal(path), ",", ps_literal(backup), ",$true)"
    )
    status <- suppressWarnings(system2(
      "powershell", c("-NoProfile", "-NonInteractive", "-Command", shQuote(command)),
      stdout = TRUE, stderr = TRUE
    ))
    if (!identical(attr(status, "status") %||% 0L, 0L) || !file.exists(path)) {
      stop("Windows atomic aggregate replacement failed: ", paste(status, collapse = " | "), call. = FALSE)
    }
  } else if (!file.rename(temporary, path)) {
    stop("Could not atomically replace file: ", path, call. = FALSE)
  }
  if (file.exists(backup) && !file.remove(backup)) warning("Atomic replacement succeeded but its rollback backup remains: ", backup)
  invisible(path)
}

mvt_save_rds_atomic <- function(x, path, lease) {
  mvt_assert_active_lease(lease, path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(x, temporary, version = 3)
  mvt_assert_active_lease(lease, path)
  mvt_replace_file_atomic(temporary, path)
  invisible(path)
}

mvt_write_csv_atomic <- function(x, path, lease) {
  owner <- mvt_assert_active_lease(lease, path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  utils::write.csv(x, temporary, row.names = FALSE, na = "")
  mvt_assert_active_lease(lease, path)
  mvt_replace_file_atomic(temporary, path)
  commit <- list(
    schema_version = 1L,
    file = basename(path), owner_role = "lease_parent", writer_pid = as.integer(owner$pid),
    lease_nonce = owner$nonce, hostname = owner$hostname,
    written_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    sha256 = mvt_sha256_file(path), rows = as.integer(nrow(x)),
    columns = names(x), schema_sha256 = mvt_hash_object(names(x)),
    bytes = as.numeric(file.info(path)$size)
  )
  mvt_save_rds_atomic(commit, paste0(path, ".commit.rds"), lease)
  mvt_save_rds_atomic(commit, paste0(path, ".ownership.rds"), lease)
  invisible(path)
}

mvt_write_lines_atomic <- function(lines, path, lease) {
  owner <- mvt_assert_active_lease(lease, path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  writeLines(lines, temporary, useBytes = TRUE)
  mvt_assert_active_lease(lease, path)
  mvt_replace_file_atomic(temporary, path)
  commit <- list(
    schema_version = 1L, file = basename(path), owner_role = "lease_parent",
    writer_pid = as.integer(owner$pid), lease_nonce = owner$nonce,
    hostname = owner$hostname, written_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    sha256 = mvt_sha256_file(path), rows = as.integer(length(lines)),
    columns = character(), schema_sha256 = mvt_hash_object(character()),
    bytes = as.numeric(file.info(path)$size)
  )
  mvt_save_rds_atomic(commit, paste0(path, ".commit.rds"), lease)
  mvt_save_rds_atomic(commit, paste0(path, ".ownership.rds"), lease)
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
    "gee_family_results.csv",
    "gee_unstructured_stress_test.csv",
    "nearest_neighbor_results.csv",
    "nearest_neighbor_paired_contrasts.csv",
    "capability_snapshot_2026-09-01.csv",
    "comparator_scope_registry.csv",
    "phase2_attempt_reconciliation.csv",
    "phase2_benchmark_audit.csv",
    "phase2_benchmark_audit.md",
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

mvt_phase2_public_output_allowlist <- function() {
  list(
    attempt_artifacts = c(
      "scenario_grid.csv", "run_metadata.csv", "package_versions.csv",
      "worker_attestations.csv", "checkpoint_rejections.csv",
      paste0(mvt_result_names(), "_by_rep.csv")
    ),
    evidence_artifacts = c(
      "gee_family_results.csv", "gee_unstructured_stress_test.csv",
      "nearest_neighbor_results.csv", "nearest_neighbor_paired_contrasts.csv",
      "capability_snapshot_2026-09-01.csv", "comparator_scope_registry.csv",
      "phase2_attempt_reconciliation.csv", "phase2_benchmark_audit.csv",
      "phase2_benchmark_audit.md"
    )
  )
}

mvt_phase2_claim_output_contract <- function() {
  list(
    schema_version = 1L,
    required_fields = c(
      "claim_id", "scenario_key", "row_key", "metric", "direction",
      "denominator", "effect_artifact", "effect_column", "mcse_column",
      "ci_lower_column", "ci_upper_column", "interval_support", "wording_strength"
    ),
    directions = c("higher", "lower", "positive", "negative", "no_direction"),
    interval_support = c("supported", "not_supported", "not_applicable"),
    wording_strength = c("exact", "directional", "cautious", "descriptive"),
    effect_artifacts = setdiff(
      mvt_phase2_public_output_allowlist()$evidence_artifacts,
      c("phase2_benchmark_audit.md")
    )
  )
}

mvt_phase2_claim_wording_matrix <- function() {
  rbind(
    expand.grid(
      direction_class = "directional", interval_support = "supported",
      wording_strength = c("exact", "directional", "cautious", "descriptive"),
      stringsAsFactors = FALSE
    ),
    expand.grid(
      direction_class = "no_direction", interval_support = "not_applicable",
      wording_strength = c("exact", "cautious", "descriptive"),
      stringsAsFactors = FALSE
    )
  )
}

mvt_add_phase2_evidence_keys <- function(x) {
  if (!is.data.frame(x) || !nrow(x)) return(x)
  scenario_fields <- intersect(
    c("scenario", "generator", "dependence", "correlation_level", "n_time", "n_subject", "total_rows", "family"),
    names(x)
  )
  row_fields <- unique(c(scenario_fields, intersect(c("method", "neighbor", "metric", "check"), names(x))))
  x$scenario_key <- if (length(scenario_fields)) {
    do.call(paste, c(lapply(x[scenario_fields], as.character), sep = "|"))
  } else {
    rep("all", nrow(x))
  }
  x$row_key <- if (length(row_fields)) {
    do.call(paste, c(lapply(x[row_fields], as.character), sep = "|"))
  } else {
    paste0("row-", seq_len(nrow(x)))
  }
  x
}

mvt_validate_phase2_claim_rows <- function(claims, evidence_dir) {
  contract <- mvt_phase2_claim_output_contract()
  missing <- setdiff(contract$required_fields, names(claims))
  if (length(missing)) stop("Phase 2 claim output contract is missing fields: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!nrow(claims) || anyDuplicated(claims$claim_id) || any(!nzchar(trimws(as.character(claims$claim_id))))) {
    stop("Phase 2 claim IDs must be nonempty and unique.", call. = FALSE)
  }
  text_fields <- setdiff(contract$required_fields, "claim_id")
  if (any(vapply(claims[text_fields], function(x) any(is.na(x) | !nzchar(trimws(as.character(x)))), logical(1L)))) {
    stop("Phase 2 claim output contract contains blank keys or evidence fields.", call. = FALSE)
  }
  if (!all(as.character(claims$direction) %in% contract$directions)) stop("Phase 2 claim direction is not registered.", call. = FALSE)
  if (!all(as.character(claims$interval_support) %in% contract$interval_support)) stop("Phase 2 claim interval support is not registered.", call. = FALSE)
  if (!all(as.character(claims$wording_strength) %in% contract$wording_strength)) stop("Phase 2 claim wording strength is not registered.", call. = FALSE)
  if (!all(as.character(claims$effect_artifact) %in% contract$effect_artifacts)) stop("Phase 2 claim effect artifact is not allowlisted.", call. = FALSE)
  for (i in seq_len(nrow(claims))) {
    claim <- claims[i, , drop = FALSE]
    path <- file.path(evidence_dir, as.character(claim$effect_artifact))
    if (!file.exists(path)) stop("Phase 2 claim effect artifact is missing: ", path, call. = FALSE)
    table <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    required_columns <- unlist(claim[c("effect_column", "denominator", "mcse_column", "ci_lower_column", "ci_upper_column")], use.names = FALSE)
    if (!all(c("scenario_key", "row_key", "metric", required_columns) %in% names(table))) {
      stop("Phase 2 claim effect artifact lacks its registered key/metric/uncertainty columns: ", basename(path), call. = FALSE)
    }
    hit <- as.character(table$scenario_key) == as.character(claim$scenario_key) &
      as.character(table$row_key) == as.character(claim$row_key) &
      as.character(table$metric) == as.character(claim$metric)
    if (sum(hit) != 1L) stop("Phase 2 claim does not resolve to exactly one evidence row: ", claim$claim_id, call. = FALSE)
    values <- suppressWarnings(as.numeric(table[hit, required_columns, drop = TRUE]))
    effect <- values[[1L]]; denominator <- values[[2L]]; mcse <- values[[3L]]
    lower <- values[[4L]]; upper <- values[[5L]]
    if (any(!is.finite(values))) stop("Phase 2 claim effect, denominator, MCSE, or CI is non-finite.", call. = FALSE)
    if (denominator <= 0 || denominator != floor(denominator)) stop("Phase 2 claim denominator must be a positive integer.", call. = FALSE)
    if (mcse < 0) stop("Phase 2 claim MCSE must be nonnegative.", call. = FALSE)
    if (lower > upper || effect < lower || effect > upper) stop("Phase 2 claim confidence interval is inverted or does not contain its effect.", call. = FALSE)
    direction <- as.character(claim$direction)
    positive <- direction %in% c("positive", "higher")
    negative <- direction %in% c("negative", "lower")
    if ((positive && effect <= 0) || (negative && effect >= 0)) stop("Phase 2 claim direction contradicts the effect sign.", call. = FALSE)
    computed_support <- if (positive) {
      if (lower > 0) "supported" else "not_supported"
    } else if (negative) {
      if (upper < 0) "supported" else "not_supported"
    } else {
      "not_applicable"
    }
    if (!identical(as.character(claim$interval_support), computed_support)) stop("Phase 2 claim interval-support label contradicts its confidence interval.", call. = FALSE)
    wording <- as.character(claim$wording_strength)
    direction_class <- if (direction == "no_direction") "no_direction" else "directional"
    allowed <- mvt_phase2_claim_wording_matrix()
    allowed_row <- allowed$direction_class == direction_class &
      allowed$interval_support == computed_support &
      allowed$wording_strength == wording
    if (!any(allowed_row)) {
      stop(
        "Phase 2 claim direction/interval/wording combination is not allowed; every directional expected direction requires a confidence interval wholly supporting it.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

mvt_validate_phase2_claim_outputs <- function(
    claims, evidence_dir,
    attestation_path = mvt_phase2_snapshot_attestation_path(),
    signature_path = mvt_phase2_snapshot_signature_path()) {
  mvt_validate_phase2_claim_evidence(evidence_dir, attestation_path, signature_path)
  mvt_validate_phase2_claim_rows(claims, evidence_dir)
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
    role = "standard_or_focal_context",
    stringsAsFactors = FALSE
  )
  comparator_table$role[comparator_table$method %in% mvt_nearest_neighbor_comparators()] <- "headline_nearest_neighbor"
  comparator_table$role[comparator_table$method == "gamCopula"] <- "legacy_alias"
  comparator_table$role[comparator_table$method == "gamCopula_vine"] <- "targeted_sensitivity"
  comparator_table$role[comparator_table$method == "glmm_slope"] <- "appendix_sensitivity"
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
      "Exactly two headline nearest-neighbor gamCopula workflows; full vine is targeted sensitivity only",
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
      "mvt_nearest_neighbor_comparators(); mvt_comparator_scope_registry(); mvt_fit_gamcopula(); mvt_fit_gamcopula_vine()",
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
      "nearest_neighbor_results.csv; nearest_neighbor_paired_contrasts.csv; comparator_scope_registry.csv; optional targeted full-vine sensitivity",
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
    capability_snapshot_2026_09_01 = mvt_capability_snapshot(),
    comparator_scope_registry = mvt_comparator_scope_registry(),
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

mvt_subprocess_identity_probe <- function() {
  namespace_path <- normalizePath(
    getNamespaceInfo(asNamespace("gamlss.longitudinal"), "path"),
    winslash = "/", mustWork = TRUE
  )
  list(
    package = "gamlss.longitudinal",
    version = as.character(utils::packageVersion("gamlss.longitudinal")),
    namespace_path = namespace_path
  )
}

mvt_parallel_equivalence_tolerances <- function() {
  c(absolute = 1e-10, relative = 1e-10)
}

mvt_run_fit_with_timeout <- function(fun_name, args = list(), timeout = Inf) {
  start <- Sys.time()
  if (is.finite(timeout) && timeout > 0 && requireNamespace("callr", quietly = TRUE)) {
    expected_execution <- getOption(
      "gamlss.longitudinal.mvt.execution_attestation",
      mvt_execution_attestation_contract(configuration_fingerprint = "standalone-subprocess")
    )
    value <- tryCatch(
      callr::r(
        function(fun_name, args, repo_root, expected_execution) {
          setwd(repo_root)
          source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))
          Sys.setenv(GAMLSS_LONGITUDINAL_MVT_SOURCE = "local")
          identity <- mvt_verify_execution_attestation(expected_execution, load = TRUE)
          out <- mvt_elapsed_do_call(get(fun_name, envir = globalenv()), args)
          out$subprocess_attestation <- identity
          out
        },
        args = list(fun_name = fun_name, args = args, repo_root = mvt_repo_root, expected_execution = expected_execution),
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
    identity <- value$subprocess_attestation
    if (!mvt_execution_attestation_matches(identity, expected_execution)) {
      return(list(
        value = simpleError("Timed fit subprocess checkout attestation failed."),
        warnings = character(), elapsed_sec = as.numeric(difftime(Sys.time(), start, units = "secs"))
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

mvt_value_converged <- function(value) {
  if (inherits(value, "error") || is.null(value)) return(FALSE)
  if (inherits(value, "merMod")) {
    messages <- value@optinfo$conv$lme4$messages
    optimizer_code <- value@optinfo$conv$opt %||% 0L
    return(!length(messages) && all(suppressWarnings(as.numeric(optimizer_code)) == 0))
  }
  if (inherits(value, "geeglm")) {
    error_code <- value$geese$error %||% NA_integer_
    return(length(error_code) == 1L && is.finite(error_code) && error_code == 0L)
  }
  if (is.list(value) || is.environment(value)) {
    if (!is.null(value[["converged"]])) return(isTRUE(value[["converged"]]))
    if (inherits(value, "mvt_gamcopula_fit") && !is.null(value[["margin_fit"]][["converged"]])) {
      return(isTRUE(value[["margin_fit"]][["converged"]]))
    }
  }
  TRUE
}

mvt_value_stop_reason <- function(value, converged, error_msg = NA_character_) {
  if (inherits(value, "error")) return("error")
  reason <- if (is.list(value) || is.environment(value)) {
    value[["stop_reason"]] %||% value[["optimizer_stop_reason"]] %||% NULL
  } else {
    NULL
  }
  if (!is.null(reason) && length(reason) && nzchar(as.character(reason[[1L]]))) return(as.character(reason[[1L]]))
  if (isTRUE(converged)) "converged" else if (!is.na(error_msg) && nzchar(error_msg)) "error" else "nonconverged"
}

mvt_method_result_row <- function(row, spec, method, fit_result, pred = NULL, extra = list()) {
  converged <- mvt_value_converged(fit_result$value)
  success <- !inherits(fit_result$value, "error") && !is.null(fit_result$value) && converged
  error_msg <- if (success) {
    NA_character_
  } else if (inherits(fit_result$value, "error")) {
    conditionMessage(fit_result$value)
  } else if (!is.null(fit_result$value) && !converged) {
    "Optimizer did not converge."
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
    attempted = TRUE,
    available = TRUE,
    success = success,
    converged = converged,
    retained = success,
    stop_reason = mvt_value_stop_reason(fit_result$value, converged, error_msg),
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
  subprocess_identity <- fit_result$subprocess_attestation %||% NULL
  out$subprocess_package_verified <- if (is.null(subprocess_identity)) NA else isTRUE(subprocess_identity$verified)
  out$subprocess_source_sha256 <- subprocess_identity$verified_source_sha256 %||% NA_character_
  out$subprocess_namespace_path <- subprocess_identity$loaded_namespace_path %||% NA_character_
  out$subprocess_dependency_fingerprint <- subprocess_identity$dependency_fingerprint %||% NA_character_
  out$subprocess_runtime_fingerprint <- subprocess_identity$runtime_fingerprint %||% NA_character_
  out$subprocess_configuration_fingerprint <- subprocess_identity$configuration_fingerprint %||% NA_character_
  out$subprocess_full_verified <- if (is.null(subprocess_identity)) NA else isTRUE(subprocess_identity$full_verified)
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
  fit_converged <- vapply(as.character(fit$results$method), function(method) {
    candidate <- fit$fits[[method]] %||% NULL
    mvt_value_converged(candidate)
  }, logical(1L))
  fit$results$converged <- fit_converged
  fit$results$success <- (fit$results$success %in% TRUE) & fit_converged
  fit$results$retained <- fit$results$success
  fit$results$stop_reason <- ifelse(fit_converged, "converged", "optimizer_nonconvergence")
  if (!"error" %in% names(fit$results)) fit$results$error <- NA_character_
  failed_convergence <- !fit_converged & (is.na(fit$results$error) | !nzchar(trimws(as.character(fit$results$error))))
  fit$results$error[failed_convergence] <- "Optimizer did not converge."
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
    mvt_registered_timeout("GAMLSS_LONGITUDINAL_MVT_GEE_UNSTRUCTURED_TIMEOUT_SEC", if (row$n_time >= 50L) 30 else 60)
  } else {
    mvt_registered_timeout("GAMLSS_LONGITUDINAL_MVT_GEE_TIMEOUT_SEC", 300)
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
  subprocess_attestation <- NULL
  value <- if (is.finite(timeout) && timeout > 0 && requireNamespace("callr", quietly = TRUE)) {
    expected_execution <- getOption(
      "gamlss.longitudinal.mvt.execution_attestation",
      mvt_execution_attestation_contract(configuration_fingerprint = "standalone-gee")
    )
    child <- tryCatch(callr::r(function(dat, spec, corstr, root, setup_path, expected_execution) {
      source(setup_path, local = .GlobalEnv)
      Sys.setenv(GAMLSS_LONGITUDINAL_MVT_SOURCE = "local")
      identity <- mvt_verify_execution_attestation(expected_execution, load = TRUE)
      result <- gamlss.longitudinal::benchmark_standard_models(
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
      list(result = result, attestation = identity)
    }, args = list(
      dat = dat, spec = spec, corstr = corstr,
      root = mvt_repo_root,
      setup_path = file.path(mvt_script_dir, "00-multivariate-setup.R"),
      expected_execution = expected_execution
    ), timeout = timeout), error = function(e) e)
    if (inherits(child, "error")) child else {
      subprocess_attestation <- child$attestation
      if (!mvt_execution_attestation_matches(subprocess_attestation, expected_execution)) {
        simpleError("Timed GEE subprocess checkout attestation failed.")
      } else child$result
    }
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
  gee_converged <- mvt_value_converged(value$fits$gee %||% NULL)
  value$results$converged <- gee_converged
  value$results$success <- (value$results$success %in% TRUE) & gee_converged
  value$results$retained <- value$results$success
  value$results$stop_reason <- if (gee_converged) "converged" else "optimizer_nonconvergence"
  if (!gee_converged) value$results$error <- "Optimizer did not converge."
  value$results$subprocess_package_verified <- if (is.null(subprocess_attestation)) NA else isTRUE(subprocess_attestation$verified)
  value$results$subprocess_source_sha256 <- subprocess_attestation$verified_source_sha256 %||% NA_character_
  value$results$subprocess_namespace_path <- subprocess_attestation$loaded_namespace_path %||% NA_character_
  value$results$subprocess_dependency_fingerprint <- subprocess_attestation$dependency_fingerprint %||% NA_character_
  value$results$subprocess_runtime_fingerprint <- subprocess_attestation$runtime_fingerprint %||% NA_character_
  value$results$subprocess_configuration_fingerprint <- subprocess_attestation$configuration_fingerprint %||% NA_character_
  value$results$subprocess_full_verified <- if (is.null(subprocess_attestation)) NA else isTRUE(subprocess_attestation$full_verified)
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
  nsim <- suppressWarnings(as.integer(nsim))
  if (length(nsim) != 1L || is.na(nsim) || nsim <= 0L) {
    return(mvt_variogram_empty(row, spec, method, "Variogram simulation skipped because nsim <= 0."))
  }
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

mvt_run_variogram_score <- function(
    dat, fit, method, row, spec, nsim = 50L, seed = NULL,
    timeout = mvt_env_num("GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_TIMEOUT_SEC", 300)) {
  nsim <- suppressWarnings(as.integer(nsim))
  if (length(nsim) != 1L || is.na(nsim) || nsim <= 0L) {
    return(mvt_variogram_score(dat, fit, method, row, spec, nsim = nsim, seed = seed))
  }
  if (is.finite(timeout) && timeout > 0 && requireNamespace("callr", quietly = TRUE)) {
    value <- tryCatch(
      callr::r(
        function(dat, fit, method, row, spec, nsim, seed, repo_root) {
          setwd(repo_root)
          source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))
          Sys.setenv(GAMLSS_LONGITUDINAL_MVT_SOURCE = "local")
          mvt_load_package()
          mvt_variogram_score(dat, fit, method, row, spec, nsim = nsim, seed = seed)
        },
        args = list(
          dat = dat, fit = fit, method = method, row = row, spec = spec,
          nsim = nsim, seed = seed, repo_root = mvt_repo_root
        ),
        timeout = timeout
      ),
      error = function(e) e
    )
    if (inherits(value, "error")) {
      out <- mvt_variogram_empty(row, spec, method, conditionMessage(value))
      out$variogram_nsim <- nsim
      return(out)
    }
    return(value)
  }
  mvt_variogram_score(dat, fit, method, row, spec, nsim = nsim, seed = seed)
}

mvt_enforce_attempt_schema <- function(results) {
  if (!is.data.frame(results) || !nrow(results)) return(results)
  if (!"attempted" %in% names(results)) results$attempted <- TRUE
  results$attempted[is.na(results$attempted)] <- TRUE
  results$attempted <- results$attempted %in% TRUE
  if (!"converged" %in% names(results)) results$converged <- results$success %in% TRUE
  results$converged <- results$converged %in% TRUE
  if (!"success" %in% names(results)) results$success <- results$converged
  results$success <- results$success %in% TRUE & results$converged
  if (!"retained" %in% names(results)) results$retained <- results$success
  results$retained <- results$retained %in% TRUE & results$converged & results$success
  if (!"stop_reason" %in% names(results)) results$stop_reason <- ifelse(results$converged, "converged", "nonconverged")
  results$stop_reason[is.na(results$stop_reason) | !nzchar(trimws(results$stop_reason))] <- ifelse(
    results$converged[is.na(results$stop_reason) | !nzchar(trimws(results$stop_reason))], "converged", "nonconverged"
  )
  if (!"error" %in% names(results)) results$error <- NA_character_
  if (!"warning" %in% names(results)) results$warning <- NA_character_
  if (!"elapsed_sec" %in% names(results)) results$elapsed_sec <- NA_real_
  if (!"status_class" %in% names(results)) results$status_class <- NA_character_
  missing_status <- is.na(results$status_class) | !nzchar(trimws(as.character(results$status_class)))
  missing_status <- missing_status |
    (!results$success & tolower(as.character(results$status_class)) %in% c("ok", "success")) |
    (results$success & tolower(as.character(results$status_class)) %in% c("error", "timeout"))
  if (any(missing_status)) results$status_class[missing_status] <- vapply(which(missing_status), function(i) {
    mvt_fit_status_class(results$success[[i]], results$error[[i]], results$warning[[i]])
  }, character(1L))
  if (!"failure_reason_short" %in% names(results)) results$failure_reason_short <- NA_character_
  missing_reason <- is.na(results$failure_reason_short) | !nzchar(trimws(as.character(results$failure_reason_short)))
  missing_reason <- missing_reason | (!results$success & results$failure_reason_short == "none")
  if (any(missing_reason)) results$failure_reason_short[missing_reason] <- mapply(
    mvt_failure_reason_short,
    results$success[missing_reason], results$status_class[missing_reason],
    results$error[missing_reason], results$warning[missing_reason], USE.NAMES = FALSE
  )
  accuracy_cols <- intersect(
    c("mae", "rmse", grep("^benchmark_", names(results), value = TRUE), "logLik", "AIC", "BIC"),
    names(results)
  )
  if (length(accuracy_cols) && any(!results$retained)) results[!results$retained, accuracy_cols] <- NA
  results
}

mvt_status_from_results <- function(results) {
  if (!is.data.frame(results) || nrow(results) == 0L) return(data.frame())
  out <- mvt_enforce_attempt_schema(results)
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
      "method", "attempted", "available", "success", "converged", "retained",
      "stop_reason", "status_class", "elapsed_sec", "nobs",
      "failure_reason_short", "warning", "error"
    ),
    names(out)
  )
  out[, keep, drop = FALSE]
}

mvt_run_case <- function(row, seed_base = 20260818L, require_gamcopula = TRUE, case_seed = NULL) {
  if (isTRUE(require_gamcopula)) {
    mvt_require_namespaces("gamCopula", strict = TRUE)
  }
  families <- mvt_family_specs(include_special = TRUE)
  deps <- mvt_dependence_specs(include_appendix = TRUE)
  spec <- families[[row$family_name]]
  dep <- mvt_resize_dependence(deps[[row$dependence_name]], row$n_time)
  active_comparators <- mvt_active_comparators()
  derived_seed <- seed_base + match(row$family_name, names(families)) * 100000L +
    match(row$dependence_name, names(deps)) * 1000L + row$n_time * 10L + row$rep
  seed <- if (is.null(case_seed)) as.integer(derived_seed) else as.integer(case_seed)
  if (length(seed) != 1L || is.na(seed) || !identical(seed, as.integer(derived_seed))) {
    stop("case_seed does not match the registered deterministic seed formula for ", row$case_id, call. = FALSE)
  }

  sim_capture <- mvt_elapsed_capture(mvt_simulate_case(row, seed = seed))
  if (inherits(sim_capture$value, "error")) {
    benchmark_results <- mvt_bind_rows_fill(lapply(
      active_comparators,
      function(method) mvt_method_result_row(row, spec, method, sim_capture)
    ))
    fit_status <- mvt_status_from_results(benchmark_results)
    return(list(
      fit_status = fit_status,
      benchmark_results = benchmark_results,
      coefficient_results = data.frame(),
      dependence_recovery = data.frame(),
      variogram_scores = data.frame(),
      runtime = benchmark_results[intersect(
        c("case_id", "scenario", "family", "rep", "method", "attempted", "success", "converged", "retained", "stop_reason", "status_class", "elapsed_sec", "failure_reason_short", "warning", "error"),
        names(benchmark_results)
      )]
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
  if (length(mvt_standard_comparators(active_comparators)) > 0L) {
    standard <- mvt_elapsed_capture(mvt_fit_standard_models(dat, spec, row, dep))
    if (!inherits(standard$value, "error")) {
      result_rows[[length(result_rows) + 1L]] <- standard$value$results
      coef_rows[[length(coef_rows) + 1L]] <- mvt_annotate_coefficients(standard$value$coefficients, spec)
      standard_fits <- standard$value$fits
      standard_fits <- standard_fits[names(standard_fits) %in% mvt_standard_comparators(active_comparators)]
      for (fit_name in names(standard_fits)) {
        if (!is.null(standard_fits[[fit_name]])) {
          vario_rows[[length(vario_rows) + 1L]] <- mvt_run_variogram_score(
            model_dat, standard_fits[[fit_name]], fit_name, row, spec,
            nsim = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM", 20L),
            seed = seed + 7000L + match(fit_name, names(standard_fits))
          )
        }
      }
    } else {
      failed_standard <- list(value = standard$value, warnings = standard$warnings, elapsed_sec = standard$elapsed_sec)
      result_rows[[length(result_rows) + 1L]] <- mvt_bind_rows_fill(lapply(
        mvt_standard_comparators(active_comparators),
        function(method) mvt_method_result_row(
          row, spec, method, failed_standard, extra = list(y = y)
        )
      ))
    }
  }

  for (corstr in intersect(dep$gee_correlations, mvt_gee_comparators(active_comparators))) {
    gee <- mvt_run_one_gee(dat, spec, row, corstr)
    result_rows[[length(result_rows) + 1L]] <- gee$results
    coef_rows[[length(coef_rows) + 1L]] <- mvt_annotate_coefficients(gee$coefficients, spec)
    if (!is.null(gee$fit)) {
      vario_rows[[length(vario_rows) + 1L]] <- mvt_run_variogram_score(
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
      vario_rows[[length(vario_rows) + 1L]] <- mvt_run_variogram_score(
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
      timeout = mvt_env_num(
        "GAMLSS_LONGITUDINAL_MVT_GAMLSS_TIMEOUT_SEC",
        mvt_env_num("GAMLSS_LONGITUDINAL_MVT_PRIMARY_MAX_ELAPSED_SEC", 180)
      )
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
      vario_rows[[length(vario_rows) + 1L]] <- mvt_run_variogram_score(
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
      vario_rows[[length(vario_rows) + 1L]] <- mvt_run_variogram_score(
        model_dat, gc_fit$value, gc_method, row, spec,
        nsim = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM", 20L),
        seed = seed + if (identical(gc_method, "gamCopula_vine")) 9700L else if (identical(gc_method, "gamCopula_vine_simplified")) 9650L else 9500L
      )
    }
  }

  benchmark_results <- mvt_enforce_attempt_schema(mvt_bind_rows_fill(result_rows))
  fit_status <- mvt_status_from_results(benchmark_results)
  coefficient_results <- mvt_bind_rows_fill(coef_rows)
  dependence_recovery <- mvt_bind_rows_fill(dep_rows)
  variogram_scores <- mvt_bind_rows_fill(vario_rows)
  retained_keys <- paste(
    benchmark_results$case_id[benchmark_results$retained %in% TRUE],
    benchmark_results$method[benchmark_results$retained %in% TRUE], sep = "\r"
  )
  retain_accuracy_rows <- function(table) {
    if (!is.data.frame(table) || !nrow(table)) return(table)
    key <- paste(table$case_id, table$method, sep = "\r")
    table[key %in% retained_keys, , drop = FALSE]
  }
  coefficient_results <- retain_accuracy_rows(coefficient_results)
  dependence_recovery <- retain_accuracy_rows(dependence_recovery)
  if (nrow(variogram_scores)) {
    variogram_key <- paste(variogram_scores$case_id, variogram_scores$method, sep = "\r")
    failed_vario <- !variogram_key %in% retained_keys
    variogram_scores$variogram_score_p05[failed_vario] <- NA_real_
    variogram_scores$variogram_score_p2[failed_vario] <- NA_real_
    variogram_scores$error[failed_vario] <- "Fit was not retained; variogram accuracy excluded."
  }
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
  runtime <- benchmark_results[, intersect(
    c("case_id", "scenario", "family", "rep", "method", "attempted", "success", "converged", "retained", "stop_reason", "status_class", "elapsed_sec", "failure_reason_short", "warning", "error"),
    names(benchmark_results)
  ), drop = FALSE]

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

mvt_checkpoint_schema_version <- 3L
mvt_identity_algorithm <- "SHA-256"
mvt_identity_version <- "mvt-identity-v3"
mvt_producer_id <- "module09-multivariate-longitudinal-benchmark"
mvt_producer_version <- "3"
mvt_active_leases <- new.env(parent = emptyenv())

mvt_hash_object <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for canonical SHA-256 identity.", call. = FALSE)
  }
  digest::digest(x, algo = "sha256", serialize = TRUE, serializeVersion = 3)
}

mvt_sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for canonical SHA-256 identity.", call. = FALSE)
  }
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

mvt_files_fingerprint <- function(paths, root = mvt_repo_root) {
  paths <- sort(unique(normalizePath(paths[file.exists(paths)], winslash = "/", mustWork = TRUE)))
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  labels <- ifelse(startsWith(paths, paste0(root, "/")), substring(paths, nchar(root) + 2L), paths)
  hashes <- vapply(paths, mvt_sha256_file, character(1L))
  manifest <- paste(labels, hashes, sep = "\t", collapse = "\n")
  mvt_hash_object(list(
    algorithm = mvt_identity_algorithm,
    identity_version = mvt_identity_version,
    sorted_manifest = manifest
  ))
}

mvt_runtime_identity <- function() {
  soft <- extSoftVersion()
  blas <- unname(soft["BLAS"] %||% NA_character_)
  if (!length(blas) || is.na(blas) || !nzchar(blas)) blas <- "default"
  lapack <- tryCatch(La_library(), error = function(e) "")
  if (!length(lapack) || is.na(lapack) || !nzchar(lapack)) lapack <- "default"
  list(
    r_version = R.version.string,
    platform = R.version$platform,
    rng_kind = RNGkind(),
    blas = blas,
    lapack = lapack
  )
}

mvt_hardware_identity <- function() {
  ram_bytes <- NA_real_
  if (.Platform$OS.type == "windows") {
    ram_text <- tryCatch(system2("powershell", c("-NoProfile", "-Command", "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory"), stdout = TRUE, stderr = FALSE), error = function(e) character())
    if (length(ram_text)) ram_bytes <- suppressWarnings(as.numeric(trimws(ram_text[[1L]])))
  }
  list(
    hostname = unname(Sys.info()[["nodename"]]),
    cpu_model = Sys.getenv("PROCESSOR_IDENTIFIER", unset = unname(Sys.info()[["machine"]])),
    logical_cores = parallel::detectCores(logical = TRUE),
    physical_cores = parallel::detectCores(logical = FALSE),
    ram_bytes = ram_bytes,
    parent_pid = as.integer(Sys.getpid())
  )
}

mvt_checkout_package_identity <- function() {
  root <- normalizePath(mvt_repo_root, winslash = "/", mustWork = TRUE)
  files <- c(
    list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE),
    file.path(root, c("DESCRIPTION", "NAMESPACE"))
  )
  files <- sort(normalizePath(files[file.exists(files)], winslash = "/", mustWork = TRUE))
  description <- read.dcf(file.path(root, "DESCRIPTION"))
  list(
    algorithm = mvt_identity_algorithm,
    identity_version = mvt_identity_version,
    package = "gamlss.longitudinal",
    version = unname(description[1L, "Version"]),
    checkout_path = root,
    source_sha256 = mvt_files_fingerprint(files, root),
    source_files = substring(files, nchar(root) + 2L),
    source_file_count = length(files)
  )
}

mvt_dependency_packages <- function() {
  c(
    "gamlss.longitudinal", "gamlss", "gamlss.dist", "gamCopula", "VineCopula",
    "mvtnorm", "geepack", "lme4", "mgcv", "callr", "scoringRules"
  )
}

mvt_dependency_source_sha256 <- function(package, namespace_path) {
  if (!nzchar(namespace_path) || !dir.exists(namespace_path)) return(NA_character_)
  if (identical(package, "gamlss.longitudinal") &&
      identical(normalizePath(namespace_path, winslash = "/", mustWork = TRUE),
        normalizePath(mvt_repo_root, winslash = "/", mustWork = TRUE))) {
    return(mvt_checkout_package_identity()$source_sha256)
  }
  candidates <- c(
    file.path(namespace_path, c("DESCRIPTION", "NAMESPACE")),
    list.files(file.path(namespace_path, "R"), full.names = TRUE),
    list.files(file.path(namespace_path, "libs"), recursive = TRUE, full.names = TRUE),
    list.files(file.path(namespace_path, "Meta"), pattern = "package[.]rds$", full.names = TRUE)
  )
  candidates <- candidates[file.exists(candidates) & !(file.info(candidates)$isdir %in% TRUE)]
  if (!length(candidates)) return(NA_character_)
  mvt_files_fingerprint(candidates, namespace_path)
}

mvt_dependency_identity_table <- function(
    packages = mvt_dependency_packages(),
    focal_identity = mvt_checkout_package_identity(),
    focal_checkout = TRUE) {
  rows <- lapply(sort(unique(packages)), function(package) {
    available <- requireNamespace(package, quietly = TRUE)
    if (identical(package, focal_identity$package) && isTRUE(focal_checkout)) {
      path <- focal_identity$checkout_path
      version <- focal_identity$version
      source_sha256 <- focal_identity$source_sha256
      available <- TRUE
    } else if (available) {
      path <- normalizePath(getNamespaceInfo(asNamespace(package), "path"), winslash = "/", mustWork = TRUE)
      version <- as.character(utils::packageVersion(package))
      source_sha256 <- mvt_dependency_source_sha256(package, path)
    } else {
      path <- NA_character_
      version <- NA_character_
      source_sha256 <- NA_character_
    }
    data.frame(
      package = package, available = available, version = version,
      namespace_path = path, source_sha256 = source_sha256,
      stringsAsFactors = FALSE
    )
  })
  mvt_canonical_dependency_identity(do.call(rbind, rows))
}

mvt_canonical_dependency_identity <- function(x) {
  required <- c("package", "available", "version", "namespace_path", "source_sha256")
  if (!is.data.frame(x) || !identical(sort(names(x)), sort(required))) {
    stop("Dependency identity table has an invalid schema.", call. = FALSE)
  }
  x <- x[, required, drop = FALSE]
  x$package <- as.character(x$package)
  x$available <- as.logical(x$available)
  x$version <- as.character(x$version)
  x$namespace_path <- as.character(x$namespace_path)
  x$source_sha256 <- as.character(x$source_sha256)
  x <- x[order(x$package), , drop = FALSE]
  rownames(x) <- NULL
  x
}

mvt_runtime_fingerprint <- function(runtime = mvt_runtime_identity()) {
  mvt_hash_object(runtime)
}

mvt_execution_attestation_contract <- function(
    fingerprints = mvt_checkpoint_fingerprints(), configuration_fingerprint = "standalone") {
  dependencies <- fingerprints$dependency_identity %||% fingerprints$package_versions
  if (is.null(dependencies) || !is.data.frame(dependencies) || !nrow(dependencies)) {
    dependencies <- mvt_dependency_identity_table()
  }
  dependencies <- mvt_canonical_dependency_identity(dependencies)
  list(
    identity_algorithm = fingerprints$algorithm %||% mvt_identity_algorithm,
    identity_version = fingerprints$identity_version %||% mvt_identity_version,
    package_identity = fingerprints$package_identity %||% mvt_checkout_package_identity(),
    dependency_identity = dependencies,
    dependency_fingerprint = mvt_hash_object(dependencies),
    runtime_identity = fingerprints$runtime_identity %||% mvt_runtime_identity(),
    runtime_fingerprint = mvt_runtime_fingerprint(fingerprints$runtime_identity %||% mvt_runtime_identity()),
    configuration_fingerprint = configuration_fingerprint
  )
}

mvt_verify_execution_attestation <- function(expected, load = TRUE) {
  focal <- mvt_verify_checkout_package(expected$package_identity, load = load)
  dependencies <- mvt_dependency_identity_table(
    packages = as.character(expected$dependency_identity$package),
    focal_identity = expected$package_identity,
    focal_checkout = TRUE
  )
  dependencies <- mvt_canonical_dependency_identity(dependencies)
  expected_dependencies <- mvt_canonical_dependency_identity(expected$dependency_identity)
  dependency_fingerprint <- mvt_hash_object(dependencies)
  expected_dependency_fingerprint <- mvt_hash_object(expected_dependencies)
  runtime <- mvt_runtime_identity()
  runtime_fingerprint <- mvt_runtime_fingerprint(runtime)
  expected_runtime_fingerprint <- mvt_runtime_fingerprint(expected$runtime_identity)
  valid <- identical(dependencies, expected_dependencies) &&
    identical(dependency_fingerprint, expected$dependency_fingerprint) &&
    identical(expected_dependency_fingerprint, expected$dependency_fingerprint) &&
    identical(runtime_fingerprint, expected$runtime_fingerprint) &&
    identical(expected_runtime_fingerprint, expected$runtime_fingerprint) &&
    is.character(expected$configuration_fingerprint) && length(expected$configuration_fingerprint) == 1L &&
    nzchar(expected$configuration_fingerprint)
  if (!valid) {
    stop("Full dependency/runtime/configuration execution attestation failed.", call. = FALSE)
  }
  c(focal, list(
    dependency_identity = dependencies,
    dependency_fingerprint = dependency_fingerprint,
    runtime_identity = runtime,
    runtime_fingerprint = runtime_fingerprint,
    configuration_fingerprint = expected$configuration_fingerprint,
    full_verified = TRUE
  ))
}

mvt_execution_attestation_matches <- function(attestation, expected) {
  is.list(attestation) && isTRUE(attestation$full_verified) && isTRUE(attestation$verified) &&
    identical(attestation$verified_source_sha256, expected$package_identity$source_sha256) &&
    identical(attestation$loaded_namespace_path, expected$package_identity$checkout_path) &&
    identical(attestation$loaded_version, expected$package_identity$version) &&
    identical(
      mvt_canonical_dependency_identity(attestation$dependency_identity),
      mvt_canonical_dependency_identity(expected$dependency_identity)
    ) &&
    identical(attestation$dependency_fingerprint, expected$dependency_fingerprint) &&
    identical(attestation$runtime_identity, expected$runtime_identity) &&
    identical(attestation$runtime_fingerprint, expected$runtime_fingerprint) &&
    identical(attestation$configuration_fingerprint, expected$configuration_fingerprint)
}

mvt_verify_checkout_package <- function(expected = mvt_checkout_package_identity(), load = TRUE) {
  if (isTRUE(load)) {
    if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required for checkout loading.", call. = FALSE)
    pkgload::load_all(expected$checkout_path, quiet = TRUE, export_all = TRUE, helpers = FALSE)
  }
  actual <- mvt_checkout_package_identity()
  namespace_path <- normalizePath(
    getNamespaceInfo(asNamespace(expected$package), "path"), winslash = "/", mustWork = TRUE
  )
  loaded_version <- as.character(utils::packageVersion(expected$package))
  valid <- identical(actual$source_sha256, expected$source_sha256) &&
    identical(actual$version, expected$version) && identical(loaded_version, expected$version) &&
    identical(namespace_path, expected$checkout_path)
  if (!valid) {
    stop(
      "Checked-out package attestation failed: expected ", expected$version, " at ", expected$checkout_path,
      " with SHA-256 ", expected$source_sha256, "; loaded ", loaded_version, " at ", namespace_path,
      " with SHA-256 ", actual$source_sha256, ".", call. = FALSE
    )
  }
  c(expected, list(
    verified = TRUE,
    loaded_namespace_path = namespace_path,
    loaded_version = loaded_version,
    verified_source_sha256 = actual$source_sha256,
    libpaths = paste(normalizePath(.libPaths(), winslash = "/", mustWork = FALSE), collapse = ";"),
    runtime_identity = mvt_runtime_identity()
  ))
}

mvt_attestation_row <- function(
    attestation, role = "worker", pid = Sys.getpid(),
    setup_path = file.path(mvt_script_dir, "00-multivariate-setup.R")) {
  runtime <- attestation$runtime_identity %||% mvt_runtime_identity()
  gc_state <- attestation$gc_state %||% gc()
  data.frame(
    role = role, pid = as.integer(pid), verified = isTRUE(attestation$verified),
    setup_path = normalizePath(setup_path, winslash = "/", mustWork = FALSE),
    package = attestation$package, package_version = attestation$loaded_version %||% attestation$version,
    checkout_path = attestation$checkout_path,
    loaded_namespace_path = attestation$loaded_namespace_path %||% NA_character_,
    source_sha256 = attestation$source_sha256,
    verified_source_sha256 = attestation$verified_source_sha256 %||% NA_character_,
    dependency_fingerprint = attestation$dependency_fingerprint %||% NA_character_,
    runtime_fingerprint = attestation$runtime_fingerprint %||% mvt_runtime_fingerprint(runtime),
    configuration_fingerprint = attestation$configuration_fingerprint %||% NA_character_,
    full_verified = isTRUE(attestation$full_verified),
    identity_algorithm = attestation$algorithm, identity_version = attestation$identity_version,
    r_version = runtime$r_version, platform = runtime$platform,
    rng_kind = paste(runtime$rng_kind, collapse = "/"), blas = runtime$blas, lapack = runtime$lapack,
    hostname = unname(Sys.info()[["nodename"]]),
    libpaths = attestation$libpaths %||% paste(.libPaths(), collapse = ";"),
    gc_peak_ncells = as.numeric(gc_state["Ncells", "max used"]),
    gc_peak_vcells = as.numeric(gc_state["Vcells", "max used"]),
    stringsAsFactors = FALSE
  )
}

mvt_checkpoint_fingerprints <- function() {
  producer_files <- file.path(
    mvt_script_dir,
    c("00-multivariate-setup.R", "01-run-pilot-grid.R", "02-run-main-grid.R")
  )
  checkout <- mvt_checkout_package_identity()
  dependencies <- mvt_dependency_identity_table(
    packages = mvt_dependency_packages(), focal_identity = checkout, focal_checkout = TRUE
  )
  dependencies <- mvt_canonical_dependency_identity(dependencies)
  list(
    algorithm = mvt_identity_algorithm,
    identity_version = mvt_identity_version,
    producer_id = mvt_producer_id,
    producer_version = mvt_producer_version,
    producer_fingerprint = mvt_hash_object(list(
      producer_id = mvt_producer_id,
      producer_version = mvt_producer_version,
      sorted_source_fingerprint = mvt_files_fingerprint(producer_files)
    )),
    code_fingerprint = checkout$source_sha256,
    package_fingerprint = mvt_hash_object(dependencies),
    package_identity = checkout,
    dependency_identity = dependencies,
    package_versions = dependencies,
    runtime_identity = mvt_runtime_identity()
  )
}

mvt_checkpoint_configuration <- function(
    seed_base, require_gamcopula, active_comparators = mvt_active_comparators()) {
  list(
    seed_base = as.integer(seed_base),
    require_gamcopula = isTRUE(require_gamcopula),
    active_comparators = as.character(active_comparators),
    source = mvt_env("GAMLSS_LONGITUDINAL_MVT_SOURCE", "installed"),
    compute_vcov = mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_COMPUTE_VCOV", TRUE),
    max_inner_iter = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_MAX_INNER_ITER", 60L),
    max_outer_iter = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_MAX_OUTER_ITER", 60L),
    repair_correlation = mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_REPAIR_CORRELATION", TRUE),
    summary_vcov = mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_SUMMARY_VCOV", TRUE),
    primary_timeout_sec = mvt_registered_timeout("GAMLSS_LONGITUDINAL_MVT_PRIMARY_MAX_ELAPSED_SEC", 180),
    gee_timeout_sec = mvt_registered_timeout("GAMLSS_LONGITUDINAL_MVT_GEE_TIMEOUT_SEC", 300),
    gee_unstructured_timeout_sec = mvt_registered_timeout("GAMLSS_LONGITUDINAL_MVT_GEE_UNSTRUCTURED_TIMEOUT_SEC", 60),
    gamcopula_timeout_sec = mvt_registered_timeout("GAMLSS_LONGITUDINAL_MVT_GAMCOPULA_TIMEOUT_SEC", 180),
    gamcopula_markov_timeout_sec = mvt_registered_timeout("GAMLSS_LONGITUDINAL_MVT_GAMCOPULA_MARKOV_TIMEOUT_SEC", 180),
    gamcopula_vine_timeout_sec = mvt_registered_timeout("GAMLSS_LONGITUDINAL_MVT_GAMCOPULA_VINE_TIMEOUT_SEC", 180),
    gamcopula_pit_seed = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_GAMCOPULA_PIT_SEED", 7349L),
    variogram_nsim = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM", 20L),
    vine_dependence_nsim = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_VINE_DEPENDENCE_NSIM", 500L),
    runtime_identity = mvt_runtime_identity()
  )
}

mvt_timeout_contract <- function(configuration, methods = configuration$active_comparators) {
  required <- c()
  if (any(grepl("^gee_", methods))) required <- c(required, "gee_timeout_sec")
  if ("gee_unstructured" %in% methods) required <- c(required, "gee_unstructured_timeout_sec")
  if ("gamlss.longitudinal" %in% methods) required <- c(required, "primary_timeout_sec")
  if ("gamCopula_markov" %in% methods) required <- c(required, "gamcopula_markov_timeout_sec")
  if (any(c("gamCopula_vine", "gamCopula_vine_simplified") %in% methods)) {
    required <- c(required, "gamcopula_vine_timeout_sec")
  }
  values <- suppressWarnings(as.numeric(unlist(configuration[required], use.names = FALSE)))
  list(
    valid = length(required) > 0L && length(values) == length(required) && all(is.finite(values) & values > 0),
    required = required,
    values = stats::setNames(values, required)
  )
}

mvt_grid_is_exact_production <- function(grid) {
  production <- mvt_phase2_production_contract()
  required <- c("case_id", "time_name", "n_time", "family_name", "dependence_name", "rep")
  if (!is.data.frame(grid) || !all(required %in% names(grid))) return(FALSE)
  expected <- expand.grid(
    time_name = production$time_names,
    dependence_name = production$dependence,
    rep = production$reps,
    family_name = production$families,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  expected$n_time <- production$n_time
  expected$case_id <- mapply(
    function(time_name, dependence_name, family_name, rep) {
      paste(time_name, dependence_name, family_name, sprintf("rep%03d", as.integer(rep)), sep = "__")
    },
    expected$time_name, expected$dependence_name, expected$family_name, expected$rep,
    USE.NAMES = FALSE
  )
  fields <- required
  actual <- grid[order(grid$case_id), fields, drop = FALSE]
  expected <- expected[order(expected$case_id), fields, drop = FALSE]
  rownames(actual) <- rownames(expected) <- NULL
  nrow(actual) == production$cases && identical(
    lapply(actual, as.character), lapply(expected, as.character)
  )
}

mvt_case_seed <- function(row, seed_base = 20260818L) {
  families <- mvt_family_specs(include_special = TRUE)
  deps <- mvt_dependence_specs(include_appendix = TRUE)
  as.integer(seed_base + match(row$family_name, names(families)) * 100000L +
    match(row$dependence_name, names(deps)) * 1000L + row$n_time * 10L + row$rep)
}

mvt_canonical_grid_order <- function(grid) {
  time_rank <- if ("time_name" %in% names(grid)) match(grid$time_name, names(mvt_time_specs())) else grid$n_time
  dependence_rank <- if ("dependence_name" %in% names(grid)) {
    match(grid$dependence_name, names(mvt_dependence_specs(include_appendix = TRUE)))
  } else {
    seq_len(nrow(grid))
  }
  family_rank <- if ("family_name" %in% names(grid)) {
    match(grid$family_name, names(mvt_family_specs(include_special = TRUE)))
  } else {
    seq_len(nrow(grid))
  }
  order(
    time_rank, dependence_rank, suppressWarnings(as.integer(grid$rep)),
    family_rank, as.character(grid$case_id),
    na.last = TRUE
  )
}

mvt_prepare_tasks <- function(
    grid, seed_base, require_gamcopula, fingerprints = mvt_checkpoint_fingerprints(),
    workers_requested = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_WORKERS", 1L),
    run_started_at = Sys.time()) {
  if (!is.data.frame(grid) || !nrow(grid) || !"case_id" %in% names(grid)) {
    stop("The multivariate task grid must contain at least one case_id.", call. = FALSE)
  }
  if (anyDuplicated(as.character(grid$case_id))) {
    stop("The multivariate task grid contains duplicate case_id values.", call. = FALSE)
  }
  configuration <- mvt_checkpoint_configuration(seed_base, require_gamcopula)
  configuration_fingerprint <- mvt_hash_object(configuration)
  execution_attestation <- mvt_execution_attestation_contract(
    fingerprints = fingerprints,
    configuration_fingerprint = configuration_fingerprint
  )
  execution_fingerprint <- mvt_hash_object(execution_attestation)
  planned_methods <- as.character(configuration$active_comparators)
  lapply(seq_len(nrow(grid)), function(i) {
    row <- grid[i, , drop = FALSE]
    case_seed <- mvt_case_seed(row, seed_base)
    case_contract <- list(
      schema_version = mvt_checkpoint_schema_version,
      identity_algorithm = fingerprints$algorithm %||% mvt_identity_algorithm,
      identity_version = fingerprints$identity_version %||% mvt_identity_version,
      producer_id = fingerprints$producer_id %||% mvt_producer_id,
      producer_version = fingerprints$producer_version %||% mvt_producer_version,
      row = lapply(row, function(x) as.character(x[[1L]])),
      case_seed = case_seed,
      configuration_fingerprint = configuration_fingerprint,
      producer_fingerprint = fingerprints$producer_fingerprint,
      code_fingerprint = fingerprints$code_fingerprint,
      package_fingerprint = fingerprints$package_fingerprint,
      execution_fingerprint = execution_fingerprint
    )
    list(
      task_index = i,
      case_id = as.character(row$case_id[[1L]]),
      row = row,
      seed_base = as.integer(seed_base),
      case_seed = case_seed,
      workers_requested = as.integer(workers_requested),
      run_started_at = mvt_iso_timestamp(run_started_at),
      planned_methods = planned_methods,
      require_gamcopula = isTRUE(require_gamcopula),
      execution_attestation = execution_attestation,
      execution_fingerprint = execution_fingerprint,
      producer_id = fingerprints$producer_id %||% mvt_producer_id,
      producer_version = fingerprints$producer_version %||% mvt_producer_version,
      configuration_fingerprint = configuration_fingerprint,
      producer_fingerprint = fingerprints$producer_fingerprint,
      code_fingerprint = fingerprints$code_fingerprint,
      package_fingerprint = fingerprints$package_fingerprint,
      contract_fingerprint = mvt_hash_object(case_contract)
    )
  })
}

mvt_case_checkpoint_dir <- function(run_dir) file.path(run_dir, "case_checkpoints")

mvt_run_lock_path <- function(run_dir) file.path(run_dir, ".mvt-run.lock")

mvt_acquire_run_lock <- function(run_dir) {
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  lock_path <- mvt_run_lock_path(run_dir)
  acquired <- dir.create(lock_path, recursive = FALSE, showWarnings = FALSE)
  if (!isTRUE(acquired)) {
    owner_path <- file.path(lock_path, "owner.rds")
    owner <- tryCatch(suppressWarnings(readRDS(owner_path)), error = function(e) NULL)
    detail <- if (is.list(owner)) paste0("pid=", owner$pid, " host=", owner$hostname, " started=", owner$started_at) else "owner unavailable"
    stop("Multivariate run directory is already leased: ", lock_path, " (", detail, ")", call. = FALSE)
  }
  owner <- list(
    pid = as.integer(Sys.getpid()), hostname = unname(Sys.info()[["nodename"]]),
    started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    nonce = mvt_hash_object(list(Sys.getpid(), Sys.time(), normalizePath(run_dir, winslash = "/", mustWork = TRUE)))
  )
  tryCatch(
    saveRDS(owner, file.path(lock_path, "owner.rds"), version = 3),
    error = function(e) {
      unlink(lock_path, recursive = TRUE, force = TRUE)
      stop("Could not initialize multivariate run-directory lease: ", conditionMessage(e), call. = FALSE)
    }
  )
  token <- new.env(parent = emptyenv())
  assign(owner$nonce, token, envir = mvt_active_leases)
  structure(
    list(
      path = normalizePath(lock_path, winslash = "/", mustWork = TRUE),
      run_dir = normalizePath(run_dir, winslash = "/", mustWork = TRUE),
      owner = owner,
      token = token
    ),
    class = "mvt_run_lease"
  )
}

mvt_assert_active_lease <- function(lease, target = NULL) {
  if (!inherits(lease, "mvt_run_lease") || !is.list(lease) || !is.environment(lease$token)) {
    stop("A valid active multivariate run lease is required.", call. = FALSE)
  }
  lock_path <- normalizePath(lease$path, winslash = "/", mustWork = FALSE)
  disk_owner <- tryCatch(readRDS(file.path(lock_path, "owner.rds")), error = function(e) NULL)
  registry_token <- if (exists(lease$owner$nonce, envir = mvt_active_leases, inherits = FALSE)) {
    get(lease$owner$nonce, envir = mvt_active_leases, inherits = FALSE)
  } else {
    NULL
  }
  valid <- dir.exists(lock_path) && is.list(disk_owner) &&
    identical(disk_owner$nonce, lease$owner$nonce) &&
    identical(as.integer(disk_owner$pid), as.integer(Sys.getpid())) &&
    identical(as.integer(lease$owner$pid), as.integer(Sys.getpid())) &&
    identical(registry_token, lease$token)
  if (!valid) stop("The multivariate run lease is inactive, forged, or owned by another process.", call. = FALSE)
  if (!is.null(target)) {
    target_norm <- normalizePath(target, winslash = "/", mustWork = FALSE)
    run_prefix <- paste0(normalizePath(lease$run_dir, winslash = "/", mustWork = TRUE), "/")
    if (!startsWith(target_norm, run_prefix)) stop("Lease does not own aggregate target: ", target, call. = FALSE)
  }
  invisible(lease$owner)
}

mvt_release_run_lock <- function(lease) {
  mvt_assert_active_lease(lease)
  lease_owner <- lease$owner
  lock_path <- normalizePath(lease$path, winslash = "/", mustWork = FALSE)
  if (!endsWith(lock_path, "/.mvt-run.lock")) stop("Refusing to release an unexpected lock path.", call. = FALSE)
  if (dir.exists(lock_path)) {
    disk_owner <- tryCatch(readRDS(file.path(lock_path, "owner.rds")), error = function(e) NULL)
    owned <- is.list(lease_owner) && is.list(disk_owner) &&
      identical(lease_owner$nonce, disk_owner$nonce) &&
      identical(as.integer(lease_owner$pid), as.integer(Sys.getpid()))
    if (!owned) stop("Refusing to release a run-directory lease owned by another process.", call. = FALSE)
    unlink(lock_path, recursive = TRUE, force = TRUE)
  }
  if (exists(lease_owner$nonce, envir = mvt_active_leases, inherits = FALSE)) {
    rm(list = lease_owner$nonce, envir = mvt_active_leases)
  }
  invisible(!dir.exists(lock_path))
}

mvt_case_checkpoint_path <- function(run_dir, case_id) {
  if (!grepl("^[A-Za-z0-9_.-]+$", case_id)) {
    stop("Unsafe case_id for checkpoint path: ", case_id, call. = FALSE)
  }
  file.path(mvt_case_checkpoint_dir(run_dir), paste0(case_id, ".rds"))
}

mvt_iso_timestamp <- function(time = Sys.time()) {
  format(as.POSIXct(time, tz = "UTC"), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
}

mvt_parse_iso_timestamp <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}[.][0-9]{3}Z$", value)) {
    return(as.POSIXct(NA, tz = "UTC"))
  }
  suppressWarnings(as.POSIXct(value, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"))
}

mvt_checkpoint_time_bounds_valid <- function(metadata) {
  fields <- c("run_started_at", "task_started_at", "task_completed_at", "created_at")
  parsed <- lapply(fields, function(field) mvt_parse_iso_timestamp(metadata[[field]]))
  if (any(vapply(parsed, is.na, logical(1L)))) return(FALSE)
  seconds <- vapply(parsed, as.numeric, numeric(1L))
  names(seconds) <- fields
  seconds[["run_started_at"]] <= seconds[["task_started_at"]] &&
    seconds[["task_started_at"]] <= seconds[["task_completed_at"]] &&
    seconds[["task_completed_at"]] == seconds[["created_at"]] &&
    seconds[["created_at"]] <= as.numeric(Sys.time()) + 300
}

mvt_checkpoint_worker_attestation_valid <- function(metadata, task) {
  attestation <- metadata$worker_attestation
  is.list(attestation) &&
    identical(as.integer(attestation$pid), as.integer(metadata$worker_pid)) &&
    identical(attestation$writer_role, metadata$writer_role) &&
    identical(attestation$execution_fingerprint, task$execution_fingerprint) &&
    identical(attestation$configuration_fingerprint, task$configuration_fingerprint) &&
    identical(attestation$contract_fingerprint, task$contract_fingerprint) &&
    identical(attestation$attestation_sha256, mvt_hash_object(attestation[setdiff(names(attestation), "attestation_sha256")]))
}

mvt_checkpoint_result_valid <- function(result, task) {
  if (!is.list(result) || !identical(sort(names(result)), sort(mvt_result_names()))) return(FALSE)
  if (!all(vapply(result, is.data.frame, logical(1L)))) return(FALSE)
  required <- result[c("fit_status", "benchmark_results", "runtime")]
  if (any(vapply(required, nrow, integer(1L)) < 1L)) return(FALSE)
  semantic_fields <- c(
    "case_id", "method", "attempted", "success", "converged", "retained",
    "stop_reason", "status_class", "elapsed_sec", "failure_reason_short", "error"
  )
  if (!all(vapply(required, function(x) all(semantic_fields %in% names(x)), logical(1L)))) return(FALSE)
  required_methods <- lapply(required, function(x) as.character(x$method))
  if (any(vapply(required_methods, anyDuplicated, integer(1L)) > 0L)) return(FALSE)
  if (!setequal(required_methods$fit_status, required_methods$benchmark_results) ||
      !setequal(required_methods$fit_status, required_methods$runtime)) return(FALSE)
  if (!setequal(required_methods$fit_status, task$planned_methods)) return(FALSE)
  semantic_ok <- vapply(required, function(x) {
    attempted <- x$attempted %in% TRUE
    success <- x$success %in% TRUE
    converged <- x$converged %in% TRUE
    retained <- x$retained %in% TRUE
    elapsed <- suppressWarnings(as.numeric(x$elapsed_sec))
    reasons <- as.character(x$failure_reason_short)
    stops <- as.character(x$stop_reason)
    statuses <- as.character(x$status_class)
    errors <- as.character(x$error)
    all(attempted) && !anyNA(x$success) && !anyNA(x$converged) && !anyNA(x$retained) &&
      all(!retained | (success & converged)) && all(success == retained) &&
      all(is.finite(elapsed) & elapsed >= 0) && all(nzchar(trimws(stops))) &&
      all(nzchar(trimws(statuses))) &&
      all(!retained | tolower(statuses) %in% c("ok", "warning", "success")) &&
      all(retained | tolower(statuses) %in% c("error", "timeout", "warning")) &&
      all(!retained | stops == "converged") &&
      all(retained | (!is.na(errors) & nzchar(trimws(errors)))) &&
      all(retained | (!is.na(reasons) & nzchar(trimws(reasons)) & reasons != "none"))
  }, logical(1L))
  if (!all(semantic_ok)) return(FALSE)
  semantic_compare <- c(
    "attempted", "success", "converged", "retained", "stop_reason",
    "status_class", "failure_reason_short", "error"
  )
  normalized_semantics <- lapply(required, function(x) {
    x <- x[order(match(as.character(x$method), task$planned_methods)), c("method", semantic_compare), drop = FALSE]
    rownames(x) <- NULL
    x
  })
  if (!identical(normalized_semantics$fit_status, normalized_semantics$benchmark_results) ||
      !identical(normalized_semantics$fit_status, normalized_semantics$runtime)) return(FALSE)
  benchmark <- result$benchmark_results
  nonretained <- !(benchmark$retained %in% TRUE)
  accuracy <- intersect(c("mae", "rmse", grep("^benchmark_", names(benchmark), value = TRUE), "logLik", "AIC", "BIC"), names(benchmark))
  if (any(nonretained) && length(accuracy) && any(vapply(benchmark[nonretained, accuracy, drop = FALSE], function(x) any(is.finite(suppressWarnings(as.numeric(x)))), logical(1L)))) return(FALSE)
  populated <- result[vapply(result, nrow, integer(1L)) > 0L]
  case_ok <- all(vapply(populated, function(x) {
    "case_id" %in% names(x) && all(as.character(x$case_id) == task$case_id)
  }, logical(1L)))
  if (!case_ok) return(FALSE)
  auxiliary_schemas <- list(
    coefficient_results = c(
      "case_id", "scenario", "family", "rep", "method", "parameter", "term",
      "estimate", "std_error", "truth", "bias", "conf.low", "conf.high"
    ),
    dependence_recovery = c(
      "case_id", "scenario", "family", "rep", "method", "dependence_scope",
      "dependence_n", "theta_mae", "theta_rmse", "tau_mae", "tau_rmse"
    ),
    variogram_scores = c(
      "case_id", "scenario", "family", "rep", "method",
      "variogram_score_p05", "variogram_score_p2", "error"
    )
  )
  for (name in names(auxiliary_schemas)) {
    table <- result[[name]]
    if (!nrow(table)) next
    if (!all(auxiliary_schemas[[name]] %in% names(table))) return(FALSE)
    if (!all(as.character(table$method) %in% task$planned_methods)) return(FALSE)
    retained_methods <- as.character(benchmark$method[benchmark$retained %in% TRUE])
    if (name != "variogram_scores" && !all(as.character(table$method) %in% retained_methods)) return(FALSE)
  }
  TRUE
}

mvt_case_checkpoint_valid <- function(checkpoint, task) {
  if (!is.list(checkpoint) || !is.list(checkpoint$metadata)) return(FALSE)
  metadata <- checkpoint$metadata
  expected <- c(
    "checkpoint_schema_version", "identity_algorithm", "identity_version", "producer_id", "producer_version", "case_id", "task_index", "case_seed",
    "workers_requested", "planned_methods",
    "configuration_fingerprint", "producer_fingerprint", "code_fingerprint",
    "package_fingerprint", "execution_fingerprint", "contract_fingerprint",
    "run_started_at", "task_started_at", "task_completed_at", "created_at",
    "worker_pid", "writer_role", "worker_attestation"
  )
  if (!all(expected %in% names(metadata))) return(FALSE)
  identical(metadata$checkpoint_schema_version, mvt_checkpoint_schema_version) &&
    identical(metadata$identity_algorithm, mvt_identity_algorithm) &&
    identical(metadata$identity_version, mvt_identity_version) &&
    identical(metadata$producer_id, task$producer_id) &&
    identical(metadata$producer_version, task$producer_version) &&
    identical(metadata$case_id, task$case_id) &&
    length(metadata$task_index) == 1L && !is.na(metadata$task_index) &&
    metadata$task_index == floor(metadata$task_index) && metadata$task_index > 0L &&
    identical(as.integer(metadata$task_index), as.integer(task$task_index)) &&
    identical(metadata$case_seed, as.integer(task$case_seed)) &&
    length(metadata$workers_requested) == 1L && is.finite(metadata$workers_requested) &&
    metadata$workers_requested >= 1L && metadata$workers_requested == floor(metadata$workers_requested) &&
    identical(metadata$planned_methods, task$planned_methods) &&
    identical(metadata$configuration_fingerprint, task$configuration_fingerprint) &&
    identical(metadata$producer_fingerprint, task$producer_fingerprint) &&
    identical(metadata$code_fingerprint, task$code_fingerprint) &&
    identical(metadata$package_fingerprint, task$package_fingerprint) &&
    identical(metadata$execution_fingerprint, task$execution_fingerprint) &&
    identical(metadata$contract_fingerprint, task$contract_fingerprint) &&
    length(metadata$worker_pid) == 1L && is.finite(metadata$worker_pid) &&
    metadata$worker_pid == floor(metadata$worker_pid) && metadata$worker_pid > 0L &&
    length(metadata$writer_role) == 1L && metadata$writer_role %in% c("psock_worker", "serial_parent") &&
    mvt_checkpoint_time_bounds_valid(metadata) &&
    mvt_checkpoint_worker_attestation_valid(metadata, task) &&
    mvt_checkpoint_result_valid(checkpoint$result, task)
}

mvt_case_checkpoint_issues <- function(checkpoint, task) {
  issues <- character()
  if (!is.list(checkpoint)) return("checkpoint_not_list")
  if (!is.list(checkpoint$metadata)) return("metadata_missing")
  metadata <- checkpoint$metadata
  comparisons <- c(
    schema = identical(metadata$checkpoint_schema_version, mvt_checkpoint_schema_version),
    identity_algorithm = identical(metadata$identity_algorithm, mvt_identity_algorithm),
    identity_version = identical(metadata$identity_version, mvt_identity_version),
    producer_id = identical(metadata$producer_id, task$producer_id),
    producer_version = identical(metadata$producer_version, task$producer_version),
    case_id = identical(metadata$case_id, task$case_id),
    task_index = length(metadata$task_index) == 1L && !is.na(metadata$task_index) &&
      metadata$task_index == floor(metadata$task_index) && metadata$task_index > 0L &&
      identical(as.integer(metadata$task_index), as.integer(task$task_index)),
    case_seed = identical(metadata$case_seed, as.integer(task$case_seed)),
    workers_requested = length(metadata$workers_requested) == 1L &&
      is.finite(metadata$workers_requested) && metadata$workers_requested >= 1L &&
      metadata$workers_requested == floor(metadata$workers_requested),
    planned_methods = identical(metadata$planned_methods, task$planned_methods),
    configuration = identical(metadata$configuration_fingerprint, task$configuration_fingerprint),
    producer = identical(metadata$producer_fingerprint, task$producer_fingerprint),
    code = identical(metadata$code_fingerprint, task$code_fingerprint),
    packages = identical(metadata$package_fingerprint, task$package_fingerprint),
    execution = identical(metadata$execution_fingerprint, task$execution_fingerprint),
    contract = identical(metadata$contract_fingerprint, task$contract_fingerprint),
    time_bounds = mvt_checkpoint_time_bounds_valid(metadata),
    worker_pid = length(metadata$worker_pid) == 1L && is.finite(metadata$worker_pid) &&
      metadata$worker_pid == floor(metadata$worker_pid) && metadata$worker_pid > 0L,
    worker_ownership = mvt_checkpoint_worker_attestation_valid(metadata, task)
  )
  issues <- c(issues, names(comparisons)[!comparisons])
  if (!mvt_checkpoint_result_valid(checkpoint$result, task)) issues <- c(issues, "semantic_result_schema")
  unique(issues)
}

mvt_read_case_checkpoint <- function(path, task) {
  if (!file.exists(path)) return(NULL)
  tryCatch({
    checkpoint <- readRDS(path)
    issues <- mvt_case_checkpoint_issues(checkpoint, task)
    if (!length(issues)) checkpoint else structure(
      list(), class = "mvt_rejected_checkpoint",
      rejection_reason = paste(issues, collapse = "|")
    )
  }, error = function(e) structure(
    list(), class = "mvt_rejected_checkpoint",
    rejection_reason = paste0("read_error:", conditionMessage(e))
  ))
}

mvt_archive_stale_checkpoint <- function(path, run_dir, reason = "invalid") {
  if (!file.exists(path)) return(invisible(FALSE))
  stale_dir <- file.path(mvt_case_checkpoint_dir(run_dir), "stale")
  dir.create(stale_dir, recursive = TRUE, showWarnings = FALSE)
  destination <- tempfile(
    paste0(basename(path), ".", gsub("[^A-Za-z0-9_-]+", "_", reason), ".", format(Sys.time(), "%Y%m%d%H%M%S"), "."),
    tmpdir = stale_dir
  )
  if (!file.rename(path, destination)) {
    stop("Could not archive stale case checkpoint: ", path, call. = FALSE)
  }
  invisible(destination)
}

mvt_write_case_checkpoint_atomic <- function(
    result, task, path,
    writer_role = if (task$workers_requested > 1L) "psock_worker" else "serial_parent",
    task_started_at = Sys.time()) {
  writer_role <- match.arg(writer_role, c("psock_worker", "serial_parent"))
  task_started_at <- as.POSIXct(task_started_at, tz = "UTC")
  completed_at <- Sys.time()
  worker_attestation <- list(
    pid = as.integer(Sys.getpid()), writer_role = writer_role,
    execution_fingerprint = task$execution_fingerprint,
    configuration_fingerprint = task$configuration_fingerprint,
    contract_fingerprint = task$contract_fingerprint
  )
  worker_attestation$attestation_sha256 <- mvt_hash_object(worker_attestation)
  checkpoint <- list(
    metadata = list(
      checkpoint_schema_version = mvt_checkpoint_schema_version,
      identity_algorithm = mvt_identity_algorithm,
      identity_version = mvt_identity_version,
      producer_id = task$producer_id,
      producer_version = task$producer_version,
      case_id = task$case_id,
      task_index = as.integer(task$task_index),
      case_seed = as.integer(task$case_seed),
      workers_requested = as.integer(task$workers_requested),
      planned_methods = task$planned_methods,
      configuration_fingerprint = task$configuration_fingerprint,
      producer_fingerprint = task$producer_fingerprint,
      code_fingerprint = task$code_fingerprint,
      package_fingerprint = task$package_fingerprint,
      execution_fingerprint = task$execution_fingerprint,
      contract_fingerprint = task$contract_fingerprint,
      run_started_at = task$run_started_at,
      task_started_at = mvt_iso_timestamp(task_started_at),
      task_completed_at = mvt_iso_timestamp(completed_at),
      created_at = mvt_iso_timestamp(completed_at),
      worker_pid = as.integer(Sys.getpid()),
      writer_role = writer_role,
      worker_attestation = worker_attestation
    ),
    result = result
  )
  if (!mvt_case_checkpoint_valid(checkpoint, task)) {
    stop("Refusing incomplete multivariate case checkpoint: ", task$case_id, call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(checkpoint, temporary, version = 3)
  if (!mvt_case_checkpoint_valid(readRDS(temporary), task)) {
    stop("Temporary multivariate case checkpoint failed validation: ", task$case_id, call. = FALSE)
  }
  if (file.exists(path)) stop("Refusing concurrent checkpoint overwrite: ", path, call. = FALSE)
  if (!file.rename(temporary, path)) stop("Could not atomically install case checkpoint: ", path, call. = FALSE)
  invisible(path)
}

mvt_run_checkpoint_task <- function(task, run_dir, run_case_fun, writer_role = "serial_parent") {
  task_started_at <- Sys.time()
  old_options <- options(gamlss.longitudinal.mvt.execution_attestation = task$execution_attestation)
  on.exit(options(old_options), add = TRUE)
  result <- run_case_fun(
    task$row,
    seed_base = task$seed_base,
    require_gamcopula = task$require_gamcopula,
    case_seed = task$case_seed
  )
  path <- mvt_case_checkpoint_path(run_dir, task$case_id)
  mvt_write_case_checkpoint_atomic(
    result, task, path, writer_role = writer_role,
    task_started_at = task_started_at
  )
  list(task_index = task$task_index, case_id = task$case_id, path = path, worker_pid = Sys.getpid())
}

mvt_order_result_rows <- function(x, grid, methods) {
  if (!is.data.frame(x) || nrow(x) < 2L) return(x)
  scenario_order <- unique(as.character(grid$scenario))
  family_order <- unique(as.character(grid$family_name))
  keys <- list(
    if ("scenario" %in% names(x)) match(as.character(x$scenario), scenario_order) else rep(1L, nrow(x)),
    if ("rep" %in% names(x)) suppressWarnings(as.integer(x$rep)) else rep(1L, nrow(x)),
    if ("method" %in% names(x)) match(as.character(x$method), methods) else rep(1L, nrow(x)),
    if ("family" %in% names(x)) match(as.character(x$family), family_order) else rep(1L, nrow(x)),
    if ("case_id" %in% names(x)) match(as.character(x$case_id), grid$case_id) else rep(1L, nrow(x))
  )
  for (column in intersect(c("term", "parameter", "dependence_scope"), names(x))) {
    keys[[length(keys) + 1L]] <- as.character(x[[column]])
  }
  keys <- lapply(keys, function(key) {
    if (is.numeric(key) || is.integer(key)) key[is.na(key)] <- Inf else key[is.na(key)] <- "\U0010ffff"
    key
  })
  out <- x[do.call(order, c(keys, list(na.last = TRUE))), , drop = FALSE]
  rownames(out) <- NULL
  out
}

mvt_collect_case_checkpoints <- function(tasks, run_dir, grid, methods) {
  checkpoints <- lapply(tasks, function(task) {
    checkpoint <- mvt_read_case_checkpoint(mvt_case_checkpoint_path(run_dir, task$case_id), task)
    if (is.null(checkpoint) || inherits(checkpoint, "mvt_rejected_checkpoint")) {
      stop("Missing or invalid completed checkpoint: ", task$case_id, call. = FALSE)
    }
    checkpoint$result
  })
  stats::setNames(lapply(mvt_result_names(), function(name) {
    mvt_order_result_rows(
      mvt_bind_rows_fill(lapply(checkpoints, `[[`, name)),
      grid = grid,
      methods = methods
    )
  }), mvt_result_names())
}

mvt_checkpoint_manifest <- function(tasks, run_dir) {
  rows <- lapply(tasks, function(task) {
    path <- mvt_case_checkpoint_path(run_dir, task$case_id)
    checkpoint <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!mvt_case_checkpoint_valid(checkpoint, task)) {
      stop("Checkpoint set failed final semantic validation: ", task$case_id, call. = FALSE)
    }
    data.frame(
      task_index = task$task_index, case_id = task$case_id,
      case_seed = task$case_seed, contract_fingerprint = task$contract_fingerprint,
      checkpoint_sha256 = mvt_sha256_file(path),
      result_sha256 = mvt_hash_object(checkpoint$result),
      writer_pid = checkpoint$metadata$worker_pid,
      writer_role = checkpoint$metadata$writer_role,
      worker_attestation_sha256 = checkpoint$metadata$worker_attestation$attestation_sha256,
      created_at = checkpoint$metadata$created_at,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

mvt_csv_serialized_sha256 <- function(x) {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  mvt_sha256_file(path)
}

mvt_validate_aggregate_file_commit <- function(
    path, expected_nonce = NULL, expected_pid = NULL,
    expected_commit_sha256 = NULL, expected_ownership_sha256 = NULL) {
  commit_path <- paste0(path, ".commit.rds")
  ownership_path <- paste0(path, ".ownership.rds")
  commit <- tryCatch(suppressWarnings(readRDS(commit_path)), error = function(e) NULL)
  ownership <- tryCatch(suppressWarnings(readRDS(ownership_path)), error = function(e) NULL)
  required <- c(
    "schema_version", "file", "owner_role", "writer_pid", "lease_nonce",
    "hostname", "written_at", "sha256", "rows", "columns", "schema_sha256", "bytes"
  )
  if (!file.exists(path) || !is.list(commit) || !is.list(ownership) ||
      !all(required %in% names(commit)) || !identical(commit, ownership)) {
    stop("Aggregate commit/ownership record is missing, incomplete, or incoherent: ", basename(path), call. = FALSE)
  }
  actual <- mvt_read_optional_csv(path)
  valid <- identical(commit$schema_version, 1L) &&
    identical(commit$file, basename(path)) && identical(commit$owner_role, "lease_parent") &&
    is.finite(commit$writer_pid) && commit$writer_pid > 0L &&
    is.character(commit$lease_nonce) && length(commit$lease_nonce) == 1L && nzchar(commit$lease_nonce) &&
    identical(commit$sha256, mvt_sha256_file(path)) &&
    identical(as.integer(commit$rows), as.integer(nrow(actual))) &&
    identical(as.character(commit$columns), names(actual)) &&
    identical(commit$schema_sha256, mvt_hash_object(names(actual))) &&
    identical(as.numeric(commit$bytes), as.numeric(file.info(path)$size))
  if (!is.null(expected_nonce)) valid <- valid && identical(commit$lease_nonce, expected_nonce)
  if (!is.null(expected_pid)) valid <- valid && identical(as.integer(commit$writer_pid), as.integer(expected_pid))
  if (!is.null(expected_commit_sha256)) valid <- valid && identical(mvt_sha256_file(commit_path), expected_commit_sha256)
  if (!is.null(expected_ownership_sha256)) valid <- valid && identical(mvt_sha256_file(ownership_path), expected_ownership_sha256)
  if (!valid) stop("Aggregate commit/ownership fields do not reconcile: ", basename(path), call. = FALSE)
  list(
    commit = commit,
    commit_sha256 = mvt_sha256_file(commit_path),
    ownership_sha256 = mvt_sha256_file(ownership_path)
  )
}

mvt_validate_lines_file_commit <- function(path) {
  commit_path <- paste0(path, ".commit.rds")
  ownership_path <- paste0(path, ".ownership.rds")
  commit <- tryCatch(readRDS(commit_path), error = function(e) NULL)
  ownership <- tryCatch(readRDS(ownership_path), error = function(e) NULL)
  lines <- if (file.exists(path)) readLines(path, warn = FALSE) else character()
  valid <- file.exists(path) && is.list(commit) && identical(commit, ownership) &&
    identical(commit$schema_version, 1L) && identical(commit$file, basename(path)) &&
    identical(commit$owner_role, "lease_parent") &&
    is.finite(commit$writer_pid) && commit$writer_pid > 0L &&
    is.character(commit$lease_nonce) && length(commit$lease_nonce) == 1L && nzchar(commit$lease_nonce) &&
    identical(commit$sha256, mvt_sha256_file(path)) &&
    identical(as.integer(commit$rows), as.integer(length(lines))) &&
    identical(commit$columns, character()) &&
    identical(commit$schema_sha256, mvt_hash_object(character())) &&
    identical(as.numeric(commit$bytes), as.numeric(file.info(path)$size))
  if (!valid) stop("Line artifact commit/ownership fields do not reconcile: ", basename(path), call. = FALSE)
  invisible(commit)
}

mvt_validate_aggregate_table_semantics <- function(grid, tables, methods) {
  required <- c("fit_status", "benchmark_results", "runtime")
  if (!all(required %in% names(tables))) stop("Committed aggregate snapshot is missing attempt tables.", call. = FALSE)
  expected_keys <- do.call(paste, c(
    merge(grid["case_id"], data.frame(method = methods, stringsAsFactors = FALSE), by = NULL, sort = FALSE),
    sep = "\r"
  ))
  for (name in required) {
    table <- tables[[name]]
    semantic <- c(
      "case_id", "method", "attempted", "success", "converged", "retained",
      "stop_reason", "status_class", "elapsed_sec", "failure_reason_short", "error"
    )
    if (!all(semantic %in% names(table))) stop("Aggregate semantic schema failure: ", name, call. = FALSE)
    keys <- paste(table$case_id, table$method, sep = "\r")
    if (nrow(table) != length(expected_keys) || anyDuplicated(keys) || !setequal(keys, expected_keys)) {
      stop("Aggregate case-method cross-product failure: ", name, call. = FALSE)
    }
    attempted <- table$attempted %in% TRUE
    success <- table$success %in% TRUE
    converged <- table$converged %in% TRUE
    retained <- table$retained %in% TRUE
    elapsed <- suppressWarnings(as.numeric(table$elapsed_sec))
    stop_reason <- as.character(table$stop_reason)
    status_class <- as.character(table$status_class)
    failure <- as.character(table$failure_reason_short)
    error <- as.character(table$error)
    if (!all(attempted) || any(success != retained) || any(retained & !converged) ||
        any(!is.finite(elapsed) | elapsed < 0) || any(!nzchar(trimws(stop_reason)) | is.na(stop_reason)) ||
        any(retained & stop_reason != "converged") || any(retained & !tolower(status_class) %in% c("ok", "warning", "success")) ||
        any(!retained & (is.na(error) | !nzchar(trimws(error)))) ||
        any(!retained & (is.na(failure) | !nzchar(trimws(failure)) | failure == "none"))) {
      stop("Aggregate attempt semantics are contradictory: ", name, call. = FALSE)
    }
  }
  semantic_compare <- c(
    "case_id", "method", "attempted", "success", "converged", "retained",
    "stop_reason", "status_class", "failure_reason_short", "error"
  )
  normalized <- lapply(tables[required], function(x) {
    x <- x[order(match(x$case_id, grid$case_id), match(x$method, methods)), semantic_compare, drop = FALSE]
    rownames(x) <- NULL
    x
  })
  if (!identical(normalized$fit_status, normalized$benchmark_results) ||
      !identical(normalized$fit_status, normalized$runtime)) {
    stop("Aggregate attempt tables disagree on semantic status fields.", call. = FALSE)
  }
  benchmark <- tables$benchmark_results
  nonretained <- !(benchmark$retained %in% TRUE)
  accuracy <- intersect(c(
    "mae", "rmse", "logLik", "AIC", "BIC", grep("^benchmark_", names(benchmark), value = TRUE)
  ), names(benchmark))
  if (any(nonretained) && length(accuracy) && any(vapply(
    benchmark[nonretained, accuracy, drop = FALSE],
    function(x) any(is.finite(suppressWarnings(as.numeric(x)))), logical(1L)
  ))) stop("Nonretained aggregate rows contribute accuracy evidence.", call. = FALSE)
  for (name in c("coefficient_results", "dependence_recovery", "variogram_scores")) {
    table <- tables[[name]]
    if (is.null(table) || !nrow(table)) next
    if (!all(c("case_id", "scenario", "family", "rep", "method") %in% names(table))) {
      stop("Aggregate truth-table binding schema failure: ", name, call. = FALSE)
    }
    if (!all(paste(table$case_id, table$method, sep = "\r") %in% expected_keys)) {
      stop("Aggregate truth table contains an unregistered case-method key: ", name, call. = FALSE)
    }
    grid_binding <- grid
    if (!"family" %in% names(grid_binding) && "family_name" %in% names(grid_binding)) grid_binding$family <- grid_binding$family_name
    binding_fields <- intersect(c("scenario", "generator", "dependence", "correlation_level", "n_time", "n_subject", "total_rows", "family", "rep"), names(grid_binding))
    lookup <- grid_binding[match(table$case_id, grid_binding$case_id), binding_fields, drop = FALSE]
    comparable <- intersect(binding_fields, names(table))
    if (length(comparable) && any(vapply(comparable, function(field) {
      !identical(as.character(table[[field]]), as.character(lookup[[field]]))
    }, logical(1L)))) stop("Aggregate truth table metadata contradicts its registered case: ", name, call. = FALSE)
  }
  if (nrow(tables$coefficient_results) &&
      !all(c("parameter", "term", "estimate", "std_error", "truth", "bias", "conf.low", "conf.high") %in%
        names(tables$coefficient_results))) {
    stop("Coefficient truth table is incomplete.", call. = FALSE)
  }
  invisible(TRUE)
}

mvt_write_aggregate_snapshot <- function(
    tasks, run_dir, grid, methods, rows, checkpoint_manifest_before,
    configuration, production_run, lease) {
  checkpoint_manifest_after <- mvt_checkpoint_manifest(tasks, run_dir)
  if (!identical(checkpoint_manifest_after, checkpoint_manifest_before)) {
    stop("Checkpoint set changed during parent aggregate commit.", call. = FALSE)
  }
  mvt_validate_aggregate_table_semantics(grid, rows, methods)
  artifact_names <- mvt_phase2_public_output_allowlist()$attempt_artifacts
  artifact_paths <- file.path(run_dir, artifact_names)
  missing <- !file.exists(artifact_paths)
  if (any(missing)) stop("Cannot commit incomplete aggregate snapshot: ", paste(artifact_names[missing], collapse = ", "), call. = FALSE)
  artifact_rows <- lapply(seq_along(artifact_paths), function(i) {
    path <- artifact_paths[[i]]
    validated <- mvt_validate_aggregate_file_commit(
      path, expected_nonce = lease$owner$nonce, expected_pid = lease$owner$pid
    )
    commit <- validated$commit
    data.frame(
      file = artifact_names[[i]], sha256 = commit$sha256,
      rows = commit$rows, schema_sha256 = commit$schema_sha256,
      bytes = commit$bytes, writer_pid = commit$writer_pid,
      owner_role = commit$owner_role, lease_nonce = commit$lease_nonce,
      commit_sha256 = validated$commit_sha256,
      ownership_sha256 = validated$ownership_sha256,
      stringsAsFactors = FALSE
    )
  })
  artifacts <- do.call(rbind, artifact_rows)
  expected_result_hashes <- vapply(names(rows), function(name) {
    mvt_csv_serialized_sha256(rows[[name]])
  }, character(1L))
  actual_result_hashes <- stats::setNames(
    artifacts$sha256[match(paste0(names(rows), "_by_rep.csv"), artifacts$file)], names(rows)
  )
  if (!identical(unname(actual_result_hashes), unname(expected_result_hashes))) {
    stop("Aggregate CSVs do not exactly serialize the validated checkpoint payloads.", call. = FALSE)
  }
  snapshot <- list(
    schema_version = 1L, status = "preparing",
    committed_at = NA_character_,
    lease_nonce = lease$owner$nonce, owner_pid = lease$owner$pid,
    producer_id = mvt_producer_id, producer_version = mvt_producer_version,
    production_run = isTRUE(production_run), configuration = configuration,
    configuration_fingerprint = tasks[[1L]]$configuration_fingerprint,
    execution_attestation = tasks[[1L]]$execution_attestation,
    grid = grid, methods = methods, tasks = tasks,
    checkpoint_manifest = checkpoint_manifest_after,
    artifacts = artifacts, expected_result_hashes = expected_result_hashes
  )
  snapshot_path <- file.path(run_dir, "aggregate_snapshot.rds")
  mvt_save_rds_atomic(snapshot, snapshot_path, lease)
  manifest_precommit <- mvt_checkpoint_manifest(tasks, run_dir)
  if (!identical(manifest_precommit, checkpoint_manifest_after) || any(vapply(
    seq_len(nrow(artifacts)), function(i) {
      !identical(mvt_sha256_file(file.path(run_dir, artifacts$file[[i]])), artifacts$sha256[[i]])
    }, logical(1L)
  ))) stop("Checkpoint or aggregate set changed immediately before snapshot commit.", call. = FALSE)
  snapshot$status <- "committed_immutable"
  snapshot$committed_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  mvt_save_rds_atomic(snapshot, snapshot_path, lease)
  manifest_postcommit <- mvt_checkpoint_manifest(tasks, run_dir)
  if (!identical(manifest_postcommit, checkpoint_manifest_after) || any(vapply(
    seq_len(nrow(artifacts)), function(i) {
      !identical(mvt_sha256_file(file.path(run_dir, artifacts$file[[i]])), artifacts$sha256[[i]])
    }, logical(1L)
  ))) stop("Checkpoint or aggregate set changed immediately after snapshot commit.", call. = FALSE)
  snapshot_commit <- list(
    schema_version = 1L, status = "committed_immutable",
    file = basename(snapshot_path), owner_role = "lease_parent",
    lease_nonce = lease$owner$nonce, writer_pid = lease$owner$pid,
    sha256 = mvt_sha256_file(snapshot_path), bytes = as.numeric(file.info(snapshot_path)$size),
    snapshot_schema_sha256 = mvt_hash_object(names(snapshot)),
    checkpoint_manifest_sha256 = mvt_hash_object(snapshot$checkpoint_manifest),
    artifacts_sha256 = mvt_hash_object(snapshot$artifacts),
    committed_at = snapshot$committed_at
  )
  mvt_save_rds_atomic(snapshot_commit, paste0(snapshot_path, ".commit.rds"), lease)
  mvt_save_rds_atomic(snapshot_commit, paste0(snapshot_path, ".ownership.rds"), lease)
  invisible(snapshot)
}

mvt_validate_committed_snapshot <- function(run_dir, require_production = FALSE) {
  snapshot_path <- file.path(run_dir, "aggregate_snapshot.rds")
  snapshot <- tryCatch(suppressWarnings(readRDS(snapshot_path)), error = function(e) NULL)
  commit <- tryCatch(suppressWarnings(readRDS(paste0(snapshot_path, ".commit.rds"))), error = function(e) NULL)
  ownership <- tryCatch(suppressWarnings(readRDS(paste0(snapshot_path, ".ownership.rds"))), error = function(e) NULL)
  root_fields <- c(
    "schema_version", "status", "file", "owner_role", "lease_nonce", "writer_pid",
    "sha256", "bytes", "snapshot_schema_sha256", "checkpoint_manifest_sha256",
    "artifacts_sha256", "committed_at"
  )
  if (!is.list(snapshot) || !identical(snapshot$status, "committed_immutable") || !is.list(commit) ||
      !is.list(ownership) || !identical(commit, ownership) || !all(root_fields %in% names(commit)) ||
      !identical(commit$status, "committed_immutable") ||
      !identical(commit$schema_version, 1L) || !identical(commit$file, basename(snapshot_path)) ||
      !identical(commit$owner_role, "lease_parent") ||
      !identical(commit$sha256, mvt_sha256_file(snapshot_path)) ||
      !identical(as.numeric(commit$bytes), as.numeric(file.info(snapshot_path)$size)) ||
      !identical(commit$snapshot_schema_sha256, mvt_hash_object(names(snapshot))) ||
      !identical(commit$checkpoint_manifest_sha256, mvt_hash_object(snapshot$checkpoint_manifest)) ||
      !identical(commit$artifacts_sha256, mvt_hash_object(snapshot$artifacts)) ||
      !identical(commit$lease_nonce, snapshot$lease_nonce) ||
      !identical(as.integer(commit$writer_pid), as.integer(snapshot$owner_pid))) {
    stop("Run directory does not contain a valid immutable aggregate snapshot.", call. = FALSE)
  }
  if (isTRUE(require_production) && !isTRUE(snapshot$production_run)) {
    stop("Aggregate snapshot is not registered as an exact production run.", call. = FALSE)
  }
  if (isTRUE(require_production)) {
    current_fingerprints <- mvt_checkpoint_fingerprints()
    current_execution <- mvt_execution_attestation_contract(
      current_fingerprints,
      configuration_fingerprint = snapshot$configuration_fingerprint
    )
    stored_task <- snapshot$tasks[[1L]]
    current_fields <- c("producer_fingerprint", "code_fingerprint", "package_fingerprint")
    stale <- any(vapply(current_fields, function(field) {
      !identical(stored_task[[field]], current_fingerprints[[field]])
    }, logical(1L))) ||
      !identical(snapshot$configuration_fingerprint, mvt_hash_object(snapshot$configuration)) ||
      !identical(stored_task$execution_fingerprint, mvt_hash_object(current_execution))
    if (stale) {
      stop("Production aggregate snapshot is stale relative to the current producer/code/package/runtime/configuration identity.", call. = FALSE)
    }
    current_attestation <- mvt_verify_execution_attestation(current_execution, load = TRUE)
    if (!mvt_execution_attestation_matches(current_attestation, current_execution)) {
      stop("Production aggregate snapshot current execution attestation failed.", call. = FALSE)
    }
  }
  for (i in seq_len(nrow(snapshot$artifacts))) {
    artifact <- snapshot$artifacts[i, , drop = FALSE]
    path <- file.path(run_dir, artifact$file)
    validated <- tryCatch(mvt_validate_aggregate_file_commit(
      path,
      expected_nonce = snapshot$lease_nonce,
      expected_pid = snapshot$owner_pid,
      expected_commit_sha256 = artifact$commit_sha256,
      expected_ownership_sha256 = artifact$ownership_sha256
    ), error = function(e) e)
    file_commit <- if (inherits(validated, "error")) NULL else validated$commit
    if (is.null(file_commit) || !identical(file_commit$sha256, artifact$sha256) ||
        !identical(as.integer(file_commit$rows), as.integer(artifact$rows)) ||
        !identical(file_commit$schema_sha256, artifact$schema_sha256) ||
        !identical(as.numeric(file_commit$bytes), as.numeric(artifact$bytes)) ||
        !identical(file_commit$owner_role, artifact$owner_role) ||
        !identical(as.integer(file_commit$writer_pid), as.integer(artifact$writer_pid)) ||
        !identical(file_commit$lease_nonce, artifact$lease_nonce)) {
      stop("Committed aggregate artifact is missing or mutable: ", artifact$file, call. = FALSE)
    }
  }
  current_manifest <- mvt_checkpoint_manifest(snapshot$tasks, run_dir)
  if (!identical(current_manifest, snapshot$checkpoint_manifest)) {
    stop("Checkpoint hashes/schema no longer reconcile with the committed snapshot.", call. = FALSE)
  }
  rows <- mvt_collect_case_checkpoints(snapshot$tasks, run_dir, snapshot$grid, snapshot$methods)
  for (name in names(rows)) {
    path <- file.path(run_dir, paste0(name, "_by_rep.csv"))
    if (!identical(mvt_csv_serialized_sha256(rows[[name]]), mvt_sha256_file(path))) {
      stop("Aggregate payload no longer reconciles with checkpoints: ", name, call. = FALSE)
    }
  }
  loaded <- stats::setNames(lapply(names(rows), function(name) {
    mvt_read_optional_csv(file.path(run_dir, paste0(name, "_by_rep.csv")))
  }), names(rows))
  mvt_validate_aggregate_table_semantics(snapshot$grid, loaded, snapshot$methods)
  invisible(c(snapshot, list(rows = loaded)))
}

mvt_phase2_snapshot_attestation_path <- function() {
  Sys.getenv("GAMLSS_LONGITUDINAL_MVT_ATTESTATION", unset = "")
}

mvt_phase2_snapshot_signature_path <- function() {
  Sys.getenv("GAMLSS_LONGITUDINAL_MVT_ATTESTATION_SIGNATURE", unset = "")
}

mvt_phase2_pinned_public_key <- function() {
  as.raw(c(66, 91, 71, 233, 21, 247, 172, 45, 215, 202, 170, 0, 64,
    43, 83, 206, 23, 50, 48, 154, 25, 217, 178, 37, 252, 59, 158, 195,
    237, 0, 31, 216))
}

mvt_phase2_parse_approval_time <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$", x)) {
    stop("Phase 2 snapshot approval timestamp is not RFC3339 UTC.", call. = FALSE)
  }
  parsed <- as.POSIXct(strptime(x, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  if (is.na(parsed) || !identical(format(parsed, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), x)) {
    stop("Phase 2 snapshot approval timestamp is not a real RFC3339 UTC instant.", call. = FALSE)
  }
  parsed
}

mvt_read_signed_snapshot_approval <- function(attestation_path, signature_path) {
  if (!requireNamespace("sodium", quietly = TRUE)) stop("sodium is required for Phase 2 snapshot approval verification.", call. = FALSE)
  if (!nzchar(attestation_path) || !nzchar(signature_path)) {
    stop("Phase 2 snapshot lacks a detached production approval signature.", call. = FALSE)
  }
  root <- paste0(tolower(normalizePath(mvt_repo_root, winslash = "/", mustWork = TRUE)), "/")
  approval_paths <- vapply(c(attestation_path, signature_path), normalizePath, character(1L),
    winslash = "/", mustWork = TRUE)
  if (any(startsWith(tolower(paste0(approval_paths, "/")), root))) {
    stop("Phase 2 snapshot approval attestation/signature must be external to the checkout.", call. = FALSE)
  }
  message_raw <- readBin(approval_paths[[1L]], "raw", n = file.info(approval_paths[[1L]])$size)
  signature_raw <- readBin(approval_paths[[2L]], "raw", n = file.info(approval_paths[[2L]])$size)
  if (!isTRUE(tryCatch(sodium::sig_verify(message_raw, signature_raw,
      mvt_phase2_pinned_public_key()), error = function(e) FALSE))) {
    stop("Phase 2 snapshot approval lacks a valid detached production signature.", call. = FALSE)
  }
  approval <- tryCatch(unserialize(message_raw), error = function(e) NULL)
  expected_names <- c(
    "schema_version", "study", "snapshot_sha256", "snapshot_schema_version",
    "producer_id", "producer_version", "configuration_fingerprint",
    "audit_sha256", "artifact_manifest_sha256", "checkpoint_manifest_sha256",
    "approved_at_utc", "approver"
  )
  if (!is.list(approval) || !identical(names(approval), expected_names) ||
      !identical(approval$schema_version, 2L) || !identical(approval$study, "multivariate-benchmark") ||
      !is.character(approval$approver) || length(approval$approver) != 1L ||
      !nzchar(trimws(approval$approver))) {
    stop("Phase 2 detached approval schema/study/approver is invalid.", call. = FALSE)
  }
  mvt_phase2_parse_approval_time(approval$approved_at_utc)
  approval
}

mvt_phase2_snapshot_trust_sha256 <- function(run_dir) {
  path <- file.path(run_dir, "aggregate_snapshot.rds")
  if (!file.exists(path)) stop("Phase 2 candidate snapshot is missing.", call. = FALSE)
  mvt_sha256_file(path)
}

mvt_write_phase2_snapshot_candidate <- function(run_dir, snapshot) {
  committed <- mvt_validate_committed_snapshot(run_dir, require_production = TRUE)
  if (!isTRUE(snapshot$production_run) ||
      !identical(mvt_hash_object(snapshot$checkpoint_manifest), mvt_hash_object(committed$checkpoint_manifest)) ||
      !identical(mvt_hash_object(snapshot$artifacts), mvt_hash_object(committed$artifacts))) {
    stop("Phase 2 candidate input does not match the reconciled exact production snapshot.", call. = FALSE)
  }
  audit <- mvt_phase2_audit_from_committed_snapshot(committed)
  if (!mvt_phase2_production_eligible(audit)) {
    stop("Phase 2 candidate emission requires every registered production audit check to pass.", call. = FALSE)
  }
  candidate_hash <- mvt_phase2_snapshot_trust_sha256(run_dir)
  candidate_dir <- file.path(mvt_output_root, "snapshot-candidates")
  dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)
  run_id <- gsub("[^A-Za-z0-9_.-]+", "_", basename(normalizePath(run_dir, winslash = "/", mustWork = TRUE)))
  path <- file.path(candidate_dir, paste0(run_id, "-", substr(candidate_hash, 1L, 16L), ".csv"))
  row <- data.frame(
    registry_version = 1L, snapshot_schema_version = snapshot$schema_version,
    profile = "full", status = "candidate_pending_independent_promotion",
    snapshot_sha256 = candidate_hash,
    producer_id = snapshot$producer_id, producer_version = snapshot$producer_version,
    configuration_fingerprint = snapshot$configuration_fingerprint,
    run_dir = normalizePath(run_dir, winslash = "/", mustWork = TRUE),
    created_at = snapshot$committed_at, stringsAsFactors = FALSE
  )
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = candidate_dir)
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  utils::write.csv(row, temporary, row.names = FALSE)
  if (file.exists(path)) {
    existing <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!identical(existing, row)) stop("Phase 2 candidate registry path already contains different content.", call. = FALSE)
  } else if (!file.rename(temporary, path)) {
    stop("Could not atomically install the external Phase 2 candidate record.", call. = FALSE)
  }
  invisible(path)
}

mvt_validate_approved_snapshot <- function(
    run_dir, attestation_path = mvt_phase2_snapshot_attestation_path(),
    require_production = TRUE,
    signature_path = mvt_phase2_snapshot_signature_path()) {
  snapshot <- mvt_validate_committed_snapshot(run_dir, require_production = require_production)
  candidate_hash <- mvt_phase2_snapshot_trust_sha256(run_dir)
  audit <- mvt_phase2_audit_from_committed_snapshot(snapshot)
  if (!mvt_phase2_production_eligible(audit)) stop("Phase 2 approved snapshot fails its production audit.", call. = FALSE)
  approval <- mvt_read_signed_snapshot_approval(attestation_path, signature_path)
  if (!identical(approval$snapshot_sha256, candidate_hash) ||
      !identical(approval$snapshot_schema_version, snapshot$schema_version) ||
      !identical(approval$producer_id, snapshot$producer_id) ||
      !identical(approval$producer_version, snapshot$producer_version) ||
      !identical(approval$configuration_fingerprint, snapshot$configuration_fingerprint) ||
      !identical(approval$audit_sha256, mvt_hash_object(audit)) ||
      !identical(approval$artifact_manifest_sha256, mvt_hash_object(snapshot$artifacts)) ||
      !identical(approval$checkpoint_manifest_sha256, mvt_hash_object(snapshot$checkpoint_manifest))) {
    stop("Phase 2 signed approval does not bind the immutable snapshot/configuration/producer/audit.", call. = FALSE)
  }
  after <- mvt_validate_committed_snapshot(run_dir, require_production = require_production)
  if (!identical(candidate_hash, mvt_phase2_snapshot_trust_sha256(run_dir)) ||
      !identical(mvt_hash_object(snapshot$artifacts), mvt_hash_object(after$artifacts)) ||
      !identical(mvt_hash_object(snapshot$checkpoint_manifest), mvt_hash_object(after$checkpoint_manifest)) ||
      !identical(mvt_hash_object(audit), mvt_hash_object(mvt_phase2_audit_from_committed_snapshot(after)))) {
    stop("Phase 2 snapshot changed after signed approval validation.", call. = FALSE)
  }
  attr(snapshot, "approved_snapshot_sha256") <- candidate_hash
  attr(snapshot, "approval") <- approval
  snapshot
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

mvt_completed_case_ids <- function(existing, planned_methods = NULL) {
  if (!is.null(planned_methods) && length(planned_methods)) {
    status <- existing[["fit_status"]]
    benchmark <- existing[["benchmark_results"]]
    if (!is.data.frame(status) || !is.data.frame(benchmark) ||
        !all(c("case_id", "method") %in% names(status)) ||
        !all(c("case_id", "method") %in% names(benchmark))) return(character())
    cases <- intersect(unique(as.character(status$case_id)), unique(as.character(benchmark$case_id)))
    complete <- vapply(cases, function(case_id) {
      status_methods <- status$method[status$case_id == case_id]
      benchmark_methods <- benchmark$method[benchmark$case_id == case_id]
      all(vapply(planned_methods, function(method) sum(status_methods == method) == 1L, logical(1L))) &&
        all(vapply(planned_methods, function(method) sum(benchmark_methods == method) == 1L, logical(1L)))
    }, logical(1L))
    return(cases[complete])
  }
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

mvt_write_preflight <- function(grid, run_dir, require_gamcopula = TRUE, resume = TRUE, lease = NULL) {
  owns_lease <- is.null(lease)
  if (owns_lease) {
    lease <- mvt_acquire_run_lock(run_dir)
    on.exit(mvt_release_run_lock(lease), add = TRUE)
  }
  active <- mvt_active_comparators()
  packages <- c("gamlss", "gamlss.dist", "mvtnorm", "VineCopula", "callr")
  if (isTRUE(require_gamcopula) || any(grepl("^gamCopula", active))) packages <- c(packages, "gamCopula")
  if (any(grepl("^gee_", active))) packages <- c(packages, "geepack")
  if (any(c("glmm", "glmm_slope") %in% active)) packages <- c(packages, "lme4")
  packages <- unique(packages)
  package_status <- mvt_require_namespaces(packages, strict = FALSE)
  existing_outputs <- file.exists(file.path(run_dir, paste0(mvt_result_names(), "_by_rep.csv"))) |
    file.exists(file.path(run_dir, paste0(mvt_result_names(), "_checkpoint.csv")))
  existing_case_checkpoints <- length(list.files(
    mvt_case_checkpoint_dir(run_dir), pattern = "[.]rds$", full.names = TRUE
  ))
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
      if (isTRUE(resume) || (!any(existing_outputs) && existing_case_checkpoints == 0L)) "pass" else "warn",
      paste(
        "resume", resume, "| existing aggregate files", sum(existing_outputs),
        "| case checkpoints", existing_case_checkpoints
      ),
      sum(existing_outputs) + existing_case_checkpoints
    ),
    mvt_preflight_check(
      "gee_unstructured_timeout",
      if (!unstructured_t50 || is.finite(unstructured_timeout)) "pass" else "warn",
      paste("T>=50 unstructured GEE", unstructured_t50, "| timeout", unstructured_timeout),
      if (is.finite(unstructured_timeout)) unstructured_timeout else NA_integer_
    )
  )
  out <- mvt_bind_rows_fill(checks)
  mvt_write_csv_atomic(out, file.path(run_dir, "preflight_checks.csv"), lease = lease)
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
  mvt_write_lines_atomic(lines, file.path(run_dir, "preflight_checks.md"), lease = lease)
  invisible(out)
}

mvt_write_run_metadata <- function(
    run_dir, grid, seed_base, require_gamcopula, resume,
    workers = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_WORKERS", 1L),
    workers_used = workers,
    fingerprints = mvt_checkpoint_fingerprints(), configuration_fingerprint = "",
    run_started = Sys.time(), run_completed = NA, run_elapsed_sec = NA_real_, lease = NULL) {
  owns_lease <- is.null(lease)
  if (owns_lease) {
    lease <- mvt_acquire_run_lock(run_dir)
    on.exit(mvt_release_run_lock(lease), add = TRUE)
  }
  runtime <- mvt_runtime_identity()
  hardware <- mvt_hardware_identity()
  gc_state <- gc()
  metadata <- data.frame(
    name = c(
      "created_at",
      "run_started_at", "run_completed_at", "run_elapsed_sec",
      "repo_root",
      "seed_base",
      "require_gamcopula",
      "resume",
      "active_comparators",
      "n_cases",
      "workers_requested",
      "workers_used",
      "checkpoint_schema_version",
      "identity_algorithm", "identity_version", "producer_id", "producer_version",
      "checkpoint_configuration_fingerprint",
      "producer_fingerprint",
      "code_fingerprint",
      "package_fingerprint",
      "r_version", "platform", "rng_kind", "blas", "lapack",
      "hostname", "cpu_model", "logical_cores", "physical_cores", "ram_bytes", "parent_pid",
      "gc_peak_ncells", "gc_peak_vcells"
    ),
    value = c(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      format(run_started, "%Y-%m-%d %H:%M:%S %Z"),
      if (length(run_completed) == 1L && !is.na(run_completed)) format(run_completed, "%Y-%m-%d %H:%M:%S %Z") else "",
      as.character(run_elapsed_sec),
      mvt_repo_root,
      as.character(seed_base),
      as.character(require_gamcopula),
      as.character(resume),
      paste(mvt_active_comparators(), collapse = ","),
      as.character(nrow(grid)),
      as.character(workers),
      as.character(workers_used),
      as.character(mvt_checkpoint_schema_version),
      fingerprints$algorithm %||% mvt_identity_algorithm,
      fingerprints$identity_version %||% mvt_identity_version,
      fingerprints$producer_id %||% mvt_producer_id,
      fingerprints$producer_version %||% mvt_producer_version,
      configuration_fingerprint,
      fingerprints$producer_fingerprint,
      fingerprints$code_fingerprint,
      fingerprints$package_fingerprint,
      runtime$r_version, runtime$platform, paste(runtime$rng_kind, collapse = "/"), runtime$blas, runtime$lapack,
      hardware$hostname, hardware$cpu_model, as.character(hardware$logical_cores), as.character(hardware$physical_cores),
      as.character(hardware$ram_bytes), as.character(hardware$parent_pid),
      as.character(gc_state["Ncells", "max used"]), as.character(gc_state["Vcells", "max used"])
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
  mvt_write_csv_atomic(metadata, file.path(run_dir, "run_metadata.csv"), lease = lease)
  recorded_versions <- fingerprints$package_versions %||% mvt_package_versions(c(
    "gamlss.longitudinal", "gamlss", "gamlss.dist", "gamCopula", "VineCopula",
    "mvtnorm", "geepack", "lme4", "mgcv", "callr", "scoringRules", "ggplot2"
  ))
  mvt_write_csv_atomic(
    recorded_versions,
    file.path(run_dir, "package_versions.csv"), lease = lease
  )
  mvt_write_lines_atomic(
    capture.output(utils::sessionInfo()), file.path(run_dir, "session_info.txt"), lease = lease
  )
  invisible(run_dir)
}

mvt_run_grid <- function(
    grid, run_dir, seed_base = 20260818L, checkpoint_every = 5L,
    require_gamcopula = TRUE,
    workers = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_WORKERS", 1L),
    run_case_fun = mvt_run_case,
    worker_load_local = TRUE,
    fingerprints = mvt_checkpoint_fingerprints()) {
  run_started <- Sys.time()
  dir.create(mvt_output_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  run_lock <- mvt_acquire_run_lock(run_dir)
  on.exit(mvt_release_run_lock(run_lock), add = TRUE)
  workers <- suppressWarnings(as.integer(workers))
  if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    stop("GAMLSS_LONGITUDINAL_MVT_WORKERS must be a positive integer.", call. = FALSE)
  }
  resume <- mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_RESUME", TRUE)
  methods <- mvt_active_comparators()
  grid <- grid[mvt_canonical_grid_order(grid), , drop = FALSE]
  rownames(grid) <- NULL
  tasks <- mvt_prepare_tasks(
    grid, seed_base, require_gamcopula, fingerprints,
    workers_requested = workers, run_started_at = run_started
  )
  production_run <- mvt_grid_is_exact_production(grid)
  if (production_run) {
    current_fingerprints <- mvt_checkpoint_fingerprints()
    fingerprint_fields <- c(
      "algorithm", "identity_version", "producer_id", "producer_version",
      "producer_fingerprint", "code_fingerprint", "package_fingerprint"
    )
    fingerprints_current <- all(vapply(fingerprint_fields, function(field) {
      identical(fingerprints[[field]] %||% NULL, current_fingerprints[[field]] %||% NULL)
    }, logical(1L)))
    if (!identical(run_case_fun, mvt_run_case) || !fingerprints_current) {
      stop("Production multivariate runs require the registered case runner and current producer/code/package fingerprints.", call. = FALSE)
    }
    timeout_check <- mvt_timeout_contract(
      mvt_checkpoint_configuration(seed_base, require_gamcopula), methods
    )
    if (!timeout_check$valid) {
      stop(
        "Production multivariate runs require finite positive registered timeouts for every GEE/callr method: ",
        paste(timeout_check$required, timeout_check$values, sep = "=", collapse = ", "),
        call. = FALSE
      )
    }
    if (!isTRUE(worker_load_local) ||
        !identical(mvt_env("GAMLSS_LONGITUDINAL_MVT_SOURCE", "installed"), "local")) {
      stop("Production multivariate runs require fully attested checked-out local source.", call. = FALSE)
    }
  }
  dir.create(mvt_case_checkpoint_dir(run_dir), recursive = TRUE, showWarnings = FALSE)

  existed <- vapply(tasks, function(task) {
    file.exists(mvt_case_checkpoint_path(run_dir, task$case_id))
  }, logical(1L))
  rejections <- list()
  checkpoints <- lapply(tasks, function(task) {
    path <- mvt_case_checkpoint_path(run_dir, task$case_id)
    checkpoint <- if (isTRUE(resume)) mvt_read_case_checkpoint(path, task) else NULL
    rejected <- inherits(checkpoint, "mvt_rejected_checkpoint")
    if (file.exists(path) && (rejected || !isTRUE(resume))) {
      reason <- if (isTRUE(resume)) attr(checkpoint, "rejection_reason") %||% "invalid" else "resume_disabled"
      quarantine <- mvt_archive_stale_checkpoint(path, run_dir, reason)
      rejections[[length(rejections) + 1L]] <<- data.frame(
        case_id = task$case_id, source_path = path,
        quarantine_path = normalizePath(quarantine, winslash = "/", mustWork = FALSE),
        rejection_reason = reason, rejected_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        parent_pid = Sys.getpid(), stringsAsFactors = FALSE
      )
      checkpoint <- NULL
    }
    checkpoint
  })
  rejection_ledger <- mvt_bind_rows_fill(rejections)
  if (!nrow(rejection_ledger)) rejection_ledger <- data.frame(
    case_id = character(), source_path = character(), quarantine_path = character(),
    rejection_reason = character(), rejected_at = character(), parent_pid = integer(),
    stringsAsFactors = FALSE
  )
  prior_rejections <- mvt_read_optional_csv(file.path(run_dir, "checkpoint_rejections.csv"))
  if (nrow(prior_rejections)) {
    rejection_ledger <- mvt_bind_rows_fill(prior_rejections, rejection_ledger)
  }
  if (nrow(rejection_ledger)) {
    rejection_ledger <- unique(rejection_ledger)
    rownames(rejection_ledger) <- NULL
  }
  mvt_write_csv_atomic(rejection_ledger, file.path(run_dir, "checkpoint_rejections.csv"), lease = run_lock)
  completed <- !vapply(checkpoints, is.null, logical(1L))
  pending <- tasks[!completed]
  workers_used <- if (length(pending)) min(workers, length(pending)) else 0L
  if (workers_used > 1L && isTRUE(worker_load_local) &&
      !identical(mvt_env("GAMLSS_LONGITUDINAL_MVT_SOURCE", "installed"), "local")) {
    stop(
      "Parallel multivariate runs require GAMLSS_LONGITUDINAL_MVT_SOURCE=local so parent and PSOCK workers use the checked-out code.",
      call. = FALSE
    )
  }
  expected_execution <- tasks[[1L]]$execution_attestation
  expected_package_identity <- expected_execution$package_identity
  parent_attestation <- if (identical(mvt_env("GAMLSS_LONGITUDINAL_MVT_SOURCE", "installed"), "local")) {
    mvt_verify_execution_attestation(expected_execution, load = identical(run_case_fun, mvt_run_case))
  } else {
    c(expected_package_identity, list(
      verified = NA, loaded_namespace_path = NA_character_, loaded_version = expected_package_identity$version,
      verified_source_sha256 = expected_package_identity$source_sha256,
      runtime_identity = mvt_runtime_identity(), libpaths = paste(.libPaths(), collapse = ";")
    ))
  }
  attestation_rows <- list(mvt_attestation_row(parent_attestation, role = "parent", pid = Sys.getpid()))

  mvt_write_csv_atomic(grid, file.path(run_dir, "scenario_grid.csv"), lease = run_lock)
  mvt_write_preflight(grid, run_dir, require_gamcopula = require_gamcopula, resume = resume, lease = run_lock)
  mvt_write_run_metadata(
    run_dir, grid, seed_base, require_gamcopula, resume,
    workers = workers, workers_used = workers_used, fingerprints = fingerprints,
    configuration_fingerprint = tasks[[1L]]$configuration_fingerprint,
    run_started = run_started, lease = run_lock
  )
  message(sprintf(
    "Multivariate checkpoint audit: existing=%d valid=%d rejected=%d pending=%d workers=%d",
    sum(existed), sum(completed), sum(existed & !completed), length(pending), workers_used
  ))

  if (length(pending) && workers_used == 1L) {
    for (i in seq_along(pending)) {
      task <- pending[[i]]
      message("[", i, "/", length(pending), "] ", task$case_id)
      mvt_run_checkpoint_task(task, run_dir, run_case_fun, writer_role = "serial_parent")
    }
  } else if (length(pending)) {
    message("Starting ", workers_used, " Windows-safe PSOCK workers for ", length(pending), " unique case(s).")
    cluster <- parallel::makePSOCKcluster(workers_used, outfile = "")
    on.exit(if (!is.null(cluster)) try(parallel::stopCluster(cluster), silent = TRUE), add = TRUE)
    setup_path <- normalizePath(
      file.path(mvt_script_dir, "00-multivariate-setup.R"),
      winslash = "/", mustWork = TRUE
    )
    worker_attestations <- parallel::clusterCall(cluster, function(path, load_local, expected_execution) {
      source(path, local = .GlobalEnv)
      if (isTRUE(load_local)) {
        Sys.setenv(GAMLSS_LONGITUDINAL_MVT_SOURCE = "local")
        identity <- mvt_verify_execution_attestation(expected_execution, load = TRUE)
      } else {
        expected_identity <- expected_execution$package_identity
        identity <- c(expected_identity, list(
          verified = NA, loaded_namespace_path = NA_character_, loaded_version = expected_identity$version,
          verified_source_sha256 = expected_identity$source_sha256, full_verified = FALSE,
          runtime_identity = mvt_runtime_identity(), libpaths = paste(.libPaths(), collapse = ";")
        ))
      }
      list(pid = Sys.getpid(), setup_path = normalizePath(path, winslash = "/", mustWork = TRUE), identity = identity)
    }, setup_path, isTRUE(worker_load_local), expected_execution)
    if (isTRUE(worker_load_local)) {
      verified <- vapply(worker_attestations, function(x) {
        mvt_execution_attestation_matches(x$identity, expected_execution)
      }, logical(1L))
      if (!all(verified)) stop("One or more PSOCK worker checkout attestations failed.", call. = FALSE)
    }
    parallel::parLapplyLB(
      cluster,
      pending,
      function(task, destination, case_runner) {
        mvt_run_checkpoint_task(task, destination, case_runner, writer_role = "psock_worker")
      },
      destination = normalizePath(run_dir, winslash = "/", mustWork = TRUE),
      case_runner = run_case_fun
    )
    post_worker_attestations <- parallel::clusterCall(cluster, function(path, load_local, expected_execution) {
      identity <- if (isTRUE(load_local)) {
        mvt_verify_execution_attestation(expected_execution, load = FALSE)
      } else {
        expected_identity <- expected_execution$package_identity
        c(expected_identity, list(
          verified = NA, loaded_namespace_path = NA_character_, loaded_version = expected_identity$version,
          verified_source_sha256 = expected_identity$source_sha256, full_verified = FALSE,
          runtime_identity = mvt_runtime_identity(), libpaths = paste(.libPaths(), collapse = ";")
        ))
      }
      identity$gc_state <- gc()
      list(pid = Sys.getpid(), setup_path = normalizePath(path, winslash = "/", mustWork = TRUE), identity = identity)
    }, setup_path, isTRUE(worker_load_local), expected_execution)
    if (isTRUE(worker_load_local)) {
      post_verified <- vapply(post_worker_attestations, function(x) {
        mvt_execution_attestation_matches(x$identity, expected_execution)
      }, logical(1L))
      if (!all(post_verified)) stop("One or more post-task PSOCK worker checkout attestations failed.", call. = FALSE)
    }
    attestation_rows <- c(attestation_rows, lapply(post_worker_attestations, function(x) {
      mvt_attestation_row(x$identity, role = "worker", pid = x$pid, setup_path = x$setup_path)
    }))
    parallel::stopCluster(cluster)
    cluster <- NULL
  }

  checkpoint_manifest_before <- mvt_checkpoint_manifest(tasks, run_dir)
  mvt_write_csv_atomic(
    mvt_bind_rows_fill(attestation_rows), file.path(run_dir, "worker_attestations.csv"), lease = run_lock
  )

  rows <- mvt_collect_case_checkpoints(tasks, run_dir, grid, methods)
  for (nm in names(rows)) {
    mvt_write_csv_atomic(rows[[nm]], file.path(run_dir, paste0(nm, "_by_rep.csv")), lease = run_lock)
  }
  run_dir_norm <- normalizePath(run_dir, winslash = "/", mustWork = FALSE)
  output_root_norm <- normalizePath(mvt_output_root, winslash = "/", mustWork = FALSE)
  if (startsWith(run_dir_norm, output_root_norm)) {
    writeLines(run_dir, file.path(mvt_output_root, "latest_run_dir.txt"), useBytes = TRUE)
  }
  run_completed <- Sys.time()
  mvt_write_run_metadata(
    run_dir, grid, seed_base, require_gamcopula, resume,
    workers = workers, workers_used = workers_used, fingerprints = fingerprints,
    configuration_fingerprint = tasks[[1L]]$configuration_fingerprint,
    run_started = run_started, run_completed = run_completed,
    run_elapsed_sec = as.numeric(difftime(run_completed, run_started, units = "secs")),
    lease = run_lock
  )
  snapshot <- mvt_write_aggregate_snapshot(
    tasks = tasks, run_dir = run_dir, grid = grid, methods = methods, rows = rows,
    checkpoint_manifest_before = checkpoint_manifest_before,
    configuration = mvt_checkpoint_configuration(seed_base, require_gamcopula),
    production_run = production_run, lease = run_lock
  )
  if (isTRUE(production_run)) mvt_write_phase2_snapshot_candidate(run_dir, snapshot)
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
  lease <- mvt_acquire_run_lock(target_dir)
  on.exit(mvt_release_run_lock(lease), add = TRUE)

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
  mvt_write_csv_atomic(grid, file.path(target_dir, "scenario_grid.csv"), lease)

  for (name in mvt_result_names()) {
    combined <- mvt_bind_rows_fill(lapply(source_dirs, read_result_table, name = name))
    if ("source_shard" %in% names(combined)) combined$source_shard <- NULL
    if (nrow(combined) > 0L) {
      key_cols <- intersect(c("case_id", "method", "parameter", "term", "dependence_scope"), names(combined))
      if (length(key_cols) > 0L) combined <- unique(combined)
      mvt_write_csv_atomic(combined, file.path(target_dir, paste0(name, "_by_rep.csv")), lease)
    }
  }

  preflight <- mvt_bind_rows_fill(lapply(source_dirs, read_table, file = "preflight_checks.csv"))
  if (nrow(preflight) > 0L) mvt_write_csv_atomic(preflight, file.path(target_dir, "preflight_checks.csv"), lease)
  metadata <- mvt_bind_rows_fill(lapply(source_dirs, read_table, file = "run_metadata.csv"))
  if (nrow(metadata) > 0L) mvt_write_csv_atomic(metadata, file.path(target_dir, "run_metadata.csv"), lease)
  versions <- unique(mvt_bind_rows_fill(lapply(source_dirs, read_table, file = "package_versions.csv")))
  if ("source_shard" %in% names(versions)) versions$source_shard <- NULL
  if (nrow(versions) > 0L) mvt_write_csv_atomic(versions, file.path(target_dir, "package_versions.csv"), lease)
  mvt_write_lines_atomic(
    c(
      "Merged multivariate longitudinal simulation shards.",
      paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      "Source shards:",
      paste0("- ", source_dirs)
    ),
    file.path(target_dir, "session_info.txt"), lease
  )
  mvt_write_lines_atomic(
    c(
      "# Merged Preflight",
      "",
      paste("Source shards:", length(source_dirs)),
      "",
      paste0("- ", basename(source_dirs))
    ),
    file.path(target_dir, "preflight_checks.md"), lease
  )
  mvt_write_lines_atomic(
    c(
      "NONPUBLICATION: merged legacy shards have no case-checkpoint hash reconciliation.",
      "Use mvt_run_grid() to produce a committed immutable aggregate snapshot."
    ),
    file.path(target_dir, "NONPUBLICATION.txt"), lease
  )
  invisible(target_dir)
}

mvt_summarise_results <- function(run_dir = mvt_read_run_dir()) {
  lease <- mvt_acquire_run_lock(run_dir)
  released <- FALSE
  on.exit(if (!released) mvt_release_run_lock(lease), add = TRUE)
  committed <- mvt_validate_committed_snapshot(run_dir, require_production = FALSE)
  results <- committed$rows$benchmark_results
  coefs <- committed$rows$coefficient_results
  if (nrow(results) == 0L) stop("benchmark_results_by_rep.csv is missing or empty.", call. = FALSE)
  if (nrow(coefs) == 0L) stop("coefficient_results_by_rep.csv is missing or empty.", call. = FALSE)
  dep <- committed$rows$dependence_recovery
  vario <- committed$rows$variogram_scores

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
  mvt_write_csv_atomic(result_summary, file.path(run_dir, "benchmark_summary.csv"), lease)
  mvt_write_csv_atomic(coef_summary, file.path(run_dir, "coefficient_summary.csv"), lease)
  mvt_write_csv_atomic(dep_summary, file.path(run_dir, "dependence_recovery_summary.csv"), lease)
  mvt_write_csv_atomic(vario_summary, file.path(run_dir, "variogram_summary.csv"), lease)
  status_path <- file.path(run_dir, "fit_status_by_rep.csv")
  if (file.exists(status_path)) {
    status <- committed$rows$fit_status
    mvt_write_csv_atomic(mvt_case_method_completion_summary(status), file.path(run_dir, "case_method_completion_summary.csv"), lease)
  }
  feasibility <- mvt_write_pilot_feasibility(run_dir, lease = lease)
  mvt_release_run_lock(lease)
  released <- TRUE
  phase2 <- if (isTRUE(committed$production_run)) {
    approval <- tryCatch(mvt_validate_approved_snapshot(run_dir), error = function(e) e)
    if (inherits(approval, "error")) {
      list(
        production_eligible = FALSE,
        reason = "candidate_pending_independent_promotion",
        detail = conditionMessage(approval),
        candidate_snapshot_sha256 = mvt_phase2_snapshot_trust_sha256(run_dir)
      )
    } else {
      mvt_write_phase2_benchmark_evidence(run_dir)
    }
  } else {
    list(production_eligible = FALSE, reason = "nonproduction aggregate snapshot")
  }
  invisible(list(
    results = result_summary,
    coefficients = coef_summary,
    dependence = dep_summary,
    variogram = vario_summary,
    feasibility = feasibility,
    phase2 = phase2
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

mvt_write_pilot_feasibility <- function(run_dir = mvt_read_run_dir(), target_reps = mvt_env_int("GAMLSS_LONGITUDINAL_MVT_TARGET_REPS", 100L), lease = NULL) {
  owns_lease <- is.null(lease)
  if (owns_lease) {
    lease <- mvt_acquire_run_lock(run_dir)
    on.exit(mvt_release_run_lock(lease), add = TRUE)
  }
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

  mvt_write_csv_atomic(method_summary, file.path(run_dir, "pilot_feasibility_by_method.csv"), lease)
  mvt_write_csv_atomic(scenario_summary, file.path(run_dir, "pilot_feasibility_by_scenario.csv"), lease)
  mvt_write_csv_atomic(overall_method, file.path(run_dir, "pilot_feasibility_overall_method.csv"), lease)

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
  mvt_write_lines_atomic(lines, file.path(run_dir, "pilot_feasibility.md"), lease)
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
      sha256 = unname(vapply(files, mvt_sha256_file, character(1L))),
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
  lines <- c(lines, "", "See `artifact_manifest.csv` for file paths, sizes, timestamps, and SHA-256 hashes.")
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
  phase2_audit <- read_optional("phase2_benchmark_audit.csv")

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
      "phase2_benchmark_audit_clean",
      if (nrow(phase2_audit) > 0L && "status" %in% names(phase2_audit) && !any(phase2_audit$status == "fail")) "pass" else "fail",
      if (nrow(phase2_audit) > 0L && "status" %in% names(phase2_audit)) {
        paste("failures", sum(phase2_audit$status == "fail"), "| checks", nrow(phase2_audit))
      } else {
        "phase2_benchmark_audit.csv missing or malformed"
      },
      if (nrow(phase2_audit) > 0L && "status" %in% names(phase2_audit)) sum(phase2_audit$status == "fail") else NA_integer_
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

mvt_wilson_interval <- function(events, attempts, level = 0.95) {
  if (!is.finite(events) || !is.finite(attempts) || attempts < 1L || events < 0L || events > attempts) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- events / attempts
  denominator <- 1 + z^2 / attempts
  centre <- (p + z^2 / (2 * attempts)) / denominator
  half_width <- z * sqrt(p * (1 - p) / attempts + z^2 / (4 * attempts^2)) / denominator
  c(lower = max(0, centre - half_width), upper = min(1, centre + half_width))
}

mvt_mean_interval <- function(x, level = 0.95) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(c(mean = NA_real_, mcse = NA_real_, lower = NA_real_, upper = NA_real_))
  estimate <- mean(x)
  mcse <- if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
  critical <- if (length(x) > 1L) stats::qt(1 - (1 - level) / 2, df = length(x) - 1L) else NA_real_
  c(
    mean = estimate,
    mcse = mcse,
    lower = if (is.finite(mcse)) estimate - critical * mcse else NA_real_,
    upper = if (is.finite(mcse)) estimate + critical * mcse else NA_real_
  )
}

mvt_phase2_metadata_value <- function(metadata, name, default = "") {
  if (!nrow(metadata) || !all(c("name", "value") %in% names(metadata))) return(default)
  value <- metadata$value[metadata$name == name]
  if (!length(value) || is.na(value[[1L]]) || !nzchar(value[[1L]])) default else as.character(value[[1L]])
}

mvt_phase2_planned_methods <- function(metadata) {
  value <- mvt_phase2_metadata_value(metadata, "active_comparators", "")
  if (!nzchar(value)) value <- mvt_phase2_metadata_value(metadata, "env.GAMLSS_LONGITUDINAL_MVT_COMPARATORS", "")
  if (!nzchar(value)) return(character())
  methods <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  unique(methods[nzchar(methods)])
}

mvt_phase2_attempt_reconciliation <- function(grid, status, benchmark, planned_methods) {
  if (!nrow(grid) || !"case_id" %in% names(grid) || !length(planned_methods)) return(data.frame())
  grid_cases <- unique(as.character(grid$case_id))
  expected <- merge(
    data.frame(case_id = grid_cases, stringsAsFactors = FALSE),
    data.frame(method = planned_methods, stringsAsFactors = FALSE),
    by = NULL
  )
  expected$planned <- TRUE
  count_keys <- function(data, count_name) {
    if (!nrow(data) || !all(c("case_id", "method") %in% names(data))) {
      out <- expected[FALSE, c("case_id", "method"), drop = FALSE]
      out[[count_name]] <- integer()
      return(out)
    }
    key <- paste(data$case_id, data$method, sep = "\r")
    tab <- table(key)
    split_key <- strsplit(names(tab), "\r", fixed = TRUE)
    out <- data.frame(
      case_id = vapply(split_key, `[[`, character(1L), 1L),
      method = vapply(split_key, `[[`, character(1L), 2L),
      stringsAsFactors = FALSE
    )
    out[[count_name]] <- as.integer(tab)
    out
  }
  status_counts <- count_keys(status, "status_rows")
  benchmark_counts <- count_keys(benchmark, "benchmark_rows")
  out <- merge(expected, status_counts, by = c("case_id", "method"), all = TRUE)
  out <- merge(out, benchmark_counts, by = c("case_id", "method"), all = TRUE)
  out$planned[is.na(out$planned)] <- FALSE
  out$status_rows[is.na(out$status_rows)] <- 0L
  out$benchmark_rows[is.na(out$benchmark_rows)] <- 0L
  out$status_exactly_once <- out$planned & out$status_rows == 1L
  out$benchmark_exactly_once <- out$planned & out$benchmark_rows == 1L
  out$unexpected <- !out$planned
  out$reconciled <- out$planned & out$status_exactly_once & out$benchmark_exactly_once
  out[order(out$case_id, out$method), , drop = FALSE]
}

mvt_phase2_penalty <- function(values) {
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]
  if (!length(values)) return(NA_real_)
  worst <- max(values)
  worst + 0.10 * max(1, abs(worst))
}

mvt_phase2_penalty_rule <- function() {
  "within each scenario/family/time cell, compute the pool from successful finite rows across the complete prespecified planned-method set; replace every failed or non-finite lower-is-better metric with the pool maximum plus 10% of max(1, abs(maximum)); fixed before summaries"
}

mvt_phase2_group_index <- function(data, group_row, group_cols) {
  idx <- rep(TRUE, nrow(data))
  for (col in group_cols) idx <- idx & data[[col]] == group_row[[col]][1L]
  idx[is.na(idx)] <- FALSE
  idx
}

mvt_phase2_method_summary <- function(status, benchmark, methods, planned_methods = methods) {
  if (!nrow(status) || !length(methods)) return(data.frame())
  all_status <- mvt_status_from_results(status)
  all_benchmark <- benchmark
  planned_methods <- unique(as.character(planned_methods[nzchar(planned_methods)]))
  if (!length(planned_methods)) planned_methods <- unique(as.character(all_status$method))
  status <- all_status
  status <- status[status$method %in% methods, , drop = FALSE]
  if (!nrow(status)) return(data.frame())
  group_cols <- intersect(
    c("scenario", "generator", "dependence", "correlation_level", "n_time", "n_subject", "total_rows", "family", "gamlss_family", "method"),
    names(status)
  )
  groups <- unique(status[group_cols])
  metric_key <- paste(benchmark$case_id, benchmark$method, sep = "\r")
  status_key <- paste(status$case_id, status$method, sep = "\r")
  metric_match <- match(status_key, metric_key)
  rmse_metric <- suppressWarnings(as.numeric(benchmark$benchmark_mean_rmse[metric_match]))
  nls_metric <- if ("benchmark_neg_log_score" %in% names(benchmark)) {
    suppressWarnings(as.numeric(benchmark$benchmark_neg_log_score[metric_match]))
  } else {
    rep(NA_real_, nrow(status))
  }
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    idx <- mvt_phase2_group_index(status, groups[i, , drop = FALSE], group_cols)
    sub <- status[idx, , drop = FALSE]
    success <- mvt_status_success(sub$success)
    rmse_values <- rmse_metric[idx]
    nls_values <- nls_metric[idx]
    elapsed <- suppressWarnings(as.numeric(sub$elapsed_sec))
    included <- success & is.finite(rmse_values)
    success_interval <- mvt_wilson_interval(sum(success), nrow(sub))
    accuracy <- mvt_mean_interval(rmse_values[included])
    nls_included <- success & is.finite(nls_values)
    nls_conditional <- mvt_mean_interval(nls_values[nls_included])
    elapsed_included <- success & is.finite(elapsed)
    elapsed_conditional <- mvt_mean_interval(elapsed[elapsed_included])
    cell_cols <- setdiff(group_cols, "method")
    cell_status_idx <- mvt_phase2_group_index(all_status, groups[i, , drop = FALSE], intersect(cell_cols, names(all_status)))
    cell_status <- all_status[cell_status_idx & all_status$method %in% planned_methods, , drop = FALSE]
    retained_keys <- unique(paste(
      cell_status$case_id[mvt_status_success(cell_status$success)],
      cell_status$method[mvt_status_success(cell_status$success)],
      sep = "\r"
    ))
    cell_idx <- mvt_phase2_group_index(all_benchmark, groups[i, , drop = FALSE], intersect(cell_cols, names(all_benchmark)))
    cell_benchmark <- all_benchmark[cell_idx & all_benchmark$method %in% planned_methods, , drop = FALSE]
    cell_benchmark_key <- paste(cell_benchmark$case_id, cell_benchmark$method, sep = "\r")
    cell_benchmark <- cell_benchmark[cell_benchmark_key %in% retained_keys, , drop = FALSE]
    rmse_penalty <- mvt_phase2_penalty(cell_benchmark$benchmark_mean_rmse)
    nls_penalty <- if ("benchmark_neg_log_score" %in% names(cell_benchmark)) mvt_phase2_penalty(cell_benchmark$benchmark_neg_log_score) else NA_real_
    elapsed_penalty <- mvt_phase2_penalty(suppressWarnings(as.numeric(
      cell_status$elapsed_sec[mvt_status_success(cell_status$success)]
    )))
    rmse_penalized <- rmse_values
    rmse_penalized[!included] <- rmse_penalty
    nls_penalized <- nls_values
    nls_penalized[!nls_included] <- nls_penalty
    elapsed_penalized <- elapsed
    elapsed_penalized[!elapsed_included] <- elapsed_penalty
    rmse_penalized_interval <- mvt_mean_interval(rmse_penalized)
    nls_penalized_interval <- mvt_mean_interval(nls_penalized)
    elapsed_penalized_interval <- mvt_mean_interval(elapsed_penalized)
    failure_detail <- if (any(!success)) {
      paste0(sub$method[!success], "=", sub$failure_reason_short[!success])
    } else {
      character()
    }
    rows[[i]] <- cbind(
      groups[i, , drop = FALSE],
      data.frame(
        evidence_role = if (identical(groups$method[[i]], "gee_unstructured") && groups$n_time[[i]] >= 20L) "high_dimensional_stress_test" else "routine_comparator",
        n_attempted = nrow(sub),
        n_converged = sum(success),
        convergence_definition = "method returned a successful finite benchmark row after method-specific convergence checks",
        n_included_accuracy = sum(included),
        n_warning = sum(sub$status_class == "warning", na.rm = TRUE),
        n_timeout = sum(sub$status_class == "timeout", na.rm = TRUE),
        n_error = sum(sub$status_class == "error", na.rm = TRUE),
        success_rate = mean(success),
        success_mcse = sqrt(mean(success) * (1 - mean(success)) / nrow(sub)),
        success_ci_lower = success_interval[["lower"]],
        success_ci_upper = success_interval[["upper"]],
        mean_rmse_conditional = accuracy[["mean"]],
        mean_rmse_mcse = accuracy[["mcse"]],
        mean_rmse_ci_lower = accuracy[["lower"]],
        mean_rmse_ci_upper = accuracy[["upper"]],
        mean_rmse_failure_penalized = rmse_penalized_interval[["mean"]],
        mean_rmse_failure_penalized_mcse = rmse_penalized_interval[["mcse"]],
        mean_rmse_failure_penalized_ci_lower = rmse_penalized_interval[["lower"]],
        mean_rmse_failure_penalized_ci_upper = rmse_penalized_interval[["upper"]],
        mean_neg_log_score_conditional = nls_conditional[["mean"]],
        mean_neg_log_score_conditional_mcse = nls_conditional[["mcse"]],
        mean_neg_log_score_conditional_ci_lower = nls_conditional[["lower"]],
        mean_neg_log_score_conditional_ci_upper = nls_conditional[["upper"]],
        mean_neg_log_score_failure_penalized = nls_penalized_interval[["mean"]],
        mean_neg_log_score_failure_penalized_mcse = nls_penalized_interval[["mcse"]],
        mean_neg_log_score_failure_penalized_ci_lower = nls_penalized_interval[["lower"]],
        mean_neg_log_score_failure_penalized_ci_upper = nls_penalized_interval[["upper"]],
        mean_elapsed_sec_conditional = elapsed_conditional[["mean"]],
        mean_elapsed_sec_conditional_mcse = elapsed_conditional[["mcse"]],
        mean_elapsed_sec_conditional_ci_lower = elapsed_conditional[["lower"]],
        mean_elapsed_sec_conditional_ci_upper = elapsed_conditional[["upper"]],
        mean_elapsed_sec_failure_penalized = elapsed_penalized_interval[["mean"]],
        mean_elapsed_sec_failure_penalized_mcse = elapsed_penalized_interval[["mcse"]],
        mean_elapsed_sec_failure_penalized_ci_lower = elapsed_penalized_interval[["lower"]],
        mean_elapsed_sec_failure_penalized_ci_upper = elapsed_penalized_interval[["upper"]],
        rmse_failure_penalty = rmse_penalty,
        neg_log_score_failure_penalty = nls_penalty,
        elapsed_failure_penalty = elapsed_penalty,
        rmse_failure_penalty_pool_n = sum(is.finite(suppressWarnings(as.numeric(cell_benchmark$benchmark_mean_rmse)))),
        neg_log_score_failure_penalty_pool_n = if ("benchmark_neg_log_score" %in% names(cell_benchmark)) sum(is.finite(suppressWarnings(as.numeric(cell_benchmark$benchmark_neg_log_score)))) else 0L,
        elapsed_failure_penalty_pool_n = sum(is.finite(suppressWarnings(as.numeric(cell_status$elapsed_sec[mvt_status_success(cell_status$success)])))),
        failure_penalty_planned_methods = paste(sort(planned_methods), collapse = ","),
        failure_penalty_definition = mvt_phase2_penalty_rule(),
        median_elapsed_sec = if (any(is.finite(elapsed))) stats::median(elapsed, na.rm = TRUE) else NA_real_,
        elapsed_iqr_sec = if (sum(is.finite(elapsed)) > 1L) stats::IQR(elapsed, na.rm = TRUE) else NA_real_,
        failure_reasons = mvt_collapse_nonempty(failure_detail, empty = "none"),
        stringsAsFactors = FALSE
      )
    )
  }
  mvt_bind_rows_fill(rows)
}

mvt_nearest_neighbor_paired_contrasts <- function(status, benchmark, scenario_grid = data.frame()) {
  neighbors <- mvt_nearest_neighbor_comparators()
  methods <- c("gamlss.longitudinal", neighbors)
  status <- mvt_status_from_results(status)
  status <- status[status$method %in% methods, , drop = FALSE]
  benchmark <- benchmark[benchmark$method %in% methods, , drop = FALSE]
  if (!nrow(scenario_grid) || !"case_id" %in% names(scenario_grid)) return(data.frame())
  if (!"family" %in% names(scenario_grid) && "family_name" %in% names(scenario_grid)) scenario_grid$family <- scenario_grid$family_name
  group_cols <- intersect(c("scenario", "generator", "dependence", "correlation_level", "n_time", "family"), names(status))
  if (!all(group_cols %in% names(scenario_grid))) return(data.frame())
  case_groups <- unique(scenario_grid[group_cols])
  metric_names <- intersect(c("benchmark_mean_rmse", "benchmark_neg_log_score", "elapsed_sec"), names(benchmark))
  rows <- list()
  k <- 1L
  for (i in seq_len(nrow(case_groups))) {
    group_idx <- mvt_phase2_group_index(status, case_groups[i, , drop = FALSE], group_cols)
    group_status <- status[group_idx, , drop = FALSE]
    benchmark_idx <- mvt_phase2_group_index(benchmark, case_groups[i, , drop = FALSE], group_cols)
    group_benchmark <- benchmark[benchmark_idx, , drop = FALSE]
    grid_idx <- mvt_phase2_group_index(scenario_grid, case_groups[i, , drop = FALSE], group_cols)
    expected_cases <- unique(as.character(scenario_grid$case_id[grid_idx]))
    focal_status <- group_status[group_status$method == "gamlss.longitudinal", , drop = FALSE]
    for (neighbor in neighbors) {
      neighbor_status <- group_status[group_status$method == neighbor, , drop = FALSE]
      planned_cases <- expected_cases
      focal_attempted <- planned_cases %in% as.character(focal_status$case_id)
      neighbor_attempted <- planned_cases %in% as.character(neighbor_status$case_id)
      focal_success <- stats::setNames(mvt_status_success(focal_status$success), focal_status$case_id)
      neighbor_success <- stats::setNames(mvt_status_success(neighbor_status$success), neighbor_status$case_id)
      for (metric in metric_names) {
        focal <- group_benchmark[group_benchmark$method == "gamlss.longitudinal", c("case_id", metric), drop = FALSE]
        neighbor_data <- group_benchmark[group_benchmark$method == neighbor, c("case_id", metric), drop = FALSE]
        names(focal)[2L] <- "focal"
        names(neighbor_data)[2L] <- "neighbor"
        paired <- merge(focal, neighbor_data, by = "case_id", all = FALSE)
        paired$focal <- suppressWarnings(as.numeric(paired$focal))
        paired$neighbor <- suppressWarnings(as.numeric(paired$neighbor))
        paired$both_success <- focal_success[paired$case_id] %in% TRUE & neighbor_success[paired$case_id] %in% TRUE
        paired <- paired[paired$both_success & is.finite(paired$focal) & is.finite(paired$neighbor), , drop = FALSE]
        difference <- paired$neighbor - paired$focal
        interval <- mvt_mean_interval(difference)
        wins <- sum(difference > 0)
        win_interval <- mvt_wilson_interval(wins, length(difference))
        rows[[k]] <- cbind(
          case_groups[i, , drop = FALSE],
          data.frame(
            neighbor = neighbor,
            metric = metric,
            direction = "neighbor_minus_gamlss.longitudinal; positive favors gamlss.longitudinal for lower-is-better metrics",
            comparison_population = "conditional on both methods succeeding with finite metric values",
            failure_handling = paste0("conditional contrast excludes failed pairs; failure-inclusive sensitivity is reported in nearest_neighbor_results.csv using rule: ", mvt_phase2_penalty_rule()),
            estimand_qualification = if (metric == "elapsed_sec") "runtime contrast is conditional on completed fits and is not an unconditional compute-cost estimand" else if (metric == "benchmark_neg_log_score") "negative-log-score contrast is conditional on completed fitted distributions" else "prediction-RMSE contrast is conditional on completed fits",
            n_expected_pairs = length(expected_cases),
            n_focal_attempted = sum(focal_attempted),
            n_neighbor_attempted = sum(neighbor_attempted),
            n_attempted_pairs = sum(focal_attempted & neighbor_attempted),
            n_both_successful = length(difference),
            n_focal_failed = sum(!(focal_success[planned_cases] %in% TRUE)),
            n_neighbor_failed = sum(!(neighbor_success[planned_cases] %in% TRUE)),
            mean_paired_difference = interval[["mean"]],
            paired_difference_mcse = interval[["mcse"]],
            paired_difference_ci_lower = interval[["lower"]],
            paired_difference_ci_upper = interval[["upper"]],
            focal_win_rate = if (length(difference)) wins / length(difference) else NA_real_,
            focal_win_rate_ci_lower = win_interval[["lower"]],
            focal_win_rate_ci_upper = win_interval[["upper"]],
            stringsAsFactors = FALSE
          )
        )
        k <- k + 1L
      }
    }
  }
  mvt_bind_rows_fill(rows)
}

mvt_phase2_production_contract <- function() {
  list(
    time_names = "t20", n_time = 20L,
    families = c("gaussian", "poisson", "gamma", "binomial"),
    dependence = c(
      "external_exchangeable", "external_ar1",
      "native_time_varying_adjacent", "native_covariate_dependent_adjacent"
    ),
    reps = 1:100,
    methods = mvt_default_comparators(),
    cases = 1600L,
    method_rows = 14400L
  )
}

mvt_phase2_audit_registry <- function() {
  c(
    "capability_snapshot_dated",
    "capability_provenance_complete",
    "headline_neighbor_count_bounded",
    "headline_neighbors_match_registry",
    "scenario_grid_unique_nonempty",
    "production_time_grid_exact",
    "production_family_set_exact",
    "production_dependence_set_exact",
    "production_replicates_exact",
    "production_case_cardinality_exact",
    "production_grid_cross_product_exact",
    "production_method_registry_exact",
    "production_attempt_cardinality_exact",
    "production_registered_timeouts_finite",
    "production_subprocess_attestations_complete",
    "planned_methods_registered",
    "attempt_cross_product_reconciles",
    "expected_cases_reconcile",
    "planned_method_sets_exact",
    "case_method_sets_exact",
    "paired_denominators_reconcile",
    "headline_neighbor_rows_present",
    "gee_results_family_specific",
    "unstructured_high_dimension_is_stress_test",
    "attempt_denominators_visible",
    "failure_penalized_uncertainty_complete",
    "comparison_uncertainty_complete",
    "production_replication_target_met"
  )
}

mvt_phase2_production_eligible <- function(audit) {
  registry <- mvt_phase2_audit_registry()
  is.data.frame(audit) && all(c("check", "status", "detail") %in% names(audit)) &&
    nrow(audit) == length(registry) && !anyDuplicated(as.character(audit$check)) &&
    identical(as.character(audit$check), registry) &&
    !anyNA(audit$status) && identical(as.character(audit$status), rep("pass", length(registry)))
}

mvt_phase2_benchmark_audit <- function(
    capabilities, scope, gee, neighbors, contrasts, grid, status, benchmark,
    planned_methods, reconciliation, configuration = NULL, execution_attestation = NULL) {
  production <- mvt_phase2_production_contract()
  retained <- scope$method[scope$headline_empirical_comparator %in% TRUE]
  required_gee_methods <- c("gee_independence", "gee_exchangeable", "gee_ar1", "gee_unstructured")
  audit_row <- function(check, ok, detail) data.frame(check = check, status = if (isTRUE(ok)) "pass" else "fail", detail = detail, stringsAsFactors = FALSE)
  expected_attempts <- nrow(grid) * length(planned_methods)
  reconciliation_ok <- nrow(reconciliation) == expected_attempts &&
    expected_attempts > 0L && all(reconciliation$planned) && all(reconciliation$reconciled) &&
    !any(reconciliation$unexpected) && all(reconciliation$status_rows == 1L) && all(reconciliation$benchmark_rows == 1L)
  contrast_uncertainty_cols <- c("paired_difference_mcse", "paired_difference_ci_lower", "paired_difference_ci_upper", "focal_win_rate_ci_lower", "focal_win_rate_ci_upper")
  contrasts_nonempty <- nrow(contrasts) > 0L && all(contrasts$n_expected_pairs > 0L)
  contrasts_uncertain <- contrasts_nonempty && all(contrast_uncertainty_cols %in% names(contrasts)) &&
    all(contrasts$n_both_successful >= 2L) &&
    all(vapply(contrast_uncertainty_cols, function(col) all(is.finite(contrasts[[col]])), logical(1L)))
  penalized_cols <- c(
    "mean_rmse_failure_penalized", "mean_rmse_failure_penalized_mcse", "mean_rmse_failure_penalized_ci_lower", "mean_rmse_failure_penalized_ci_upper",
    "mean_neg_log_score_failure_penalized", "mean_neg_log_score_failure_penalized_mcse", "mean_neg_log_score_failure_penalized_ci_lower", "mean_neg_log_score_failure_penalized_ci_upper",
    "mean_elapsed_sec_failure_penalized", "mean_elapsed_sec_failure_penalized_mcse", "mean_elapsed_sec_failure_penalized_ci_lower", "mean_elapsed_sec_failure_penalized_ci_upper"
  )
  summary_rows <- mvt_bind_rows_fill(gee, neighbors)
  penalty_pool_cols <- c(
    "rmse_failure_penalty_pool_n", "neg_log_score_failure_penalty_pool_n",
    "elapsed_failure_penalty_pool_n", "failure_penalty_planned_methods"
  )
  expected_penalty_methods <- paste(sort(planned_methods), collapse = ",")
  penalized_uncertain <- nrow(summary_rows) > 0L && all(penalized_cols %in% names(summary_rows)) &&
    all(penalty_pool_cols %in% names(summary_rows)) &&
    all(vapply(penalized_cols, function(col) all(is.finite(summary_rows[[col]])), logical(1L))) &&
    all(summary_rows$rmse_failure_penalty_pool_n > 0L) &&
    all(summary_rows$neg_log_score_failure_penalty_pool_n > 0L) &&
    all(summary_rows$elapsed_failure_penalty_pool_n > 0L) &&
    nzchar(expected_penalty_methods) &&
    all(summary_rows$failure_penalty_planned_methods == expected_penalty_methods)
  grid_family <- if ("family_name" %in% names(grid)) grid$family_name else if ("family" %in% names(grid)) grid$family else rep("", nrow(grid))
  grid_for_audit <- grid
  grid_for_audit$family <- as.character(grid_family)
  gee_cell_cols <- intersect(
    c("scenario", "generator", "dependence", "correlation_level", "n_time", "n_subject", "total_rows", "family"),
    intersect(names(grid_for_audit), names(gee))
  )
  row_key <- function(data, cols) {
    if (!nrow(data) || !length(cols)) return(character())
    do.call(paste, c(lapply(data[cols], as.character), sep = "\r"))
  }
  expected_gee_cells <- if (length(gee_cell_cols)) unique(grid_for_audit[gee_cell_cols]) else data.frame()
  expected_gee <- if (nrow(expected_gee_cells)) merge(
    expected_gee_cells,
    data.frame(method = required_gee_methods, stringsAsFactors = FALSE),
    by = NULL
  ) else data.frame()
  gee_key_cols <- c(gee_cell_cols, "method")
  expected_gee_keys <- row_key(expected_gee, gee_key_cols)
  actual_gee_keys <- row_key(gee, gee_key_cols)
  expected_families <- sort(unique(as.character(grid_family[nzchar(grid_family)])))
  actual_families <- if ("family" %in% names(gee)) sort(unique(as.character(gee$family[nzchar(gee$family)]))) else character()
  gee_reconciles <- nrow(expected_gee) > 0L && nrow(gee) == nrow(expected_gee) &&
    !anyDuplicated(actual_gee_keys) && setequal(actual_gee_keys, expected_gee_keys) &&
    setequal(actual_families, expected_families)
  grid_n_time <- if ("n_time" %in% names(grid_for_audit)) suppressWarnings(as.numeric(grid_for_audit$n_time)) else rep(NA_real_, nrow(grid_for_audit))
  high_grid <- grid_for_audit[is.finite(grid_n_time) & grid_n_time >= 20L, , drop = FALSE]
  expected_stress <- if (nrow(high_grid) && length(gee_cell_cols)) unique(high_grid[gee_cell_cols]) else data.frame()
  gee_n_time <- if ("n_time" %in% names(gee)) suppressWarnings(as.numeric(gee$n_time)) else rep(NA_real_, nrow(gee))
  gee_method <- if ("method" %in% names(gee)) as.character(gee$method) else rep("", nrow(gee))
  actual_stress <- gee[gee_method == "gee_unstructured" & is.finite(gee_n_time) & gee_n_time >= 20L, , drop = FALSE]
  expected_stress_keys <- row_key(expected_stress, gee_cell_cols)
  actual_stress_keys <- row_key(actual_stress, gee_cell_cols)
  stress_reconciles <- nrow(expected_stress) > 0L && nrow(actual_stress) == nrow(expected_stress) &&
    !anyDuplicated(actual_stress_keys) && setequal(actual_stress_keys, expected_stress_keys) &&
    all(actual_stress$evidence_role == "high_dimensional_stress_test")
  rep_groups <- if (nrow(grid)) interaction(grid$scenario, grid_family, drop = TRUE) else factor()
  reps_per_cell <- if (length(rep_groups)) vapply(split(grid$rep, rep_groups), function(x) length(unique(x)), integer(1L)) else integer()
  minimum_reps <- if (length(reps_per_cell)) min(reps_per_cell) else 0L
  capability_validation <- mvt_capability_provenance_validation(capabilities)
  status_cases <- if (all(c("case_id", "method") %in% names(status))) unique(as.character(status$case_id)) else character()
  status_methods <- if (all(c("case_id", "method") %in% names(status))) unique(as.character(status$method)) else character()
  benchmark_cases <- if (all(c("case_id", "method") %in% names(benchmark))) unique(as.character(benchmark$case_id)) else character()
  benchmark_methods <- if (all(c("case_id", "method") %in% names(benchmark))) unique(as.character(benchmark$method)) else character()
  grid_cases <- if ("case_id" %in% names(grid)) unique(as.character(grid$case_id)) else character()
  paired_cols <- c("n_expected_pairs", "n_focal_attempted", "n_neighbor_attempted", "n_attempted_pairs", "n_both_successful")
  paired_denominators_ok <- contrasts_nonempty && all(paired_cols %in% names(contrasts)) &&
    all(contrasts$n_focal_attempted == contrasts$n_expected_pairs) &&
    all(contrasts$n_neighbor_attempted == contrasts$n_expected_pairs) &&
    all(contrasts$n_attempted_pairs == contrasts$n_expected_pairs) &&
    all(contrasts$n_both_successful <= contrasts$n_attempted_pairs)
  exact_time <- nrow(grid) == production$cases && "time_name" %in% names(grid) &&
    identical(sort(unique(as.character(grid$time_name))), production$time_names) &&
    all(suppressWarnings(as.integer(grid$n_time)) == production$n_time)
  exact_families <- setequal(unique(as.character(grid_family)), production$families) &&
    length(unique(as.character(grid_family))) == length(production$families)
  exact_dependence <- "dependence_name" %in% names(grid) &&
    setequal(unique(as.character(grid$dependence_name)), production$dependence) &&
    length(unique(as.character(grid$dependence_name))) == length(production$dependence)
  exact_reps <- setequal(unique(as.integer(grid$rep)), production$reps) &&
    all(reps_per_cell == length(production$reps))
  exact_methods <- identical(as.character(planned_methods), production$methods) &&
    setequal(status_methods, production$methods) && setequal(benchmark_methods, production$methods)
  exact_cardinality <- nrow(grid) == production$cases && nrow(status) == production$method_rows &&
    nrow(benchmark) == production$method_rows && nrow(reconciliation) == production$method_rows
  timeout_check <- if (is.list(configuration)) mvt_timeout_contract(configuration, planned_methods) else list(valid = FALSE, required = character(), values = numeric())
  subprocess_methods <- intersect(
    c("gee_independence", "gee_exchangeable", "gee_ar1", "gee_unstructured",
      "gamlss.longitudinal", "gamCopula_markov", "gamCopula_vine_simplified"),
    planned_methods
  )
  attestation_cols <- c(
    "subprocess_full_verified", "subprocess_dependency_fingerprint",
    "subprocess_runtime_fingerprint", "subprocess_configuration_fingerprint"
  )
  attested <- benchmark$method %in% subprocess_methods
  subprocess_ok <- length(subprocess_methods) > 0L && all(attestation_cols %in% names(benchmark)) &&
    sum(attested) == nrow(grid) * length(subprocess_methods) &&
    all(benchmark$subprocess_full_verified[attested] %in% TRUE) &&
    is.list(execution_attestation) &&
    all(benchmark$subprocess_dependency_fingerprint[attested] == execution_attestation$dependency_fingerprint) &&
    all(benchmark$subprocess_runtime_fingerprint[attested] == execution_attestation$runtime_fingerprint) &&
    all(benchmark$subprocess_configuration_fingerprint[attested] == execution_attestation$configuration_fingerprint)
  rows <- list(
    audit_row("capability_snapshot_dated", nrow(capabilities) > 0L && all(capabilities$retrieved_on == "2026-09-01"), paste("rows", nrow(capabilities))),
    audit_row(
      "capability_provenance_complete",
      capability_validation$valid,
      if (capability_validation$valid) paste("unique atomic claims", nrow(capabilities)) else paste(capability_validation$problems, collapse = "; ")
    ),
    audit_row("headline_neighbor_count_bounded", length(retained) >= 1L && length(retained) <= 2L, paste("headline comparators", length(retained))),
    audit_row("headline_neighbors_match_registry", setequal(retained, mvt_nearest_neighbor_comparators()), paste(retained, collapse = ",")),
    audit_row("scenario_grid_unique_nonempty", nrow(grid) > 0L && "case_id" %in% names(grid) && !anyDuplicated(grid$case_id), paste("grid cases", nrow(grid))),
    audit_row("production_time_grid_exact", exact_time, paste("required t20/T=20; actual", paste(unique(grid$time_name %||% "missing"), collapse = ","))),
    audit_row("production_family_set_exact", exact_families, paste("required", paste(production$families, collapse = ","), "actual", paste(sort(unique(grid_family)), collapse = ","))),
    audit_row("production_dependence_set_exact", exact_dependence, paste("required", paste(production$dependence, collapse = ","), "actual", paste(sort(unique(grid$dependence_name %||% "missing")), collapse = ","))),
    audit_row("production_replicates_exact", exact_reps, paste("required reps 1:100; actual unique", length(unique(grid$rep)))),
    audit_row("production_case_cardinality_exact", nrow(grid) == production$cases, paste("required", production$cases, "actual", nrow(grid))),
    audit_row("production_grid_cross_product_exact", mvt_grid_is_exact_production(grid), "canonical t20 x dependence x rep x family case IDs and fields"),
    audit_row("production_method_registry_exact", exact_methods, paste("required ordered", paste(production$methods, collapse = ","), "actual", paste(planned_methods, collapse = ","))),
    audit_row("production_attempt_cardinality_exact", exact_cardinality, paste("required", production$method_rows, "status", nrow(status), "benchmark", nrow(benchmark), "reconciliation", nrow(reconciliation))),
    audit_row("production_registered_timeouts_finite", timeout_check$valid, paste(timeout_check$required, timeout_check$values, sep = "=", collapse = ",")),
    audit_row("production_subprocess_attestations_complete", subprocess_ok, paste("required methods", paste(subprocess_methods, collapse = ","), "attested rows", sum(attested))),
    audit_row("planned_methods_registered", length(planned_methods) > 0L && all(c("gamlss.longitudinal", retained, required_gee_methods) %in% planned_methods), paste(planned_methods, collapse = ",")),
    audit_row("attempt_cross_product_reconciles", reconciliation_ok, paste("expected", expected_attempts, "reconciliation rows", nrow(reconciliation), "status rows", nrow(status), "benchmark rows", nrow(benchmark))),
    audit_row("expected_cases_reconcile", length(grid_cases) == nrow(grid) && setequal(status_cases, grid_cases) && setequal(benchmark_cases, grid_cases), paste("expected", length(grid_cases), "status", length(status_cases), "benchmark", length(benchmark_cases))),
    audit_row("planned_method_sets_exact", length(planned_methods) > 0L && setequal(status_methods, planned_methods) && setequal(benchmark_methods, planned_methods), paste("expected", paste(sort(planned_methods), collapse = ","), "status", paste(sort(status_methods), collapse = ","), "benchmark", paste(sort(benchmark_methods), collapse = ","))),
    audit_row("case_method_sets_exact", reconciliation_ok && length(status_cases) == nrow(grid) && setequal(status_methods, planned_methods), paste("cases", length(status_cases), "methods", paste(sort(status_methods), collapse = ","))),
    audit_row("paired_denominators_reconcile", paired_denominators_ok, paste("contrast rows", nrow(contrasts), "minimum attempted pairs", if (nrow(contrasts)) min(contrasts$n_attempted_pairs) else 0L)),
    audit_row("headline_neighbor_rows_present", nrow(neighbors) > 0L && all(c("gamlss.longitudinal", retained) %in% unique(neighbors$method)), paste("methods present", paste(sort(unique(neighbors$method)), collapse = ","))),
    audit_row("gee_results_family_specific", gee_reconciles, paste("expected rows", nrow(expected_gee), "actual rows", nrow(gee), "expected families", paste(expected_families, collapse = ","), "actual families", paste(actual_families, collapse = ","))),
    audit_row("unstructured_high_dimension_is_stress_test", stress_reconciles, paste("required T>=20 cells", nrow(expected_stress), "actual stress rows", nrow(actual_stress))),
    audit_row("attempt_denominators_visible", all(c("n_attempted", "n_converged", "n_included_accuracy") %in% names(gee)) && all(c("n_attempted", "n_converged") %in% names(neighbors)), paste("GEE rows", nrow(gee), "neighbor rows", nrow(neighbors))),
    audit_row("failure_penalized_uncertainty_complete", penalized_uncertain, paste("summary rows", nrow(summary_rows), "rule", mvt_phase2_penalty_rule())),
    audit_row("comparison_uncertainty_complete", contrasts_uncertain, paste("paired contrast rows", nrow(contrasts), "minimum both-successful", if (nrow(contrasts)) min(contrasts$n_both_successful) else 0L)),
    audit_row("production_replication_target_met", minimum_reps >= 100L, paste("minimum unique reps per scenario/family", minimum_reps, "required 100"))
  )
  out <- mvt_bind_rows_fill(rows)
  registry <- mvt_phase2_audit_registry()
  if (!identical(as.character(out$check), registry) || anyDuplicated(out$check)) {
    stop("Internal Phase 2 audit registry/order mismatch.", call. = FALSE)
  }
  out
}

mvt_phase2_audit_from_committed_snapshot <- function(committed) {
  if (!is.list(committed) || !isTRUE(committed$production_run) ||
      !is.data.frame(committed$grid) || !is.list(committed$rows)) {
    stop("Phase 2 production audit requires a reconciled exact production snapshot.", call. = FALSE)
  }
  grid <- committed$grid
  status <- committed$rows$fit_status
  benchmark <- committed$rows$benchmark_results
  planned_methods <- committed$methods
  reconciliation <- mvt_phase2_attempt_reconciliation(grid, status, benchmark, planned_methods)
  capabilities <- mvt_capability_snapshot()
  scope <- mvt_comparator_scope_registry()
  gee <- mvt_phase2_method_summary(
    status, benchmark,
    c("gee_independence", "gee_exchangeable", "gee_ar1", "gee_unstructured"),
    planned_methods = planned_methods
  )
  neighbors <- mvt_phase2_method_summary(
    status, benchmark,
    c("gamlss.longitudinal", mvt_nearest_neighbor_comparators()),
    planned_methods = planned_methods
  )
  contrasts <- mvt_nearest_neighbor_paired_contrasts(status, benchmark, scenario_grid = grid)
  mvt_phase2_benchmark_audit(
    capabilities, scope, gee, neighbors, contrasts,
    grid = grid, status = status, benchmark = benchmark,
    planned_methods = planned_methods, reconciliation = reconciliation,
    configuration = committed$configuration,
    execution_attestation = committed$execution_attestation
  )
}

mvt_phase2_expected_evidence <- function(committed) {
  grid <- committed$grid
  status <- committed$rows$fit_status
  benchmark <- committed$rows$benchmark_results
  methods <- committed$methods
  capabilities <- mvt_capability_snapshot()
  scope <- mvt_comparator_scope_registry()
  gee <- mvt_phase2_method_summary(
    status, benchmark,
    c("gee_independence", "gee_exchangeable", "gee_ar1", "gee_unstructured"),
    planned_methods = methods
  )
  stress <- gee[gee$method == "gee_unstructured" & gee$n_time >= 20L, , drop = FALSE]
  neighbors <- mvt_phase2_method_summary(
    status, benchmark,
    c("gamlss.longitudinal", mvt_nearest_neighbor_comparators()),
    planned_methods = methods
  )
  contrasts <- mvt_nearest_neighbor_paired_contrasts(status, benchmark, scenario_grid = grid)
  reconciliation <- mvt_phase2_attempt_reconciliation(grid, status, benchmark, methods)
  audit <- mvt_phase2_audit_from_committed_snapshot(committed)
  keyed <- lapply(
    list(gee = gee, stress = stress, neighbors = neighbors, contrasts = contrasts,
         reconciliation = reconciliation, audit = audit),
    mvt_add_phase2_evidence_keys
  )
  tables <- c(
    list(
      capability_snapshot_2026.09.01.csv = capabilities,
      comparator_scope_registry.csv = scope
    ),
    stats::setNames(
      keyed,
      c(
        "gee_family_results.csv", "gee_unstructured_stress_test.csv",
        "nearest_neighbor_results.csv", "nearest_neighbor_paired_contrasts.csv",
        "phase2_attempt_reconciliation.csv", "phase2_benchmark_audit.csv"
      )
    )
  )
  names(tables)[[1L]] <- "capability_snapshot_2026-09-01.csv"
  audit_md <- c(
    "# Phase 2 benchmark evidence audit", "",
    "Capability snapshot date: 2026-09-01",
    paste("Strict production eligibility:", if (mvt_phase2_production_eligible(audit)) "eligible" else "INELIGIBLE"),
    paste("Evidence readiness:", if (all(audit$status == "pass")) "complete" else "INCOMPLETE - do not use for publication claims"),
    "",
    "The only headline nearest-neighbor empirical workflows are `gamCopula_markov` and `gamCopula_vine_simplified`. Comparisons are paired by simulation case and remain conditional on both methods succeeding; failure-inclusive method summaries are reported separately.",
    "",
    "High-dimensional unstructured GEE rows (`T >= 20`) are stress-test feasibility evidence, not routine comparator evidence.",
    "Conditional paired contrasts exclude any pair with a failed or non-finite method. Failure-inclusive sensitivity summaries use the prespecified penalty rule recorded in the CSV outputs.",
    "",
    paste0("- ", audit$check, ": ", audit$status, " (", audit$detail, ")")
  )
  list(tables = tables, audit = audit, audit_md = audit_md)
}

mvt_validate_phase2_claim_evidence <- function(
    run_dir, attestation_path = mvt_phase2_snapshot_attestation_path(),
    signature_path = mvt_phase2_snapshot_signature_path()) {
  committed <- mvt_validate_approved_snapshot(
    run_dir, attestation_path = attestation_path, require_production = TRUE,
    signature_path = signature_path
  )
  expected <- mvt_phase2_expected_evidence(committed)
  if (!mvt_phase2_production_eligible(expected$audit)) {
    stop("Approved Phase 2 snapshot no longer passes the exact production audit registry.", call. = FALSE)
  }
  for (name in names(expected$tables)) {
    path <- file.path(run_dir, name)
    mvt_validate_aggregate_file_commit(path)
    if (!identical(mvt_csv_serialized_sha256(expected$tables[[name]]), mvt_sha256_file(path))) {
      stop("Phase 2 claim evidence does not reconcile with the approved snapshot: ", name, call. = FALSE)
    }
  }
  md_path <- file.path(run_dir, "phase2_benchmark_audit.md")
  mvt_validate_lines_file_commit(md_path)
  if (!identical(readLines(md_path, warn = FALSE), expected$audit_md)) {
    stop("Phase 2 audit narrative does not reconcile with the approved snapshot.", call. = FALSE)
  }
  identity <- list(
    snapshot_sha256 = attr(committed, "approved_snapshot_sha256"),
    producer_id = committed$producer_id,
    producer_version = committed$producer_version,
    configuration_fingerprint = committed$configuration_fingerprint,
    audit_registry_sha256 = mvt_hash_object(mvt_phase2_audit_registry()),
    evidence_sha256 = vapply(
      mvt_phase2_public_output_allowlist()$evidence_artifacts,
      function(name) mvt_sha256_file(file.path(run_dir, name)), character(1L)
    )
  )
  invisible(identity)
}

mvt_integrate_approved_phase2_snapshot <- function(
    run_dir, destination_dir,
    attestation_path = mvt_phase2_snapshot_attestation_path(),
    signature_path = mvt_phase2_snapshot_signature_path()) {
  identity <- mvt_validate_phase2_claim_evidence(run_dir, attestation_path, signature_path)
  committed <- mvt_validate_approved_snapshot(run_dir, attestation_path,
    require_production = TRUE, signature_path = signature_path)
  destination <- normalizePath(destination_dir, winslash = "/", mustWork = FALSE)
  if (file.exists(destination) || dir.exists(destination)) {
    stop("Approved Phase 2 integration destination must not already exist.", call. = FALSE)
  }
  parent <- dirname(destination)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  staging <- tempfile(paste0(".", basename(destination), "-staging-"), tmpdir = parent)
  dir.create(staging, recursive = FALSE, showWarnings = FALSE)
  installed <- FALSE
  on.exit(if (!installed && dir.exists(staging)) unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)
  snapshot_files <- c("aggregate_snapshot.rds", "aggregate_snapshot.rds.commit.rds", "aggregate_snapshot.rds.ownership.rds")
  aggregate_files <- unlist(lapply(committed$artifacts$file, function(name) {
    c(name, paste0(name, ".commit.rds"), paste0(name, ".ownership.rds"))
  }), use.names = FALSE)
  checkpoint_files <- file.path(
    "case_checkpoints",
    paste0(vapply(committed$tasks, `[[`, character(1L), "case_id"), ".rds")
  )
  evidence_files <- unlist(lapply(
    mvt_phase2_public_output_allowlist()$evidence_artifacts,
    function(name) c(name, paste0(name, ".commit.rds"), paste0(name, ".ownership.rds"))
  ), use.names = FALSE)
  relative <- unique(c(snapshot_files, aggregate_files, checkpoint_files, evidence_files))
  for (name in relative) {
    source <- file.path(run_dir, name)
    target <- file.path(staging, name)
    if (!file.exists(source)) stop("Approved Phase 2 integration source is incomplete: ", name, call. = FALSE)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(source, target, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE)) {
      stop("Could not stage approved Phase 2 artifact: ", name, call. = FALSE)
    }
  }
  staged_identity <- mvt_validate_phase2_claim_evidence(staging, attestation_path, signature_path)
  if (!identical(identity, staged_identity)) stop("Staged Phase 2 snapshot identity changed during integration.", call. = FALSE)
  if (!file.rename(staging, destination)) stop("Could not atomically install approved Phase 2 snapshot integration.", call. = FALSE)
  installed <- TRUE
  final_identity <- mvt_validate_phase2_claim_evidence(destination, attestation_path, signature_path)
  if (!identical(identity, final_identity)) stop("Installed Phase 2 snapshot identity failed post-copy validation.", call. = FALSE)
  invisible(c(final_identity, list(integration_dir = destination, files = relative)))
}

mvt_quarantine_phase2_publication_evidence <- function(run_dir, quarantine_dir, lease) {
  mvt_assert_active_lease(lease, file.path(run_dir, "publication-evidence"))
  dir.create(quarantine_dir, recursive = TRUE, showWarnings = FALSE)
  artifacts <- mvt_phase2_public_output_allowlist()$evidence_artifacts
  for (artifact in artifacts) {
    source <- file.path(run_dir, artifact)
    related <- c(source, paste0(source, ".commit.rds"), paste0(source, ".ownership.rds"))
    for (path in related[file.exists(related)]) {
      destination <- file.path(
        quarantine_dir,
        paste0(format(Sys.time(), "%Y%m%d%H%M%OS3"), "-", basename(path))
      )
      if (!file.rename(path, destination)) stop("Could not quarantine stale publication evidence: ", path, call. = FALSE)
    }
  }
  invisible(TRUE)
}

mvt_write_phase2_benchmark_evidence <- function(run_dir = mvt_read_run_dir()) {
  lease <- mvt_acquire_run_lock(run_dir)
  on.exit(mvt_release_run_lock(lease), add = TRUE)
  committed <- mvt_validate_approved_snapshot(run_dir)
  grid <- committed$grid
  status <- committed$rows$fit_status
  benchmark <- committed$rows$benchmark_results
  planned_methods <- committed$methods
  reconciliation <- mvt_phase2_attempt_reconciliation(grid, status, benchmark, planned_methods)
  capabilities <- mvt_capability_snapshot()
  scope <- mvt_comparator_scope_registry()
  gee_methods <- c("gee_independence", "gee_exchangeable", "gee_ar1", "gee_unstructured")
  gee <- mvt_phase2_method_summary(status, benchmark, gee_methods, planned_methods = planned_methods)
  stress <- gee[gee$method == "gee_unstructured" & gee$n_time >= 20L, , drop = FALSE]
  neighbors <- mvt_phase2_method_summary(
    status,
    benchmark,
    c("gamlss.longitudinal", mvt_nearest_neighbor_comparators()),
    planned_methods = planned_methods
  )
  contrasts <- mvt_nearest_neighbor_paired_contrasts(status, benchmark, scenario_grid = grid)
  audit <- mvt_phase2_audit_from_committed_snapshot(committed)
  if (!mvt_phase2_production_eligible(audit)) {
    quarantine_dir <- file.path(run_dir, "nonpublication_evidence")
    mvt_quarantine_phase2_publication_evidence(run_dir, quarantine_dir, lease)
    dir.create(quarantine_dir, recursive = TRUE, showWarnings = FALSE)
    mvt_write_csv_atomic(
      audit, file.path(quarantine_dir, "phase2_benchmark_audit.csv"), lease = lease
    )
    mvt_write_lines_atomic(
      c(
        "NONPUBLICATION: Phase 2 audit failed.",
        paste0(audit$check, "=", audit$status, ": ", audit$detail)
      ),
      file.path(quarantine_dir, "NONPUBLICATION.txt"), lease = lease
    )
    stop(
      "Phase 2 aggregate snapshot is ineligible; publication evidence was not emitted. See nonpublication_evidence/.",
      call. = FALSE
    )
  }
  gee <- mvt_add_phase2_evidence_keys(gee)
  stress <- mvt_add_phase2_evidence_keys(stress)
  neighbors <- mvt_add_phase2_evidence_keys(neighbors)
  contrasts <- mvt_add_phase2_evidence_keys(contrasts)
  reconciliation <- mvt_add_phase2_evidence_keys(reconciliation)
  audit <- mvt_add_phase2_evidence_keys(audit)
  paths <- c(
    capability = file.path(run_dir, "capability_snapshot_2026-09-01.csv"),
    scope = file.path(run_dir, "comparator_scope_registry.csv"),
    gee = file.path(run_dir, "gee_family_results.csv"),
    stress = file.path(run_dir, "gee_unstructured_stress_test.csv"),
    neighbors = file.path(run_dir, "nearest_neighbor_results.csv"),
    contrasts = file.path(run_dir, "nearest_neighbor_paired_contrasts.csv"),
    reconciliation = file.path(run_dir, "phase2_attempt_reconciliation.csv"),
    audit = file.path(run_dir, "phase2_benchmark_audit.csv")
  )
  mvt_write_csv_atomic(capabilities, paths[["capability"]], lease = lease)
  mvt_write_csv_atomic(scope, paths[["scope"]], lease = lease)
  mvt_write_csv_atomic(gee, paths[["gee"]], lease = lease)
  mvt_write_csv_atomic(stress, paths[["stress"]], lease = lease)
  mvt_write_csv_atomic(neighbors, paths[["neighbors"]], lease = lease)
  mvt_write_csv_atomic(contrasts, paths[["contrasts"]], lease = lease)
  mvt_write_csv_atomic(reconciliation, paths[["reconciliation"]], lease = lease)
  mvt_write_csv_atomic(audit, paths[["audit"]], lease = lease)
  md_path <- file.path(run_dir, "phase2_benchmark_audit.md")
  mvt_write_lines_atomic(
    c(
      "# Phase 2 benchmark evidence audit",
      "",
      "Capability snapshot date: 2026-09-01",
      paste("Strict production eligibility:", if (mvt_phase2_production_eligible(audit)) "eligible" else "INELIGIBLE"),
      paste("Evidence readiness:", if (all(audit$status == "pass")) "complete" else "INCOMPLETE - do not use for publication claims"),
      "",
      "The only headline nearest-neighbor empirical workflows are `gamCopula_markov` and `gamCopula_vine_simplified`. Comparisons are paired by simulation case and remain conditional on both methods succeeding; failure-inclusive method summaries are reported separately.",
      "",
      "High-dimensional unstructured GEE rows (`T >= 20`) are stress-test feasibility evidence, not routine comparator evidence.",
      "Conditional paired contrasts exclude any pair with a failed or non-finite method. Failure-inclusive sensitivity summaries use the prespecified penalty rule recorded in the CSV outputs.",
      "",
      paste0("- ", audit$check, ": ", audit$status, " (", audit$detail, ")")
    ),
    md_path, lease = lease
  )
  invisible(list(
    paths = c(paths, audit_md = md_path), audit = audit,
    production_eligible = mvt_phase2_production_eligible(audit)
  ))
}
