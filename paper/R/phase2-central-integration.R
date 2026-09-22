# Central Phase 2 integration gates. Production modules retain their own
# numerical contracts; this file composes those contracts without weakening
# their trust or artifact boundaries.

jss_phase2_normalize_output_path <- function(path, settings) {
  root <- paste0(normalizePath(settings$out_dir, winslash = "/", mustWork = TRUE), "/")
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!startsWith(tolower(normalized), tolower(root))) {
    stop("Phase 2 module returned an output outside the replication directory: ", path, call. = FALSE)
  }
  substring(normalized, nchar(root) + 1L)
}

jss_phase2_artifact_identity <- function(paths) {
  paths <- sort(unique(as.character(paths)))
  if (!length(paths) || any(!file.exists(paths)) || any(file.info(paths)$isdir)) {
    stop("Phase 2 artifact identity requires an exact set of existing files.", call. = FALSE)
  }
  entries <- paste(basename(paths), vapply(paths, jss_file_sha256, character(1L)), sep = "\t")
  digest::digest(paste(entries, collapse = "\n"), algo = "sha256", serialize = FALSE)
}

jss_phase2_central_gate_registry <- function() {
  c(
    "all_phase2_modules_current",
    "all_module_production_gates_pass",
    "phase2_manifest_exact_bidirectional",
    "no_registered_stale_claims",
    "verified_claim_evidence_valid",
    "claim_schema_and_scope_registered",
    "design_contract_and_tex_assets_verified",
    "canonical_manuscript_artifacts_registered",
    "main_recovery_code_gate",
    "main_recovery_production_evidence_gate"
  )
}

jss_phase2_artifact_contract <- function(root = getwd()) {
  path <- file.path(root, "paper", "phase2-artifact-contract.csv")
  if (!file.exists(path)) stop("Canonical Phase 2 artifact contract is missing.", call. = FALSE)
  contract <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(
    "contract_version", "artifact_id", "artifact_type", "generated_path",
    "producer", "profiles", "input_bundle", "access", "publication_status", "verification"
  )
  if (!identical(names(contract), required) || nrow(contract) != 84L ||
      anyDuplicated(contract$artifact_id) || anyDuplicated(contract$generated_path) ||
      !all(contract$contract_version == "phase2-artifacts-2026-09-03.2-lean")) {
    stop("Versioned Phase 2 artifact contract is malformed.", call. = FALSE)
  }
  contract
}

jss_phase2_module_output_allowlist <- function(settings, modules, claim_evidence = NULL,
    gate_path = file.path(settings$tables_dir, "phase2-gate-audit.csv")) {
  required_modules <- c(
    "main_recovery", "optimizer", "missingness", "copula_selection", "fit_scaling"
  )
  if (!identical(names(modules), required_modules)) {
    stop("Phase 2 module registry is missing, reordered, or contains an unknown module.", call. = FALSE)
  }
  returned <- lapply(modules, function(module) {
    if (!is.list(module) || !all(c("module_id", "status", "data", "tables", "figures") %in% names(module))) {
      stop("Phase 2 module output contract is incomplete.", call. = FALSE)
    }
    c(module$data, module$tables, module$figures)
  })
  contract <- jss_phase2_artifact_contract(settings$root)
  producers <- c(
    main_recovery = "module_phase2_main_recovery",
    optimizer = "module_03_joint_vs_separate",
    missingness = "module_04_missingness",
    copula_selection = "module_07_copula_misspecification",
    fit_scaling = "module_10_fit_scaling"
  )
  for (name in names(returned)) {
    actual <- vapply(returned[[name]][nzchar(returned[[name]])],
      jss_phase2_normalize_output_path, character(1L), settings = settings)
    expected <- contract$generated_path[contract$producer == producers[[name]]]
    if (anyDuplicated(actual) || !setequal(actual, expected) || length(actual) != length(expected)) {
      stop("Phase 2 module output differs from the fixed artifact contract: ", name, call. = FALSE)
    }
  }
  expected_claim <- contract$generated_path[contract$producer == "phase2_claim_evidence"]
  expected_gate <- contract$generated_path[contract$producer == "phase2_gate_audit"]
  actual_claim <- if (is.null(claim_evidence)) character() else
    vapply(claim_evidence, jss_phase2_normalize_output_path, character(1L), settings = settings)
  actual_gate <- jss_phase2_normalize_output_path(gate_path, settings)
  if (!setequal(actual_claim, expected_claim) || length(actual_claim) != length(expected_claim) ||
      !identical(actual_gate, expected_gate)) {
    stop("Phase 2 claim/gate paths differ from the fixed artifact contract.", call. = FALSE)
  }
  contract$generated_path
}

jss_phase2_validate_manifest_allowlist <- function(settings, modules, claim_evidence = NULL,
    manifest_path = file.path(settings$root, "paper", "manifest.csv")) {
  expected <- jss_phase2_module_output_allowlist(settings, modules, claim_evidence)
  contract <- jss_phase2_artifact_contract(settings$root)
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- setdiff(names(contract), "contract_version")
  if (!all(required %in% names(manifest)) || anyDuplicated(manifest$artifact_id)) {
    stop("Phase 2 manifest schema or artifact IDs are invalid.", call. = FALSE)
  }
  active <- manifest$access == "public" & manifest$publication_status == "active" &
    grepl("^p2_", manifest$artifact_id)
  registered <- manifest[active, required, drop = FALSE]
  registered$generated_path <- gsub("\\\\", "/", registered$generated_path)
  canonical <- contract[required]
  rownames(registered) <- NULL; rownames(canonical) <- NULL
  if (!identical(registered, canonical) || !identical(registered$generated_path, expected)) {
    stop("Phase 2 manifest is missing, reordered, mutated, or contains reverse/unproduced rows relative to the fixed versioned contract.", call. = FALSE)
  }
  invisible(contract)
}

jss_phase2_validate_production_modules <- function(settings, modules) {
  expected <- c(
    main_recovery = "current", optimizer = "current", missingness = "regenerated",
    copula_selection = "regenerated", fit_scaling = "current"
  )
  if (!identical(names(modules), names(expected))) {
    stop("Phase 2 production module registry is incomplete or reordered.", call. = FALSE)
  }
  actual <- vapply(modules, `[[`, character(1L), "status")
  if (!identical(unname(actual), unname(expected))) {
    stop("Phase 2 production modules are not all publication-current: ", paste(names(actual), actual, sep = "=", collapse = "; "), call. = FALSE)
  }

  jss_submission_validate_main_recovery(file.path(settings$data_dir, "main-recovery"), settings$root)
  jss_submission_validate_optimizer(settings$data_dir, settings$tables_dir,
    settings$figures_dir, settings$root)
  jss_missingness_production_gate(file.path(settings$data_dir, "missingness"), settings$root)

  scaling <- list(
    design = utils::read.csv(file.path(settings$data_dir, "10-fit-scaling-scenarios.csv"), stringsAsFactors = FALSE),
    attempts = utils::read.csv(file.path(settings$data_dir, "10-fit-scaling-attempts.csv"), stringsAsFactors = FALSE),
    summary = utils::read.csv(file.path(settings$tables_dir, "10-fit-scaling-summary.csv"), stringsAsFactors = FALSE),
    contrasts = utils::read.csv(file.path(settings$tables_dir, "10-fit-scaling-contrasts.csv"), stringsAsFactors = FALSE),
    hardware = utils::read.csv(file.path(settings$tables_dir, "10-fit-scaling-hardware.csv"), stringsAsFactors = FALSE)
  )
  jss_scaling_assert_publication_eligible(
    scaling$design, scaling$attempts, scaling$summary, scaling$contrasts,
    jss_scaling_precision_rule("full"), hardware = scaling$hardware
  )
  jss_submission_validate_scaling(
    file.path(settings$data_dir, "10-fit-scaling-scenarios.csv"),
    file.path(settings$data_dir, "10-fit-scaling-attempts.csv"),
    file.path(settings$tables_dir, "10-fit-scaling-summary.csv"),
    file.path(settings$tables_dir, "10-fit-scaling-contrasts.csv"),
    file.path(settings$tables_dir, "10-fit-scaling-hardware.csv"), settings$root)
  jss_submission_validate_copula(
    file.path(settings$data_dir, "07-gamma-copula-misspecification-results.csv"), settings$root)
  invisible(TRUE)
}

jss_phase2_validate_claim_sources <- function(settings) {
  recovery <- jss_submission_validate_main_recovery(
    file.path(settings$data_dir, "main-recovery"), settings$root)
  optimizer <- jss_submission_validate_optimizer(settings$data_dir, settings$tables_dir,
    settings$figures_dir, settings$root)
  missingness_dir <- file.path(settings$data_dir, "missingness")
  jss_missingness_production_gate(missingness_dir, settings$root)
  missingness <- jss_submission_validate_missingness(missingness_dir, settings$root)

  module07 <- jss_submission_validate_copula(
    file.path(settings$public_data_dir, "copula-misspecification", "results.csv"), settings$root)
  staged07_path <- file.path(settings$data_dir, "07-gamma-copula-misspecification-results.csv")
  if (!file.exists(staged07_path)) stop("Staged Module 07 claim source is missing.", call. = FALSE)
  expected07 <- jss_misspec_add_deltas(module07)
  staged07 <- utils::read.csv(staged07_path, stringsAsFactors = FALSE, check.names = FALSE)
  jss_phase2_compare_frame(staged07, expected07, "fit_id", "Staged Module 07 attempt source")

  scaling_paths <- c(
    design = file.path(settings$data_dir, "10-fit-scaling-scenarios.csv"),
    attempts = file.path(settings$data_dir, "10-fit-scaling-attempts.csv"),
    summary = file.path(settings$tables_dir, "10-fit-scaling-summary.csv"),
    contrasts = file.path(settings$tables_dir, "10-fit-scaling-contrasts.csv"),
    hardware = file.path(settings$tables_dir, "10-fit-scaling-hardware.csv")
  )
  if (any(!file.exists(scaling_paths))) stop("Fit-scaling claim evidence is incomplete.", call. = FALSE)
  scaling <- lapply(scaling_paths, utils::read.csv, stringsAsFactors = FALSE)
  jss_scaling_assert_publication_eligible(
    scaling$design, scaling$attempts, scaling$summary, scaling$contrasts,
    jss_scaling_precision_rule("full"), hardware = scaling$hardware
  )
  scaling_identity <- jss_submission_validate_scaling(
    scaling_paths[["design"]], scaling_paths[["attempts"]], scaling_paths[["summary"]],
    scaling_paths[["contrasts"]], scaling_paths[["hardware"]], settings$root)
  invisible(list(
    `main-recovery` = recovery,
    `optimizer-benchmark` = optimizer,
    missingness = missingness,
    `copula-misspecification` = attr(module07, "submission_schema"),
    `fit-scaling` = scaling_identity
  ))
}

jss_phase2_validate_claim_evidence_rows <- function(claims, registered_claims, settings = NULL,
    expected = NULL) {
  if (!is.null(expected)) {
    stop(
      "Caller-supplied expected claim evidence is not a trusted input; exact source reconstruction is mandatory.",
      call. = FALSE
    )
  }
  if (is.null(settings) || !is.list(settings) || is.null(names(settings))) {
    stop(
      "Phase 2 claim validation requires a trusted settings context for current-source reconstruction; the two-argument form is disabled.",
      call. = FALSE
    )
  }
  expected <- jss_build_phase2_claim_evidence(settings)
  required <- c(
    "claim_id", "estimate", "mcse", "lower", "upper", "denominator",
    "attempted", "converged", "retained", "failed", "failure_reasons",
    "source_schema_version", "status"
  )
  if (!all(required %in% names(claims)) ||
      !all(c("claim_id", "expected_direction") %in% names(registered_claims)) ||
      anyDuplicated(claims$claim_id) || anyDuplicated(registered_claims$claim_id)) {
    stop("Phase 2 claim evidence/register schema is malformed.", call. = FALSE)
  }
  joined <- merge(
    registered_claims[c("claim_id", "expected_direction")], claims,
    by = "claim_id", all = TRUE, sort = FALSE
  )
  if (nrow(joined) != nrow(registered_claims) || anyNA(joined$expected_direction) || anyNA(joined$estimate)) {
    stop("Phase 2 claim evidence is missing, stale, duplicated, or unregistered.", call. = FALSE)
  }
  direction_supported <-
    joined$expected_direction == "no_direction" |
    (joined$expected_direction == "positive" & joined$lower > 0) |
    (joined$expected_direction == "negative" & joined$upper < 0)
  failure_reason_ok <- ifelse(joined$failed > 0,
    !is.na(joined$failure_reasons) & nzchar(trimws(joined$failure_reasons)) & joined$failure_reasons != "none",
    !is.na(joined$failure_reasons) & joined$failure_reasons == "none")
  count_fields <- c("denominator", "attempted", "converged", "retained", "failed")
  valid <- all(is.finite(joined$estimate)) && all(is.finite(joined$mcse) & joined$mcse >= 0) &&
    all(is.finite(joined$lower) & is.finite(joined$upper) &
      joined$lower <= joined$estimate & joined$estimate <= joined$upper) &&
    all(vapply(joined[count_fields], function(x) all(is.finite(x) & x >= 0 & x == floor(x)), logical(1L))) &&
    all(joined$denominator > 0 & joined$denominator == joined$attempted) &&
    all(joined$attempted == joined$retained + joined$failed) &&
    all(joined$converged >= joined$retained) &&
    all(failure_reason_ok) &&
    all(direction_supported) &&
    all(joined$status %in% c("verified_exact_design_count", "verified_exact_evidence_count", "verified_directional_effect")) &&
    all(joined$source_schema_version == jss_submission_schema_version())
  if (!valid) stop("Phase 2 claim evidence is stale, internally contradictory, or lacks interval-supported provenance.", call. = FALSE)
  if (!identical(names(claims), names(expected)) || nrow(claims) != nrow(expected)) {
    stop("Phase 2 claim evidence schema/count differs from exact source reconstruction.", call. = FALSE)
  }
  actual <- claims[match(expected$claim_id, claims$claim_id), names(expected), drop = FALSE]
  rownames(actual) <- NULL; rownames(expected) <- NULL
  for (name in names(expected)) {
    same <- if (is.numeric(expected[[name]])) {
      isTRUE(all.equal(as.numeric(actual[[name]]), as.numeric(expected[[name]]), tolerance = 0, check.attributes = FALSE))
    } else identical(as.character(actual[[name]]), as.character(expected[[name]]))
    if (!same) {
      stop("Phase 2 claim evidence disagrees with current source in field: ", name, call. = FALSE)
    }
  }
  invisible(joined)
}
