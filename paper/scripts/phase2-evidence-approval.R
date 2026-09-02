#!/usr/bin/env Rscript

# Explicit, checkout-external Phase 2 evidence key ceremony and signing helper.
# This script never prints or stores private-key bytes inside the checkout.

phase2_pinned_public_key <- function() {
  as.raw(c(66, 91, 71, 233, 21, 247, 172, 45, 215, 202, 170, 0, 64,
    43, 83, 206, 23, 50, 48, 154, 25, 217, 178, 37, 252, 59, 158, 195,
    237, 0, 31, 216))
}

phase2_refuse_checkout_path <- function(path, root) {
  root <- paste0(tolower(normalizePath(root, winslash = "/", mustWork = TRUE)), "/")
  path <- paste0(tolower(normalizePath(path, winslash = "/", mustWork = FALSE)), "/")
  if (startsWith(path, root)) stop("Approval keys and signatures must be outside the checkout.", call. = FALSE)
}

phase2_lock_private_key <- function(path) {
  Sys.chmod(path, mode = "0600")
  if (.Platform$OS.type == "windows" && nzchar(Sys.getenv("USERNAME"))) {
    status <- suppressWarnings(system2("icacls", c(shQuote(normalizePath(path, winslash = "\\", mustWork = TRUE)),
      "/inheritance:r", "/grant:r", shQuote(paste0(Sys.getenv("USERNAME"), ":(R,W)"))),
      stdout = FALSE, stderr = FALSE))
    if (!identical(status, 0L)) warning("Could not restrict the private-key ACL automatically; reviewer must restrict it manually.")
  }
  invisible(path)
}

phase2_write_signed_attestation <- function(attestation, private_key, public_key,
    attestation_path, signature_path) {
  message_raw <- serialize(attestation, NULL, version = 3L)
  signature_raw <- sodium::sig_sign(message_raw, private_key)
  if (!isTRUE(sodium::sig_verify(message_raw, signature_raw, public_key))) {
    stop("Detached signature self-check failed.", call. = FALSE)
  }
  writeBin(message_raw, attestation_path)
  writeBin(signature_raw, signature_path)
  invisible(c(
    attestation = normalizePath(attestation_path, winslash = "/", mustWork = TRUE),
    signature = normalizePath(signature_path, winslash = "/", mustWork = TRUE)
  ))
}

phase2_generate_key <- function(external_dir, root) {
  if (!requireNamespace("sodium", quietly = TRUE) || !requireNamespace("digest", quietly = TRUE)) {
    stop("sodium and digest are required.", call. = FALSE)
  }
  phase2_refuse_checkout_path(external_dir, root)
  dir.create(external_dir, recursive = TRUE, showWarnings = FALSE)
  private_path <- file.path(external_dir, "phase2-production-ed25519-private.key")
  public_path <- file.path(external_dir, "phase2-production-ed25519-public.key")
  if (file.exists(private_path) || file.exists(public_path)) stop("Refusing to overwrite an approval key.", call. = FALSE)
  private_key <- sodium::sig_keygen()
  public_key <- sodium::sig_pubkey(private_key)
  writeBin(private_key, private_path)
  phase2_lock_private_key(private_path)
  writeBin(public_key, public_path)
  Sys.chmod(public_path, mode = "0644")
  message("Public-key SHA-256: ", digest::digest(public_key, "sha256", serialize = FALSE))
  message("Private key written outside checkout: ", normalizePath(private_path, winslash = "/"))
  message("Do not sign until the public key has been independently reviewed and pinned in all five validators.")
  invisible(list(private_path = private_path, public_path = public_path))
}

phase2_sign <- function(study, bundle_dir, root, approver, private_path,
    attestation_path, signature_path) {
  if (!requireNamespace("sodium", quietly = TRUE) || !requireNamespace("digest", quietly = TRUE)) {
    stop("sodium and digest are required.", call. = FALSE)
  }
  if (!nzchar(approver)) stop("A named independent approver is required.", call. = FALSE)
  for (path in c(private_path, attestation_path, signature_path)) phase2_refuse_checkout_path(path, root)
  private_key <- readBin(private_path, "raw", n = file.info(private_path)$size)
  public_key <- sodium::sig_pubkey(private_key)
  if (!identical(public_key, phase2_pinned_public_key())) {
    stop("Private key does not correspond to the independently pinned production public key.", call. = FALSE)
  }
  if (identical(study, "main-recovery")) {
    source(file.path(root, "paper", "R", "main-recovery-evidence.R"), local = .GlobalEnv)
    jss_recovery_require_clean_checkout(jss_recovery_git_identity(root))
    original <- file.path(bundle_dir, jss_main_recovery_artifact_names())
    before <- unname(jss_recovery_sha256(original))
    source_before <- jss_recovery_control_source_identity(root)
    package_before <- jss_recovery_package_source_sha256(root)
    producer_before <- jss_recovery_producer_sha256(root)
    snapshot <- tempfile("main-recovery-approval-snapshot-"); dir.create(snapshot)
    on.exit(unlink(snapshot, recursive = TRUE, force = TRUE), add = TRUE)
    required <- file.path(snapshot, basename(original))
    if (!all(file.copy(original, required, overwrite = FALSE)) || !identical(before, unname(jss_recovery_sha256(original))) ||
        !identical(before, unname(jss_recovery_sha256(required)))) stop("Main-recovery bundle changed while creating the immutable approval snapshot.", call. = FALSE)
    jss_main_recovery_validate_public_bundle(snapshot, repo_root = root, require_attestation = FALSE)
    identity <- jss_recovery_bundle_identity(snapshot)
    approved_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    jss_recovery_validate_approval_order(approved_at, identity$execution_manifest$execution_completed_at_utc)
    attestation <- c(identity, list(approved_at_utc = approved_at, approver = approver))
    jss_main_recovery_validate_public_bundle(bundle_dir, repo_root = root, require_attestation = FALSE)
    jss_main_recovery_validate_public_bundle(snapshot, repo_root = root, require_attestation = FALSE)
    if (!identical(before, unname(jss_recovery_sha256(original))) || !identical(before, unname(jss_recovery_sha256(required))) ||
        !identical(source_before, jss_recovery_control_source_identity(root)) || !identical(package_before, jss_recovery_package_source_sha256(root)) ||
        !identical(producer_before, jss_recovery_producer_sha256(root)) || !identical(identity, jss_recovery_bundle_identity(snapshot))) {
      stop("Main-recovery bundle or source identity changed during independent approval/re-attestation.", call. = FALSE)
    }
    return(phase2_write_signed_attestation(attestation, private_key, public_key, attestation_path, signature_path))
  }
  if (identical(study, "copula-misspecification")) {
    source(file.path(root, "paper", "R", "07-gamma-copula-misspecification.R"), local = .GlobalEnv)
    results_path <- file.path(bundle_dir, "results.csv")
    source_paths <- jss_misspec_evidence_paths(results_path)
    if (any(!file.exists(source_paths))) stop("Module 07 bundle is incomplete.", call. = FALSE)
    before <- vapply(source_paths, jss_misspec_sha256_file, character(1L))
    snapshot_dir <- tempfile("copula-approval-snapshot-"); dir.create(snapshot_dir)
    on.exit(unlink(snapshot_dir, recursive = TRUE, force = TRUE), add = TRUE)
    snapshot_paths <- file.path(snapshot_dir, basename(source_paths)); names(snapshot_paths) <- names(source_paths)
    if (!all(file.copy(source_paths, snapshot_paths, overwrite = FALSE)) ||
        !identical(before, vapply(source_paths, jss_misspec_sha256_file, character(1L))) ||
        !identical(before, vapply(snapshot_paths, jss_misspec_sha256_file, character(1L)))) {
      stop("Module 07 bundle changed while creating the approval snapshot.", call. = FALSE)
    }
    snapshot <- snapshot_paths[["results"]]
    config <- jss_misspec_config(list(root = root, profile = "full", seed = 20260528L), stage = "full")
    identity <- jss_misspec_candidate_identity(snapshot, config)
    approved_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    attestation <- list(
      schema_version = 3L, study = study,
      results_sha256 = identity$results_sha256[[1L]], results_rows = identity$results_rows[[1L]],
      bundle_sha256 = identity$bundle_sha256[[1L]],
      execution_manifest_sha256 = identity$execution_manifest_sha256[[1L]],
      configuration_sha256 = identity$configuration_sha256[[1L]],
      producer_sha256 = identity$producer_sha256[[1L]],
      package_source_sha256 = identity$package_source_sha256[[1L]],
      approved_at_utc = approved_at, approver = approver
    )
    if (!identical(before, vapply(source_paths, jss_misspec_sha256_file, character(1L))) ||
        !identical(before, vapply(snapshot_paths, jss_misspec_sha256_file, character(1L))) ||
        !identical(identity, jss_misspec_candidate_identity(snapshot, config))) {
      stop("Module 07 bundle changed during independent approval.", call. = FALSE)
    }
    return(phase2_write_signed_attestation(attestation, private_key, public_key,
      attestation_path, signature_path))
  }
  if (identical(study, "multivariate-benchmark")) {
    source(file.path(root, "paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"), local = .GlobalEnv)
    snapshot <- mvt_validate_committed_snapshot(bundle_dir, require_production = TRUE)
    audit <- mvt_phase2_audit_from_committed_snapshot(snapshot)
    if (!mvt_phase2_production_eligible(audit)) stop("Module 09 candidate fails its exact production audit.", call. = FALSE)
    snapshot_sha <- mvt_phase2_snapshot_trust_sha256(bundle_dir)
    attestation <- list(
      schema_version = 2L, study = study, snapshot_sha256 = snapshot_sha,
      snapshot_schema_version = snapshot$schema_version,
      producer_id = snapshot$producer_id, producer_version = snapshot$producer_version,
      configuration_fingerprint = snapshot$configuration_fingerprint,
      audit_sha256 = mvt_hash_object(audit),
      artifact_manifest_sha256 = mvt_hash_object(snapshot$artifacts),
      checkpoint_manifest_sha256 = mvt_hash_object(snapshot$checkpoint_manifest),
      approved_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      approver = approver
    )
    after <- mvt_validate_committed_snapshot(bundle_dir, require_production = TRUE)
    if (!identical(snapshot_sha, mvt_phase2_snapshot_trust_sha256(bundle_dir)) ||
        !identical(mvt_hash_object(snapshot$artifacts), mvt_hash_object(after$artifacts)) ||
        !identical(mvt_hash_object(snapshot$checkpoint_manifest), mvt_hash_object(after$checkpoint_manifest))) {
      stop("Module 09 snapshot changed during independent approval.", call. = FALSE)
    }
    return(phase2_write_signed_attestation(attestation, private_key, public_key,
      attestation_path, signature_path))
  }
  if (identical(study, "missingness")) {
    source(file.path(root, "paper", "R", "missingness-study-helpers.R"), local = .GlobalEnv)
    source(file.path(root, "paper", "R", "04-missingness-dropout-sensitivity.R"), local = .GlobalEnv)
    original <- file.path(bundle_dir, jss_missingness_public_artifacts())
    before <- vapply(original, jss_missing_sha256_file, character(1))
    snapshot <- tempfile("missingness-approval-snapshot-"); dir.create(snapshot)
    on.exit(unlink(snapshot, recursive = TRUE, force = TRUE), add = TRUE)
    required <- file.path(snapshot, basename(original))
    if (!all(file.copy(original, required, overwrite = FALSE)) ||
        !identical(before, vapply(original, jss_missing_sha256_file, character(1))) ||
        !identical(unname(before), unname(vapply(required, jss_missing_sha256_file, character(1))))) {
      stop("Missingness bundle changed while creating the approval snapshot.", call. = FALSE)
    }
    jss_missingness_validate_candidate_bundle(snapshot, root = root, require_promotion = FALSE)
    checkpoints <- utils::read.csv(file.path(snapshot, "missingness_checkpoint_status.csv"), stringsAsFactors = FALSE)
    package_sha <- jss_missing_checkout_identity(root)$source_sha256
    producer_sha <- jss_missing_producer_sha256(c(file.path(root, "paper", "R", "missingness-study-helpers.R"),
      file.path(root, "paper", "scripts", "final-simulations", "missingness", "run_missingness_study.R")))
    file_sha <- jss_missing_sha256_file
  } else if (identical(study, "optimizer-benchmark")) {
    source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = .GlobalEnv)
    source(file.path(root, "paper", "R", "phase2-paper-evidence.R"), local = .GlobalEnv)
    original_dirs <- file.path(bundle_dir, c("data", "tables", "figures"))
    files <- jss_phase2_optimizer_files()
    original <- c(file.path(original_dirs[[1L]], files$data), file.path(original_dirs[[2L]], files$tables),
      file.path(original_dirs[[3L]], files$figures))
    before <- vapply(original, jss_joint_sha256_file, character(1))
    snapshot <- tempfile("optimizer-approval-snapshot-"); dir.create(snapshot)
    on.exit(unlink(snapshot, recursive = TRUE, force = TRUE), add = TRUE)
    snapshot_dirs <- file.path(snapshot, c("data", "tables", "figures"))
    invisible(lapply(snapshot_dirs, dir.create, recursive = TRUE))
    required <- c(file.path(snapshot_dirs[[1L]], files$data), file.path(snapshot_dirs[[2L]], files$tables),
      file.path(snapshot_dirs[[3L]], files$figures))
    if (!all(file.copy(original, required, overwrite = FALSE)) ||
        !identical(before, vapply(original, jss_joint_sha256_file, character(1))) ||
        !identical(unname(before), unname(vapply(required, jss_joint_sha256_file, character(1))))) {
      stop("Optimizer bundle changed while creating the approval snapshot.", call. = FALSE)
    }
    data_dir <- snapshot_dirs[[1L]]; tables_dir <- snapshot_dirs[[2L]]
    figures_dir <- snapshot_dirs[[3L]]
    jss_phase2_optimizer_validate_candidate(data_dir, tables_dir, figures_dir, root = root,
      require_promotion = FALSE)
    checkpoints <- utils::read.csv(file.path(data_dir, files$data[[3L]]), stringsAsFactors = FALSE)
    package_sha <- jss_joint_checkout_package_identity(list(root = root))$source_sha256
    producer_sha <- jss_joint_producer_fingerprint(list(root = root))
    file_sha <- jss_joint_sha256_file
  } else stop("study must be main-recovery, missingness, optimizer-benchmark, copula-misspecification, or multivariate-benchmark.", call. = FALSE)
  manifest <- paste(sort(paste(basename(required), vapply(required, file_sha, character(1)), sep = "\t")), collapse = "\n")
  approved_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  if (identical(study, "missingness")) {
    approved <- jss_missing_parse_rfc3339_utc(approved_at, "Missingness approval timestamp")
    executed <- jss_missing_parse_rfc3339_utc(checkpoints$timestamp_utc, "Missingness checkpoint timestamp")
    if (approved <= max(executed)) stop("Approval must occur after every checkpoint execution.", call. = FALSE)
  } else {
    jss_phase2_validate_approval_order(approved_at, checkpoints$checkpoint_timestamp_utc)
  }
  attestation <- list(schema_version = 1L, study = study,
    bundle_sha256 = digest::digest(manifest, "sha256", serialize = FALSE),
    package_source_sha256 = package_sha, producer_sha256 = producer_sha,
    approved_at_utc = approved_at,
    approver = approver, checkpoint_manifest = checkpoints)
  final_manifest <- paste(sort(paste(basename(required), vapply(required, file_sha, character(1)), sep = "\t")), collapse = "\n")
  if (!identical(final_manifest, manifest)) stop("Approval snapshot changed during review.", call. = FALSE)
  if (identical(study, "missingness")) {
    if (!identical(jss_missing_checkout_identity(root)$source_sha256, package_sha) ||
        !identical(jss_missing_producer_sha256(c(file.path(root, "paper", "R", "missingness-study-helpers.R"),
          file.path(root, "paper", "scripts", "final-simulations", "missingness", "run_missingness_study.R"))), producer_sha)) {
      stop("Missingness source identity changed during approval.", call. = FALSE)
    }
  } else if (!identical(jss_joint_checkout_package_identity(list(root = root))$source_sha256, package_sha) ||
      !identical(jss_joint_producer_fingerprint(list(root = root)), producer_sha)) {
    stop("Optimizer source identity changed during approval.", call. = FALSE)
  }
  phase2_write_signed_attestation(attestation, private_key, public_key,
    attestation_path, signature_path)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (length(args) == 2L && identical(args[[1L]], "generate-key")) {
    phase2_generate_key(args[[2L]], root)
  } else if (length(args) == 7L && identical(args[[1L]], "sign")) {
    phase2_sign(args[[2L]], args[[3L]], root, args[[4L]], args[[5L]], args[[6L]], args[[7L]])
  } else {
    stop("Usage: generate-key EXTERNAL_DIR; or sign STUDY BUNDLE_DIR APPROVER PRIVATE_KEY ATTESTATION SIGNATURE", call. = FALSE)
  }
}
