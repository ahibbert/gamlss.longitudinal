local_jss_repo_root <- function() {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  skip_if_not(file.exists(file.path(root, "paper", "manifest.csv")), "paper sources excluded")
  root
}

test_that("public profiles and isolated stores are declared", {
  root <- local_jss_repo_root()
  replicate <- readLines(file.path(root, "paper", "replicate.R"), warn = FALSE)
  expect_true(any(grepl('c\\("smoke", "paper", "full"\\)', replicate)))
  expect_true(any(grepl('file.path\\("paper", "_targets", profile\\)', replicate)))
  expect_true(any(grepl("expanded.*deprecated", replicate, ignore.case = TRUE)))
})

test_that("manifest classifies every active artifact without stubs", {
  root <- local_jss_repo_root()
  x <- utils::read.csv(file.path(root, "paper", "manifest.csv"), stringsAsFactors = FALSE)
  required <- c("manuscript_path", "manuscript_label", "producer", "profiles", "input_bundle", "access", "verification", "approved_sha256", "publication_status")
  expect_true(all(required %in% names(x)))
  expect_false(any(x$publication_status == "stub"))
  expect_true(all(x$access %in% c("public", "private", "non-data-static")))
  public <- x$access == "public" & x$publication_status == "active"
  expect_true(all(nzchar(x$producer[public])))
  expect_true(all(nzchar(x$generated_path[public])))
})

test_that("private modules are absent from the public target graph", {
  root <- local_jss_repo_root()
  graph <- paste(readLines(file.path(root, "paper", "_targets.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl("06-application-rand", graph, fixed = TRUE))
})

test_that("checkpoint resumption excludes completed fits", {
  root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "07-gamma-copula-misspecification.R"), local = env)
  settings <- list(profile = "full", seed = 20260528L)
  config <- env$jss_misspec_config(settings, stage = "smoke")
  grid <- env$jss_misspec_grid(config)
  checkpoint_dir <- tempfile("checkpoints-"); dir.create(checkpoint_dir)
  result <- cbind(
    grid[1L, , drop = FALSE],
    data.frame(
      stage = config$stage, correctly_specified = TRUE, attempted = TRUE,
      success = TRUE, converged = TRUE, retained = TRUE,
      stop_reason = "converged", failure_type = "none", error = NA_character_,
      warnings = NA_character_, elapsed_sec = 0.1, stringsAsFactors = FALSE
    )
  )
  env$jss_misspec_write_csv_atomic(result, env$jss_misspec_checkpoint_path(checkpoint_dir, grid$fit_id[[1]]))
  expect_equal(nrow(env$jss_misspec_pending_grid(grid, checkpoint_dir, config)), nrow(grid) - 1L)
})

test_that("copula misspecification checkpoints are atomic and worker-aware", {
  root <- local_jss_repo_root()
  code <- paste(readLines(file.path(root, "paper", "R", "07-gamma-copula-misspecification.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("file.rename(temporary, path)", code, fixed = TRUE))
  expect_true(grepl("parallel::parLapplyLB", code, fixed = TRUE))
  expect_true(grepl("workers = settings$workers", code, fixed = TRUE))
})

test_that("copula misspecification paper and full outputs pass a fatal binding gate", {
  root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  module_path <- file.path(root, "paper", "R", "07-gamma-copula-misspecification.R")
  public_path <- file.path(root, "paper", "R", "public-paper-producers.R")
  source(module_path, local = env)
  module_code <- paste(readLines(module_path, warn = FALSE), collapse = "\n")
  public_code <- paste(readLines(public_path, warn = FALSE), collapse = "\n")
  expect_true(grepl("jss_misspec_binding_review_gate(installed_results, grid, paths", module_code, fixed = TRUE))
  expect_true(grepl('jss_misspec_config(settings, stage = "full")', public_code, fixed = TRUE))
  expect_true(grepl("jss_misspec_validate_public_full_bundle(results, evidence_grid, evidence_config)", public_code, fixed = TRUE))
  expect_true(grepl("installed_results, grid = evidence_grid, paths = paths", public_code, fixed = TRUE))
  expect_false(grepl("grid = NULL", public_code, fixed = TRUE))

  output_dir <- tempfile("m07-binding-")
  dir.create(output_dir)
  paths <- list(
    results = file.path(output_dir, "results.csv"),
    heatmap = file.path(output_dir, "heatmap.png"),
    convergence = file.path(output_dir, "convergence.png")
  )
  settings <- list(profile = "full", seed = 20260528L)
  config <- env$jss_misspec_config(settings, stage = "smoke")
  grid <- env$jss_misspec_grid(config)[1L, , drop = FALSE]
  results <- cbind(
    grid,
    data.frame(
      stage = config$stage, correctly_specified = TRUE, attempted = TRUE,
      success = TRUE, converged = TRUE, retained = TRUE,
      stop_reason = "converged", failure_type = "none", error = NA_character_,
      warnings = NA_character_, elapsed_sec = 0.1, stringsAsFactors = FALSE
    )
  )
  env$jss_misspec_write_csv_atomic(results, paths$results)
  expect_false(env$jss_misspec_is_png(paths$results))
  review <- env$jss_misspec_binding_review_gate(results, grid, paths, "test-full", config)
  expect_true(all(review$status == "pass"))

  for (field in c("generating_copula", "fitted_copula", "rep", "seed", "dataset_seed")) {
    forged <- results
    forged[[field]] <- if (field %in% c("generating_copula", "fitted_copula")) "C" else forged[[field]] + 1L
    env$jss_misspec_write_csv_atomic(forged, paths$results)
    expect_error(
      env$jss_misspec_binding_review_gate(forged, grid, paths, paste0("forged-", field), config),
      "Fatal Module 07 binding review gate failed", info = field
    )
  }
  forged <- results; forged$target_tau <- 0.99
  env$jss_misspec_write_csv_atomic(forged, paths$results)
  expect_error(env$jss_misspec_binding_review_gate(forged, grid, paths, "forged-tau", config), "Fatal Module 07")
  forged <- results; forged$stage <- "full"
  env$jss_misspec_write_csv_atomic(forged, paths$results)
  expect_error(env$jss_misspec_binding_review_gate(forged, grid, paths, "forged-stage", config), "Fatal Module 07")
  forged <- results
  forged$failure_type <- "timeout"; forged$stop_reason <- "timeout"; forged$error <- "timed out"
  env$jss_misspec_write_csv_atomic(forged, paths$results)
  expect_error(env$jss_misspec_binding_review_gate(forged, grid, paths, "forged-success-timeout", config), "Fatal Module 07")
  env$jss_misspec_write_csv_atomic(results, paths$results)

  writeBin(as.raw(c(137, 80, 78, 71, 13, 10, 26, 10)), paths$results)
  expect_true(env$jss_misspec_is_png(paths$results))
  expect_error(
    env$jss_misspec_binding_review_gate(results, grid, paths, "test-mutated", config),
    "Fatal Module 07 binding review gate failed"
  )
  paths$results <- paths$heatmap
  expect_error(
    env$jss_misspec_binding_review_gate(results, grid, paths, "test-alias", config),
    "png_paths_are_distinct"
  )
})

test_that("Module 07 paper/full fits require finite killable subprocess timeouts", {
  root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "07-gamma-copula-misspecification.R"), local = env)
  settings <- list(profile = "full", seed = 20260528L, root = root)
  withr::local_envvar(GAMLSS_LONGITUDINAL_JSS_MISSPEC_TIMEOUT_SEC = "Inf")
  expect_error(env$jss_misspec_config(settings, stage = "full"), "finite positive")
  withr::local_envvar(GAMLSS_LONGITUDINAL_JSS_MISSPEC_TIMEOUT_SEC = "0.25")
  config <- env$jss_misspec_config(settings, stage = "full")
  start <- Sys.time()
  killed <- env$jss_misspec_killable_call(
    "jss_misspec_timeout_probe", list(seconds = 5), config
  )
  expect_s3_class(killed, "error")
  expect_lt(as.numeric(difftime(Sys.time(), start, units = "secs")), 5)
})

test_that("Module 07 full design is immutable R=100 and public input is full-only", {
  root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "07-gamma-copula-misspecification.R"), local = env)
  settings <- list(profile = "full", seed = 20260528L, root = root)

  withr::local_envvar(GAMLSS_LONGITUDINAL_JSS_MISSPEC_REPS = "99")
  expect_error(env$jss_misspec_config(settings, stage = "full"), "exactly 100")
  withr::local_envvar(GAMLSS_LONGITUDINAL_JSS_MISSPEC_REPS = "100")
  expect_identical(env$jss_misspec_config(settings, stage = "full")$reps, 100L)
  withr::local_envvar(GAMLSS_LONGITUDINAL_JSS_MISSPEC_REPS = "2")
  expect_error(env$jss_misspec_config(settings), "explicit smoke or pilot")
  expect_error(env$jss_run_07_gamma_copula_misspecification(settings), "explicit smoke or pilot")
  expect_identical(env$jss_misspec_config(settings, stage = "pilot")$reps, 2L)

  withr::local_envvar(GAMLSS_LONGITUDINAL_JSS_MISSPEC_REPS = NA_character_)
  config <- env$jss_misspec_config(settings, stage = "full")
  grid <- env$jss_misspec_grid(config)
  expect_equal(nrow(grid), 2400L)
  registered_routes <- gamlss.longitudinal::longitudinal_capabilities("routes")
  expect_identical(
    env$jss_misspec_copulas(),
    as.character(registered_routes$copula[registered_routes$margin_family == "GA"])
  )
  expect_identical(config$capability_registry_version, "2026.2")
  expect_identical(length(unique(grid$fit_id)), 2400L)
  expect_identical(length(unique(grid$seed)), 2400L)
  expect_identical(length(unique(grid$dataset_seed)), 1200L)
  expect_true(all(table(grid$dataset_seed) == 2L))
  expect_identical(sum(grid$generating_copula == grid$fitted_copula), 1200L)
  expect_identical(sum(grid$generating_copula != grid$fitted_copula), 1200L)

  paths <- env$jss_misspec_paths(list(
    data_dir = tempfile("m07-data-"), tables_dir = tempfile("m07-tables-"),
    figures_dir = tempfile("m07-figures-")
  ))
  expect_identical(basename(paths$checkpoints), "07-gamma-copula-misspecification-checkpoints-ga-nc-v2")

  unsupported_config <- config
  unsupported_config$copulas <- c(unsupported_config$copulas, "F")
  expect_error(
    env$jss_misspec_validate_config(unsupported_config),
    "diverges from the locked Gamma capability contract"
  )

  drift_root <- tempfile("m07-capability-drift-")
  dir.create(file.path(drift_root, "R"), recursive = TRUE)
  registry_lines <- readLines(file.path(root, "R", "capability-registry.R"), warn = FALSE)
  registry_lines <- sub(
    'compatible_copulas = c("N", "C")',
    'compatible_copulas = c("N", "C", "F")',
    registry_lines,
    fixed = TRUE
  )
  writeLines(registry_lines, file.path(drift_root, "R", "capability-registry.R"))
  expect_error(
    env$jss_misspec_capability_contract(drift_root),
    "no longer matches the locked Phase 1 capability registry"
  )
  results <- cbind(
    grid,
    data.frame(
      stage = "full", correctly_specified = grid$generating_copula == grid$fitted_copula,
      attempted = TRUE, success = TRUE, converged = TRUE, retained = TRUE,
      stop_reason = "converged", failure_type = "none", error = NA_character_,
      warnings = NA_character_, elapsed_sec = 0.1, stringsAsFactors = FALSE
    )
  )
  expect_true(env$jss_misspec_validate_public_full_bundle(results, grid, config))
  unsupported_result <- results
  unsupported_result$fitted_copula[[1L]] <- "F"
  expect_error(
    env$jss_misspec_validate_public_full_bundle(unsupported_result, grid, config),
    "outside the registered grid"
  )
  unsupported_checkpoint_dir <- tempfile("m07-unsupported-checkpoint-")
  dir.create(unsupported_checkpoint_dir)
  unsupported_checkpoint <- env$jss_misspec_checkpoint_path(
    unsupported_checkpoint_dir, grid$fit_id[[1L]]
  )
  utils::write.csv(unsupported_result[1L, , drop = FALSE], unsupported_checkpoint, row.names = FALSE)
  expect_identical(
    nrow(env$jss_misspec_pending_grid(grid[1L, , drop = FALSE], unsupported_checkpoint_dir, config)),
    1L
  )
  checkpoint_rejections <- utils::read.csv(
    file.path(unsupported_checkpoint_dir, "checkpoint-rejections-ledger.csv"),
    stringsAsFactors = FALSE
  )
  expect_match(checkpoint_rejections$reason, "outside the registered grid")
  expect_false(file.exists(unsupported_checkpoint))

  stale_checkpoint_dir <- tempfile("m07-stale-source-checkpoint-")
  dir.create(stale_checkpoint_dir)
  stale_checkpoint <- env$jss_misspec_checkpoint_path(stale_checkpoint_dir, grid$fit_id[[1L]])
  stale_result <- results[1L, , drop = FALSE]
  stale_result$configuration_sha256 <- paste(rep("b", 64L), collapse = "")
  utils::write.csv(stale_result, stale_checkpoint, row.names = FALSE)
  expect_identical(
    nrow(env$jss_misspec_pending_grid(grid[1L, , drop = FALSE], stale_checkpoint_dir, config)),
    1L
  )
  stale_rejections <- utils::read.csv(
    file.path(stale_checkpoint_dir, "checkpoint-rejections-ledger.csv"),
    stringsAsFactors = FALSE
  )
  expect_match(stale_rejections$reason, "fingerprint is stale")
  expect_false(file.exists(stale_checkpoint))
  mixed <- results; mixed$stage[[1L]] <- "pilot"
  expect_error(env$jss_misspec_validate_public_full_bundle(mixed, grid, config), "mixed-stage")
  expect_error(env$jss_misspec_validate_public_full_bundle(results[-1L, ], grid, config), "exact complete")

  approved_dir <- tempfile("m07-approved-full-"); dir.create(approved_dir)
  approved_results <- file.path(approved_dir, "results.csv")
  candidate_results <- transform(
    results,
    joint_loglik = ifelse(correctly_specified, -10, -11),
    copula_loglik = ifelse(correctly_specified, -2, -3),
    aic_joint = ifelse(correctly_specified, 24, 26),
    bic_joint = ifelse(correctly_specified, 30, 32),
    margin_param_rmse = ifelse(correctly_specified, 0.1, 0.2),
    tau_abs_error = ifelse(correctly_specified, 0.02, 0.04)
  )
  candidate_results <- env$jss_misspec_add_deltas(candidate_results)
  candidate_paths <- env$jss_misspec_evidence_paths(approved_results)
  env$jss_misspec_write_csv_atomic(candidate_results, candidate_paths[["results"]])
  env$jss_misspec_write_csv_atomic(env$jss_misspec_paired_effects(candidate_results), candidate_paths[["paired_effects"]])
  env$jss_misspec_write_csv_atomic(env$jss_misspec_warning_audit(candidate_results), candidate_paths[["warning_audit"]])
  env$jss_misspec_write_csv_atomic(env$jss_misspec_selection_attempts(candidate_results), candidate_paths[["selection_attempts"]])
  env$jss_misspec_write_csv_atomic(env$jss_misspec_selection(candidate_results), candidate_paths[["selection"]])
  current_sha <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)[[1L]]
  original_git_identity <- env$jss_misspec_git_identity
  env$jss_misspec_git_identity <- function(root) list(sha = current_sha, state = "clean")
  manifest_start <- list(
    started_at_utc = "2026-09-02T00:00:00Z", git_sha = current_sha, git_state = "clean",
    r_version = R.version.string, platform = R.version$platform,
    os = paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "),
    dependency_versions = env$jss_misspec_dependency_versions(), workers_requested = 1L
  )
  manifest <- env$jss_misspec_execution_manifest(
    manifest_start, config, grid, candidate_results,
    candidate_paths[c("results", "paired_effects", "warning_audit", "selection_attempts", "selection")]
  )
  env$jss_misspec_write_csv_atomic(manifest, candidate_paths[["execution_manifest"]])
  attestation <- tempfile("m07-attestation-", fileext = ".bin")
  signature <- tempfile("m07-signature-", fileext = ".sig")
  expect_error(
    env$jss_misspec_validate_approved_public_bundle(approved_results, grid, config, "", ""),
    "detached production approval signature"
  )
  identity <- env$jss_misspec_candidate_identity(approved_results, config)

  split_root <- tempfile("m07-split-layout-")
  dir.create(file.path(split_root, "data"), recursive = TRUE)
  dir.create(file.path(split_root, "tables"), recursive = TRUE)
  split_paths <- c(
    results = file.path(split_root, "data", "07-gamma-copula-misspecification-results.csv"),
    paired_effects = file.path(split_root, "tables", "07-gamma-copula-misspecification-paired-effects.csv"),
    warning_audit = file.path(split_root, "tables", "07-gamma-copula-misspecification-warning-audit.csv"),
    execution_manifest = file.path(split_root, "data", "07-gamma-copula-misspecification-execution-manifest.csv"),
    selection_attempts = file.path(split_root, "data", "07-gamma-copula-misspecification-selection-attempts.csv"),
    selection = file.path(split_root, "tables", "07-gamma-copula-misspecification-selection.csv")
  )
  expect_true(all(file.copy(candidate_paths, split_paths)))
  expect_silent(env$jss_misspec_validate_evidence_bundle(
    split_paths[["results"]], config, split_paths
  ))
  expect_error(env$jss_misspec_validate_evidence_bundle(
    split_paths[["results"]], config, split_paths[rev(names(split_paths))]
  ), "six-artifact contract")

  post_production_config <- config
  post_production_config$producer_sha256 <- paste(rep("a", 64L), collapse = "")
  post_production_config$package_source_sha256 <- paste(rep("b", 64L), collapse = "")
  post_production_config$configuration_sha256 <- env$jss_misspec_sha256_object(
    env$jss_misspec_registered_config_identity(post_production_config)
  )
  expect_silent(env$jss_misspec_validate_evidence_bundle(
    approved_results, post_production_config
  ))
  expect_identical(
    env$jss_misspec_candidate_identity(approved_results, post_production_config),
    identity
  )

  forged_manifest <- manifest
  forged_manifest$workers_requested <- 0L
  env$jss_misspec_write_csv_atomic(forged_manifest, candidate_paths[["execution_manifest"]])
  expect_error(env$jss_misspec_candidate_identity(approved_results, config), "execution manifest")
  env$jss_misspec_write_csv_atomic(manifest, candidate_paths[["execution_manifest"]])
  forged <- list(
    schema_version = 3L, study = "copula-misspecification",
    results_sha256 = identity$results_sha256[[1L]], results_rows = identity$results_rows[[1L]],
    bundle_sha256 = identity$bundle_sha256[[1L]],
    execution_manifest_sha256 = identity$execution_manifest_sha256[[1L]],
    configuration_sha256 = identity$configuration_sha256[[1L]],
    producer_sha256 = identity$producer_sha256[[1L]],
    package_source_sha256 = identity$package_source_sha256[[1L]],
    approved_at_utc = "2026-09-02T00:00:00Z", approver = "self-appointed"
  )
  message <- serialize(forged, NULL, version = 3L)
  fake_key <- sodium::sig_keygen()
  writeBin(message, attestation); writeBin(sodium::sig_sign(message, fake_key), signature)
  expect_error(
    env$jss_misspec_validate_approved_public_bundle(approved_results, grid, config, attestation, signature),
    "valid detached production signature"
  )
  original_verify <- env$jss_misspec_verify_signed_approval
  env$jss_misspec_verify_signed_approval <- function(identity, config, attestation_path, signature_path)
    list(approved_at_utc = "2026-09-02T00:00:00Z", approver = "test-only-mock")
  approval <- env$jss_misspec_validate_approved_public_bundle(approved_results, grid, config, attestation, signature)
  mutated <- candidate_results; mutated$warnings[[1L]] <- "post-validation mutation"
  env$jss_misspec_write_csv_atomic(mutated, approved_results)
  expect_error(env$jss_misspec_revalidate_approved_source(approved_results, approval, config), "changed after immutable validation")
  env$jss_misspec_verify_signed_approval <- original_verify
  env$jss_misspec_git_identity <- original_git_identity

  tracked_path <- file.path(
    root, "paper", "data", "public-derived", "copula-misspecification", "results.csv"
  )
  tracked <- utils::read.csv(tracked_path, stringsAsFactors = FALSE)
  tracked <- env$jss_misspec_upgrade_result_contract(tracked)
  expect_identical(nrow(tracked), 2400L)
  expect_identical(unique(as.character(tracked$stage)), "full")
  env$jss_misspec_git_identity <- function(...) list(sha = "test-clean-checkout", state = "clean")
  tracked_bundle <- expect_silent(env$jss_misspec_validate_evidence_bundle(tracked_path, config))
  expect_identical(nrow(tracked_bundle$results), 2400L)
  expect_identical(unique(as.character(tracked_bundle$results$stage)), "full")
  env$jss_misspec_git_identity <- original_git_identity

  selection_results <- transform(
    grid,
    success = TRUE, converged = TRUE,
    aic_joint = ifelse(generating_copula == fitted_copula, 10, 12),
    bic_joint = ifelse(generating_copula == fitted_copula, 11, 13),
    joint_loglik = ifelse(generating_copula == fitted_copula, -5, -6),
    copula_loglik = ifelse(generating_copula == fitted_copula, -1, -2),
    elapsed_sec = 1, mu_bias = 0, sigma_bias = 0,
    margin_param_rmse = ifelse(generating_copula == fitted_copula, 0.1, 0.2),
    tau_abs_error = ifelse(generating_copula == fitted_copula, 0.01, 0.02),
    benchmark_neg_log_score = 1, benchmark_pit_mean_abs_error = 0.1,
    failure_type = "none"
  )
  selection_attempts <- env$jss_misspec_selection_attempts(selection_results)
  fit_summary <- env$jss_misspec_summary(env$jss_misspec_add_deltas(selection_results))
  expect_identical(nrow(fit_summary), 24L)
  expect_identical(nrow(selection_attempts), 1200L)
  expect_identical(nrow(env$jss_misspec_selection(selection_results)), 12L)
  expect_identical(nrow(env$jss_misspec_selection_confusion(selection_attempts)), 36L)
  expect_identical(nrow(env$jss_misspec_selection_failures(selection_attempts)), 12L)
})

test_that("malformed Module 07 resume checkpoints are quarantined and rerun", {
  root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "07-gamma-copula-misspecification.R"), local = env)
  config <- env$jss_misspec_config(list(profile = "full", seed = 20260528L), stage = "smoke")
  grid <- env$jss_misspec_grid(config)
  checkpoint_dir <- tempfile("m07-invalid-resume-"); dir.create(checkpoint_dir)
  bad_path <- env$jss_misspec_checkpoint_path(checkpoint_dir, grid$fit_id[[1L]])
  utils::write.csv(data.frame(fit_id = grid$fit_id[[1L]], success = TRUE), bad_path, row.names = FALSE)
  pending <- env$jss_misspec_pending_grid(grid, checkpoint_dir, config)
  expect_equal(nrow(pending), nrow(grid))
  expect_false(file.exists(bad_path))
  ledger <- utils::read.csv(file.path(checkpoint_dir, "checkpoint-rejections-ledger.csv"), stringsAsFactors = FALSE)
  expect_equal(ledger$fit_id, grid$fit_id[[1L]])
  expect_match(ledger$reason, "missing fields")
  expect_true(file.exists(ledger$quarantined_path))
})

test_that("Module 07 PSOCK workers replace installed-package precedence with attested checkout source", {
  root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "07-gamma-copula-misspecification.R"), local = env)
  config <- env$jss_misspec_config(
    list(profile = "full", seed = 20260528L, root = root), stage = "smoke"
  )
  withr::local_envvar(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  cl <- parallel::makePSOCKcluster(2L, outfile = tempfile("m07-worker-attestation-"))
  on.exit(parallel::stopCluster(cl), add = TRUE)

  installed_paths <- unlist(parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(gamlss.longitudinal))
    normalizePath(
      getNamespaceInfo(asNamespace("gamlss.longitudinal"), "path"),
      winslash = "/", mustWork = TRUE
    )
  }), use.names = FALSE)
  expect_true(all(installed_paths != normalizePath(root, winslash = "/", mustWork = TRUE)))

  attestations <- env$jss_misspec_prepare_workers(cl, config, source_env = env)
  expected <- c(
    namespace_path = normalizePath(root, winslash = "/", mustWork = TRUE),
    package_source_sha256 = config$package_source_sha256
  )
  expect_length(attestations, 2L)
  expect_true(all(vapply(attestations, identical, logical(1L), y = expected)))
  post_paths <- unlist(parallel::clusterCall(cl, function() {
    normalizePath(
      getNamespaceInfo(asNamespace("gamlss.longitudinal"), "path"),
      winslash = "/", mustWork = TRUE
    )
  }), use.names = FALSE)
  expect_true(all(post_paths == expected[["namespace_path"]]))
})

test_that("copula selection confusion is rebuilt from complete attempt rows", {
  root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "07-gamma-copula-misspecification.R"), local = env)
  candidates <- env$jss_misspec_copulas()
  results <- do.call(rbind, lapply(1:2, function(rep_id) {
    data.frame(
      generating_copula = "N",
      fitted_copula = candidates,
      tau_label = "moderate",
      target_tau = 0.25,
      n_subject = 50L,
      rep = rep_id,
      success = TRUE,
      converged = TRUE,
      aic_joint = seq_along(candidates) + 9 + rep_id,
      bic_joint = seq_along(candidates) + 11 + rep_id,
      joint_loglik = -seq_along(candidates) - 4,
      failure_type = "none",
      stringsAsFactors = FALSE
    )
  }))
  failed_candidate <- tail(candidates, 1L)
  results$success[results$rep == 2L & results$fitted_copula == failed_candidate] <- FALSE
  results$converged[results$rep == 2L & results$fitted_copula == failed_candidate] <- FALSE
  results$failure_type[results$rep == 2L & results$fitted_copula == failed_candidate] <- "timeout"

  attempts <- env$jss_misspec_selection_attempts(results)
  confusion <- env$jss_misspec_selection_confusion(attempts)
  summary <- env$jss_misspec_selection(results)
  failures <- env$jss_misspec_selection_failures(attempts)

  expect_equal(nrow(attempts), 2L)
  expect_equal(sum(attempts$selection_eligible), 1L)
  expect_identical(attempts$failure_reasons[attempts$rep == 1L], "")
  expect_false(any(attempts$failure_reasons == "="))
  expect_equal(sum(confusion$n_selected), 2L)
  expect_equal(confusion$n_selected[confusion$selected_copula == "<no_selection>"], 1L)
  expect_equal(summary$aic_correct_selection_rate, 0.5)
  expect_true(is.finite(summary$aic_correct_selection_mcse))
  expect_true(is.finite(summary$close_aic_mcse))
  expect_true(is.finite(summary$close_aic_ci_lower))
  expect_true(is.finite(summary$median_aic_gap_bootstrap_se))
  expect_true(is.finite(summary$median_aic_gap_ci_lower))
  expect_equal(failures$n_no_selection, 1L)
  expect_match(failures$failure_reasons, paste0(failed_candidate, "=timeout"), fixed = TRUE)
})

test_that("copula selection names missing and duplicate candidate failures", {
  root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "07-gamma-copula-misspecification.R"), local = env)
  candidates <- env$jss_misspec_copulas()
  results <- data.frame(
    generating_copula = "N",
    fitted_copula = candidates,
    tau_label = "moderate",
    target_tau = 0.25,
    n_subject = 50L,
    rep = 1L,
    success = TRUE,
    converged = TRUE,
    aic_joint = seq_along(candidates),
    bic_joint = seq_along(candidates) + 1,
    joint_loglik = -seq_along(candidates),
    failure_type = "none",
    stringsAsFactors = FALSE
  )
  missing_candidate <- tail(candidates, 1L)
  results$fitted_copula[results$fitted_copula == missing_candidate] <- candidates[[1L]]

  attempts <- env$jss_misspec_selection_attempts(results)
  failures <- env$jss_misspec_selection_failures(attempts)

  expect_false(attempts$selection_eligible)
  expect_identical(attempts$missing_candidates, missing_candidate)
  expect_identical(attempts$duplicate_candidates, candidates[[1L]])
  expect_match(attempts$failure_reasons, paste0("missing_candidate=", missing_candidate), fixed = TRUE)
  expect_match(attempts$failure_reasons, paste0("duplicate_candidate=", candidates[[1L]], "(n=2)"), fixed = TRUE)
  expect_match(failures$failure_reasons, paste0("missing_candidate=", missing_candidate), fixed = TRUE)
  expect_match(failures$failure_reasons, paste0("duplicate_candidate=", candidates[[1L]], "(n=2)"), fixed = TRUE)
})

test_that("correlation benchmark checkpoints are replace-safe and resumable", {
  root <- local_jss_repo_root()
  scripts <- file.path(
    root, "paper", "R", "08-simulation-sensitivity-correlation-misspecification",
    "standard-model-benchmarking"
  )
  setup <- paste(readLines(file.path(scripts, "00-benchmark-setup.R"), warn = FALSE), collapse = "\n")
  runner <- paste(readLines(file.path(scripts, "02-run-rs-joint-standard-model-grid.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("file.rename(temporary, path)", setup, fixed = TRUE))
  expect_true(grepl("active_run_dir.txt", runner, fixed = TRUE))
  expect_true(grepl("GAMLSS_LONGITUDINAL_BENCHMARK_RESUME_RUN_DIR", runner, fixed = TRUE))
  expect_true(grepl("skipping completed", runner, fixed = TRUE))

  env <- new.env(parent = globalenv())
  source(file.path(scripts, "00-benchmark-setup.R"), local = env)
  output <- tempfile(fileext = ".csv")
  env$bmk_write_csv(data.frame(value = 1), output)
  env$bmk_write_csv(data.frame(value = 2), output)
  expect_equal(utils::read.csv(output)$value, 2)
  expect_false(file.exists(paste0(output, ".bak")))
  expect_length(Sys.glob(paste0(output, ".tmp-*")), 0L)
})

test_that("promoted Monte Carlo runners validate and atomically replace checkpoints", {
  root <- local_jss_repo_root()
  contains_call <- function(expr, name) {
    if (!is.call(expr) && !is.expression(expr) && !is.pairlist(expr)) return(FALSE)
    if (is.call(expr) && identical(as.character(expr[[1L]]), name)) return(TRUE)
    any(vapply(as.list(expr), contains_call, logical(1L), name = name))
  }
  protected_read <- function(expr) {
    if (!is.call(expr) && !is.expression(expr) && !is.pairlist(expr)) return(FALSE)
    if (is.call(expr) && identical(as.character(expr[[1L]]), "tryCatch") && contains_call(expr, "readRDS")) return(TRUE)
    any(vapply(as.list(expr), protected_read, logical(1L)))
  }
  scripts <- file.path(root, "paper", "scripts", "final-simulations", c(
    "missingness/run_missingness_study.R",
    "bcpe-t/simulation_bcpe_t_gamlss_comparison.R"
  ))
  for (script in scripts) {
    parsed <- parse(file = script, keep.source = FALSE)
    code <- paste(readLines(script, warn = FALSE), collapse = "\n")
    expect_true(protected_read(parsed))
    expect_true(grepl("checkpoint_result_issues", code, fixed = TRUE) || grepl("jss_missing_checkpoint_valid", code, fixed = TRUE))
    expect_true(grepl("saveRDS(result, temporary_path)", code, fixed = TRUE))
    expect_true(grepl("file.rename(temporary_path, final_path)", code, fixed = TRUE))
    expect_true(grepl("max_elapsed_sec", code, fixed = TRUE))
  }
})

test_that("joint-versus-separate paper inputs are per-replicate deltas", {
  root <- local_jss_repo_root()
  base <- file.path(root, "paper", "data", "public-derived", "joint-vs-separate")
  bundles <- c("normal", "gamma", "nbi")
  expect_true(all(file.exists(file.path(base, bundles, "data", "03-joint-vs-separate-optimization-deltas.csv"))))
  expect_false(any(file.exists(file.path(base, c(
    "normal-joint-vs-separate-six-case-median-iqr-table.tex",
    "gamma-joint-vs-separate-six-case-median-iqr-table.tex",
    "negative-binomial-joint-vs-separate-six-case-median-iqr-table.tex"
  )))))
})

test_that("publisher is allowlist-only and never commits or edits TeX", {
  root <- local_jss_repo_root()
  code <- paste(readLines(file.path(root, "paper", "publish-assets.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl('access == "public"', code, fixed = TRUE))
  expect_false(grepl("main.tex", code, fixed = TRUE))
  expect_false(grepl("git commit", code, fixed = TRUE))
  expect_false(grepl("git push", code, fixed = TRUE))
})

test_that("generated reviewer manifest is public-only and portable", {
  root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "replication-helpers.R"), local = env)
  out <- tempfile("reviewer-manifest-"); dir.create(out)
  settings <- list(root = root, out_dir = out, profile = "smoke")
  path <- env$jss_write_manifest(settings)
  manifest <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_true(all(manifest$access == "public"))
  expect_false(any(grepl("rand", manifest$artifact_id, ignore.case = TRUE)))
  expect_false(any(grepl("^([A-Za-z]:|/)", manifest$output_path)))
})

test_that("nonconverged smoke fits are structured events", {
  root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "replication-helpers.R"), local = env)
  logs <- tempfile("fit-events-"); dir.create(logs)
  results <- data.frame(
    family = "NO", copula = "N", design = "intercept", method = c("a", "b"),
    success = TRUE, converged = c(TRUE, FALSE), failure_type = "ok",
    stringsAsFactors = FALSE
  )
  path <- env$jss_write_fit_event_audit(results, list(logs_dir = logs))
  audit <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_equal(nrow(audit), 1L)
  expect_identical(audit$event_type, "optimizer_nonconvergence")
})

test_that("full-profile tolerance comparisons enforce registered bounds", {
  root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "replication-helpers.R"), local = env)
  tolerances <- utils::read.csv(file.path(root, "paper", "tolerances.csv"), stringsAsFactors = FALSE)
  reference <- data.frame(
    study = "fixture", key = "case=1", metric = "bias",
    tolerance_group = "fixed_effect_bias", statistic = "mean",
    value = 1, n = 100L, stringsAsFactors = FALSE
  )
  actual <- reference
  actual$value <- 1.01
  pass <- env$jss_compare_full_metrics(actual, reference, tolerances)
  expect_identical(pass$status, "pass")
  actual$value <- 1.2
  fail <- env$jss_compare_full_metrics(actual, reference, tolerances)
  expect_identical(fail$status, "outside_tolerance")
})

test_that("full joint-versus-separate bundles define reproducible cases and checkpoints", {
  root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "public-paper-producers.R"), local = env)
  settings <- list(public_data_dir = file.path(root, "paper", "data", "public-derived"))
  for (bundle in c("normal", "gamma", "nbi")) {
    cases <- env$jss_jvs_full_cases(settings, bundle)
    expect_identical(as.character(cases$case_id), sprintf("JVS%02d", 1:6))
    checkpoint <- tempfile(fileext = ".csv")
    rows <- data.frame(
      case_id = "JVS01", joint_review_rep = 1L,
      method = c("rs_separate", "rs_joint"), value = c(1, 2),
      audit_note = rep(paste(rep("checkpoint", 12), collapse = "-"), 2)
    )
    env$jss_jvs_full_write_checkpoint(rows, checkpoint)
    expect_true(env$jss_jvs_full_checkpoint_complete(checkpoint, "JVS01", 1L))
    expect_false(env$jss_jvs_full_checkpoint_complete(checkpoint, "JVS01", 2L))
  }
})
