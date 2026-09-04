local_main_recovery_root <- function() {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  skip_if_not(file.exists(file.path(root, "paper", "R", "main-recovery-evidence.R")), "paper sources excluded")
  root
}

local_main_recovery_env <- function() {
  root <- local_main_recovery_root()
  env <- new.env(parent = baseenv())
  sys.source(file.path(root, "paper", "R", "main-recovery-evidence.R"), envir = env)
  env
}

local_repin_main_recovery_bundle <- function(env, path) {
  input <- utils::read.csv(file.path(path, "input_provenance.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  for (i in seq_len(nrow(input))) {
    raw <- file.path(path, input$path[[i]])
    input$sha256[[i]] <- unname(env$jss_recovery_sha256(raw)); input$bytes[[i]] <- file.info(raw)$size
  }
  utils::write.csv(input, file.path(path, "input_provenance.csv"), row.names = FALSE, na = "")
  manifest <- utils::read.csv(file.path(path, "output_manifest.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  for (i in seq_len(nrow(manifest))) {
    raw <- file.path(path, manifest$artifact[[i]])
    manifest$sha256[[i]] <- unname(env$jss_recovery_sha256(raw)); manifest$bytes[[i]] <- file.info(raw)$size
  }
  utils::write.csv(manifest, file.path(path, "output_manifest.csv"), row.names = FALSE, na = "")
  checkpoint <- utils::read.csv(file.path(path, "bundle_checkpoint.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  checkpoint$input_provenance_sha256 <- digest::digest(paste(input$input_id, input$sha256, sep = "=", collapse = "\n"), algo = "sha256", serialize = FALSE)
  checkpoint$output_manifest_sha256 <- unname(env$jss_recovery_sha256(file.path(path, "output_manifest.csv")))
  utils::write.csv(checkpoint, file.path(path, "bundle_checkpoint.csv"), row.names = FALSE, na = "")
  invisible(path)
}

local_main_recovery_fixture <- function(env) {
  root <- local_main_recovery_root()
  actual_git_identity <- env$jss_recovery_git_identity(root)
  env$jss_recovery_git_identity <- local({ identity <- actual_git_identity; function(repo_root) {
    identity[["git_state"]] <- "clean"; identity
  }})
  env$jss_recovery_source_git_state <- function(repo_root) "clean"
  bcpe_fields <- env$jss_recovery_bcpe_settings_fields()
  bcpe_settings <- as.data.frame(setNames(as.list(rep("fixture_registered_value", length(bcpe_fields))), bcpe_fields), stringsAsFactors = FALSE)
  bcpe_settings$runner_contract_version <- "bcpe-t-main-recovery-2026-09-02.7"
  bcpe_settings$phase1_contract_version <- env$jss_recovery_phase1_contract_version()
  bcpe_settings$n_fits <- "100"; bcpe_settings$rep_ids <- paste(1:100, collapse = ",")
  bcpe_settings$max_attempts_per_fit <- "3"; bcpe_settings$n_subject <- "500"; bcpe_settings$n_time <- "4"
  nbi_settings <- data.frame(runner_contract_version = "nbi-clayton-main-recovery-2026-09-02.5",
    phase1_contract_version = env$jss_recovery_phase1_contract_version(), n_subject = 500L, times = "0|1|2|3",
    reps = 100L, base_seed = 20260401L, sigma_signal_multiplier = 2, engines = "gamlss|ours_rs_joint",
    max_elapsed_sec = 180, max_outer_iter = 1000L, max_inner_iter = 100L, start_step_size = 0.5,
    step_adjustment_env = NA_character_, lambda_start = 1, rs_update_lambda = TRUE, warm_start_joint_iter = 5L,
    compute_se = TRUE, vcov_method_longitudinal = "analytical", compute_predictive_scores = TRUE,
    predictive_nsim = 100L, variogram_p_values = "0.5|2", theta_intercept = 0.7, theta_time_coef = 0.2,
    max_attempts_per_fit = 3L, runtime_n_cores = 1L, runtime_backend = "sequential",
    rscript_path = normalizePath(file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"), winslash = "/"),
    rscript_sha256 = unname(env$jss_recovery_sha256(file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"))),
    rscript_version = R.version.string)
  bcpe_signature <- env$jss_recovery_settings_signature(bcpe_settings, "BCPE")
  nbi_signature <- env$jss_recovery_settings_signature(nbi_settings, "NBI")
  git_identity <- env$jss_recovery_git_identity(root)
  bcpe_settings$git_sha <- nbi_settings$git_sha <- unname(git_identity[["git_sha"]])
  bcpe_settings$git_state <- nbi_settings$git_state <- unname(git_identity[["git_state"]])
  package_sha <- env$jss_recovery_package_source_sha256(root)
  bcpe_runner_sha <- unname(env$jss_recovery_sha256(file.path(root, "paper", "scripts", "final-simulations", "bcpe-t", "simulation_bcpe_t_gamlss_comparison.R")))
  nbi_runner_sha <- unname(env$jss_recovery_sha256(file.path(root, "paper", "scripts", "final-simulations", "nbi-clayton", "compare_gamlss_ours_nbi_sigma_smooth.R")))
  bcpe <- data.frame(
    evidence_status = "post_phase1_production", study_id = "bcpe_t_main_recovery",
    margin_family = "BCPE", copula_family = "t", copula_code = "t",
    runner_contract_version = "bcpe-t-main-recovery-2026-09-02.7", phase1_contract_version = env$jss_recovery_phase1_contract_version(), runner_settings_signature = bcpe_signature,
    runner_settings_sha256 = digest::digest(bcpe_signature, algo = "sha256", serialize = FALSE), runner_sha256 = bcpe_runner_sha, package_source_path = ".", package_version = "0.1.0", package_source_sha256 = package_sha,
    scenario = "n500_d4", n = 500L, d = 4L,
    target_replicates = 100L, rep = 1:100, model = "gamlss.longitudinal",
    success = TRUE, converged = TRUE, elapsed_sec = seq_len(100), error = NA_character_, seed = 110001:110100,
    execution_completed_at_utc = sprintf("2026-09-01T00:%02d:%02dZ", ((0:99) %/% 60), (0:99) %% 60)
  )
  bcpe <- rbind(bcpe, transform(bcpe, model = "gamlss2", elapsed_sec = elapsed_sec + 10))
  nbi <- data.frame(
    evidence_status = "post_phase1_production", study_id = "nbi_clayton_main_recovery",
    margin_family = "NBI", copula_family = "Clayton", copula_code = "C",
    runner_contract_version = "nbi-clayton-main-recovery-2026-09-02.5", phase1_contract_version = env$jss_recovery_phase1_contract_version(), runner_settings_signature = nbi_signature,
    runner_settings_sha256 = digest::digest(nbi_signature, algo = "sha256", serialize = FALSE), runner_sha256 = nbi_runner_sha, package_source_path = ".", package_version = "0.1.0", package_source_sha256 = package_sha,
    scenario = "n500_d4_nbi_signal2", n = 500L, d = 4L,
    target_replicates = 100L, rep = 1:100, engine = "ours_rs_joint",
    success = TRUE, converged = TRUE, elapsed_sec = 101:200, error = NA_character_, seed = 210001:210100,
    execution_completed_at_utc = sprintf("2026-09-01T01:%02d:%02dZ", ((0:99) %/% 60), (0:99) %% 60)
  )
  nbi <- rbind(nbi, transform(nbi, engine = "gamlss", elapsed_sec = elapsed_sec + 10))
  add_raw_convergence <- function(x, margin, method_column) {
    method <- as.character(x[[method_column]])
    contract <- lapply(method, function(value) env$jss_recovery_convergence_method_contract(margin, value))
    x$raw_convergence_schema <- "raw-convergence-2026-09-01.1"
    x$raw_convergence_api <- vapply(contract, `[[`, character(1), "api")
    x$raw_convergence_basis <- vapply(contract, `[[`, character(1), "basis")
    x$raw_convergence_status <- ifelse(method == "gamlss2", "finite_fit_below_outer_iteration_cap", "explicit_optimizer_convergence")
    x$raw_convergence_indicator_name <- vapply(contract, `[[`, character(1), "indicator")
    x$raw_convergence_indicator_value <- ifelse(method == "gamlss2", NA, TRUE)
    x$raw_convergence_loglik <- -1; x$raw_convergence_deviance <- 2
    x$raw_convergence_coefficient_count <- 4L; x$raw_convergence_coefficient_nonfinite_count <- 0L
    x$raw_convergence_fitted_count <- 4L; x$raw_convergence_fitted_nonfinite_count <- 0L
    x$raw_convergence_iteration_count <- 1L; x$raw_convergence_iteration_cap <- 20L
    x$raw_convergence_hit_outer_limit <- FALSE; x$raw_convergence_hit_max_stall <- FALSE
    x$raw_convergence_hit_raw_loglik_deterioration <- FALSE
    x
  }
  bcpe <- add_raw_convergence(bcpe, "BCPE", "model")
  nbi <- add_raw_convergence(nbi, "NBI", "engine")
  ledger <- env$jss_recovery_attempt_ledger(bcpe, nbi)
  outputs <- unique(ledger[c("study_id", "scenario_id", "method", "replicate", "margin_family")])
  pieces <- lapply(seq_len(nrow(outputs)), function(i) {
    identity <- outputs[i, c("study_id", "scenario_id", "method", "replicate")]
    contract <- env$jss_recovery_metric_contract(outputs$margin_family[i], outputs$method[i])
    fixed_keys <- strsplit(contract$fixed, "\r", fixed = TRUE)
    list(
      fixed = data.frame(identity[rep(1, length(fixed_keys)), ], parameter = vapply(fixed_keys, `[`, character(1), 1), term = vapply(fixed_keys, `[`, character(1), 2), estimate = 1, std_error = 0.1, true_value = 1),
      smooth = data.frame(identity[rep(1, length(contract$smooth)), ], parameter = contract$smooth, bias_abs_integrated = 0.1, irmse = 0.2),
      predictive = data.frame(identity[rep(1, length(contract$variogram_p)), ], test_log_score_joint = -2000,
        test_log_score_marginal = -1800, test_log_score_copula = -200, test_log_score_per_obs = -1,
        predictive_nsim = 100L, variogram_score = 0.5, variogram_p = contract$variogram_p),
      diagnostic = transform(identity, logLik = -1, rosenblatt_ks = 0.1, rosenblatt_cvm = 0.2,
        abs_rosenblatt_lag1_cor = 0.05, abs_rosenblatt_normal_lag1_cor = 0.04,
        rosenblatt_mean_abs_time_mean = 0.03, rosenblatt_normal_mean_abs_time_mean = 0.02)
    )
  })
  fixed <- do.call(rbind, lapply(pieces, `[[`, "fixed"))
  smooth <- do.call(rbind, lapply(pieces, `[[`, "smooth"))
  predictive <- do.call(rbind, lapply(pieces, `[[`, "predictive"))
  diagnostic <- do.call(rbind, lapply(pieces, `[[`, "diagnostic"))
  ledger <- env$jss_recovery_reconcile_retention(ledger, fixed, smooth, predictive, diagnostic)
  list(ledger = ledger, fixed = fixed, smooth = smooth, predictive = predictive, diagnostic = diagnostic,
    bcpe_settings = bcpe_settings, nbi_settings = nbi_settings)
}

test_that("attempt ledger keeps failures and unknown convergence visible", {
  env <- local_main_recovery_env()
  bcpe <- data.frame(
    scenario = "n500_d4", n = 500L, d = 4L, rep = 1:3, model = "gamlss.longitudinal",
    success = c(TRUE, TRUE, FALSE), converged = c(TRUE, FALSE, NA), elapsed_sec = c(1, 2, 3),
    stop_reason = c("tolerance", "outer_limit", NA), error = c(NA, NA, "fit error")
  )
  nbi <- data.frame(
    rep = 1:3, engine = "ours_rs_joint", success = TRUE, converged = c(TRUE, NA, TRUE),
    elapsed_sec = c(4, 5, 6), error = NA_character_
  )

  ledger <- env$jss_recovery_attempt_ledger(bcpe, nbi)
  counts <- env$jss_recovery_counts(ledger)
  failures <- env$jss_recovery_failure_summary(ledger)

  bcpe_id <- unique(ledger$study_id[ledger$margin_family == "BCPE"])
  nbi_id <- unique(ledger$study_id[ledger$margin_family == "NBI"])
  bcpe_count <- counts[counts$study_id == bcpe_id, ]
  expect_equal(unname(unlist(bcpe_count[c("attempted", "successful", "converged", "retained")])), c(3, 2, 1, 1))
  expect_true(any(grepl("optimizer_nonconvergence: outer_limit", failures$failure_reason, fixed = TRUE)))
  expect_true(any(grepl("execution_error: fit error", failures$failure_reason, fixed = TRUE)))
  expect_identical(ledger$convergence_status[ledger$study_id == nbi_id & ledger$replicate == 2L], "not_reported")
})

test_that("design table is generated only from attempt metadata", {
  env <- local_main_recovery_env()
  bcpe <- data.frame(
    scenario = "metadata_scenario", n = 500L, d = 4L, rep = 1:100,
    model = "method_a", success = TRUE, converged = TRUE, elapsed_sec = 1, error = NA_character_,
    seed = 9001:9100, target_replicates = 100L
  )
  nbi <- data.frame(
    scenario = "metadata_nbi", n = 500L, d = 4L, rep = 1:100,
    engine = "method_b", success = TRUE, converged = TRUE, elapsed_sec = 1, error = NA_character_,
    seed = 12001:12100, target_replicates = 100L
  )

  design <- env$jss_recovery_design_table(env$jss_recovery_attempt_ledger(bcpe, nbi))

  expect_equal(design$n_subjects, c(500L, 500L))
  expect_equal(design$n_time, c(4L, 4L))
  expect_equal(design$attempted_replicates, c(100L, 100L))
  expect_equal(design$data_seed_min, c(9001L, 12001L))
  expect_true(all(design$seed_source == "runner_metadata"))
})

test_that("legacy retry rows receive unique attempt identifiers", {
  env <- local_main_recovery_env()
  bcpe <- data.frame(
    scenario = "n500_d4", n = 500L, d = 4L, rep = 1L, model = "method_a",
    success = TRUE, converged = TRUE, elapsed_sec = 1, error = NA_character_
  )
  nbi <- data.frame(
    rep = c(1L, 1L), engine = "method_b", success = c(FALSE, TRUE), converged = c(NA, TRUE),
    elapsed_sec = c(2, 3), error = c("first attempt", NA)
  )

  ledger <- env$jss_recovery_attempt_ledger(bcpe, nbi)
  retry <- ledger[ledger$margin_family == "NBI", ]

  expect_equal(retry$retry_index, 1:2)
  expect_equal(anyDuplicated(retry$attempt_id), 0L)
})

test_that("recovery summaries include estimands, denominators, MCSEs and intervals", {
  env <- local_main_recovery_env()
  bcpe_log <- data.frame(
    scenario = "n500_d4", n = 500L, d = 4L, rep = 1:4, model = "gamlss.longitudinal",
    success = TRUE, converged = TRUE, elapsed_sec = 1:4, error = NA_character_
  )
  nbi_log <- data.frame(rep = 1:4, engine = "ours_rs_joint", success = TRUE, converged = TRUE, elapsed_sec = 5:8, error = NA_character_)
  ledger <- env$jss_recovery_attempt_ledger(bcpe_log, nbi_log)
  fixed <- data.frame(
    study_id = "bcpe_t_legacy_reconciliation", scenario_id = "n500_d4", method = "gamlss.longitudinal", replicate = 1:4,
    parameter = "mu", term = "x1", estimate = c(0.8, 0.9, 1.1, 1.2), std_error = 0.1, true_value = 1
  )

  out <- env$jss_recovery_fixed_summary(fixed, ledger)

  expect_equal(out$bias, 0, tolerance = 1e-12)
  expect_equal(out$rmse, sqrt(mean(c(-0.2, -0.1, 0.1, 0.2)^2)))
  expect_equal(out$empirical_sd, stats::sd(c(0.8, 0.9, 1.1, 1.2)))
  expect_equal(out$mean_se, 0.1)
  expect_equal(out$coverage, 0.5)
  expect_equal(out[c("attempted", "converged", "retained", "n_metric", "n_se")], data.frame(attempted = 4L, converged = 4L, retained = 4L, n_metric = 4L, n_se = 4L))
  expect_true(all(is.finite(as.matrix(out[c("bias_mcse", "bias_ci_low", "bias_ci_high", "rmse_mcse", "empirical_sd_mcse", "empirical_sd_ci_low", "empirical_sd_ci_high", "coverage_ci_low", "coverage_ci_high")]))))
})

test_that("weak Student t shape recovery is a separate natural-scale output", {
  env <- local_main_recovery_env()
  bcpe_log <- data.frame(
    scenario = "n500_d4", n = 500L, d = 4L, rep = 1:3, model = "gamlss.longitudinal",
    success = TRUE, converged = TRUE, elapsed_sec = 1, error = NA_character_
  )
  nbi_log <- data.frame(rep = 1:3, engine = "ours_rs_joint", success = TRUE, converged = TRUE, elapsed_sec = 1, error = NA_character_)
  ledger <- env$jss_recovery_attempt_ledger(bcpe_log, nbi_log)
  fixed <- rbind(
    data.frame(study_id = "bcpe_t_legacy_reconciliation", scenario_id = "n500_d4", method = "gamlss.longitudinal", replicate = 1:3, parameter = "theta", term = "intercept", estimate = atanh(0.5), std_error = 0.1, true_value = atanh(0.5)),
    data.frame(study_id = "bcpe_t_legacy_reconciliation", scenario_id = "n500_d4", method = "gamlss.longitudinal", replicate = 1:3, parameter = "zeta", term = "intercept", estimate = log(3), std_error = 0.2, true_value = log(3))
  )

  out <- env$jss_recovery_t_shape_summary(fixed, ledger)

  expect_setequal(out$quantity, c("zeta_link", "degrees_of_freedom", "rho", "lower_tail_dependence"))
  expect_equal(out$true_value_mean[out$quantity == "degrees_of_freedom"], 5)
  expect_equal(out$true_value_mean[out$quantity == "rho"], 0.5)
  expect_true(all(out$reference_profile == "intercept predictor: t=x1=x2=0; fitted smooth centering retained"))
  expect_true(all(c("mean_zeta_interval_width_mcse", "mean_df_interval_width_ci_high", "median_zeta_interval_width_mcse", "median_df_interval_width_ci_low") %in% names(out)))
})

test_that("staged n=500 T=4 R=100 evidence reconciles and resumes", {
  env <- local_main_recovery_env()
  root <- local_main_recovery_root()
  input <- file.path(root, "paper", "data", "public-derived")
  output <- tempfile("main-recovery-evidence-")

  first <- env$jss_build_main_recovery_evidence(input, output, resume = TRUE, repo_root = root)
  second <- env$jss_build_main_recovery_evidence(input, output, resume = TRUE, repo_root = root)
  design <- utils::read.csv(file.path(output, "design_table.csv"), stringsAsFactors = FALSE)
  status <- utils::read.csv(file.path(output, "attempt_status_summary.csv"), stringsAsFactors = FALSE)
  fixed <- utils::read.csv(file.path(output, "fixed_parameter_recovery.csv"), stringsAsFactors = FALSE)
  provenance <- utils::read.csv(file.path(output, "input_provenance.csv"), stringsAsFactors = FALSE)
  bundle_status <- utils::read.csv(file.path(output, "bundle_status.csv"), stringsAsFactors = FALSE)
  production_provenance <- utils::read.csv(file.path(output, "production_provenance.csv"), stringsAsFactors = FALSE)

  expect_false(isTRUE(attr(first, "resumed")))
  expect_true(isTRUE(attr(second, "resumed")))
  expect_equal(nrow(design), 2L)
  expect_true(all(design$n_subjects == 500L & design$n_time == 4L & design$attempted_replicates == 100L))
  expect_true(all(c("attempted", "successful", "converged", "convergence_not_reported", "retained") %in% names(status)))
  expect_true(all(c("bias", "rmse", "empirical_sd", "mean_se", "coverage", "bias_mcse", "rmse_ci_low", "coverage_ci_high") %in% names(fixed)))
  expect_false(any(grepl("Aydin|Users|OneDrive", provenance$path, ignore.case = TRUE)))
  expect_identical(bundle_status$status, "legacy_reconciliation_not_authoritative")
  expect_false(bundle_status$publication_eligible)
  expect_true(any(design$margin_family == "NBI" & design$copula_family == "Gaussian" & design$copula_code == "N"))
  expect_true(all(c("attempt_metadata", "fixed_by_attempt", "smooth_by_attempt", "predictive_by_attempt", "diagnostic_by_attempt") %in% provenance$input_id))
  expect_true(all(c("phase1_contract_version", "package_version", "git_sha", "source_tree_sha256", "platform") %in% production_provenance$key))
  expect_false("tracked_worktree_sha256" %in% production_provenance$key)

  protected_output <- file.path(output, "fixed_parameter_recovery.csv")
  expected_hash <- unname(env$jss_recovery_sha256(protected_output))
  writeLines("forced corruption", protected_output)
  repaired <- env$jss_build_main_recovery_evidence(input, output, resume = TRUE, repo_root = root)
  expect_false(isTRUE(attr(repaired, "resumed")))
  expect_identical(unname(env$jss_recovery_sha256(protected_output)), expected_hash)
  checkpoint <- readRDS(file.path(output, ".main-recovery-checkpoint.rds"))
  expect_setequal(checkpoint$output_files$file, basename(repaired))
  expect_false(any(grepl("[.]staging-|[.]backup-", list.files(dirname(output)), fixed = FALSE)))
})

test_that("validation rejects mixed metadata, bad retries, cardinality errors, and diagnostic orphans", {
  env <- local_main_recovery_env()
  fixture <- local_main_recovery_fixture(env)
  check <- function(id, ledger = fixture$ledger, fixed = fixture$fixed, smooth = fixture$smooth,
                    predictive = fixture$predictive, diagnostic = fixture$diagnostic) {
    x <- env$jss_recovery_validate(ledger, fixed, smooth, predictive, diagnostic)
    isTRUE(x$pass[x$check_id == id])
  }
  expect_true(all(env$jss_recovery_validate(fixture$ledger, fixture$fixed, fixture$smooth, fixture$predictive, fixture$diagnostic)$pass))

  mixed <- fixture$ledger
  mixed$n_subjects[[1]] <- 499L
  expect_false(check("scenario_metadata_consistent", ledger = mixed))
  mixed <- fixture$ledger
  mixed$n_time[[1]] <- 3L
  mixed$copula_family[[1]] <- "Gaussian"
  expect_false(check("scenario_metadata_consistent", ledger = mixed))
  wrong_family <- fixture$ledger
  wrong_family$copula_family[wrong_family$margin_family == "NBI"] <- "Gaussian"
  expect_false(check("copula_code_family_consistent", ledger = wrong_family))
  stale_runner <- fixture$ledger
  stale_runner$runner_contract_version[stale_runner$margin_family == "NBI"] <- "legacy-runner"
  expect_false(check("production_runner_contract_current", ledger = stale_runner))
  expect_identical(env$jss_recovery_bundle_status(stale_runner), "legacy_reconciliation_not_authoritative")
  missing_signature <- fixture$ledger
  missing_signature$runner_settings_signature[missing_signature$margin_family == "BCPE"] <- "legacy_unrecorded"
  expect_false(check("production_settings_signature_recorded", ledger = missing_signature))

  duplicate <- rbind(fixture$ledger, fixture$ledger[1, ])
  expect_false(check("unique_attempt_ids", ledger = duplicate))
  expect_false(check("retry_indices_contiguous", ledger = duplicate))

  missing_rep <- fixture$ledger[-which(fixture$ledger$margin_family == "NBI" & fixture$ledger$replicate == 100L), ]
  expect_false(check("method_replicate_cardinality", ledger = missing_rep))

  orphan <- fixture$diagnostic[1, ]
  orphan$replicate <- 999L
  expect_false(check("diagnostic_rows_reference_attempts", diagnostic = rbind(fixture$diagnostic, orphan)))
  expect_false(check("one_diagnostic_row_per_method_replicate", diagnostic = rbind(fixture$diagnostic, fixture$diagnostic[1, ])))
  expect_false(check("predictive_exact_schema_cardinality", predictive = fixture$predictive[-1, ]))
})

test_that("failure sensitivity and runtime quantiles include Monte Carlo uncertainty", {
  env <- local_main_recovery_env()
  fixture <- local_main_recovery_fixture(env)
  failed <- fixture$ledger$margin_family == "NBI" & fixture$ledger$replicate == 100L
  fixture$ledger$retained[failed] <- FALSE
  fixture$ledger$success[failed] <- FALSE
  keys <- c("study_id", "scenario_id", "method", "replicate")
  keep <- !do.call(paste, c(fixture$fixed[keys], sep = "\r")) %in% do.call(paste, c(fixture$ledger[failed, keys], sep = "\r"))
  fixture$fixed <- fixture$fixed[keep, ]
  fixture$smooth <- fixture$smooth[keep, ]
  fixture$predictive <- fixture$predictive[keep, ]
  sensitivity <- env$jss_recovery_failure_sensitivity(fixture$ledger, fixture$fixed, fixture$smooth, fixture$predictive)
  nbi <- sensitivity[sensitivity$study_id == "nbi_clayton_main_recovery", ]
  expect_true(all(nbi$attempted_replicates == 100L & nbi$observed_replicates == 99L))
  expect_true(all(nbi$failure_penalized_mean >= nbi$conditional_mean))
  expect_true(all(c("conditional_mcse", "conditional_ci_low", "failure_rate_mcse", "failure_rate_ci_high") %in% names(nbi)))
  runtime <- env$jss_recovery_runtime_summary(fixture$ledger)
  expect_true(all(c("median_runtime_mcse", "median_runtime_ci_low", "runtime_q25_mcse", "runtime_q75_ci_high") %in% names(runtime)))
  expect_true(all(is.finite(runtime$median_runtime_mcse)))
})

test_that("failure sensitivity retains an all-failed method", {
  env <- local_main_recovery_env()
  fixture <- local_main_recovery_fixture(env)
  failed <- fixture$ledger[fixture$ledger$margin_family == "NBI", ]
  failed$method <- "all_failed_method"
  failed$success <- FALSE
  failed$converged <- NA
  failed$retained <- FALSE
  failed$attempt_id <- sub("ours_rs_joint", "all_failed_method", failed$attempt_id, fixed = TRUE)
  failed$failure_reason <- "execution_error: forced"
  ledger <- rbind(fixture$ledger, failed)
  sensitivity <- env$jss_recovery_failure_sensitivity(ledger, fixture$fixed, fixture$smooth, fixture$predictive)
  all_failed <- sensitivity[sensitivity$method == "all_failed_method", ]
  expect_setequal(all_failed$metric, unique(sensitivity$metric))
  expect_true(all(all_failed$attempted_replicates == 100L & all_failed$observed_replicates == 0L & all_failed$failed_or_missing_replicates == 100L))
  expect_true(all(is.finite(all_failed$failure_penalty) & is.finite(all_failed$failure_penalized_mean)))
  expect_true(all(all_failed$failure_penalty_source == "same_scenario_metric_maximum_observed"))
})

test_that("production runners emit native metadata and atomic checkpoints", {
  root <- local_main_recovery_root()
  bcpe <- paste(readLines(file.path(root, "paper", "scripts", "final-simulations", "bcpe-t", "simulation_bcpe_t_gamlss_comparison.R"), warn = FALSE), collapse = "\n")
  nbi <- paste(readLines(file.path(root, "paper", "scripts", "final-simulations", "nbi-clayton", "compare_gamlss_ours_nbi_sigma_smooth.R"), warn = FALSE), collapse = "\n")

  expect_match(bcpe, 'study_id = "bcpe_t_main_recovery"', fixed = TRUE)
  expect_match(bcpe, 'seed_source = "runner_metadata"', fixed = TRUE)
  expect_match(bcpe, "extract_gamlss2_convergence", fixed = TRUE)
  expect_match(nbi, "NBI_COMPARE_MAX_ATTEMPTS_PER_FIT", fixed = TRUE)
  expect_match(nbi, 'seed_source = "runner_metadata"', fixed = TRUE)
  expect_match(nbi, "write_checkpoint_csv", fixed = TRUE)
  expect_match(nbi, "file.rename(temporary, path)", fixed = TRUE)
  expect_match(nbi, 'nbi_copula_code <- "C"', fixed = TRUE)
  expect_match(nbi, 'copula_dist = nbi_copula_code', fixed = TRUE)
  expect_match(nbi, 'theta_inverse_link = "exp"', fixed = TRUE)
  expect_false(grepl('copula_dist = "N"', nbi, fixed = TRUE))
  expect_match(nbi, "training_covariate_seed", fixed = TRUE)
  expect_match(nbi, "test_response_seed", fixed = TRUE)
  expect_match(nbi, "valid_attempt_checkpoint", fixed = TRUE)
  expect_match(bcpe, "checkpoint_result_issues", fixed = TRUE)
  expect_match(bcpe, "read_task_checkpoint", fixed = TRUE)
  expect_match(bcpe, "runner_settings_signature", fixed = TRUE)
  expect_match(bcpe, "training_covariate_seed", fixed = TRUE)
  expect_match(bcpe, "with_preserved_seed", fixed = TRUE)
})

test_that("NBI generation and fitting use the Clayton family and positive theta scale", {
  skip_if_not_installed("gamlss.longitudinal")
  skip_if_not_installed("gamlss.dist")
  dat <- gamlss.longitudinal::simulate_longitudinal_dataset(
    n = 24, times = 0:3, margin_dist = gamlss.dist::NBI(), copula_dist = "C",
    margin_params = list(mu = 3, sigma = 0.4), copula_params = list(theta = exp(0.4)), seed = 42
  )
  expect_true(all(dat$true_theta[is.finite(dat$true_theta)] > 0))
  fit <- suppressWarnings(gamlss.longitudinal::gamlss_longitudinal(
    dataset = dat, margin_dist = gamlss.dist::NBI(), copula_dist = "C",
    time_var = "time", subject_var = "subject", mu.formula = "response ~ 1",
    sigma.formula = "~ 1", theta.formula = "~ 1", method = "RS",
    max_outer_iter = 3, max_inner_iter = 5, compute_vcov = FALSE, verbose = 0
  ))
  expect_identical(fit$copula_dist, "C")
  expect_equal(gamlss.longitudinal:::.copula_family_number(fit$copula_dist), 3)
  expect_gt(exp(unname(fit$par[["theta.intercept"]])), 0)
})

test_that("retention requires convergence and complete post-fit outputs while sensitivity keeps failures", {
  env <- local_main_recovery_env()
  fixture <- local_main_recovery_fixture(env)
  target <- fixture$ledger$study_id == "nbi_clayton_main_recovery" & fixture$ledger$method == "ours_rs_joint" & fixture$ledger$replicate %in% c(98L, 99L)
  fixture$ledger$converged[target & fixture$ledger$replicate == 98L] <- FALSE
  fixture$ledger$convergence_eligible[target & fixture$ledger$replicate == 98L] <- FALSE
  fixture$predictive <- fixture$predictive[!(fixture$predictive$study_id == "nbi_clayton_main_recovery" & fixture$predictive$method == "ours_rs_joint" & fixture$predictive$replicate == 99L), ]
  ledger <- env$jss_recovery_reconcile_retention(fixture$ledger, fixture$fixed, fixture$smooth, fixture$predictive, fixture$diagnostic)
  failed <- ledger[target, ]

  expect_false(any(failed$retained))
  expect_match(failed$failure_reason[failed$replicate == 98L], "optimizer_nonconvergence")
  expect_match(failed$failure_reason[failed$replicate == 99L], "incomplete_postfit_outputs")
  ordinary <- env$jss_recovery_fixed_summary(fixture$fixed, ledger)
  row <- ordinary[ordinary$study_id == "nbi_clayton_main_recovery" & ordinary$method == "ours_rs_joint", ]
  expect_true(all(row$retained == 98L & row$n_metric == 98L))
  sensitivity <- env$jss_recovery_failure_sensitivity(ledger, fixture$fixed, fixture$smooth, fixture$predictive)
  row <- sensitivity[sensitivity$study_id == "nbi_clayton_main_recovery" & sensitivity$method == "ours_rs_joint", ]
  expect_true(all(row$failed_or_missing_replicates >= 2L))
})

test_that("rate summaries are Wilson, retry-aware, and distinguish planned cells from attempt rows", {
  env <- local_main_recovery_env()
  fixture <- local_main_recovery_fixture(env)
  retry <- fixture$ledger[fixture$ledger$study_id == "bcpe_t_main_recovery" & fixture$ledger$method == "gamlss.longitudinal" & fixture$ledger$replicate == 1L, ]
  retry$retry_index <- 2L
  retry$attempt_id <- sub("try01$", "try02", retry$attempt_id)
  retry$execution_success <- FALSE; retry$success <- FALSE; retry$converged <- NA; retry$convergence_eligible <- FALSE; retry$retained <- FALSE
  retry$failure_reason <- "execution_error: adversarial retry"; retry$retention_reason <- retry$failure_reason
  ledger <- rbind(fixture$ledger, retry)
  counts <- env$jss_recovery_counts(ledger)
  row <- counts[counts$study_id == "bcpe_t_main_recovery" & counts$method == "gamlss.longitudinal", ]
  expect_equal(row$planned_method_replicate_cells, 100L)
  expect_equal(row$actual_attempt_rows, 101L)
  expect_equal(row$attempted, 100L)
  expect_equal(row$retry_rows, 1L)
  expect_true(all(c("attempted_rate_mcse", "convergence_rate_ci_low", "retention_rate_ci_high", "failure_rate_mcse", "rate_denominator") %in% names(row)))
  failures <- env$jss_recovery_failure_summary(ledger)
  expect_true(all(c("attempt_row_denominator", "rate_mcse", "rate_ci_low", "rate_ci_high") %in% names(failures)))
  fixed_base <- env$jss_recovery_attach_attempt_identity(fixture$fixed, ledger)
  retry_fixed <- fixed_base[fixed_base$study_id == retry$study_id & fixed_base$method == retry$method & fixed_base$replicate == retry$replicate, ]
  retry_fixed$retry_index <- 2L; retry_fixed$attempt_id <- retry$attempt_id; retry_fixed$estimate <- 999
  sensitivity <- env$jss_recovery_failure_sensitivity(ledger, rbind(fixed_base, retry_fixed), fixture$smooth, fixture$predictive)
  fixed_loss <- sensitivity$conditional_mean[sensitivity$study_id == retry$study_id & sensitivity$method == retry$method & sensitivity$metric == "fixed_effect_mse"]
  expect_equal(fixed_loss, 0)
})

test_that("weak t-shape output retains infinite widths and registers the near-Gaussian threshold", {
  env <- local_main_recovery_env()
  fixture <- local_main_recovery_fixture(env)
  fixed <- fixture$fixed
  zeta <- fixed$study_id == "bcpe_t_main_recovery" & fixed$method == "gamlss.longitudinal" & fixed$parameter == "zeta" & fixed$term == "intercept"
  fixed$estimate[zeta & fixed$replicate == 1L] <- 1000
  fixed$std_error[zeta & fixed$replicate == 2L] <- Inf
  fixed$estimate[zeta & fixed$replicate == 3L] <- NA_real_
  out <- env$jss_recovery_t_shape_summary(fixed, fixture$ledger)
  expect_true(all(out$near_gaussian_zeta_threshold == 7.5))
  expect_true(all(out$df_overflow_count >= 1L))
  expect_true(all(out$df_interval_width_infinite_count + out$df_interval_width_overflow_count >= 1L))
  expect_true(all(out$zeta_missing_count >= 1L))
  expect_match(unique(out$near_gaussian_threshold_contract), "retained, not discarded", fixed = TRUE)
})

test_that("all-NA substantive output is a named incomplete attempt and never retained", {
  env <- local_main_recovery_env()
  fixture <- local_main_recovery_fixture(env)
  target_id <- fixture$ledger$attempt_id[[1]]
  fixed <- env$jss_recovery_attach_attempt_identity(fixture$fixed, fixture$ledger)
  smooth <- env$jss_recovery_attach_attempt_identity(fixture$smooth, fixture$ledger)
  predictive <- env$jss_recovery_attach_attempt_identity(fixture$predictive, fixture$ledger)
  diagnostic <- env$jss_recovery_attach_attempt_identity(fixture$diagnostic, fixture$ledger)
  fixed$estimate[fixed$attempt_id == target_id] <- NA_real_
  smooth$irmse[smooth$attempt_id == target_id] <- NA_real_
  predictive$variogram_score[predictive$attempt_id == target_id] <- NA_real_
  diagnostic$rosenblatt_ks[diagnostic$attempt_id == target_id] <- NA_real_
  ledger <- env$jss_recovery_reconcile_retention(fixture$ledger, fixed, smooth, predictive, diagnostic)
  attacked <- ledger$attempt_id == target_id
  expect_false(ledger$retained[attacked])
  expect_false(ledger$descriptive_outputs_complete[attacked])
  expect_match(ledger$failure_reason[attacked], "^incomplete_postfit_outputs:")
})

test_that("BCPE gamlss2 is marginal-only and non-common paired estimands are explicit", {
  env <- local_main_recovery_env()
  fixture <- local_main_recovery_fixture(env)
  contract <- env$jss_recovery_metric_contract("BCPE", "gamlss2")
  expect_false(any(grepl("^(theta|zeta)\\r", contract$fixed)))
  gamlss2 <- fixture$fixed$study_id == "bcpe_t_main_recovery" & fixture$fixed$method == "gamlss2"
  expect_setequal(unique(fixture$fixed$parameter[gamlss2]), c("mu", "sigma", "nu", "tau"))
  paired <- env$jss_recovery_paired_method_differences(fixture$ledger, fixture$fixed, fixture$smooth, fixture$predictive, fixture$diagnostic)
  noncommon <- paired$study_id == "bcpe_t_main_recovery" & grepl("^(theta|zeta):", paired$stratum)
  expect_true(any(noncommon))
  expect_true(all(paired$analysis[noncommon] == "not_applicable_noncommon_estimand"))
  expect_true(all(paired$paired_denominator[noncommon] == 0L))
  expect_true(all(is.na(paired$estimate[noncommon])))
})

test_that("BCPE checkpoints are append-only per-method retries with a registered cap", {
  root <- local_main_recovery_root()
  source <- paste(readLines(file.path(root, "paper", "scripts", "final-simulations", "bcpe-t", "simulation_bcpe_t_gamlss_comparison.R"), warn = FALSE), collapse = "\n")
  expect_match(source, "attempt_checkpoints", fixed = TRUE)
  expect_match(source, "BCPE_MAX_ATTEMPTS_PER_FIT", fixed = TRUE)
  expect_match(source, "publication_candidate %in% TRUE", fixed = TRUE)
  expect_match(source, "cell$retry_index <- nrow(prior) + 1L", fixed = TRUE)
  expect_match(source, "Append-only checkpoint already exists", fixed = TRUE)
})

test_that("paired fixed coverage excludes unavailable inference instead of coding it as zero", {
  env <- local_main_recovery_env()
  fixture <- local_main_recovery_fixture(env)
  fixed <- fixture$fixed
  fixed$inference_status <- "available"
  fixed$inference_denominator <- 1L
  unavailable <- fixed$study_id == "bcpe_t_main_recovery" & fixed$method == "gamlss2" &
    fixed$replicate == 1L & fixed$parameter == "mu" & fixed$term == "intercept"
  expect_equal(sum(unavailable), 1L)
  fixed$std_error[unavailable] <- NA_real_
  fixed$inference_status[unavailable] <- "unavailable_computation_failed"
  fixed$inference_denominator[unavailable] <- 0L
  paired <- env$jss_recovery_paired_method_differences(fixture$ledger, fixed, fixture$smooth, fixture$predictive, fixture$diagnostic)
  row <- paired$study_id == "bcpe_t_main_recovery" & paired$metric == "fixed_coverage" & paired$stratum == "mu:intercept" & paired$analysis == "complete_pair"
  expect_equal(sum(row), 1L)
  expect_equal(paired$paired_denominator[row], 99L)
  expect_equal(paired$estimate[row], 0)
})

test_that("fixed inference availability contract accepts named unavailability only", {
  env <- local_main_recovery_env()
  contract <- env$jss_recovery_metric_contract("BCPE", "gamlss2")
  keys <- strsplit(contract$fixed, "\r", fixed = TRUE)
  fixed <- data.frame(parameter = vapply(keys, `[`, character(1), 1), term = vapply(keys, `[`, character(1), 2),
    estimate = 1, true_value = 1, std_error = 0.1, inference_status = "available", inference_denominator = 1L)
  smooth <- data.frame(parameter = contract$smooth, bias_abs_integrated = 0.1, irmse = 0.2)
  predictive <- data.frame(variogram_p = contract$variogram_p, test_log_score_joint = -1, test_log_score_marginal = -1,
    test_log_score_copula = -1, test_log_score_per_obs = -1, variogram_score = 0.1, predictive_nsim = 100L)
  diagnostic <- data.frame(logLik = -1, rosenblatt_ks = 0.1, rosenblatt_cvm = 0.1, abs_rosenblatt_lag1_cor = 0.1,
    abs_rosenblatt_normal_lag1_cor = 0.1, rosenblatt_mean_abs_time_mean = 0.1, rosenblatt_normal_mean_abs_time_mean = 0.1)
  expect_true(all(env$jss_recovery_attempt_output_contract(fixed, smooth, predictive, diagnostic, contract)))
  fixed$std_error[[1]] <- NA_real_; fixed$inference_status[[1]] <- "unavailable_computation_failed"; fixed$inference_denominator[[1]] <- 0L
  expect_true(env$jss_recovery_attempt_output_contract(fixed, smooth, predictive, diagnostic, contract)[["fixed"]])
  bad <- fixed; bad$inference_status[[1]] <- NA_character_
  expect_false(env$jss_recovery_attempt_output_contract(bad, smooth, predictive, diagnostic, contract)[["fixed"]])
  bad <- fixed; bad$std_error[[1]] <- -0.1; bad$inference_status[[1]] <- "available"; bad$inference_denominator[[1]] <- 1L
  expect_false(env$jss_recovery_attempt_output_contract(bad, smooth, predictive, diagnostic, contract)[["fixed"]])
})

test_that("checkout source identity is clone-portable and rejects paths outside the checkout", {
  env <- local_main_recovery_env()
  root <- local_main_recovery_root()
  expect_identical(env$jss_recovery_portable_source_identity(c(".", "./", ".\\")), rep(".", 3L))
  expect_true(all(is.na(env$jss_recovery_portable_source_identity(c("../outside", "C:/outside", "/outside", "subdir")))))
  clone <- tempfile("portable-package-clone-"); dir.create(clone); dir.create(file.path(clone, "R"))
  source_files <- env$jss_recovery_package_source_files(root)
  destinations <- file.path(clone, ifelse(basename(dirname(source_files)) == "R", file.path("R", basename(source_files)), basename(source_files)))
  expect_true(all(file.copy(source_files, destinations)))
  expect_identical(env$jss_recovery_package_source_sha256(clone), env$jss_recovery_package_source_sha256(root))
})

test_that("producer identity binds every publication control file by exact bytes", {
  env <- local_main_recovery_env(); root <- local_main_recovery_root()
  identity <- env$jss_recovery_control_source_identity(root)
  required <- c("paper/R/main-recovery-evidence.R",
    "paper/scripts/final-simulations/bcpe-t/simulation_bcpe_t_gamlss_comparison.R",
    "paper/scripts/final-simulations/nbi-clayton/compare_gamlss_ours_nbi_sigma_smooth.R",
    "paper/scripts/final-simulations/main-recovery/run_main_recovery_evidence.R",
    "paper/R/phase2-paper-evidence.R", "paper/R/phase2-evidence-contracts.R",
    "paper/R/phase2-central-integration.R", "paper/R/public-paper-producers.R",
    "paper/R/replication-helpers.R", "paper/_targets.R", "paper/manifest.csv",
    "paper/phase2-artifact-contract.csv", "paper/phase2-claims.csv",
    "paper/scripts/phase2-evidence-approval.R", "paper/notes/phase2-evidence-signing.md")
  expect_true(all(required %in% identity$source_path))
  expect_true(all(grepl("^[0-9a-f]{64}$", identity$sha256)))
  expect_true(all(identity$bytes > 0))
  expected <- digest::digest(paste(identity$source_path, identity$sha256, identity$bytes, sep = "\t", collapse = "\n"), algo = "sha256", serialize = FALSE)
  expect_identical(env$jss_recovery_producer_sha256(root), expected)
})

test_that("publication checkout policy rejects an added untracked file", {
  skip_if(Sys.which("git") == "")
  env <- local_main_recovery_env()
  repo <- tempfile("clean-checkout-policy-"); dir.create(repo)
  system2("git", c("-C", shQuote(repo), "init", "--quiet"))
  writeLines("registered", file.path(repo, "registered.txt"))
  system2("git", c("-C", shQuote(repo), "add", "registered.txt"))
  system2("git", c("-C", shQuote(repo), "-c", "user.name=Evidence_Test", "-c", "user.email=evidence@example.invalid", "commit", "--quiet", "-m", "fixture"))
  clean <- env$jss_recovery_git_identity(repo)
  expect_identical(unname(clean[["git_state"]]), "clean")
  expect_silent(env$jss_recovery_require_clean_checkout(clean))
  writeLines("untracked", file.path(repo, "untracked-attack.txt"))
  dirty <- env$jss_recovery_git_identity(repo)
  expect_identical(unname(dirty[["git_state"]]), "dirty")
  expect_error(env$jss_recovery_require_clean_checkout(dirty), "no modified or untracked files")
})

test_that("main-recovery promotion binds stable registered sources rather than evidence files", {
  skip_if(Sys.which("git") == "")
  env <- local_main_recovery_env()
  root <- local_main_recovery_root()
  identity <- env$jss_recovery_source_tree_identity(root)
  expect_false(any(grepl("^paper/data/public-derived/", identity$source_path)))
  expect_true(all(grepl("^[0-9a-f]{64}$", identity$sha256)))
  provenance <- env$jss_recovery_production_provenance(root, local_main_recovery_fixture(env)$ledger)
  expect_identical(provenance$key[provenance$key == "source_tree_sha256"], "source_tree_sha256")
  expect_identical(
    provenance$value[provenance$key == "source_tree_sha256"],
    env$jss_recovery_source_tree_sha256(root)
  )

  recorded_sha <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)[[1L]]
  expect_silent(env$jss_recovery_validate_recorded_git(recorded_sha, "clean", root, "fixture"))
  expect_error(
    env$jss_recovery_validate_recorded_git(paste(rep("f", 40L), collapse = ""), "clean", root, "fixture"),
    "clean ancestor"
  )
  env$jss_recovery_source_git_state <- function(repo_root) "clean"
  expect_silent(env$jss_recovery_require_clean_source_checkout(root))
  env$jss_recovery_source_git_state <- function(repo_root) "dirty"
  expect_error(env$jss_recovery_require_clean_source_checkout(root), "unchanged registered source tree")
})

test_that("BCPE canonical configuration covers every registered DGP, fit, and checkpoint control", {
  env <- local_main_recovery_env()
  fixture <- local_main_recovery_fixture(env)
  fields <- env$jss_recovery_bcpe_settings_fields()
  expect_identical(names(fixture$bcpe_settings)[seq_along(fields)], fields)
  baseline <- env$jss_recovery_settings_signature(fixture$bcpe_settings, "BCPE")
  for (field in fields) {
    mutated <- fixture$bcpe_settings
    mutated[[field]] <- paste0(mutated[[field]], "_mutation")
    expect_false(identical(env$jss_recovery_settings_signature(mutated, "BCPE"), baseline), info = field)
  }
  source <- paste(readLines(file.path(local_main_recovery_root(), "paper", "scripts", "final-simulations", "bcpe-t", "simulation_bcpe_t_gamlss_comparison.R"), warn = FALSE), collapse = "\n")
  expect_true(all(vapply(fields, function(field) grepl(paste0(field, " ="), source, fixed = TRUE), logical(1))))
  expected_env <- c("N_FITS", "REP_IDS", "BCPE_MAX_ATTEMPTS_PER_FIT", "N_CORES", "R_SCRIPT", "VERBOSE_FITS", "BCPE_N", "BCPE_T", "SMOOTH_K", "THETA_EFFECT_MODE", "THETA_BINARY_EFFECT", "THETA_BINARY_PROB", "COMPUTE_SE", "SAVE_FITS", "VCOV_METHOD_LONGITUDINAL", "INCLUDE_DLCOPDPAR", "OPT_METHOD", "MAX_OUTER_ITER", "MAX_INNER_ITER", "MAX_ELAPSED_SEC", "OUTER_STOP_CRIT", "INNER_STOP_CRIT", "USE_BACKTRACKING", "BACKTRACKING_MAX_HALVES", "START_STEP_SIZE", "STEP_ADJUSTMENT", "MAX_STEPS", "LAMBDA_START", "WARM_START_JOINT", "CG_MAX_DELTA", "CG_ARMIJO_C1", "CG_MAX_STALL", "CG_UPDATE_LAMBDA", "CG_LINE_SEARCH", "CG_MAX_LINE_SEARCH_EVALS", "CG_GRADIENT_METHOD", "CG_ZETA_HESSIAN", "CG_LAMBDA_UPDATE_EVERY", "CG_MAX_LAMBDA_UPDATES", "CG_RAW_LOGLIK_DROP_TOL", "COMPUTE_PREDICTIVE_SCORES", "PREDICTIVE_NSIM", "VARIOGRAM_P", "VARIOGRAM_P_VALUES")
  expect_true(all(vapply(expected_env, function(name) grepl(paste0('Sys.getenv("', name, '"'), source, fixed = TRUE), logical(1))))
})

test_that("recovery manuscript guards are local and do not capture unrelated figures or sections", {
  lines <- readLines(file.path(local_main_recovery_root(), "paper", "manuscript", "main.tex"), warn = FALSE)
  opens <- grep("^\\\\IfFileExists\\{charts/paper_simulation_", lines)
  expect_length(opens, 4L)
  for (index in opens) {
    close <- grep("\\\\end\\{figure\\}\\}\\{\\}", lines[(index + 1L):min(length(lines), index + 8L)])
    expect_true(length(close) >= 1L, info = paste("guard line", index))
    block <- lines[index:(index + close[[1]])]
    expect_false(any(grepl("^\\\\section|^\\\\subsection", block)))
  }
  expect_false(any(grepl("IfFileExists.*intro - copula examples|IfFileExists.*software - plot_copula", lines)))
})

test_that("public bundle validation recomputes evidence and rejects a one-cell mutation of every artifact", {
  env <- local_main_recovery_env()
  fixture <- local_main_recovery_fixture(env)
  ledger <- fixture$ledger
  validation <- env$jss_recovery_validate(ledger, fixture$fixed, fixture$smooth, fixture$predictive, fixture$diagnostic)
  expect_true(all(validation$pass))
  provenance <- data.frame(schema_version = env$jss_recovery_schema_version(), producer_sha256 = env$jss_recovery_producer_sha256(local_main_recovery_root()), analysis_config = "conf.level=0.95", input_id = "fixture", path = "fixture.csv", sha256 = paste(rep("a", 64), collapse = ""), bytes = 1, stringsAsFactors = FALSE)
  outputs <- list(
    attempt_metadata = ledger, design_table = env$jss_recovery_design_table(ledger), attempt_status_summary = env$jss_recovery_counts(ledger),
    failure_reason_summary = env$jss_recovery_failure_summary(ledger), fixed_parameter_recovery = env$jss_recovery_fixed_summary(fixture$fixed, ledger),
    smooth_recovery = env$jss_recovery_smooth_summary(fixture$smooth, ledger), runtime_summary = env$jss_recovery_runtime_summary(ledger),
    predictive_metrics = env$jss_recovery_predictive_summary(fixture$predictive, ledger), diagnostic_metrics = env$jss_recovery_long_metric_summary(fixture$diagnostic, ledger, c("rosenblatt_ks", "rosenblatt_cvm", "abs_rosenblatt_lag1_cor", "abs_rosenblatt_normal_lag1_cor", "rosenblatt_mean_abs_time_mean", "rosenblatt_normal_mean_abs_time_mean")),
    weak_t_copula_shape_recovery = env$jss_recovery_t_shape_summary(fixture$fixed, ledger), failure_inclusive_sensitivity = env$jss_recovery_failure_sensitivity(ledger, fixture$fixed, fixture$smooth, fixture$predictive),
    paired_method_differences = env$jss_recovery_paired_method_differences(ledger, fixture$fixed, fixture$smooth, fixture$predictive, fixture$diagnostic),
    fixed_by_attempt = env$jss_recovery_attach_attempt_identity(fixture$fixed, ledger), smooth_by_attempt = env$jss_recovery_attach_attempt_identity(fixture$smooth, ledger),
    predictive_by_attempt = env$jss_recovery_attach_attempt_identity(fixture$predictive, ledger), diagnostic_by_attempt = env$jss_recovery_attach_attempt_identity(fixture$diagnostic, ledger),
    runner_settings_identity = unique(ledger[c("margin_family", "study_id", "runner_settings_signature", "runner_settings_sha256", "runner_sha256", "package_source_path", "package_version", "package_source_sha256")]),
    runner_settings_bcpe = fixture$bcpe_settings, runner_settings_nbi = fixture$nbi_settings,
    control_source_manifest = env$jss_recovery_control_source_identity(local_main_recovery_root()),
    bundle_status = data.frame(schema_version = env$jss_recovery_schema_version(), status = "candidate_post_phase1_unapproved", publication_eligible = FALSE, reason = "candidate fixture requires external approval"),
    production_provenance = env$jss_recovery_production_provenance(local_main_recovery_root(), ledger), evidence_validation = validation
  )
  outputs$input_provenance <- env$jss_recovery_frame_identity(outputs[c("attempt_metadata", "fixed_by_attempt", "smooth_by_attempt", "predictive_by_attempt", "diagnostic_by_attempt", "runner_settings_identity", "runner_settings_bcpe", "runner_settings_nbi")])
  pristine <- tempfile("authoritative-main-recovery-")
  env$jss_recovery_install_bundle(outputs, pristine, provenance, local_main_recovery_root())
  expect_silent(env$jss_main_recovery_validate_public_bundle(pristine, repo_root = local_main_recovery_root(), require_attestation = FALSE))
  self_approved_csv <- tempfile("self-approved-main-recovery-", fileext = ".csv")
  utils::write.csv(data.frame(status = "approved", bundle_sha256 = env$jss_recovery_bundle_trust_hash(pristine)), self_approved_csv, row.names = FALSE)
  expect_error(env$jss_main_recovery_validate_public_bundle(pristine, repo_root = local_main_recovery_root(), trust_registry = self_approved_csv), "unused argument")
  expect_false("public_key" %in% names(formals(env$jss_main_recovery_validate_public_bundle)))
  expect_identical(digest::digest(env$jss_recovery_pinned_public_key(), "sha256", serialize = FALSE), "cb73c05cede55bfd56357b1780b90c1bf254413d82765c7a682fc9db3a0d8587")
  fake_attestation <- tempfile(fileext = ".bin"); fake_signature <- tempfile(fileext = ".sig")
  fake_key <- sodium::sig_keygen(); fake_identity <- env$jss_recovery_bundle_identity(pristine)
  fake_message <- serialize(c(fake_identity, list(approved_at_utc = "2026-09-02T00:00:00Z", approver = "self-appointed")), NULL, version = 3L)
  writeBin(fake_message, fake_attestation); writeBin(sodium::sig_sign(fake_message, fake_key), fake_signature)
  expect_error(env$jss_main_recovery_validate_public_bundle(pristine, repo_root = local_main_recovery_root(), attestation_path = fake_attestation, signature_path = fake_signature), "pinned key")
  expect_error(env$jss_recovery_parse_rfc3339_utc("2026-02-30T00:00:00Z"), "not a real RFC3339")
  expect_error(env$jss_recovery_validate_approval_order("2026-08-31T23:59:59Z", ledger$execution_completed_at_utc), "after every execution")
  for (attack in c("named_status", "basis", "extra_field")) {
    forged_convergence <- tempfile(paste0("forged-convergence-", attack, "-")); dir.create(forged_convergence)
    expect_true(all(file.copy(list.files(pristine, full.names = TRUE, all.files = TRUE, no.. = TRUE), forged_convergence, recursive = TRUE)))
    attempt_attack <- utils::read.csv(file.path(forged_convergence, "attempt_metadata.csv"), stringsAsFactors = FALSE, check.names = FALSE)
    target <- which(attempt_attack$method == "gamlss2")[[1]]
    if (attack == "named_status") attempt_attack$raw_convergence_status[target] <- "explicit_optimizer_convergence"
    if (attack == "basis") attempt_attack$raw_convergence_basis[target] <- "forged_but_named_basis"
    if (attack == "extra_field") attempt_attack$raw_convergence_forged <- "coherent-looking"
    utils::write.csv(attempt_attack, file.path(forged_convergence, "attempt_metadata.csv"), row.names = FALSE, na = "")
    local_repin_main_recovery_bundle(env, forged_convergence)
    expect_error(env$jss_main_recovery_validate_public_bundle(forged_convergence, repo_root = local_main_recovery_root(), require_attestation = FALSE),
      "named status contradicts|API/basis/indicator contradiction|exactly match the registered schema", info = attack)
  }
  forged <- tempfile("coherently-rehashed-main-recovery-"); dir.create(forged)
  expect_true(all(file.copy(list.files(pristine, full.names = TRUE, all.files = TRUE, no.. = TRUE), forged, recursive = TRUE)))
  attempt <- utils::read.csv(file.path(forged, "attempt_metadata.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  identity <- utils::read.csv(file.path(forged, "runner_settings_identity.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  attempt$runner_sha256[attempt$margin_family == "BCPE"] <- paste(rep("0", 64), collapse = "")
  identity$runner_sha256[identity$margin_family == "BCPE"] <- paste(rep("0", 64), collapse = "")
  utils::write.csv(attempt, file.path(forged, "attempt_metadata.csv"), row.names = FALSE, na = "")
  utils::write.csv(identity, file.path(forged, "runner_settings_identity.csv"), row.names = FALSE, na = "")
  input_hashes <- utils::read.csv(file.path(forged, "input_provenance.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  for (name in c("attempt_metadata", "runner_settings_identity")) {
    index <- input_hashes$input_id == name; raw <- file.path(forged, paste0(name, ".csv"))
    input_hashes$sha256[index] <- unname(env$jss_recovery_sha256(raw)); input_hashes$bytes[index] <- file.info(raw)$size
  }
  utils::write.csv(input_hashes, file.path(forged, "input_provenance.csv"), row.names = FALSE, na = "")
  manifest <- utils::read.csv(file.path(forged, "output_manifest.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  for (artifact in c("attempt_metadata.csv", "runner_settings_identity.csv", "input_provenance.csv")) {
    index <- manifest$artifact == artifact; raw <- file.path(forged, artifact)
    manifest$sha256[index] <- unname(env$jss_recovery_sha256(raw)); manifest$bytes[index] <- file.info(raw)$size
  }
  utils::write.csv(manifest, file.path(forged, "output_manifest.csv"), row.names = FALSE, na = "")
  checkpoint <- utils::read.csv(file.path(forged, "bundle_checkpoint.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  checkpoint$input_provenance_sha256 <- digest::digest(paste(input_hashes$input_id, input_hashes$sha256, sep = "=", collapse = "\n"), algo = "sha256", serialize = FALSE)
  checkpoint$output_manifest_sha256 <- unname(env$jss_recovery_sha256(file.path(forged, "output_manifest.csv")))
  utils::write.csv(checkpoint, file.path(forged, "bundle_checkpoint.csv"), row.names = FALSE, na = "")
  expect_error(env$jss_main_recovery_validate_public_bundle(forged, repo_root = local_main_recovery_root(), require_attestation = FALSE), "runner source SHA mismatch")
  dropped_input <- tempfile("dropped-input-provenance-"); dir.create(dropped_input)
  expect_true(all(file.copy(list.files(pristine, full.names = TRUE, all.files = TRUE, no.. = TRUE), dropped_input, recursive = TRUE)))
  dropped_provenance <- utils::read.csv(file.path(dropped_input, "input_provenance.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  dropped_provenance <- dropped_provenance[dropped_provenance$input_id != "fixed_by_attempt", , drop = FALSE]
  utils::write.csv(dropped_provenance, file.path(dropped_input, "input_provenance.csv"), row.names = FALSE, na = "")
  local_repin_main_recovery_bundle(env, dropped_input)
  expect_error(env$jss_main_recovery_validate_public_bundle(dropped_input, repo_root = local_main_recovery_root(), require_attestation = FALSE), "canonical normalized-input allowlist")
  dropped_fixed <- tempfile("dropped-fixed-row-"); dir.create(dropped_fixed)
  expect_true(all(file.copy(list.files(pristine, full.names = TRUE, all.files = TRUE, no.. = TRUE), dropped_fixed, recursive = TRUE)))
  fixed_attack <- utils::read.csv(file.path(dropped_fixed, "fixed_by_attempt.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  utils::write.csv(fixed_attack[-1L, , drop = FALSE], file.path(dropped_fixed, "fixed_by_attempt.csv"), row.names = FALSE, na = "")
  local_repin_main_recovery_bundle(env, dropped_fixed)
  expect_error(env$jss_main_recovery_validate_public_bundle(dropped_fixed, repo_root = local_main_recovery_root(), require_attestation = FALSE), "recomputed validation failed|recomputation mismatch")
  forged_git <- tempfile("forged-git-"); dir.create(forged_git)
  expect_true(all(file.copy(list.files(pristine, full.names = TRUE, all.files = TRUE, no.. = TRUE), forged_git, recursive = TRUE)))
  for (artifact in c("runner_settings_bcpe.csv", "runner_settings_nbi.csv")) {
    settings_attack <- utils::read.csv(file.path(forged_git, artifact), stringsAsFactors = FALSE, check.names = FALSE)
    settings_attack$git_sha <- paste(rep("f", 40), collapse = ""); settings_attack$git_state <- "clean"
    utils::write.csv(settings_attack, file.path(forged_git, artifact), row.names = FALSE, na = "")
  }
  provenance_attack <- utils::read.csv(file.path(forged_git, "production_provenance.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  provenance_attack$value[provenance_attack$key == "git_sha"] <- paste(rep("f", 40), collapse = "")
  provenance_attack$value[provenance_attack$key == "git_state"] <- "clean"
  utils::write.csv(provenance_attack, file.path(forged_git, "production_provenance.csv"), row.names = FALSE, na = "")
  local_repin_main_recovery_bundle(env, forged_git)
  expect_error(env$jss_main_recovery_validate_public_bundle(forged_git, repo_root = local_main_recovery_root(), require_attestation = FALSE), "Git provenance|clean ancestor")
  public_root <- tempfile("paper-profile-public-"); dir.create(public_root)
  tracked <- file.path(public_root, "main-recovery"); dir.create(tracked)
  expect_true(all(file.copy(list.files(pristine, full.names = TRUE, all.files = TRUE, no.. = TRUE), tracked, recursive = TRUE)))
  paper_out <- tempfile("paper-profile-output-")
  original_attestation_verifier <- env$jss_recovery_require_external_attestation
  attestation_calls <- 0L
  env$jss_recovery_require_external_attestation <- function(identity, root, attestation_path, signature_path) { attestation_calls <<- attestation_calls + 1L; list(approver = "test-only") }
  mutation_source <- tempfile("post-validation-mutation-"); dir.create(mutation_source)
  expect_true(all(file.copy(list.files(pristine, full.names = TRUE, all.files = TRUE, no.. = TRUE), mutation_source, recursive = TRUE)))
  original_copy_exact <- env$jss_recovery_copy_exact
  env$jss_recovery_copy_exact <- function(source, destination) {
    result <- original_copy_exact(source, destination)
    write("post-validation mutation", file.path(mutation_source, "mutation-marker.txt"))
    bundle_status <- utils::read.csv(file.path(mutation_source, "bundle_status.csv"), stringsAsFactors = FALSE)
    bundle_status$reason[[1L]] <- "mutated after first validation"
    utils::write.csv(bundle_status, file.path(mutation_source, "bundle_status.csv"), row.names = FALSE)
    result
  }
  expect_error(env$jss_recovery_stage_approved_bundle(mutation_source, tempfile("mutation-destination-"), local_main_recovery_root(), fake_attestation, fake_signature), "source changed")
  env$jss_recovery_copy_exact <- original_copy_exact
  settings <- list(profile = "paper", public_data_dir = public_root, data_dir = file.path(paper_out, "data"),
    out_dir = paper_out, figures_dir = file.path(paper_out, "figures"), tables_dir = file.path(paper_out, "tables"), root = local_main_recovery_root(), workers = 1L,
    main_recovery_attestation = fake_attestation, main_recovery_attestation_signature = fake_signature)
  dir.create(settings$data_dir, recursive = TRUE); dir.create(settings$figures_dir); dir.create(settings$tables_dir)
  paper_result <- env$jss_run_phase2_main_recovery(settings)
  expect_identical(paper_result$status, "current")
  expect_gte(attestation_calls, 3L)
  expect_true(all(file.exists(c(paper_result$data, paper_result$tables, paper_result$figures))))
  for (artifact in env$jss_main_recovery_artifact_names()) {
    mutated <- tempfile("mutated-main-recovery-")
    dir.create(mutated)
    expect_true(all(file.copy(list.files(pristine, full.names = TRUE, all.files = TRUE, no.. = TRUE), mutated, recursive = TRUE)))
    path <- file.path(mutated, artifact)
    x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    x[[1]][[1]] <- paste0(as.character(x[[1]][[1]]), "-mutated")
    utils::write.csv(x, path, row.names = FALSE, na = "")
    expect_error(env$jss_main_recovery_validate_public_bundle(mutated, repo_root = local_main_recovery_root(), require_attestation = FALSE), "hash mismatch|checkpoint is invalid|allowlist|invalid schema", info = artifact)
  }
  env$jss_recovery_require_external_attestation <- original_attestation_verifier
})

test_that("runners ignore a stale installed copy and PSOCK workers verify checkout identity", {
  skip_if_not_installed("pkgload")
  skip_if_not_installed("digest")
  root <- local_main_recovery_root()
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  fake_lib <- tempfile("fake-installed-first-")
  fake_package <- file.path(fake_lib, "gamlss.longitudinal")
  dir.create(file.path(fake_package, "R"), recursive = TRUE)
  writeLines(c("Package: gamlss.longitudinal", "Version: 99.99.99", "Title: stale adversarial copy", "Description: must never load", "License: MIT"), file.path(fake_package, "DESCRIPTION"))
  writeLines("exportPattern('^[^.]')", file.path(fake_package, "NAMESPACE"))
  out <- tempfile("nbi-source-audit-"); dir.create(out)
  script <- file.path(root, "paper", "scripts", "final-simulations", "nbi-clayton", "compare_gamlss_ours_nbi_sigma_smooth.R")
  env_values <- c(R_LIBS_USER = paste(c(fake_lib, .libPaths()), collapse = .Platform$path.sep), NBI_SOURCE_IDENTITY_AUDIT_ONLY = "1", NBI_COMPARE_OUT_DIR = out, GAMLSS_LONGITUDINAL_SOURCE_ROOT = root)
  old_values <- Sys.getenv(names(env_values), unset = NA_character_)
  on.exit(for (i in seq_along(env_values)) if (is.na(old_values[[i]])) Sys.unsetenv(names(env_values)[[i]]) else do.call(Sys.setenv, stats::setNames(list(old_values[[i]]), names(env_values)[[i]])), add = TRUE)
  for (i in seq_along(env_values)) do.call(Sys.setenv, stats::setNames(list(env_values[[i]]), names(env_values)[[i]]))
  result <- system2(rscript, c("--vanilla", shQuote(script)), stdout = TRUE, stderr = TRUE)
  process_status <- attr(result, "status"); if (is.null(process_status)) process_status <- 0L
  expect_equal(process_status, 0L, info = paste(result, collapse = "\n"))
  identity <- utils::read.csv(file.path(out, "source_identity_audit.csv"), stringsAsFactors = FALSE)
  expect_equal(identity$package_version, unname(read.dcf(file.path(root, "DESCRIPTION"))[1, "Version"]))
  expect_equal(normalizePath(identity$loaded_namespace_path, winslash = "/"), normalizePath(root, winslash = "/"))
  expect_false(grepl("99.99.99", identity$runner_settings_signature, fixed = TRUE))

  bcpe <- paste(readLines(file.path(root, "paper", "scripts", "final-simulations", "bcpe-t", "simulation_bcpe_t_gamlss_comparison.R"), warn = FALSE), collapse = "\n")
  expect_match(bcpe, "PSOCK worker package source identity mismatch", fixed = TRUE)
  expect_match(bcpe, "worker_source_identity_audit.csv", fixed = TRUE)
  expect_match(bcpe, "pkgload::load_all(source_path", fixed = TRUE)
  worker_out <- tempfile("bcpe-worker-audit-"); dir.create(worker_out)
  worker_env <- c(OUT_DIR = worker_out, N_FITS = "2", N_CORES = "2", PARALLEL_SETUP_ONLY = "1")
  worker_old <- Sys.getenv(names(worker_env), unset = NA_character_)
  on.exit(for (i in seq_along(worker_env)) if (is.na(worker_old[[i]])) Sys.unsetenv(names(worker_env)[[i]]) else do.call(Sys.setenv, stats::setNames(list(worker_old[[i]]), names(worker_env)[[i]])), add = TRUE)
  for (i in seq_along(worker_env)) do.call(Sys.setenv, stats::setNames(list(worker_env[[i]]), names(worker_env)[[i]]))
  worker_script <- file.path(root, "paper", "scripts", "final-simulations", "bcpe-t", "simulation_bcpe_t_gamlss_comparison.R")
  worker_result <- system2(rscript, c("--vanilla", shQuote(worker_script)), stdout = TRUE, stderr = TRUE)
  worker_status <- attr(worker_result, "status"); if (is.null(worker_status)) worker_status <- 0L
  expect_equal(worker_status, 0L, info = paste(worker_result, collapse = "\n"))
  workers <- utils::read.csv(file.path(worker_out, "worker_source_identity_audit.csv"), stringsAsFactors = FALSE)
  expect_equal(nrow(workers), 2L)
  expect_equal(length(unique(workers$pid)), 2L)
  expect_true(all(normalizePath(workers$package_source_path, winslash = "/") == normalizePath(root, winslash = "/")))
  expect_true(all(workers$package_source_sha256 == identity$package_source_sha256))
})

test_that("actual gamlss2 fit emits only applicable BCPE marginal parameters", {
  skip_if_not_installed("gamlss2"); skip_if_not_installed("gamlss.dist"); skip_if_not_installed("VineCopula")
  root <- local_main_recovery_root()
  out <- tempfile("gamlss2-applicability-"); dir.create(out)
  values <- c(OUT_DIR = out, N_FITS = "1", REP_IDS = "1", N_CORES = "1", BCPE_N = "80", BCPE_T = "4",
    SMOOTH_K = "5", COMPUTE_SE = "0", BCPE_GAMLSS2_APPLICABILITY_AUDIT_ONLY = "1", GAMLSS_LONGITUDINAL_SOURCE_ROOT = root)
  old <- Sys.getenv(names(values), unset = NA_character_)
  on.exit(for (i in seq_along(values)) if (is.na(old[[i]])) Sys.unsetenv(names(values)[[i]]) else do.call(Sys.setenv, stats::setNames(list(old[[i]]), names(values)[[i]])), add = TRUE)
  for (i in seq_along(values)) do.call(Sys.setenv, stats::setNames(list(values[[i]]), names(values)[[i]]))
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  script <- file.path(root, "paper", "scripts", "final-simulations", "bcpe-t", "simulation_bcpe_t_gamlss_comparison.R")
  result <- system2(rscript, c("--vanilla", shQuote(script)), stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status"); if (is.null(status)) status <- 0L
  expect_equal(status, 0L, info = paste(result, collapse = "\n"))
  fixed <- utils::read.csv(file.path(out, "gamlss2_applicability_audit.csv"), stringsAsFactors = FALSE)
  expect_setequal(unique(fixed$parameter), c("mu", "sigma", "nu", "tau"))
  expect_false(any(fixed$parameter %in% c("theta", "zeta")))
  expect_equal(nrow(fixed), 16L)
  convergence <- utils::read.csv(file.path(out, "gamlss2_convergence_audit.csv"), stringsAsFactors = FALSE)
  expect_true(convergence$converged[convergence$case == "actual_complete_fit"])
  expect_true(convergence$publication_candidate[convergence$case == "actual_complete_fit"])
  expect_false(any(convergence$converged[convergence$case != "actual_complete_fit"]))
})

test_that("one target owns full production and paper consumes only the frozen bundle", {
  root <- local_main_recovery_root()
  targets <- paste(readLines(file.path(root, "paper", "_targets.R"), warn = FALSE), collapse = "\n")
  evidence <- paste(readLines(file.path(root, "paper", "R", "main-recovery-evidence.R"), warn = FALSE), collapse = "\n")
  gate <- paste(readLines(file.path(root, "paper", "R", "phase2-paper-evidence.R"), warn = FALSE), collapse = "\n")
  expect_equal(lengths(regmatches(targets, gregexpr("tar_target\\(module_phase2_main_recovery", targets))), 1L)
  expect_false(grepl("tar_target(module_01_bcpe_t", targets, fixed = TRUE))
  expect_false(grepl("tar_target(module_02_nbi_clayton", targets, fixed = TRUE))
  expect_match(evidence, 'identical(settings$profile, "paper")', fixed = TRUE)
  expect_match(evidence, "jss_recovery_stage_approved_bundle(tracked", fixed = TRUE)
  expect_match(evidence, "simulation_bcpe_t_gamlss_comparison.R", fixed = TRUE)
  expect_match(evidence, "compare_gamlss_ours_nbi_sigma_smooth.R", fixed = TRUE)
  expect_match(evidence, "jss_recovery_write_manuscript_assets(output_dir, settings)", fixed = TRUE)
  expect_match(evidence, 'status = "candidate_pending_external_approval"', fixed = TRUE)
  expect_match(evidence, "detached Ed25519 attestation", fixed = TRUE)
  expect_match(gate, "jss_main_recovery_validate_public_bundle(tracked_bundle", fixed = TRUE)
  expect_match(gate, "require_attestation = TRUE", fixed = TRUE)
  expect_false(grepl("bundle_status$publication_eligible", gate, fixed = TRUE))
})
