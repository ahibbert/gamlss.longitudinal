local_phase2_evidence_env <- function() {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "phase2-evidence-contracts.R"), local = env)
  env
}

local_phase2_attempt_fixture <- function() {
  expand <- expand.grid(
    scenario_id = c("base", "more-visits"),
    method = c("rs_separate", "rs_joint"),
    attempt_id = 1:4,
    stringsAsFactors = FALSE
  )
  expand$study_id <- "optimizer"
  expand$seed <- 1000L + expand$attempt_id
  expand$n_subjects <- 100L
  expand$n_visits <- ifelse(expand$scenario_id == "base", 4L, 8L)
  expand$attempted <- TRUE
  expand$converged <- !(expand$scenario_id == "more-visits" & expand$method == "rs_joint" & expand$attempt_id == 4L)
  expand$retained <- expand$converged
  expand$failure_reason <- ifelse(expand$retained, "none", "optimizer_nonconvergence")
  expand
}

test_that("Phase 2 attempt and denominator contracts expose every failure", {
  env <- local_phase2_evidence_env()
  attempts <- local_phase2_attempt_fixture()
  expect_true(env$jss_phase2_validate_attempts(attempts))
  denominator <- env$jss_phase2_denominators(attempts)
  expect_equal(sum(denominator$attempted), nrow(attempts))
  expect_equal(sum(denominator$failed), 1L)
  expect_match(denominator$failure_reasons[denominator$failed == 1L], "optimizer_nonconvergence:1")

  invalid <- attempts
  invalid$failure_reason[!invalid$retained] <- "none"
  expect_error(env$jss_phase2_validate_attempts(invalid), "named failure reason")
})

test_that("Phase 2 scenario metadata and summary denominators reconcile", {
  env <- local_phase2_evidence_env()
  attempts <- local_phase2_attempt_fixture()
  scenarios <- unique(attempts[c("scenario_id", "n_subjects", "n_visits")])
  expect_true(env$jss_phase2_reconcile_scenarios(attempts, scenarios, c("n_subjects", "n_visits")))
  summary <- env$jss_phase2_denominators(attempts)
  expect_true(env$jss_phase2_reconcile_denominators(attempts, summary))

  scenarios$n_visits[scenarios$scenario_id == "base"] <- 5L
  expect_error(env$jss_phase2_reconcile_scenarios(attempts, scenarios, c("n_subjects", "n_visits")), "mismatch")
})

test_that("Phase 2 paired effects retain uncertainty and drive claim directions", {
  env <- local_phase2_evidence_env()
  attempts <- local_phase2_attempt_fixture()
  metrics <- attempts[attempts$retained, c("study_id", "scenario_id", "method", "attempt_id", "seed")]
  metrics$metric <- "runtime_sec"
  metrics$value <- ifelse(metrics$method == "rs_joint", 2, 1) + metrics$attempt_id / 100
  effects <- env$jss_phase2_paired_difference(metrics, "rs_joint", "rs_separate")
  expect_true(all(c("paired_attempts", "estimate", "mcse", "lower", "upper", "sign_probability") %in% names(effects)))
  expect_true(all(effects$estimate > 0))

  claims <- effects[c("study_id", "scenario_id", "metric", "lhs_method", "rhs_method")]
  claims$claim_id <- paste0("claim-", seq_len(nrow(claims)))
  claims$expected_direction <- "positive"
  claims$attempt_source <- "optimizer-attempts.csv"
  claims <- claims[c("claim_id", "study_id", "scenario_id", "metric", "lhs_method", "rhs_method", "expected_direction", "attempt_source")]
  verified <- env$jss_phase2_validate_claims(claims, effects)
  expect_equal(nrow(verified), nrow(claims))

  claims$expected_direction[[1L]] <- "negative"
  expect_error(env$jss_phase2_validate_claims(claims, effects), "confidence-interval support")

  crossing <- effects
  crossing$lower[[1L]] <- -0.1
  expect_error(env$jss_phase2_validate_claims(transform(claims, expected_direction = "positive"), crossing), "confidence-interval support")
  invalid_mcse <- effects; invalid_mcse$mcse[[1L]] <- -0.1
  expect_error(env$jss_phase2_validate_claims(transform(claims, expected_direction = "positive"), invalid_mcse), "Monte Carlo uncertainty")
})

test_that("Phase 2 proportion uncertainty is explicit", {
  env <- local_phase2_evidence_env()
  interval <- env$jss_phase2_wilson_interval(c(50, 90), c(100, 100))
  expect_equal(interval$estimate, c(0.5, 0.9))
  expect_true(all(interval$lower < interval$estimate & interval$upper > interval$estimate))
  expect_true(all(interval$mcse > 0))
})

test_that("Phase 2 central target and seed registrations are explicit", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  target_code <- paste(readLines(file.path(root, "paper", "_targets.R"), warn = FALSE), collapse = "\n")
  expect_true(all(vapply(c(
    "phase2-evidence-contracts.R", "main-recovery-evidence.R",
    "missingness-study-helpers.R", "phase2-central-integration.R", "10-fit-scaling.R",
    "09-simulation-multivariate-longitudinal"
  ), grepl, logical(1), x = target_code, fixed = TRUE)))
  expect_match(target_code, "value$target_integrated <- TRUE", fixed = TRUE)
  expect_match(target_code, "module_09_multivariate_benchmark", fixed = TRUE)
  expect_false(grepl("module_08_correlation_misspecification", target_code, fixed = TRUE))
  producers <- paste(readLines(file.path(root, "paper", "R", "public-paper-producers.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl("(?m)^jss_run_03_joint_vs_separate <-", producers, perl = TRUE))
  expect_false(grepl("(?m)^jss_run_08_public <-", producers, perl = TRUE))

  seeds <- utils::read.csv(file.path(root, "paper", "seeds.csv"), stringsAsFactors = FALSE)
  expect_true(all(c(
    "main-recovery-bcpe-t", "main-recovery-nbi-clayton",
    "main-recovery-summary-bootstrap", "optimizer-benchmark",
    "missingness-data", "missingness-mechanism", "copula-misspecification",
    "multivariate-benchmark", "fit-scaling"
  ) %in% seeds$study))
  expect_identical(seeds$base_seed[seeds$study == "copula-misspecification"], 20960528L)
  expect_match(seeds$replicate_rule[seeds$study == "optimizer-benchmark"], "3000 + 100*rep", fixed = TRUE)
})

test_that("Phase 2 manifest retires superseded evidence and registers current artifacts", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  manifest <- utils::read.csv(file.path(root, "paper", "manifest.csv"), stringsAsFactors = FALSE)
  expect_identical(anyDuplicated(manifest$artifact_id), 0L)
  retired <- c(
    "tab_jvs_normal", "tab_jvs_gamma", "tab_jvs_nbi", "inline_jvs_design",
    "tab_corr_design", "tab_corr_external", "tab_corr_flexible", "inline_sim_design"
  )
  expect_true(all(manifest$publication_status[match(retired, manifest$artifact_id)] == "retired"))
  active <- manifest$artifact_id[manifest$publication_status == "active" & manifest$access == "public"]
  expect_true(all(c(
    "p2_recovery_attempts", "p2_optimizer_attempts", "p2_optimizer_uncertainty",
    "p2_missing_attempts", "p2_copula_selection_attempts",
    "p2_optimizer_checkpoint_archive", "p2_optimizer_checkpoint_manifest",
    "p2_missing_checkpoint_archive", "p2_missing_checkpoint_manifest",
    "p2_missing_warning_events", "p2_missing_warning_audit",
    "p2_mvt_fit_status", "p2_mvt_benchmark_results", "p2_mvt_audit_csv",
    "p2_scaling_attempts", "p2_scaling_hardware"
  ) %in% active))
  claims <- utils::read.csv(file.path(root, "paper", "phase2-claims.csv"), stringsAsFactors = FALSE)
  expect_true(all(claims$attempt_artifact_id %in% active))
  expect_true(all(claims$effect_artifact_id %in% active))
  copula_claim <- claims[claims$claim_id == "p2-copula-selection-attempts", , drop = FALSE]
  expect_identical(nrow(copula_claim), 1L)
  expect_identical(copula_claim$wording_strength, "exact")
  expect_identical(copula_claim$expected_value, 2400L)
  expect_false(any(claims$expected_value == 21600L, na.rm = TRUE))
})

test_that("Phase 2 manifest is an exact bidirectional module allowlist", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "phase2-central-integration.R"), local = env)
  manifest_path <- file.path(root, "paper", "manifest.csv")
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  active <- manifest$access == "public" & manifest$publication_status == "active" & grepl("^p2_", manifest$artifact_id)
  p2 <- manifest[active, , drop = FALSE]
  out <- tempfile("phase2-manifest-out-"); dir.create(out)
  settings <- list(root = root, out_dir = out, data_dir = file.path(out, "data"), tables_dir = file.path(out, "tables"), figures_dir = file.path(out, "figures"))
  producer_map <- c(
    main_recovery = "module_phase2_main_recovery", optimizer = "module_03_joint_vs_separate",
    missingness = "module_04_missingness", copula_selection = "module_07_copula_misspecification",
    multivariate_benchmark = "module_09_multivariate_benchmark", fit_scaling = "module_10_fit_scaling"
  )
  modules <- lapply(producer_map, function(producer) {
    paths <- p2$generated_path[p2$producer == producer]
    paths <- file.path(out, paths)
    list(module_id = producer, status = "current", data = paths, tables = character(), figures = character())
  })
  claim <- file.path(out, p2$generated_path[p2$artifact_id == "p2_claim_evidence"])
  expect_silent(env$jss_phase2_validate_manifest_allowlist(settings, modules, claim, manifest_path))

  missing_manifest <- tempfile(fileext = ".csv")
  utils::write.csv(manifest[manifest$artifact_id != "p2_mvt_fit_status", ], missing_manifest, row.names = FALSE, na = "")
  expect_error(env$jss_phase2_validate_manifest_allowlist(settings, modules, claim, missing_manifest), "missing")
  forged <- manifest[1L, , drop = FALSE]
  forged$artifact_id <- "p2_reverse_unproduced"
  forged$generated_path <- "data/arbitrary-14400-row-bypass.csv"
  forged$access <- "public"; forged$publication_status <- "active"; forged$profiles <- "paper|full"
  reverse_manifest <- tempfile(fileext = ".csv")
  utils::write.csv(rbind(manifest, forged), reverse_manifest, row.names = FALSE, na = "")
  expect_error(env$jss_phase2_validate_manifest_allowlist(settings, modules, claim, reverse_manifest), "reverse/unproduced")
  forged_modules <- modules
  forged_modules$multivariate_benchmark$data <- c(
    forged_modules$multivariate_benchmark$data,
    file.path(out, "data", "arbitrary-14400-row-bypass.csv")
  )
  expect_error(env$jss_phase2_validate_manifest_allowlist(
    settings, forged_modules, claim, reverse_manifest
  ), "fixed artifact contract")
  p2_index <- which(active)
  reordered <- manifest
  reordered[p2_index[1:2], ] <- reordered[p2_index[2:1], ]
  reordered_manifest <- tempfile(fileext = ".csv")
  utils::write.csv(reordered, reordered_manifest, row.names = FALSE, na = "")
  expect_error(env$jss_phase2_validate_manifest_allowlist(
    settings, modules, claim, reordered_manifest
  ), "reordered")
})

test_that("central claim evidence rejects stale hashes and unsupported directions", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "phase2-central-integration.R"), local = env)
  registered <- data.frame(claim_id = "claim-a", expected_direction = "positive", stringsAsFactors = FALSE)
  evidence <- data.frame(
    claim_id = "claim-a", scenario_id = "all", source_row_key = "all",
    estimate = 0.2, mcse = 0.03, lower = 0.1, upper = 0.3,
    denominator = 100L, attempted = 100L, converged = 95L, retained = 94L,
    failed = 6L, failure_reasons = "optimizer_nonconvergence:6",
    source_sha256 = paste(rep("a", 64L), collapse = ""),
    effect_sha256 = paste(rep("a", 64L), collapse = ""),
    source_identity_sha256 = paste(rep("b", 64L), collapse = ""),
    status = "verified_directional_effect", stringsAsFactors = FALSE
  )
  trusted_settings <- list(root = root)
  env$jss_build_phase2_claim_evidence <- function(settings) {
    if (!identical(settings, trusted_settings)) stop("untrusted settings fixture", call. = FALSE)
    evidence
  }
  expect_error(
    env$jss_phase2_validate_claim_evidence_rows(evidence, registered),
    "trusted settings context.*two-argument form is disabled"
  )
  expect_error(
    env$jss_phase2_validate_claim_evidence_rows(
      evidence, registered, settings = trusted_settings, expected = evidence
    ),
    "Caller-supplied expected claim evidence is not a trusted input"
  )
  expect_silent(env$jss_phase2_validate_claim_evidence_rows(
    evidence, registered, settings = trusted_settings
  ))
  stale <- evidence; stale$source_sha256 <- "forged"
  expect_error(env$jss_phase2_validate_claim_evidence_rows(
    stale, registered, settings = trusted_settings
  ), "stale")
  crossing <- evidence; crossing$lower <- -0.1
  expect_error(env$jss_phase2_validate_claim_evidence_rows(
    crossing, registered, settings = trusted_settings
  ), "interval-supported")
  missing <- evidence[0L, ]
  expect_error(env$jss_phase2_validate_claim_evidence_rows(
    missing, registered, settings = trusted_settings
  ), "missing, stale")
  forged_estimate <- evidence; forged_estimate$estimate <- forged_estimate$lower <- forged_estimate$upper <- 999
  expect_error(env$jss_phase2_validate_claim_evidence_rows(
    forged_estimate, registered, settings = trusted_settings
  ), "current approved source")
  hidden_failure <- evidence; hidden_failure$failure_reasons <- "none"
  expect_error(env$jss_phase2_validate_claim_evidence_rows(
    hidden_failure, registered, settings = trusted_settings
  ), "stale")
  wrong_scenario <- evidence; wrong_scenario$scenario_id <- "wrong-scenario"
  expect_error(env$jss_phase2_validate_claim_evidence_rows(
    wrong_scenario, registered, settings = trusted_settings
  ), "current approved source")
  wrong_selector <- evidence; wrong_selector$source_row_key <- "scenario_id=wrong"
  expect_error(env$jss_phase2_validate_claim_evidence_rows(
    wrong_selector, registered, settings = trusted_settings
  ), "current approved source")
  fabricated_hashes <- evidence
  fabricated_hashes$source_sha256 <- paste(rep("c", 64L), collapse = "")
  fabricated_hashes$effect_sha256 <- paste(rep("c", 64L), collapse = "")
  fabricated_hashes$source_identity_sha256 <- paste(rep("d", 64L), collapse = "")
  expect_error(env$jss_phase2_validate_claim_evidence_rows(
    fabricated_hashes, registered, settings = trusted_settings
  ), "current approved source")
})

test_that("protected manuscript contains no retired Phase 2 counts or table inputs", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  manuscript <- paste(readLines(file.path(root, "paper", "manuscript", "main.tex"), warn = FALSE), collapse = " ")
  stale <- c(
    "run 1,000 simulations", "N=1000", "N=5000", "approximately 20x slower",
    "58/60 Clayton", "Gaussian 1/10 simulations", "Clayton 3/10",
    "Frank 4/10", "Gumbel 2/10",
    "simulate all available copulas", "fit with all available copulas",
    "every available copula",
    "tables/normal-joint-vs-separate-six-case-median-iqr-table",
    "tables/gamma-joint-vs-separate-six-case-median-iqr-table",
    "tables/negative-binomial-joint-vs-separate-six-case-median-iqr-table",
    "tables/08-simulation-sensitivity-correlation-misspecification-scenario-design",
    "tables/08-simulation-sensitivity-correlation-misspecification-external-correlation",
    "tables/08-simulation-sensitivity-correlation-misspecification-flexible-correlation"
  )
  hits <- stale[vapply(stale, grepl, logical(1), x = manuscript, fixed = TRUE)]
  expect_identical(hits, character())
  expect_match(manuscript, "subject-level monotone dropout", fixed = TRUE)
  expect_match(manuscript, "time-dependent intermittent MAR", fixed = TRUE)
  expect_match(manuscript, "Gaussian and Clayton copulas", fixed = TRUE)
  expect_match(manuscript, "100 replicates per generating scenario", fixed = TRUE)
})

test_that("replication output override is explicit and scoped", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "replication-helpers.R"), local = env)
  old_profile <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_PROFILE", unset = NA_character_)
  old_output <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_OUT_DIR", unset = NA_character_)
  on.exit({
    if (is.na(old_profile)) Sys.unsetenv("GAMLSS_LONGITUDINAL_JSS_PROFILE") else Sys.setenv(GAMLSS_LONGITUDINAL_JSS_PROFILE = old_profile)
    if (is.na(old_output)) Sys.unsetenv("GAMLSS_LONGITUDINAL_JSS_OUT_DIR") else Sys.setenv(GAMLSS_LONGITUDINAL_JSS_OUT_DIR = old_output)
  }, add = TRUE)
  Sys.setenv(
    GAMLSS_LONGITUDINAL_JSS_PROFILE = "full",
    GAMLSS_LONGITUDINAL_JSS_OUT_DIR = "results/jss-replication/phase2-test-output"
  )
  settings <- env$jss_settings(create = FALSE)
  expect_match(settings$out_dir, "results/jss-replication/phase2-test-output$", perl = TRUE)
})

test_that("full-profile tolerances use current Phase 2 evidence only", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  helper <- paste(readLines(file.path(root, "paper", "R", "replication-helpers.R"), warn = FALSE), collapse = "\n")
  expect_match(helper, "optimizer_benchmark", fixed = TRUE)
  expect_match(helper, "multivariate_benchmark_status", fixed = TRUE)
  expect_match(helper, "fit_scaling", fixed = TRUE)
  expect_match(helper, "recovery = file.path(public, \"main-recovery\")", fixed = TRUE)
  expect_match(helper, "recovery = file.path(settings$data_dir, \"main-recovery\")", fixed = TRUE)
  expect_false(grepl("bcpe-t-full", helper, fixed = TRUE))
  expect_false(grepl("nbi-clayton-full", helper, fixed = TRUE))
  expect_false(grepl("table_2_external_correlation.csv", helper, fixed = TRUE))
  expect_false(grepl("joint-vs-separate-optimization-deltas.csv", helper, fixed = TRUE))
  tolerances <- utils::read.csv(file.path(root, "paper", "tolerances.csv"), stringsAsFactors = FALSE)
  expect_true(all(c("paired_effect", "runtime_scaling") %in% tolerances$artifact_group))
})
