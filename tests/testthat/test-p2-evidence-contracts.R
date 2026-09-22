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

local_phase2_claim_env <- function() {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "replication-helpers.R"), local = env)
  source(file.path(root, "paper", "R", "submission-replication.R"), local = env)
  source(file.path(root, "paper", "R", "phase2-paper-evidence.R"), local = env)
  env
}

local_phase2_design_fixture <- function(root, out) {
  manifest <- utils::read.csv(file.path(root, "paper", "manifest.csv"), stringsAsFactors = FALSE)
  write_artifact <- function(id, x) {
    rel <- manifest$generated_path[manifest$artifact_id == id]
    path <- file.path(out, rel)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(x, path, row.names = FALSE, na = "")
  }
  recovery <- expand.grid(study_id = c("bcpe_t", "nbi_clayton"),
    method = c("gamlss.longitudinal", "gamlss2"), replicate = 1:100,
    stringsAsFactors = FALSE)
  recovery$scenario_id <- ifelse(recovery$study_id == "bcpe_t", "bcpe", "nbi")
  recovery$retry_index <- 1L; recovery$n_subjects <- 500L; recovery$n_time <- 4L
  recovery$target_replicates <- 100L
  recovery$margin_family <- ifelse(recovery$study_id == "bcpe_t", "BCPE", "NBI")
  recovery$copula_code <- ifelse(recovery$study_id == "bcpe_t", "t", "C")
  write_artifact("p2_recovery_attempts", recovery)
  write_artifact("p2_optimizer_design", data.frame(
    case_id = sprintf("JVS%02d", 1:3),
    changed_factor = c("none", "theta_strength", "sigma_strength")))
  write_artifact("p2_optimizer_attempts", expand.grid(case_id = sprintf("JVS%02d", 1:3),
    joint_review_rep = 1:3, method = c("rs_joint", "rs_separate"), stringsAsFactors = FALSE))
  missing <- expand.grid(missing_mechanism = c("monotone_dropout", "time_dependent_intermittent_mar"),
    target_missing_rate = seq(0, .5, .1), rep = 1:20,
    model = c("gamlss.longitudinal", "gamlss2"), stringsAsFactors = FALSE)
  missing$n <- 500L; missing$d <- 4L
  write_artifact("p2_missing_attempts", missing)
  copulas <- c("N", "C", "F", "G", "J", "t")
  copula <- expand.grid(generating_copula = copulas, fitted_copula = copulas,
    tau_label = "moderate", n_subject = 150L, rep = 1:100,
    stringsAsFactors = FALSE)
  copula$n_time <- 4L
  write_artifact("p2_copula_results", copula)
  scaling <- data.frame(
    scenario_id = c("continuous-base", "continuous-subjects", "continuous-visits", "continuous-smooth_basis",
      "discrete-base", "discrete-subjects", "discrete-visits", "discrete-smooth_basis"),
    model_class = rep(c("continuous", "discrete"), each = 4L),
    n_subjects = rep(c(100L, 250L, 100L, 100L), 2L),
    n_visits = rep(c(4L, 4L, 8L, 4L), 2L), smooth_k = rep(c(5L, 5L, 5L, 10L), 2L))
  write_artifact("p2_scaling_scenarios", scaling)
}

test_that("Phase 2 claim markers are unique, scoped, and mutation-sensitive", {
  env <- local_phase2_claim_env()
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  claims <- env$jss_phase2_claim_register(root)
  expect_identical(anyDuplicated(claims$claim_marker), 0L)
  expect_true(any(claims$expected_direction == "negative"))
  expect_true(any(claims$expected_direction == "positive"))
  expect_setequal(unique(claims$publication_scope), c("public_manuscript", "public_supplemental"))
  expect_silent(env$jss_phase2_validate_claim_markers("", claims, allow_pending = TRUE))
  expect_error(env$jss_phase2_validate_claim_markers("", claims, allow_pending = FALSE), "pending Overleaf")

  present <- claims[1L, , drop = FALSE]
  present$manuscript_state <- "present"
  expect_silent(env$jss_phase2_validate_claim_markers(present$claim_marker, present, allow_pending = FALSE))
  expect_error(env$jss_phase2_validate_claim_markers("", present, allow_pending = FALSE), "exactly once")
  expect_error(env$jss_phase2_validate_claim_markers(paste(present$claim_marker, present$claim_marker), present), "exactly once")
  expect_error(env$jss_phase2_validate_claim_markers("\\PhaseTwoClaim{unregistered}", claims), "unregistered")
  supplemental <- claims[claims$publication_scope == "public_supplemental", , drop = FALSE][1L, ]
  expect_error(env$jss_phase2_validate_claim_markers(supplemental$claim_marker, supplemental), "leaked")
  interval <- claims[claims$render_kind == "estimate_interval" & claims$evidence_state == "verified_source", , drop = FALSE][1L, ]
  interval$manuscript_state <- "present"
  bare <- sub("PhaseTwoClaimInterval", "PhaseTwoClaim", interval$claim_marker, fixed = TRUE)
  expect_error(env$jss_phase2_validate_claim_markers(bare, interval, allow_pending = FALSE), "unregistered")
  pending <- claims[claims$evidence_state == "pending_specification", , drop = FALSE]
  expect_setequal(unique(pending$study_id), c(
    "main-recovery", "optimizer-benchmark", "missingness", "fit-scaling"
  ))
  expect_true(all(nzchar(pending$pending_reason)))
  expect_true(all(pending$manuscript_state == "pending_evidence"))
})

test_that("claim lifecycle excludes warm-ups without hiding genuine failures", {
  env <- local_phase2_claim_env()
  attempts <- data.frame(
    attempted = TRUE,
    converged = c(TRUE, TRUE, FALSE),
    retained = c(FALSE, TRUE, FALSE),
    exclusion_reason = c("warmup", "none", "none"),
    failure_reason = c("none", "none", "optimizer_failure"),
    stringsAsFactors = FALSE
  )
  lifecycle <- env$jss_phase2_claim_lifecycle(attempts, "fit-scaling")
  expect_identical(lifecycle$attempted, 2L)
  expect_identical(lifecycle$converged, 1L)
  expect_identical(lifecycle$retained, 1L)
  expect_identical(lifecycle$failed, 1L)
  expect_identical(lifecycle$failure_reasons, "optimizer_failure:1")
})

test_that("Phase 2 empirical prose lint rejects unregistered numbers and directions", {
  env <- local_phase2_claim_env()
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  manuscript <- paste(readLines(file.path(root, "paper", "manuscript", "main.tex"), warn = FALSE), collapse = "\n")
  expect_silent(env$jss_phase2_lint_empirical_prose(root, manuscript, allow_pending = TRUE))
  numeric_mutation <- sub("\\label{sec:sim_recovery}",
    "\\label{sec:sim_recovery} The method had 77.7 percent success and 321 failures.", manuscript, fixed = TRUE)
  expect_error(env$jss_phase2_lint_empirical_prose(root, numeric_mutation), "Unregistered empirical")
  direction_mutation <- sub("\\label{sec:sim_recovery}",
    "\\label{sec:sim_recovery} The method substantially outperforms its comparator.", manuscript, fixed = TRUE)
  expect_error(env$jss_phase2_lint_empirical_prose(root, direction_mutation), "Unregistered empirical")
})

test_that("Phase 2 exact design registry reconstructs inputs and renders TeX-ready assets", {
  env <- local_phase2_claim_env()
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  out <- tempfile("phase2-design-"); dir.create(out)
  local_phase2_design_fixture(root, out)
  settings <- list(root = root, out_dir = out, tables_dir = file.path(out, "tables"))
  design <- env$jss_phase2_validate_design_contract(settings)
  expect_equal(nrow(design), nrow(utils::read.csv(file.path(root, "paper", "phase2-design-contract.csv"))))
  expect_true(all(grepl("^[0-9a-f]{64}$", design$source_sha256)))
  expect_equal(nrow(env$jss_phase2_validate_joint_design_contract(settings)), 5L)

  claims <- env$jss_phase2_claim_register(root)
  verified <- claims[claims$evidence_state == "verified_source", , drop = FALSE]
  evidence <- data.frame(claim_id = verified$claim_id, estimate = seq_len(nrow(verified)),
    lower = seq_len(nrow(verified)) - .1, upper = seq_len(nrow(verified)) + .1,
    attempted = 100L, retained = 99L, stringsAsFactors = FALSE)
  assets <- env$jss_write_phase2_claim_assets(settings, evidence, claims)
  expect_true(all(file.exists(assets[-1L])))
  macros <- paste(readLines(assets[["claim_macros"]], warn = FALSE), collapse = "\n")
  expect_match(macros, "PhaseTwoClaim@p2-recovery-planned-cells", fixed = TRUE)
  scenarios <- paste(readLines(assets[["scenario_design"]], warn = FALSE), collapse = "\n")
  expect_match(scenarios, "monotone\\_dropout", fixed = TRUE)
  expect_false(grepl("Measured subject counts", scenarios, fixed = TRUE))

  manifest <- utils::read.csv(file.path(root, "paper", "manifest.csv"), stringsAsFactors = FALSE)
  bad <- file.path(out, manifest$generated_path[manifest$artifact_id == "p2_scaling_scenarios"])
  x <- utils::read.csv(bad, stringsAsFactors = FALSE); x$n_subjects <- 999L; utils::write.csv(x, bad, row.names = FALSE)
  expect_error(env$jss_phase2_validate_design_contract(settings), "design values disagree")

  local_phase2_design_fixture(root, out)
  missing_path <- file.path(out, manifest$generated_path[manifest$artifact_id == "p2_missing_attempts"])
  missing <- utils::read.csv(missing_path, stringsAsFactors = FALSE)
  missing[1L, c("missing_mechanism", "target_missing_rate", "rep", "model")] <-
    missing[2L, c("missing_mechanism", "target_missing_rate", "rep", "model")]
  utils::write.csv(missing, missing_path, row.names = FALSE)
  expect_error(env$jss_phase2_validate_joint_design_contract(settings), "duplicated or empty")

})

test_that("canonical manuscript paths, labels, and private scope are explicit", {
  env <- local_phase2_claim_env()
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  expect_silent(env$jss_phase2_validate_manuscript_artifact_registry(root))
  assets <- utils::read.csv(file.path(root, "paper", "phase2-manuscript-artifacts.csv"), stringsAsFactors = FALSE)
  present <- assets$canonical_state == "present"
  manuscript <- paste(c("data cannot be redistributed",
    paste0("\\label{", utils::read.csv(file.path(root, "paper", "phase2-manuscript-scope.csv"), stringsAsFactors = FALSE)$tex_anchor, "}"),
    unlist(lapply(which(present), function(i) if (assets$usage_kind[[i]] == "includegraphics")
      c("\\begin{figure}", assets$canonical_path[[i]], paste0("\\label{", assets$canonical_label[[i]], "}"), "\\end{figure}") else
      c(paste0("\\input{", assets$canonical_path[[i]], "}"), paste0("\\ref{", assets$canonical_label[[i]], "}"))))
  ), collapse = "\n")
  expect_silent(env$jss_phase2_validate_manuscript_assets(root, manuscript, allow_pending = TRUE))
  expect_silent(env$jss_phase2_validate_manuscript_scope(root, manuscript, allow_pending = TRUE))
  expect_error(env$jss_phase2_validate_manuscript_assets(root, manuscript, allow_pending = FALSE), "pending")
  expect_error(env$jss_phase2_validate_manuscript_assets(root, sub(assets$canonical_path[[1L]], "wrong.png", manuscript, fixed = TRUE)),
    "missing registered path")
  scrambled <- sub(paste0("\\label{", assets$canonical_label[[1L]], "}"), "\\label{temporary-swap}", manuscript, fixed = TRUE)
  scrambled <- sub(paste0("\\label{", assets$canonical_label[[2L]], "}"), paste0("\\label{", assets$canonical_label[[1L]], "}"), scrambled, fixed = TRUE)
  scrambled <- sub("\\label{temporary-swap}", paste0("\\label{", assets$canonical_label[[2L]], "}"), scrambled, fixed = TRUE)
  expect_error(env$jss_phase2_validate_manuscript_assets(root, scrambled), "missing registered path/label")
  expect_error(env$jss_phase2_validate_manuscript_scope(root,
    sub("data cannot be redistributed", "private data are discussed", manuscript, fixed = TRUE)), "positive nonredistribution")

  tmp <- tempfile("phase2-root-"); dir.create(file.path(tmp, "paper"), recursive = TRUE)
  file.copy(file.path(root, "paper", c("manifest.csv", "phase2-manuscript-artifacts.csv")), file.path(tmp, "paper"))
  forged <- utils::read.csv(file.path(tmp, "paper", "manifest.csv"), stringsAsFactors = FALSE)
  forged$manuscript_label[forged$artifact_id == "p2_copula_heatmap"] <- "fig:wrong"
  utils::write.csv(forged, file.path(tmp, "paper", "manifest.csv"), row.names = FALSE, na = "")
  expect_error(env$jss_phase2_validate_manuscript_artifact_registry(tmp), "path/label disagrees")
})

test_that("audited canonical manuscript exposes every registered stale Phase 2 block", {
  canonical <- "C:/w/paper-canonical-audit"
  skip_if_not(dir.exists(canonical))
  env <- local_phase2_claim_env()
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  expect_silent(env$jss_phase2_validate_canonical_source(root, canonical))
  audit <- env$jss_phase2_known_canonical_issues(root, file.path(canonical, "main.tex"))
  expect_true(all(audit$excerpt_present[audit$status == "fail"]))
  expect_identical(sum(audit$status == "fail"), 8L)
})

test_that("Phase 2 central target and seed registrations are explicit", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  target_code <- paste(readLines(file.path(root, "paper", "_targets.R"), warn = FALSE), collapse = "\n")
  expect_true(all(vapply(c(
    "phase2-evidence-contracts.R", "main-recovery-evidence.R",
    "missingness-study-helpers.R", "phase2-central-integration.R", "10-fit-scaling.R",
    "submission-replication.R"
  ), grepl, logical(1), x = target_code, fixed = TRUE)))
  expect_match(target_code, "value$target_integrated <- TRUE", fixed = TRUE)
  expect_false(grepl("module_09_multivariate_benchmark", target_code, fixed = TRUE))
  expect_false(grepl("module_08_correlation_misspecification", target_code, fixed = TRUE))
  producers <- paste(readLines(file.path(root, "paper", "R", "public-paper-producers.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl("(?m)^jss_run_03_joint_vs_separate <-", producers, perl = TRUE))
  expect_false(grepl("(?m)^jss_run_08_public <-", producers, perl = TRUE))

  seeds <- utils::read.csv(file.path(root, "paper", "seeds.csv"), stringsAsFactors = FALSE)
  expect_true(all(c(
    "main-recovery-bcpe-t", "main-recovery-nbi-clayton",
    "main-recovery-summary-bootstrap", "optimizer-benchmark",
    "missingness-data", "missingness-mechanism", "copula-misspecification",
    "fit-scaling"
  ) %in% seeds$study))
  expect_identical(seeds$base_seed[seeds$study == "copula-misspecification"], 20960528L)
  expect_match(seeds$replicate_rule[seeds$study == "optimizer-benchmark"], "3000 + 100*rep", fixed = TRUE)
})

test_that("central gate validates frozen evidence while reporting Phase 3 claim work", {
  env <- local_phase2_claim_env()
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  source(file.path(root, "paper", "R", "phase2-central-integration.R"), local = env)
  out <- tempfile("phase2-central-gate-")
  settings <- list(root = root, out_dir = out, data_dir = file.path(out, "data"),
    tables_dir = file.path(out, "tables"), public_data_dir = file.path(out, "public"))
  dir.create(settings$data_dir, recursive = TRUE)
  dir.create(settings$tables_dir, recursive = TRUE)
  dir.create(settings$public_data_dir, recursive = TRUE)

  ledger <- env$jss_phase2_claim_register(root)
  verified <- ledger[ledger$evidence_state == "verified_source", , drop = FALSE]
  claim_path <- file.path(settings$tables_dir, "phase2-claim-evidence.csv")
  utils::write.csv(data.frame(claim_id = verified$claim_id), claim_path, row.names = FALSE)
  for (name in c("phase2-claim-macros.tex", "phase2-scenario-design.tex", "phase2-manuscript-checklist.csv"))
    writeLines("fixture", file.path(settings$tables_dir, name))
  recovery <- file.path(settings$data_dir, "main-recovery")
  dir.create(recovery)
  utils::write.csv(data.frame(pass = TRUE), file.path(recovery, "evidence_validation.csv"), row.names = FALSE)

  env$jss_phase2_validate_production_modules <- function(...) TRUE
  env$jss_phase2_validate_manifest_allowlist <- function(...) TRUE
  env$jss_phase2_module_output_allowlist <- function(...) "fixture"
  env$jss_phase2_validate_claim_evidence_rows <- function(...) TRUE
  env$jss_phase2_validate_design_contract <- function(...) data.frame()
  env$jss_phase2_validate_joint_design_contract <- function(...) data.frame()
  env$jss_phase2_validate_manuscript_artifact_registry <- function(...) TRUE
  env$jss_main_recovery_validate_public_bundle <- function(...) TRUE
  env$jss_submission_validate_main_recovery <- function(...) TRUE
  modules <- stats::setNames(rep(list(list(status = "current")), 5L),
    c("main_recovery", "optimizer", "missingness", "copula_selection", "fit_scaling"))

  expect_silent(env$jss_write_phase2_gate_audit(settings, modules, claim_path))
  audit <- utils::read.csv(file.path(settings$tables_dir, "phase2-gate-audit.csv"), stringsAsFactors = FALSE)
  value <- function(id) audit$pass[audit$check_id == id]
  expect_true(value("no_registered_stale_claims"))
  expect_true(value("claim_schema_and_scope_registered"))
  expect_true(value("verified_claim_evidence_valid"))
  expect_identical(
    audit$detail[audit$check_id == "verified_claim_evidence_valid"],
    paste0("verified_claims=", nrow(verified), "; pending_phase3_claims=20")
  )
  expect_true(all(audit$pass))
  expect_identical(
    audit$source[audit$check_id == "main_recovery_code_gate"],
    "data/main-recovery/evidence_validation.csv"
  )
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
    "p2_recovery_attempts", "p2_optimizer_attempts", "p2_optimizer_pairs", "p2_optimizer_uncertainty",
    "p2_missing_attempts", "p2_copula_selection_attempts",
    "p2_missing_configuration",
    "p2_missing_warning_events", "p2_missing_warning_audit",
    "p2_scaling_attempts", "p2_scaling_hardware"
  ) %in% active))
  expect_false(any(grepl("^p2_mvt_", manifest$artifact_id)))
  claims <- utils::read.csv(file.path(root, "paper", "phase2-claims.csv"), stringsAsFactors = FALSE)
  pending_claims <- utils::read.csv(file.path(root, "paper", "phase2-pending-claims.csv"), stringsAsFactors = FALSE)
  expect_true(all(claims$attempt_artifact_id %in% active))
  expect_true(all(claims$effect_artifact_id %in% active))
  expect_true(all(pending_claims$attempt_artifact_id %in% active))
  expect_true(all(pending_claims$effect_artifact_id %in% active))
  copula_claim <- claims[claims$claim_id == "p2-copula-selection-attempts", , drop = FALSE]
  expect_identical(nrow(copula_claim), 1L)
  expect_identical(copula_claim$wording_strength, "exact")
  expect_identical(copula_claim$expected_value, 3600L)
  expect_false(any(claims$expected_value == 21600L, na.rm = TRUE))
})

test_that("Phase 2 manifest is an exact bidirectional module allowlist", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "submission-replication.R"), local = env)
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
    fit_scaling = "module_10_fit_scaling"
  )
  modules <- lapply(producer_map, function(producer) {
    paths <- p2$generated_path[p2$producer == producer]
    paths <- file.path(out, paths)
    list(module_id = producer, status = "current", data = paths, tables = character(), figures = character())
  })
  claim <- file.path(out, p2$generated_path[p2$producer == "phase2_claim_evidence"])
  expect_silent(env$jss_phase2_validate_manifest_allowlist(settings, modules, claim, manifest_path))

  missing_manifest <- tempfile(fileext = ".csv")
  utils::write.csv(manifest[manifest$artifact_id != "p2_scaling_attempts", ], missing_manifest, row.names = FALSE, na = "")
  expect_error(env$jss_phase2_validate_manifest_allowlist(settings, modules, claim, missing_manifest), "missing")
  forged <- manifest[1L, , drop = FALSE]
  forged$artifact_id <- "p2_reverse_unproduced"
  forged$generated_path <- "data/arbitrary-14400-row-bypass.csv"
  forged$access <- "public"; forged$publication_status <- "active"; forged$profiles <- "paper|full"
  reverse_manifest <- tempfile(fileext = ".csv")
  utils::write.csv(rbind(manifest, forged), reverse_manifest, row.names = FALSE, na = "")
  expect_error(env$jss_phase2_validate_manifest_allowlist(settings, modules, claim, reverse_manifest), "reverse/unproduced")
  forged_modules <- modules
  forged_modules$fit_scaling$data <- c(
    forged_modules$fit_scaling$data,
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

test_that("central claim evidence rejects stale schemas and unsupported directions", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "submission-replication.R"), local = env)
  source(file.path(root, "paper", "R", "phase2-central-integration.R"), local = env)
  registered <- data.frame(claim_id = "claim-a", expected_direction = "positive", stringsAsFactors = FALSE)
  evidence <- data.frame(
    claim_id = "claim-a", scenario_id = "all", source_row_key = "all",
    estimate = 0.2, mcse = 0.03, lower = 0.1, upper = 0.3,
    denominator = 100L, attempted = 100L, converged = 95L, retained = 94L,
    failed = 6L, failure_reasons = "optimizer_nonconvergence:6",
    source_schema_version = "jss-submission-inputs-2",
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
  stale <- evidence; stale$source_schema_version <- "unsupported"
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
  ), "current source")
  hidden_failure <- evidence; hidden_failure$failure_reasons <- "none"
  expect_error(env$jss_phase2_validate_claim_evidence_rows(
    hidden_failure, registered, settings = trusted_settings
  ), "stale")
  wrong_scenario <- evidence; wrong_scenario$scenario_id <- "wrong-scenario"
  expect_error(env$jss_phase2_validate_claim_evidence_rows(
    wrong_scenario, registered, settings = trusted_settings
  ), "current source")
  wrong_selector <- evidence; wrong_selector$source_row_key <- "scenario_id=wrong"
  expect_error(env$jss_phase2_validate_claim_evidence_rows(
    wrong_selector, registered, settings = trusted_settings
  ), "current source")
  fabricated_hashes <- evidence
  fabricated_hashes$source_schema_version <- "jss-submission-inputs-forged"
  expect_error(env$jss_phase2_validate_claim_evidence_rows(
    fabricated_hashes, registered, settings = trusted_settings
  ), "stale")
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
  expect_match(manuscript, "six supported copula families", fixed = TRUE)
  expect_match(manuscript, "100 replicates per generating family", fixed = TRUE)
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
  expect_false(grepl("multivariate_benchmark_status", helper, fixed = TRUE))
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
