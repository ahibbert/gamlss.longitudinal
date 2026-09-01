local_phase2_scaling_env <- function() {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "10-fit-scaling.R"), local = env)
  env
}

test_that("fit-scaling scenarios vary exactly one registered factor", {
  env <- local_phase2_scaling_env()
  design <- env$jss_scaling_design("full")
  expect_true(env$jss_scaling_validate_design(design))
  expect_equal(nrow(design), 8L)
  expect_setequal(design$model_class, c("continuous", "discrete"))
  expect_setequal(design$changed_factor, c("none", "n_subjects", "n_visits", "smooth_k"))

  broken <- design
  row <- which(broken$changed_factor == "n_subjects")[[1L]]
  broken$n_visits[[row]] <- broken$n_visits[[row]] + 1L
  expect_error(env$jss_scaling_validate_design(broken), "does not change exactly")

  missing_cell <- design[-1L, , drop = FALSE]
  expect_error(env$jss_scaling_validate_design(missing_cell), "exactly the registered 8 cells")
  changed_family <- design
  changed_family$family[changed_family$model_class == "continuous"] <- "TF"
  expect_error(env$jss_scaling_validate_design(changed_family), "registered 8-cell design")
  changed_copula <- design
  changed_copula$copula[changed_copula$model_class == "discrete"] <- "N"
  expect_error(env$jss_scaling_validate_design(changed_copula), "registered 8-cell design")
})

test_that("fit-scaling precision and summary expose uncertainty and failures", {
  env <- local_phase2_scaling_env()
  rule <- list(min_attempts = 4L, max_attempts = 8L, target_relative_mcse = 0.50, bootstrap_draws = 199L)
  attempts <- data.frame(
    scenario_id = rep(c("continuous-base", "discrete-base"), each = 5L),
    attempt_id = rep(1:5, 2L),
    model_class = rep(c("continuous", "discrete"), each = 5L),
    family = rep(c("NO", "NBI"), each = 5L),
    copula = rep(c("N", "C"), each = 5L),
    changed_factor = "none",
    n_subjects = 100L,
    n_visits = 4L,
    smooth_k = 5L,
    n_observations = 400L,
    attempted = TRUE,
    converged = c(rep(TRUE, 5L), rep(TRUE, 4L), FALSE),
    retained = c(rep(TRUE, 5L), rep(TRUE, 4L), FALSE),
    elapsed_sec = c(1.0, 1.1, 0.9, 1.0, 1.05, 2.0, 2.2, 1.8, 2.1, 2.4),
    failure_reason = c(rep("none", 9L), "optimizer_nonconvergence"),
    stringsAsFactors = FALSE
  )
  summary <- env$jss_scaling_summary(attempts, rule, seed = 41L)
  reordered_summary <- env$jss_scaling_summary(attempts[nrow(attempts):1L, ], rule, seed = 41L)
  expect_equal(summary$attempted, c(5L, 5L))
  expect_equal(summary$failed, c(0L, 1L))
  expect_true(all(c("median_elapsed_sec", "q1_elapsed_sec", "q3_elapsed_sec", "median_mcse_sec", "median_ci_lower_sec", "median_ci_upper_sec", "precision_met") %in% names(summary)))
  expect_match(summary$failure_reasons[[2L]], "optimizer_nonconvergence:1")
  expect_equal(summary, reordered_summary, ignore_attr = TRUE)

  contrasts <- env$jss_scaling_contrasts(attempts, draws = 199L, seed = 41L)
  expect_equal(nrow(contrasts), 0L)
})

test_that("fit-scaling contrasts include uncertainty for every stated difference", {
  env <- local_phase2_scaling_env()
  attempts <- expand.grid(
    model_class = c("continuous", "discrete"),
    changed_factor = c("none", "n_subjects", "n_visits", "smooth_k"),
    attempt_id = 1:6,
    stringsAsFactors = FALSE
  )
  attempts$scenario_id <- paste(attempts$model_class, ifelse(attempts$changed_factor == "none", "base", attempts$changed_factor), sep = "-")
  attempts$retained <- TRUE
  attempts$elapsed_sec <- 1 + attempts$attempt_id / 100 +
    ifelse(attempts$model_class == "discrete", 1, 0) +
    ifelse(attempts$changed_factor == "none", 0, 0.5)
  contrasts <- env$jss_scaling_contrasts(attempts, draws = 199L, seed = 19L)
  reordered <- env$jss_scaling_contrasts(attempts[nrow(attempts):1L, ], draws = 199L, seed = 19L)
  expect_equal(nrow(contrasts), 6L)
  expect_true(all(c("median_difference_sec", "difference_mcse_sec", "difference_ci_lower_sec", "difference_ci_upper_sec", "median_runtime_ratio", "ratio_ci_lower", "ratio_ci_upper") %in% names(contrasts)))
  expect_true(all(contrasts$median_difference_sec > 0))
  expect_equal(contrasts, reordered, ignore_attr = TRUE)
})

test_that("fit-scaling checkpoints replace safely", {
  env <- local_phase2_scaling_env()
  path <- tempfile(fileext = ".csv")
  env$jss_scaling_atomic_write_csv(data.frame(value = 1), path)
  env$jss_scaling_atomic_write_csv(data.frame(value = 2), path)
  expect_equal(utils::read.csv(path)$value, 2)
  expect_length(Sys.glob(paste0(path, ".tmp-*")), 0L)
})

local_scaling_stub_attempt <- function(
    scenario, attempt_id, seed, warmup = FALSE, checkpoint_spec = NULL,
    execution_order = NA_integer_) {
  row <- data.frame(
    scenario_id = scenario$scenario_id[[1L]], attempt_id = as.integer(attempt_id),
    execution_order = as.integer(execution_order), seed = as.integer(seed),
    model_class = scenario$model_class[[1L]],
    family = scenario$family[[1L]], copula = scenario$copula[[1L]],
    changed_factor = scenario$changed_factor[[1L]], n_subjects = scenario$n_subjects[[1L]],
    n_visits = scenario$n_visits[[1L]], smooth_k = scenario$smooth_k[[1L]],
    n_observations = scenario$n_subjects[[1L]] * scenario$n_visits[[1L]], attempted = TRUE,
    warmup = isTRUE(warmup), fit_returned = TRUE, converged = TRUE,
    retained = !isTRUE(warmup), elapsed_sec = 1, fit_elapsed_sec = 1,
    simulation_elapsed_sec = 0.1, setup_elapsed_sec = 0.01,
    total_elapsed_sec = 1.11, timing_scope = "fit_only",
    exclusion_reason = if (isTRUE(warmup)) "warmup" else "none",
    failure_reason = "none", warning_count = 0L, warning_messages = "",
    error_message = NA_character_, stringsAsFactors = FALSE
  )
  if (!is.null(checkpoint_spec)) for (field in names(checkpoint_spec)) row[[field]] <- checkpoint_spec[[field]]
  row
}

test_that("fit-scaling seeds and checkpoints are invariant to design row order", {
  env <- local_phase2_scaling_env()
  env$jss_scaling_fit_attempt <- local_scaling_stub_attempt
  design <- env$jss_scaling_design("smoke")
  forward <- env$jss_scaling_run(design, tempfile("scaling-forward-"), "smoke", 812L)
  reverse <- env$jss_scaling_run(design[nrow(design):1L, ], tempfile("scaling-reverse-"), "smoke", 812L)
  fields <- c("scenario_id", "attempt_id", "execution_order", "seed", "warmup")
  order_rows <- function(x) x[order(x$scenario_id, x$attempt_id), fields, drop = FALSE]
  expect_equal(order_rows(forward), order_rows(reverse), ignore_attr = TRUE)
  expect_true(all(table(forward$scenario_id[forward$warmup %in% TRUE]) == 1L))
  expect_false(any(forward$retained[forward$warmup %in% TRUE] %in% TRUE))
  expect_identical(forward$execution_order, seq_len(nrow(forward)))
  expected_warmup <- env$jss_scaling_attempt_order(design$scenario_id, 0L, 812L, "warmup-order")
  expect_identical(forward$scenario_id[forward$warmup %in% TRUE], expected_warmup)
  for (round_id in 1:2) {
    expected_round <- env$jss_scaling_attempt_order(design$scenario_id, round_id, 812L, "round-order")
    actual_round <- forward$scenario_id[!forward$warmup & forward$attempt_id == round_id]
    expect_identical(actual_round, expected_round)
  }
})

test_that("fit-scaling checkpoints reject schema, configuration, design, code, package, and seed drift", {
  env <- local_phase2_scaling_env()
  env$jss_scaling_fit_attempt <- local_scaling_stub_attempt
  design <- env$jss_scaling_design("smoke")
  checkpoint_dir <- tempfile("scaling-stale-")
  env$jss_scaling_run(design, checkpoint_dir, "smoke", 91L)
  path <- file.path(checkpoint_dir, paste0(design$scenario_id[[1L]], ".csv"))
  stale <- utils::read.csv(path, stringsAsFactors = FALSE)
  valid <- stale
  bound_fields <- c(
    "checkpoint_schema_version", "design_fingerprint", "configuration_fingerprint",
    "code_fingerprint", "package_fingerprint", "timing_environment_fingerprint",
    "seed_scheme_version", "base_seed", "checkpoint_fingerprint",
    "timing_git_sha", "timing_git_dirty", "timing_r_version", "timing_platform",
    "timing_cpu_model", "timing_execution_context", "timing_dependency_versions"
  )
  expect_true(all(bound_fields %in% names(stale)))
  stale$configuration_fingerprint <- "stale"
  env$jss_scaling_atomic_write_csv(stale, path)
  expect_error(env$jss_scaling_run(design, checkpoint_dir, "smoke", 91L), "Stale or incompatible")

  env$jss_scaling_atomic_write_csv(valid, path)
  expect_error(env$jss_scaling_run(design, checkpoint_dir, "smoke", 92L), "Stale or incompatible")

  env$jss_scaling_atomic_write_csv(valid, path)
  expect_error(env$jss_scaling_run(design, checkpoint_dir, "full", 91L), "Stale or incompatible")

  env$jss_scaling_atomic_write_csv(valid, path)
  expect_error(env$jss_scaling_run(
    design, checkpoint_dir, "smoke", 91L,
    execution_context = "targets"
  ), "Stale or incompatible")

  original_environment <- env$jss_scaling_timing_environment
  env$jss_scaling_timing_environment <- function(...) {
    value <- original_environment(...)
    value$timing_platform <- "adversarial-platform"
    value
  }
  env$jss_scaling_atomic_write_csv(valid, path)
  expect_error(env$jss_scaling_run(design, checkpoint_dir, "smoke", 91L), "Stale or incompatible")
  env$jss_scaling_timing_environment <- original_environment

  original_simulate <- env$jss_scaling_simulate
  env$jss_scaling_simulate <- function(...) stop("changed implementation")
  env$jss_scaling_atomic_write_csv(valid, path)
  expect_error(env$jss_scaling_run(design, checkpoint_dir, "smoke", 91L), "Stale or incompatible")
  env$jss_scaling_simulate <- original_simulate

  original_versions <- env$jss_scaling_package_versions
  env$jss_scaling_package_versions <- function() c(original_versions(), adversarial = "1.0")
  env$jss_scaling_atomic_write_csv(valid, path)
  expect_error(env$jss_scaling_run(design, checkpoint_dir, "smoke", 91L), "Stale or incompatible")
  env$jss_scaling_package_versions <- original_versions

  original_schema <- env$jss_scaling_schema_version
  env$jss_scaling_schema_version <- function() 999L
  env$jss_scaling_atomic_write_csv(valid, path)
  expect_error(env$jss_scaling_run(design, checkpoint_dir, "smoke", 91L), "Stale or incompatible")
  env$jss_scaling_schema_version <- original_schema
})

test_that("simulation and setup failures become structured attempts", {
  env <- local_phase2_scaling_env()
  env$jss_scaling_simulate <- function(...) stop("adversarial simulation failure")
  row <- env$jss_scaling_fit_attempt(env$jss_scaling_design("smoke")[1L, ], 1L, 42L)
  expect_true(row$attempted)
  expect_false(row$fit_returned)
  expect_false(row$retained)
  expect_identical(row$failure_reason, "simulation_error")
  expect_match(row$error_message, "adversarial simulation failure")

  env <- local_phase2_scaling_env()
  env$jss_scaling_simulate <- function(...) data.frame(
    response = 1, x = 0, time_scaled = 0, s1 = 0, time = 1, subject = 1
  )
  scenario <- env$jss_scaling_design("smoke")[1L, ]
  scenario$family <- "NOT_A_FAMILY"
  row <- env$jss_scaling_fit_attempt(scenario, 1L, 42L)
  expect_true(row$attempted)
  expect_false(row$fit_returned)
  expect_identical(row$failure_reason, "setup_error")
  expect_match(row$error_message, "NOT_A_FAMILY")

  scenario$family <- "NO"
  scenario$copula <- "NOT_A_COPULA"
  row <- env$jss_scaling_fit_attempt(scenario, 1L, 42L)
  expect_true(row$attempted)
  expect_false(row$fit_returned)
  expect_identical(row$failure_reason, "fit_error")
  expect_match(row$error_message, "NOT_A_COPULA")
})

test_that("elapsed scaling metric is fit-only with setup phases separated", {
  env <- local_phase2_scaling_env()
  row <- env$jss_scaling_fit_attempt(env$jss_scaling_design("smoke")[1L, ], 1L, 42L)
  expect_identical(row$timing_scope, "fit_only")
  expect_equal(row$elapsed_sec, row$fit_elapsed_sec, tolerance = 1e-12)
  expect_true(is.finite(row$simulation_elapsed_sec))
  expect_true(is.finite(row$setup_elapsed_sec))
  expect_gte(row$total_elapsed_sec + 1e-12, row$fit_elapsed_sec)
})

test_that("warmups are excluded from timing summaries and failure counts", {
  env <- local_phase2_scaling_env()
  attempts <- data.frame(
    scenario_id = "continuous-base", attempt_id = c(0L, 1:3),
    model_class = "continuous", family = "NO", copula = "N", changed_factor = "none",
    n_subjects = 20L, n_visits = 3L, smooth_k = 4L, n_observations = 60L,
    attempted = TRUE, warmup = c(TRUE, FALSE, FALSE, FALSE), converged = TRUE,
    retained = c(FALSE, TRUE, TRUE, TRUE), elapsed_sec = c(999, 1, 1, 1),
    failure_reason = "none", stringsAsFactors = FALSE
  )
  summary <- env$jss_scaling_summary(attempts, env$jss_scaling_precision_rule("smoke"), 1L)
  expect_equal(summary$warmup_attempted, 1L)
  expect_equal(summary$attempted, 3L)
  expect_equal(summary$failed, 0L)
  expect_equal(summary$median_elapsed_sec, 1)
})

test_that("publication gate rejects all-failed, incomplete, and imprecise cells", {
  env <- local_phase2_scaling_env()
  design <- env$jss_scaling_design("full")
  attempts <- do.call(rbind, lapply(seq_len(nrow(design)), function(i) {
    scenario <- design[i, ]
    data.frame(scenario, attempt_id = 0L, seed = i, n_observations = NA_integer_,
      attempted = TRUE, warmup = TRUE, fit_returned = FALSE, converged = FALSE,
      retained = FALSE, elapsed_sec = 0.01, exclusion_reason = "warmup",
      failure_reason = "fit_error", warning_count = 0L, warning_messages = "",
      error_message = "failed", stringsAsFactors = FALSE)
  }))
  summary <- data.frame(
    scenario_id = design$scenario_id, retained = 0L, precision_met = FALSE,
    cell_eligible = FALSE, stringsAsFactors = FALSE
  )
  status <- env$jss_scaling_publication_status(
    design, attempts, summary, data.frame(), env$jss_scaling_precision_rule("full")
  )
  expect_false(status$publication_eligible)
  expect_match(status$reasons, "cell_precision_or_success_ineligible")
  expect_match(status$reasons, "contrast_cells_incomplete")
  expect_error(env$jss_scaling_assert_publication_eligible(
    design, attempts, summary, data.frame(), env$jss_scaling_precision_rule("full")
  ), "not publication eligible")
})

test_that("precision cap is explicit and unmet precision stays ineligible", {
  env <- local_phase2_scaling_env()
  env$jss_scaling_fit_attempt <- local_scaling_stub_attempt
  env$jss_scaling_precision_rule <- function(profile = "smoke") list(
    min_attempts = 2L, initial_max_attempts = 2L, extension_attempts = 1L,
    max_attempts = 3L, target_relative_mcse = -1, bootstrap_draws = 19L,
    warmup_attempts = 1L
  )
  design <- env$jss_scaling_design("smoke")
  attempts <- env$jss_scaling_run(design, tempfile("scaling-cap-"), "smoke", 10L)
  timed_counts <- table(attempts$scenario_id[!attempts$warmup])
  expect_true(all(timed_counts == 3L), info = paste(names(timed_counts), timed_counts, collapse = "; "))
  summary <- env$jss_scaling_summary(attempts, env$jss_scaling_precision_rule(), 10L)
  expect_false(any(summary$precision_met))
  expect_true(all(summary$extension_used), info = paste(summary$scenario_id, summary$attempted, collapse = "; "))
  expect_false(env$jss_scaling_publication_status(
    design, attempts, summary, data.frame(), env$jss_scaling_precision_rule()
  )$publication_eligible)
})

test_that("tracked fit-scaling evidence reconciles with its attempts", {
  env <- local_phase2_scaling_env()
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  bundle <- file.path(root, "paper", "data", "public-derived", "fit-scaling")
  expected <- c("scenarios.csv", "attempts.csv", "summary.csv", "contrasts.csv", "hardware.csv", "README.md")
  expect_true(all(file.exists(file.path(bundle, expected))))
  scenarios <- utils::read.csv(file.path(bundle, "scenarios.csv"), stringsAsFactors = FALSE)
  attempts <- utils::read.csv(file.path(bundle, "attempts.csv"), stringsAsFactors = FALSE)
  summary <- utils::read.csv(file.path(bundle, "summary.csv"), stringsAsFactors = FALSE)
  contrasts <- utils::read.csv(file.path(bundle, "contrasts.csv"), stringsAsFactors = FALSE)
  hardware <- utils::read.csv(file.path(bundle, "hardware.csv"), stringsAsFactors = FALSE)
  expect_true(env$jss_scaling_validate_design(scenarios))
  expect_setequal(unique(attempts$scenario_id), scenarios$scenario_id)
  expect_equal(sum(summary$attempted), sum(!(attempts$warmup %in% TRUE)))
  expect_equal(sum(summary$warmup_attempted), sum(attempts$warmup %in% TRUE))
  expect_equal(sum(summary$failed), sum(!(attempts$retained %in% TRUE) & !(attempts$warmup %in% TRUE)))
  expect_equal(nrow(contrasts), 6L)
  expect_true(all(is.finite(contrasts$difference_mcse_sec)))
  expect_true(all(is.finite(contrasts$ratio_ci_lower) & is.finite(contrasts$ratio_ci_upper)))
  expect_true(env$jss_scaling_publication_status(
    scenarios, attempts, summary, contrasts, env$jss_scaling_precision_rule("full"),
    hardware = hardware
  )$publication_eligible)
  expect_true(all(c("cpu_model", "dependency_versions", "checkpoint_schema_version", "seed_scheme_version", "execution_context") %in% names(hardware)))
  provenance <- env$jss_scaling_provenance_status(hardware, attempts)
  expect_false(provenance$production_eligible)
  expect_match(provenance$reasons, "worktree_not_clean")
  expect_match(provenance$reasons, "not_target_integrated")
})

test_that("publication gate recomputes every reported statistic from exact attempts", {
  env <- local_phase2_scaling_env()
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  bundle <- file.path(root, "paper", "data", "public-derived", "fit-scaling")
  design <- utils::read.csv(file.path(bundle, "scenarios.csv"), stringsAsFactors = FALSE)
  attempts <- utils::read.csv(file.path(bundle, "attempts.csv"), stringsAsFactors = FALSE)
  summary <- utils::read.csv(file.path(bundle, "summary.csv"), stringsAsFactors = FALSE)
  contrasts <- utils::read.csv(file.path(bundle, "contrasts.csv"), stringsAsFactors = FALSE)
  hardware <- utils::read.csv(file.path(bundle, "hardware.csv"), stringsAsFactors = FALSE)
  rule <- env$jss_scaling_precision_rule("full")

  all_failed <- attempts
  timed <- !all_failed$warmup %in% TRUE
  all_failed$converged[timed] <- FALSE
  all_failed$retained[timed] <- FALSE
  all_failed$failure_reason[timed] <- "optimizer_nonconvergence"
  status <- env$jss_scaling_publication_status(
    design, all_failed, summary, contrasts, rule, hardware
  )
  expect_false(status$publication_eligible)
  expect_match(status$reasons, "summary_reconciliation_failed")
  expect_match(status$reasons, "contrast_reconciliation_failed")
  expect_match(status$reasons, "cell_precision_or_success_ineligible")

  fabricated_summary <- summary
  fabricated_summary$attempted[[1L]] <- fabricated_summary$attempted[[1L]] + 1L
  fabricated_summary$median_elapsed_sec[[2L]] <- fabricated_summary$median_elapsed_sec[[2L]] + 10
  status <- env$jss_scaling_publication_status(
    design, attempts, fabricated_summary, contrasts, rule, hardware
  )
  expect_false(status$publication_eligible)
  expect_match(status$reasons, "summary_reconciliation_failed")

  fabricated_contrasts <- contrasts
  fabricated_contrasts$median_runtime_ratio[[1L]] <- fabricated_contrasts$median_runtime_ratio[[1L]] + 1
  fabricated_contrasts$scenario_id[[2L]] <- "fabricated-scenario"
  status <- env$jss_scaling_publication_status(
    design, attempts, summary, fabricated_contrasts, rule, hardware
  )
  expect_false(status$publication_eligible)
  expect_match(status$reasons, "contrast_reconciliation_failed")
  expect_match(status$reasons, "contrast_cells_incomplete")

  duplicated_attempt <- rbind(attempts, attempts[1L, , drop = FALSE])
  duplicated_attempt$execution_order[[nrow(duplicated_attempt)]] <- nrow(duplicated_attempt)
  status <- env$jss_scaling_publication_status(
    design, duplicated_attempt, summary, contrasts, rule, hardware
  )
  expect_false(status$publication_eligible)
  expect_match(status$reasons, "checkpoint_contract_invalid")
})

test_that("clean target hardware cannot relabel dirty standalone attempts", {
  env <- local_phase2_scaling_env()
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  bundle <- file.path(root, "paper", "data", "public-derived", "fit-scaling")
  design <- utils::read.csv(file.path(bundle, "scenarios.csv"), stringsAsFactors = FALSE)
  attempts <- utils::read.csv(file.path(bundle, "attempts.csv"), stringsAsFactors = FALSE)
  summary <- utils::read.csv(file.path(bundle, "summary.csv"), stringsAsFactors = FALSE)
  contrasts <- utils::read.csv(file.path(bundle, "contrasts.csv"), stringsAsFactors = FALSE)
  hardware <- utils::read.csv(file.path(bundle, "hardware.csv"), stringsAsFactors = FALSE)
  substituted <- hardware
  substituted$git_dirty <- "clean"
  substituted$execution_context <- "targets"
  status <- env$jss_scaling_publication_status(
    design, attempts, summary, contrasts, env$jss_scaling_precision_rule("full"),
    hardware = substituted
  )
  expect_false(status$publication_eligible)
  expect_match(status$reasons, "timing_provenance_reconciliation_failed")
  provenance <- env$jss_scaling_provenance_status(substituted, attempts)
  expect_false(provenance$production_eligible)
  expect_match(provenance$reasons, "attempt_hardware_provenance_mismatch")
  expect_match(provenance$reasons, "attempt_worktree_not_clean")
  expect_match(provenance$reasons, "attempt_not_target_integrated")
})
