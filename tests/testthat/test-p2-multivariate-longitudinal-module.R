local_mvt_repo_root <- function() {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  skip_if_not(
    file.exists(file.path(repo_root, "paper", "manifest.csv")),
    "paper replication sources are excluded from source-package checks"
  )
  repo_root
}

local_mvt_env <- function(repo_root = local_mvt_repo_root()) {
  env <- new.env(parent = globalenv())
  source(
    file.path(repo_root, "paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"),
    local = env
  )
  env
}

test_that("multivariate scenario grid expands time, family, dependence, and reps", {
  env <- local_mvt_env()

  grid <- env$mvt_expand_grid(
    time_names = c("t5", "t20"),
    family_names = c("gaussian", "gamma"),
    dependence_names = c("external_exchangeable", "native_covariate_dependent_adjacent"),
    reps = 3L
  )

  expect_equal(nrow(grid), 2L * 2L * 2L * 3L)
  expect_setequal(grid$n_time, c(5L, 20L))
  expect_equal(unique(grid$total_rows[grid$n_time == 5L]), 2000L)
  expect_equal(unique(grid$total_rows[grid$n_time == 20L]), 4000L)
  expect_true(all(grepl("__rep", grid$case_id)))
})

test_that("multivariate scenario grid accepts explicit replicate ids for shards", {
  env <- local_mvt_env()

  grid <- env$mvt_expand_grid(
    time_names = "t20",
    family_names = "gaussian",
    dependence_names = "external_exchangeable",
    reps = c(2L, 5L, 9L)
  )

  expect_equal(grid$rep, c(2L, 5L, 9L))
  expect_true(any(grepl("rep005", grid$case_id, fixed = TRUE)))
})

test_that("all multivariate workflow scripts parse", {
  repo_root <- local_mvt_repo_root()
  script_dir <- file.path(repo_root, "paper", "R", "09-simulation-multivariate-longitudinal")
  files <- list.files(script_dir, pattern = "[.]R$", full.names = TRUE)

  expect_true(length(files) >= 10L)
  for (file in files) {
    expect_silent(parse(file))
  }
})

test_that("multivariate grid supports smoke-size subject override", {
  env <- local_mvt_env()

  withr::local_envvar(GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = "10")
  grid <- env$mvt_expand_grid(
    time_names = "t5",
    family_names = "gaussian",
    dependence_names = "external_exchangeable",
    reps = 1L
  )

  expect_equal(grid$n_subject, 10L)
  expect_equal(grid$total_rows, 50L)
})

test_that("main-grid time defaults separate core and appendix scopes", {
  env <- local_mvt_env()

  expect_equal(env$mvt_default_main_time_names("core"), "t20")
  expect_equal(env$mvt_default_main_time_names("appendix"), c("t5", "t20", "t50"))
  expect_equal(env$mvt_default_main_time_names("all"), c("t5", "t20", "t50"))
})

test_that("comparator filters split standard, GEE, copula, and appendix methods", {
  env <- local_mvt_env()

  withr::local_envvar(GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm,gee,gamCopula,glmm_slope")

  active <- env$mvt_active_comparators()
  expect_setequal(
    active,
    c("glm", "gee_independence", "gee_exchangeable", "gee_ar1", "gee_unstructured", "gamCopula_markov", "glmm_slope")
  )
  expect_equal(env$mvt_standard_comparators(active), "glm")
  expect_setequal(env$mvt_gee_comparators(active), c("independence", "exchangeable", "ar1", "unstructured"))
  expect_false("gamlss.longitudinal" %in% active)
})

test_that("Phase 2 benchmark scope is dated and bounded to two nearest neighbors", {
  env <- local_mvt_env()
  capability <- env$mvt_capability_snapshot()
  scope <- env$mvt_comparator_scope_registry()

  expect_true(all(capability$as_of == "2026-09-01"))
  expect_true(all(c(
    "claim_id", "source_version", "installed_version", "retrieved_on",
    "documentation_topic", "documentation_url", "empirical_evidence_link", "empirical_evidence"
  ) %in% names(capability)))
  expect_true(all(nzchar(capability$documentation_url)))
  expect_true(all(nzchar(capability$empirical_evidence_link)))
  expect_equal(anyDuplicated(capability$claim_id), 0L)
  expect_true(env$mvt_capability_provenance_validation(capability)$valid)
  expect_equal(
    scope$method[scope$headline_empirical_comparator],
    env$mvt_nearest_neighbor_comparators()
  )
  expect_lte(length(env$mvt_nearest_neighbor_comparators()), 2L)
  expect_match(scope$limitation_or_infeasibility[scope$method == "GJRM"], "changes dimension|not equivalent")
})

test_that("capability provenance rejects duplicate claims unknown versions and non-primary URLs", {
  env <- local_mvt_env()
  capability <- env$mvt_capability_snapshot()

  duplicate <- capability
  duplicate$claim_id[2L] <- duplicate$claim_id[1L]
  expect_false(env$mvt_capability_provenance_validation(duplicate)$valid)
  expect_match(
    paste(env$mvt_capability_provenance_validation(duplicate)$problems, collapse = ";"),
    "duplicate claim_id",
    fixed = TRUE
  )

  unknown <- capability
  unknown$source_version[1L] <- "unknown"
  expect_false(env$mvt_capability_provenance_validation(unknown)$valid)

  non_primary <- capability
  non_primary$documentation_url[1L] <- "https://rdrr.io/cran/geepack/man/geeglm.html"
  expect_false(env$mvt_capability_provenance_validation(non_primary)$valid)

  generic_manual <- capability
  generic_manual$documentation_url[1L] <- "https://cran.r-project.org/web/packages/geepack/geepack.pdf"
  expect_false(env$mvt_capability_provenance_validation(generic_manual)$valid)

  nonexistent_official <- capability
  nonexistent_official$documentation_url[1L] <- "https://search.r-project.org/CRAN/refmans/geepack/html/does-not-exist.html"
  expect_false(env$mvt_capability_provenance_validation(nonexistent_official)$valid)

  missing_as_of <- capability[setdiff(names(capability), "as_of")]
  expect_false(env$mvt_capability_provenance_validation(missing_as_of)$valid)
})

test_that("family-specific GEE evidence exposes failures and stress-test status", {
  env <- local_mvt_env()
  methods <- c(
    "gee_independence", "gee_exchangeable", "gee_ar1", "gee_unstructured",
    "gamlss.longitudinal", "gamCopula_markov", "gamCopula_vine_simplified"
  )
  cases <- expand.grid(family = c("gaussian", "gamma"), rep = 1:2, stringsAsFactors = FALSE)
  status <- do.call(rbind, lapply(seq_len(nrow(cases)), function(i) {
    data.frame(
      case_id = paste(cases$family[[i]], cases$rep[[i]], sep = "-"),
      scenario = "external_exchangeable_t20",
      generator = "external",
      dependence = "exchangeable",
      correlation_level = "moderate",
      n_time = 20L,
      n_subject = 20L,
      total_rows = 400L,
      family = cases$family[[i]],
      gamlss_family = if (cases$family[[i]] == "gaussian") "NO" else "GA",
      rep = cases$rep[[i]],
      method = methods,
      success = TRUE,
      status_class = "ok",
      elapsed_sec = seq_along(methods),
      warning = NA_character_,
      error = NA_character_,
      failure_reason_short = "",
      stringsAsFactors = FALSE
    )
  }))
  failed <- status$family == "gamma" & status$rep == 2L & status$method == "gee_unstructured"
  status$success[failed] <- FALSE
  status$status_class[failed] <- "timeout"
  status$error[failed] <- "timed out"
  status$failure_reason_short[failed] <- "timeout"
  benchmark <- status
  benchmark$benchmark_mean_rmse <- 0.1 + seq_len(nrow(benchmark)) / 100
  benchmark$benchmark_neg_log_score <- 1 + seq_len(nrow(benchmark)) / 100
  benchmark$benchmark_mean_rmse[failed] <- NA_real_
  grid <- unique(status[c("case_id", "scenario", "generator", "dependence", "correlation_level", "n_time", "n_subject", "total_rows", "family", "rep")])
  names(grid)[names(grid) == "family"] <- "family_name"

  gee <- env$mvt_phase2_method_summary(status, benchmark, c("gee_exchangeable", "gee_unstructured"))
  contrasts <- env$mvt_nearest_neighbor_paired_contrasts(status, benchmark, scenario_grid = grid)

  expect_setequal(gee$family, c("gaussian", "gamma"))
  expect_true(all(gee$evidence_role[gee$method == "gee_unstructured"] == "high_dimensional_stress_test"))
  gamma_unstructured <- gee[gee$family == "gamma" & gee$method == "gee_unstructured", ]
  expect_equal(gamma_unstructured$n_attempted, 2L)
  expect_equal(gamma_unstructured$n_converged, 1L)
  expect_equal(gamma_unstructured$n_timeout, 1L)
  expect_true(is.finite(gamma_unstructured$success_ci_lower))
  expect_true(nrow(contrasts) > 0L)
  expect_true(all(c(
    "n_expected_pairs", "n_focal_attempted", "n_neighbor_attempted", "n_attempted_pairs",
    "paired_difference_mcse", "paired_difference_ci_lower"
  ) %in% names(contrasts)))
  expect_true(all(contrasts$n_attempted_pairs == contrasts$n_expected_pairs))

  run_dir <- file.path(tempdir(), paste0("mvt-phase2-", sample.int(1e6, 1)))
  dir.create(run_dir, recursive = TRUE)
  env$mvt_write_csv(grid, file.path(run_dir, "scenario_grid.csv"))
  env$mvt_write_csv(status, file.path(run_dir, "fit_status_by_rep.csv"))
  env$mvt_write_csv(benchmark, file.path(run_dir, "benchmark_results_by_rep.csv"))
  env$mvt_write_csv(
    data.frame(name = "active_comparators", value = paste(methods, collapse = ",")),
    file.path(run_dir, "run_metadata.csv")
  )
  expect_error(
    env$mvt_write_phase2_benchmark_evidence(run_dir),
    "valid immutable aggregate snapshot"
  )
})

test_that("failure penalty uses successful rows from the complete planned-method cell", {
  env <- local_mvt_env()
  methods <- c("gamlss.longitudinal", "glm", "gee_exchangeable")
  status <- expand.grid(
    case_id = c("case-1", "case-2"), method = methods,
    stringsAsFactors = FALSE
  )
  status$scenario <- "external_exchangeable_t20"
  status$generator <- "external"
  status$dependence <- "exchangeable"
  status$correlation_level <- "moderate"
  status$n_time <- 20L
  status$n_subject <- 10L
  status$total_rows <- 200L
  status$family <- "gaussian"
  status$gamlss_family <- "NO"
  status$success <- TRUE
  status$status_class <- "ok"
  status$elapsed_sec <- c(1, 2, 3, 4, 5, 6)
  status$warning <- NA_character_
  status$error <- NA_character_
  status$failure_reason_short <- ""
  failed <- status$case_id == "case-2" & status$method == "gamlss.longitudinal"
  status$success[failed] <- FALSE
  status$status_class[failed] <- "error"
  status$elapsed_sec[failed] <- 999

  benchmark <- status
  benchmark$benchmark_mean_rmse <- c(1, 999, 4, 3, 2, 2.5)
  benchmark$benchmark_neg_log_score <- c(1, 999, 4, 3, 2, 2.5)
  summary <- env$mvt_phase2_method_summary(
    status, benchmark, "gamlss.longitudinal", planned_methods = methods
  )

  expect_equal(summary$rmse_failure_penalty, 4.4)
  expect_equal(summary$neg_log_score_failure_penalty, 4.4)
  expect_equal(summary$mean_rmse_failure_penalized, mean(c(1, 4.4)))
  expect_lt(summary$elapsed_failure_penalty, 999)
  expect_equal(summary$rmse_failure_penalty_pool_n, 5L)
  expect_equal(summary$elapsed_failure_penalty_pool_n, 5L)
  expect_identical(summary$failure_penalty_planned_methods, paste(sort(methods), collapse = ","))
})

test_that("Phase 2 reconciliation fails adversarial missing, duplicate, and blank comparisons", {
  env <- local_mvt_env()
  methods <- c("gee_exchangeable", "gamlss.longitudinal", "gamCopula_markov", "gamCopula_vine_simplified")
  grid <- data.frame(
    case_id = c("case-1", "case-2"), scenario = "s", generator = "external",
    dependence = "exchangeable", correlation_level = "moderate", n_time = 20L,
    n_subject = 10L, total_rows = 200L, family_name = "gaussian", rep = 1:2,
    stringsAsFactors = FALSE
  )
  expected <- merge(grid["case_id"], data.frame(method = methods), by = NULL)
  status <- merge(expected, transform(grid, family = family_name)[setdiff(names(transform(grid, family = family_name)), "case_id")], by = NULL)
  status <- status[!duplicated(paste(status$case_id, status$method)), , drop = FALSE]
  status$success <- TRUE
  status$status_class <- "ok"
  status$elapsed_sec <- 1
  status$warning <- NA_character_
  status$error <- NA_character_
  status$failure_reason_short <- ""
  benchmark <- status
  benchmark$benchmark_mean_rmse <- 1
  benchmark$benchmark_neg_log_score <- 1
  bad_status <- status[status$method != "gamCopula_vine_simplified", , drop = FALSE]
  bad_status <- rbind(bad_status, bad_status[1, , drop = FALSE])
  reconciliation <- env$mvt_phase2_attempt_reconciliation(grid, bad_status, benchmark, methods)

  expect_true(any(reconciliation$status_rows != 1L))
  run_dir <- file.path(tempdir(), paste0("mvt-phase2-adversarial-", sample.int(1e6, 1)))
  dir.create(run_dir, recursive = TRUE)
  env$mvt_write_csv(grid, file.path(run_dir, "scenario_grid.csv"))
  env$mvt_write_csv(bad_status, file.path(run_dir, "fit_status_by_rep.csv"))
  env$mvt_write_csv(benchmark, file.path(run_dir, "benchmark_results_by_rep.csv"))
  env$mvt_write_csv(data.frame(name = "active_comparators", value = paste(methods, collapse = ",")), file.path(run_dir, "run_metadata.csv"))
  expect_error(
    env$mvt_write_phase2_benchmark_evidence(run_dir),
    "valid immutable aggregate snapshot"
  )
})

test_that("external correlation matrices are positive definite or repaired", {
  env <- local_mvt_env()
  deps <- env$mvt_dependence_specs(include_appendix = TRUE)

  for (dep_name in c("external_exchangeable", "external_ar1", "external_block")) {
    dep <- env$mvt_resize_dependence(deps[[dep_name]], 20L)
    R <- env$mvt_external_correlation_matrix(20L, dep)
    expect_equal(dim(R), c(20L, 20L))
    expect_true(isSymmetric(R))
    expect_true(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) > 0)
    expect_equal(diag(R), rep(1, 20L))
  }
})

test_that("native simulator emits required truth columns", {
  env <- local_mvt_env()
  skip_if_not_installed("gamlss.dist")

  withr::local_envvar(GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = "8")
  grid <- env$mvt_expand_grid(
    time_names = "t5",
    family_names = "gaussian",
    dependence_names = "native_covariate_dependent_adjacent",
    reps = 1L
  )
  dat <- env$mvt_simulate_case(grid[1, , drop = FALSE], seed = 123)

  expect_equal(nrow(dat), 40L)
  expect_true(all(c("true_mu", "true_beta_time", "true_theta", "true_tau", "u") %in% names(dat)))
  expect_true(any(is.finite(dat$true_theta)))
})

test_that("gamCopula runtime dependency check fails clearly when unavailable", {
  env <- local_mvt_env()

  if (requireNamespace("gamCopula", quietly = TRUE)) {
    expect_s3_class(env$mvt_require_namespaces("gamCopula", strict = TRUE), "data.frame")
  } else {
    expect_error(env$mvt_require_namespaces("gamCopula", strict = TRUE), "gamCopula")
  }
})

test_that("row binder tolerates failed and partial result rows", {
  env <- local_mvt_env()

  a <- data.frame(method = "a", success = TRUE, metric_a = 1)
  b <- data.frame(method = "b", error = "failed")
  out <- env$mvt_bind_rows_fill(list(a, b))
  out_mixed <- env$mvt_bind_rows_fill(a, list(b))

  expect_equal(nrow(out), 2L)
  expect_equal(nrow(out_mixed), 2L)
  expect_true(all(c("method", "success", "metric_a", "error") %in% names(out)))
  expect_true(is.na(out$metric_a[out$method == "b"]))
})

test_that("fit status rows classify all benchmark result rows", {
  env <- local_mvt_env()

  results <- data.frame(
    case_id = c("case1", "case1", "case1"),
    method = c("glm", "gee_unstructured", "gamlss.longitudinal"),
    success = c(TRUE, FALSE, TRUE),
    status_class = c("", NA, "warning"),
    warning = c("", "", "not converged"),
    error = c("", "timed out after 1 sec", ""),
    elapsed_sec = c(0.1, 1, 2),
    stringsAsFactors = FALSE
  )

  status <- env$mvt_status_from_results(results)

  expect_equal(status$status_class[status$method == "glm"], "ok")
  expect_equal(status$status_class[status$method == "gee_unstructured"], "timeout")
  expect_equal(status$status_class[status$method == "gamlss.longitudinal"], "warning")
})

test_that("distribution metrics are finite at the known Gaussian truth", {
  env <- local_mvt_env()
  skip_if_not_installed("gamlss.dist")

  spec <- env$mvt_family_specs()[["gaussian"]]
  dat <- env$mvt_design_long(n_subject = 6L, n_time = 5L)
  dat <- env$mvt_add_truth_columns(dat, spec)
  dat$response <- dat$true_mu

  metrics <- env$mvt_distribution_metrics(dat, env$mvt_truth_params(dat), spec)

  expect_equal(unname(metrics[["benchmark_mean_rmse"]]), 0)
  expect_equal(unname(metrics[["benchmark_q90_mae"]]), 0)
  expect_true(is.finite(metrics[["benchmark_neg_log_score"]]))
  expect_true(is.finite(metrics[["benchmark_interval_coverage_95"]]))
})

test_that("random-slope GLMM appendix sensitivity emits status and coefficient rows", {
  env <- local_mvt_env()
  skip_if_not_installed("gamlss.dist")
  skip_if_not_installed("mvtnorm")
  skip_if_not_installed("lme4")

  withr::local_envvar(
    GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = "8",
    GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glmm_slope"
  )
  grid <- env$mvt_expand_grid(
    time_names = "t5",
    family_names = "gaussian",
    dependence_names = "external_exchangeable",
    reps = 1L
  )

  out <- env$mvt_run_case(grid[1, , drop = FALSE], seed_base = 123, require_gamcopula = FALSE)

  expect_true("glmm_slope" %in% out$fit_status$method)
  expect_true("glmm_slope" %in% out$benchmark_results$method)
  if (isTRUE(out$fit_status$retained[out$fit_status$method == "glmm_slope"])) {
    expect_true(all(c("intercept", "time", "x", "z") %in% out$coefficient_results$term))
  } else {
    expect_equal(nrow(out$coefficient_results), 0L)
    expect_false(any(is.finite(out$benchmark_results$benchmark_mean_rmse)))
  }
})

test_that("resume helpers read existing outputs and identify completed cases", {
  env <- local_mvt_env()
  run_dir <- file.path(tempdir(), paste0("mvt-resume-", sample.int(1e6, 1)))
  dir.create(run_dir, recursive = TRUE)

  env$mvt_write_csv(
    data.frame(case_id = "case-a", method = "glm", success = TRUE),
    file.path(run_dir, "benchmark_results_by_rep.csv")
  )

  existing <- env$mvt_read_existing_results(run_dir)
  completed <- env$mvt_completed_case_ids(existing)

  expect_true("benchmark_results" %in% names(existing))
  expect_equal(completed, "case-a")
})

test_that("resume requires one attempt for every planned method", {
  env <- local_mvt_env()
  existing <- stats::setNames(
    replicate(length(env$mvt_result_names()), data.frame(), simplify = FALSE),
    env$mvt_result_names()
  )
  planned <- c("glm", "gee_unstructured")
  partial <- data.frame(
    case_id = "case-a",
    method = "glm",
    success = TRUE,
    stringsAsFactors = FALSE
  )
  existing$fit_status <- partial
  existing$benchmark_results <- partial

  expect_length(env$mvt_completed_case_ids(existing, planned_methods = planned), 0L)

  complete <- rbind(
    partial,
    transform(partial, method = "gee_unstructured")
  )
  existing$fit_status <- complete
  existing$benchmark_results <- complete
  expect_equal(
    env$mvt_completed_case_ids(existing, planned_methods = planned),
    "case-a"
  )

  existing$benchmark_results <- rbind(complete, complete[2L, , drop = FALSE])
  expect_length(env$mvt_completed_case_ids(existing, planned_methods = planned), 0L)
})

test_that("case-checkpoint runner is serial-parallel equivalent and canonically ordered", {
  env <- local_mvt_env()
  withr::local_envvar(c(
    GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm,gamlss.longitudinal",
    GAMLSS_LONGITUDINAL_MVT_RESUME = "true"
  ))
  grid <- env$mvt_expand_grid(
    time_names = "t5",
    family_names = c("gaussian", "gamma"),
    dependence_names = "external_exchangeable",
    reps = 1:2
  )
  mock_run_case <- function(row, seed_base, require_gamcopula, case_seed) {
    methods <- if (row$rep[[1L]] %% 2L) c("gamlss.longitudinal", "glm") else c("glm", "gamlss.longitudinal")
    value <- seed_base + sum(utf8ToInt(as.character(row$case_id[[1L]])))
    benchmark <- data.frame(
      case_id = row$case_id[[1L]], scenario = row$scenario[[1L]],
      family = row$family_name[[1L]], rep = row$rep[[1L]], method = methods,
      attempted = TRUE, success = TRUE, converged = TRUE, retained = TRUE,
      stop_reason = "converged", status_class = "success", elapsed_sec = 0,
      failure_reason_short = "none", error = NA_character_,
      metric = value + seq_along(methods), stringsAsFactors = FALSE
    )
    list(
      fit_status = benchmark,
      benchmark_results = benchmark,
      coefficient_results = data.frame(),
      dependence_recovery = data.frame(),
      variogram_scores = data.frame(),
      runtime = benchmark
    )
  }
  serial_dir <- tempfile("mvt-serial-")
  parallel_dir <- tempfile("mvt-parallel-")
  fingerprints <- list(
    producer_fingerprint = "test-producer",
    code_fingerprint = "test-code",
    package_fingerprint = "test-packages"
  )
  env$mvt_run_grid(
    grid, serial_dir, seed_base = 811L, require_gamcopula = FALSE,
    workers = 1L, run_case_fun = mock_run_case, fingerprints = fingerprints
  )
  env$mvt_run_grid(
    grid, parallel_dir, seed_base = 811L, require_gamcopula = FALSE,
    workers = 2L, run_case_fun = mock_run_case, worker_load_local = FALSE,
    fingerprints = fingerprints
  )

  for (name in env$mvt_result_names()) {
    serial <- readLines(file.path(serial_dir, paste0(name, "_by_rep.csv")), warn = FALSE)
    parallel <- readLines(file.path(parallel_dir, paste0(name, "_by_rep.csv")), warn = FALSE)
    expect_identical(parallel, serial, info = name)
  }
  benchmark <- utils::read.csv(file.path(parallel_dir, "benchmark_results_by_rep.csv"), stringsAsFactors = FALSE)
  expected <- benchmark[order(
    match(benchmark$scenario, unique(grid$scenario)), benchmark$rep,
    match(benchmark$method, c("glm", "gamlss.longitudinal")),
    match(benchmark$family, c("gaussian", "gamma"))
  ), , drop = FALSE]
  rownames(expected) <- NULL
  expect_identical(benchmark, expected)
})

test_that("valid case checkpoints resume without duplicate execution", {
  env <- local_mvt_env()
  withr::local_envvar(c(
    GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm",
    GAMLSS_LONGITUDINAL_MVT_RESUME = "true"
  ))
  grid <- env$mvt_expand_grid(
    time_names = "t5", family_names = "gaussian",
    dependence_names = "external_exchangeable", reps = 1:2
  )
  run_dir <- tempfile("mvt-resume-case-")
  marker <- tempfile("mvt-executions-")
  fingerprints <- list(
    producer_fingerprint = "test-producer",
    code_fingerprint = "test-code",
    package_fingerprint = "test-packages"
  )
  mock_run_case <- function(row, seed_base, require_gamcopula, case_seed) {
    cat(as.character(row$case_id[[1L]]), "\n", file = marker, append = TRUE)
    result <- data.frame(
      case_id = row$case_id[[1L]], scenario = row$scenario[[1L]],
      family = row$family_name[[1L]], rep = row$rep[[1L]], method = "glm",
      attempted = TRUE, success = TRUE, converged = TRUE, retained = TRUE,
      stop_reason = "converged", status_class = "success", elapsed_sec = 0,
      failure_reason_short = "none", error = NA_character_, stringsAsFactors = FALSE
    )
    list(
      fit_status = result, benchmark_results = result,
      coefficient_results = data.frame(), dependence_recovery = data.frame(),
      variogram_scores = data.frame(), runtime = result
    )
  }
  env$mvt_run_grid(
    grid, run_dir, seed_base = 812L, require_gamcopula = FALSE,
    workers = 1L, run_case_fun = mock_run_case, fingerprints = fingerprints
  )
  expect_length(readLines(marker), 2L)
  env$mvt_run_grid(
    grid, run_dir, seed_base = 812L, require_gamcopula = FALSE,
    workers = 2L, run_case_fun = function(...) stop("resume reran a completed case"),
    worker_load_local = FALSE, fingerprints = fingerprints
  )
  expect_length(readLines(marker), 2L)
  expect_equal(length(list.files(env$mvt_case_checkpoint_dir(run_dir), pattern = "[.]rds$")), 2L)
  metadata <- utils::read.csv(file.path(run_dir, "run_metadata.csv"), stringsAsFactors = FALSE)
  expect_equal(metadata$value[metadata$name == "workers_requested"], "2")
  expect_equal(metadata$value[metadata$name == "workers_used"], "0")
  checkpoint <- readRDS(env$mvt_case_checkpoint_path(run_dir, grid$case_id[[1L]]))
  expect_equal(checkpoint$metadata$workers_requested, 1L)
})

test_that("case checkpoints reject stale producer and numerical configuration", {
  env <- local_mvt_env()
  withr::local_envvar(GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm")
  grid <- env$mvt_expand_grid(
    time_names = "t5", family_names = "gaussian",
    dependence_names = "external_exchangeable", reps = 1L
  )
  fingerprints <- list(
    producer_fingerprint = "producer-a",
    code_fingerprint = "code-a",
    package_fingerprint = "packages-a"
  )
  task <- env$mvt_prepare_tasks(grid, 813L, FALSE, fingerprints)[[1L]]
  result <- data.frame(
    case_id = task$case_id, scenario = task$row$scenario,
    family = task$row$family_name, rep = task$row$rep,
    method = "glm", attempted = TRUE, success = TRUE, converged = TRUE, retained = TRUE,
    stop_reason = "converged", status_class = "success", elapsed_sec = 0,
    failure_reason_short = "none", error = NA_character_, stringsAsFactors = FALSE
  )
  payload <- list(
    fit_status = result, benchmark_results = result,
    coefficient_results = data.frame(), dependence_recovery = data.frame(),
    variogram_scores = data.frame(), runtime = result
  )
  path <- tempfile(fileext = ".rds")
  env$mvt_write_case_checkpoint_atomic(payload, task, path)
  expect_true(env$mvt_case_checkpoint_valid(readRDS(path), task))

  stale_code <- task
  stale_code$producer_fingerprint <- "producer-b"
  expect_false(env$mvt_case_checkpoint_valid(readRDS(path), stale_code))
  stale_package_code <- task
  stale_package_code$code_fingerprint <- "code-b"
  expect_false(env$mvt_case_checkpoint_valid(readRDS(path), stale_package_code))
  stale_dependencies <- task
  stale_dependencies$package_fingerprint <- "packages-b"
  expect_false(env$mvt_case_checkpoint_valid(readRDS(path), stale_dependencies))
  stale_runtime <- task
  stale_runtime$execution_fingerprint <- "runtime-b"
  expect_false(env$mvt_case_checkpoint_valid(readRDS(path), stale_runtime))

  withr::local_envvar(GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM = "999")
  stale_configuration <- env$mvt_prepare_tasks(grid, 813L, FALSE, fingerprints)[[1L]]
  expect_false(env$mvt_case_checkpoint_valid(readRDS(path), stale_configuration))
})

test_that("production multivariate grid has exact case and method cardinalities", {
  env <- local_mvt_env()
  withr::local_envvar(GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = NA)
  grid <- env$mvt_expand_grid(
    time_names = "t20",
    family_names = c("gaussian", "poisson", "gamma", "binomial"),
    dependence_names = c(
      "external_exchangeable", "external_ar1",
      "native_time_varying_adjacent", "native_covariate_dependent_adjacent"
    ),
    reps = 100L
  )
  expect_equal(nrow(grid), 1600L)
  expect_equal(anyDuplicated(grid$case_id), 0L)
  contract <- env$mvt_phase2_production_contract()
  expect_identical(contract$time_names, "t20")
  expect_identical(sort(unique(grid$time_name)), "t20")
  expect_setequal(unique(grid$family_name), contract$families)
  expect_setequal(unique(grid$dependence_name), contract$dependence)
  expect_setequal(unique(grid$rep), contract$reps)
  expect_identical(env$mvt_default_comparators(), contract$methods)
  expect_length(contract$methods, 9L)
  expect_equal(nrow(grid) * length(contract$methods), 14400L)
  expect_true(env$mvt_grid_is_exact_production(grid))
  forged_grid <- grid
  forged_grid$case_id[[1L]] <- "forged-but-unique"
  expect_false(env$mvt_grid_is_exact_production(forged_grid))

  tasks <- env$mvt_prepare_tasks(
    grid, seed_base = 20260818L, require_gamcopula = TRUE,
    fingerprints = list(
      algorithm = env$mvt_identity_algorithm,
      identity_version = env$mvt_identity_version,
      producer_fingerprint = "producer", code_fingerprint = "code",
      package_fingerprint = "packages"
    )
  )
  seeds <- vapply(tasks, `[[`, integer(1L), "case_seed")
  expected <- vapply(seq_len(nrow(grid)), function(i) {
    env$mvt_case_seed(grid[i, , drop = FALSE], 20260818L)
  }, integer(1L))
  expect_identical(seeds, expected)
  expect_equal(anyDuplicated(seeds), 0L)
  expect_length(seeds, 1600L)
})

test_that("production eligibility rejects any filtered Phase 2 subset", {
  env <- local_mvt_env()
  grid <- env$mvt_expand_grid(
    time_names = "t20",
    family_names = c("gaussian", "poisson", "gamma", "binomial"),
    dependence_names = c(
      "external_exchangeable", "external_ar1",
      "native_time_varying_adjacent", "native_covariate_dependent_adjacent"
    ),
    reps = 1:100
  )
  methods <- env$mvt_default_comparators()
  pairs <- merge(
    grid[c("case_id")], data.frame(method = methods, stringsAsFactors = FALSE),
    by = NULL, sort = FALSE
  )
  reconciliation <- transform(
    pairs, planned = TRUE, reconciled = TRUE, unexpected = FALSE,
    status_rows = 1L, benchmark_rows = 1L
  )
  audit_for <- function(use_grid, use_pairs, use_reconciliation) {
    env$mvt_phase2_benchmark_audit(
      env$mvt_capability_snapshot(), env$mvt_comparator_scope_registry(),
      data.frame(), data.frame(), data.frame(),
      grid = use_grid, status = use_pairs, benchmark = use_pairs,
      planned_methods = methods, reconciliation = use_reconciliation
    )
  }
  full <- audit_for(grid, pairs, reconciliation)
  expect_false(env$mvt_phase2_production_eligible(full))
  strict <- grep("^production_", full$check, value = TRUE)
  cardinality <- setdiff(strict, c("production_registered_timeouts_finite", "production_subprocess_attestations_complete"))
  expect_true(all(full$status[match(cardinality, full$check)] == "pass"))
  registry <- env$mvt_phase2_audit_registry()
  forged_eight <- data.frame(
    check = registry[1:8], status = "pass", detail = "forged partial audit",
    stringsAsFactors = FALSE
  )
  expect_false(env$mvt_phase2_production_eligible(forged_eight))
  complete_pass <- data.frame(check = registry, status = "pass", detail = "fixture", stringsAsFactors = FALSE)
  expect_true(env$mvt_phase2_production_eligible(complete_pass))
  duplicate <- complete_pass; duplicate$check[[length(registry)]] <- duplicate$check[[1L]]
  expect_false(env$mvt_phase2_production_eligible(duplicate))
  unknown <- complete_pass; unknown$check[[1L]] <- "unknown_check"
  expect_false(env$mvt_phase2_production_eligible(unknown))
  reordered <- complete_pass[rev(seq_len(nrow(complete_pass))), , drop = FALSE]
  expect_false(env$mvt_phase2_production_eligible(reordered))
  failed_complete <- complete_pass; failed_complete$status[[length(registry)]] <- "fail"
  expect_false(env$mvt_phase2_production_eligible(failed_complete))

  removed <- grid$case_id[[1L]]
  filtered_grid <- grid[grid$case_id != removed, , drop = FALSE]
  filtered_pairs <- pairs[pairs$case_id != removed, , drop = FALSE]
  filtered_reconciliation <- reconciliation[reconciliation$case_id != removed, , drop = FALSE]
  filtered <- audit_for(filtered_grid, filtered_pairs, filtered_reconciliation)
  expect_false(env$mvt_phase2_production_eligible(filtered))
  expect_equal(
    filtered$status[filtered$check == "production_attempt_cardinality_exact"],
    "fail"
  )
})

test_that("public integration contract exposes only allowlisted Phase 2 artifacts", {
  env <- local_mvt_env()
  allowlist <- env$mvt_phase2_public_output_allowlist()
  expect_identical(allowlist$attempt_artifacts, c(
    "scenario_grid.csv", "run_metadata.csv", "package_versions.csv",
    "worker_attestations.csv", "checkpoint_rejections.csv",
    paste0(env$mvt_result_names(), "_by_rep.csv")
  ))
  expect_identical(allowlist$evidence_artifacts, c(
    "gee_family_results.csv", "gee_unstructured_stress_test.csv",
    "nearest_neighbor_results.csv", "nearest_neighbor_paired_contrasts.csv",
    "capability_snapshot_2026-09-01.csv", "comparator_scope_registry.csv",
    "phase2_attempt_reconciliation.csv", "phase2_benchmark_audit.csv",
    "phase2_benchmark_audit.md"
  ))
  expect_false(any(grepl("checkpoint|sensitivity|summary", unlist(allowlist), ignore.case = TRUE) &
    !unlist(allowlist) %in% "checkpoint_rejections.csv"))
})

test_that("checkpoint validation enforces full attempt semantics and excludes nonconvergence", {
  env <- local_mvt_env()
  withr::local_envvar(GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm")
  grid <- env$mvt_expand_grid(
    time_names = "t5", family_names = "gaussian",
    dependence_names = "external_exchangeable", reps = 1L
  )
  task <- env$mvt_prepare_tasks(
    grid, 814L, FALSE,
    list(producer_fingerprint = "producer", code_fingerprint = "code", package_fingerprint = "packages")
  )[[1L]]
  forged <- data.frame(case_id = task$case_id, method = "glm", stringsAsFactors = FALSE)
  forged_payload <- stats::setNames(
    replicate(length(env$mvt_result_names()), data.frame(), simplify = FALSE),
    env$mvt_result_names()
  )
  forged_payload$fit_status <- forged
  forged_payload$benchmark_results <- forged
  forged_payload$runtime <- forged
  expect_false(env$mvt_checkpoint_result_valid(forged_payload, task))

  nonconverged <- data.frame(
    case_id = task$case_id, method = "glm", attempted = TRUE,
    success = FALSE, converged = FALSE, retained = FALSE,
    stop_reason = "optimizer_nonconvergence", status_class = "error",
    elapsed_sec = 0.1, failure_reason_short = "optimizer_nonconvergence",
    error = "optimizer did not converge", benchmark_mean_rmse = 1,
    stringsAsFactors = FALSE
  )
  payload <- forged_payload
  payload$fit_status <- nonconverged
  payload$benchmark_results <- nonconverged
  payload$runtime <- nonconverged
  expect_false(env$mvt_checkpoint_result_valid(payload, task))
  payload$fit_status$benchmark_mean_rmse <- NA_real_
  payload$benchmark_results$benchmark_mean_rmse <- NA_real_
  payload$runtime$benchmark_mean_rmse <- NA_real_
  expect_true(env$mvt_checkpoint_result_valid(payload, task))
})

test_that("a concurrent process lease prevents duplicate run invocation", {
  env <- local_mvt_env()
  repo_root <- local_mvt_repo_root()
  setup_path <- file.path(repo_root, "paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R")
  run_dir <- tempfile("mvt-concurrent-")
  dir.create(run_dir)
  ready <- tempfile("mvt-lock-ready-")
  release <- tempfile("mvt-lock-release-")
  process <- callr::r_bg(function(setup_path, run_dir, ready, release) {
    source(setup_path, local = .GlobalEnv)
    lock <- mvt_acquire_run_lock(run_dir)
    on.exit(mvt_release_run_lock(lock), add = TRUE)
    writeLines("ready", ready)
    deadline <- Sys.time() + 10
    while (!file.exists(release) && Sys.time() < deadline) Sys.sleep(0.05)
  }, args = list(setup_path, run_dir, ready, release), supervise = TRUE)
  withr::defer({
    if (process$is_alive()) process$kill()
    unlink(c(ready, release), force = TRUE)
  }, teardown_env())
  deadline <- Sys.time() + 5
  while (!file.exists(ready) && process$is_alive() && Sys.time() < deadline) Sys.sleep(0.05)
  expect_true(file.exists(ready))
  expect_error(env$mvt_acquire_run_lock(run_dir), "already leased")
  writeLines("release", release)
  process$wait(timeout = 5000)
  expect_false(process$is_alive())
  expect_false(dir.exists(env$mvt_run_lock_path(run_dir)))
})

test_that("aggregate replacement is atomic, rollback-safe, and parent-owned", {
  env <- local_mvt_env()
  run_dir <- tempfile("mvt-aggregate-"); dir.create(run_dir)
  path <- file.path(run_dir, "aggregate.csv")
  lease <- env$mvt_acquire_run_lock(run_dir)
  withr::defer(env$mvt_release_run_lock(lease), teardown_env())
  env$mvt_write_csv_atomic(data.frame(value = 1), path, lease)
  expect_error(env$mvt_replace_file_atomic(paste0(path, ".missing"), path), "does not exist")
  expect_equal(utils::read.csv(path)$value, 1)
  env$mvt_write_csv_atomic(data.frame(value = 2), path, lease)
  expect_equal(utils::read.csv(path)$value, 2)
  expect_false(file.exists(paste0(path, ".previous")))
  owner <- readRDS(paste0(path, ".ownership.rds"))
  expect_identical(owner$owner_role, "lease_parent")
  expect_identical(owner$writer_pid, Sys.getpid())
  code <- paste(deparse(body(env$mvt_write_csv_atomic)), collapse = "\n")
  expect_false(grepl("file.remove(path)", code, fixed = TRUE))
  expect_false(grepl("file.rename(path, backup)", code, fixed = TRUE))
})

test_that("aggregate lease identity cannot be forged or reused by a child", {
  env <- local_mvt_env()
  setup_path <- file.path(local_mvt_repo_root(), "paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R")
  run_dir <- tempfile("mvt-lease-adversarial-"); dir.create(run_dir)
  path <- file.path(run_dir, "aggregate.csv")
  lease <- env$mvt_acquire_run_lock(run_dir)
  withr::defer(env$mvt_release_run_lock(lease), teardown_env())
  env$mvt_write_csv_atomic(data.frame(value = 1L), path, lease)

  forged <- lease
  forged$token <- new.env(parent = emptyenv())
  expect_error(env$mvt_write_csv_atomic(data.frame(value = 2L), path, forged), "forged")
  expect_error(
    env$mvt_write_csv_atomic(data.frame(value = 2L), path, lease, owner_role = "parent"),
    "unused argument"
  )
  child <- callr::r(function(setup_path, path, lease) {
    source(setup_path, local = .GlobalEnv)
    tryCatch({
      mvt_write_csv_atomic(data.frame(value = 99L), path, lease)
      "overwrote"
    }, error = function(e) conditionMessage(e))
  }, args = list(setup_path = setup_path, path = path, lease = lease))
  expect_match(child, "inactive, forged, or owned by another process")
  expect_equal(utils::read.csv(path)$value, 1L)
})

test_that("full execution identity rejects dependency, runtime, and configuration mutations", {
  env <- local_mvt_env()
  expected <- env$mvt_execution_attestation_contract(configuration_fingerprint = "registered-config")
  actual <- env$mvt_verify_execution_attestation(expected, load = TRUE)
  expect_true(env$mvt_execution_attestation_matches(actual, expected))

  bad_dependency <- expected
  bad_dependency$dependency_identity$version[[1L]] <- "99.0.0"
  expect_error(env$mvt_verify_execution_attestation(bad_dependency, load = FALSE), "Full dependency")
  bad_runtime <- expected
  bad_runtime$runtime_identity$rng_kind[[1L]] <- "Super-Duper"
  expect_error(env$mvt_verify_execution_attestation(bad_runtime, load = FALSE), "Full dependency")
  bad_config <- expected
  bad_config$configuration_fingerprint <- ""
  expect_error(env$mvt_verify_execution_attestation(bad_config, load = FALSE), "Full dependency")

  mutated_actual <- actual
  mutated_actual$dependency_identity$namespace_path[[1L]] <- "C:/stale/library"
  expect_false(env$mvt_execution_attestation_matches(mutated_actual, expected))
})

test_that("production timeout and aggregate semantic gates reject adversarial inputs", {
  env <- local_mvt_env()
  withr::local_envvar(GAMLSS_LONGITUDINAL_MVT_GEE_TIMEOUT_SEC = "Inf")
  configuration <- env$mvt_checkpoint_configuration(
    seed_base = 1L, require_gamcopula = TRUE,
    active_comparators = env$mvt_default_comparators()
  )
  expect_false(env$mvt_timeout_contract(configuration, env$mvt_default_comparators())$valid)

  grid <- data.frame(case_id = "case-a", scenario = "s", family_name = "gaussian", rep = 1L, stringsAsFactors = FALSE)
  attempt <- data.frame(
    case_id = "case-a", method = "glm", attempted = TRUE, success = TRUE,
    converged = TRUE, retained = TRUE, stop_reason = "converged",
    status_class = "success", elapsed_sec = 0.1, failure_reason_short = "none",
    error = NA_character_, stringsAsFactors = FALSE
  )
  tables <- list(fit_status = attempt, benchmark_results = attempt, runtime = attempt,
    coefficient_results = data.frame(), dependence_recovery = data.frame(), variogram_scores = data.frame())
  expect_true(env$mvt_validate_aggregate_table_semantics(grid, tables, "glm"))
  tables$benchmark_results$converged <- FALSE
  expect_error(env$mvt_validate_aggregate_table_semantics(grid, tables, "glm"), "contradictory|disagree")
})

test_that("claim integration contract resolves exact keyed effects and uncertainty", {
  env <- local_mvt_env()
  evidence_dir <- tempfile("mvt-claim-contract-"); dir.create(evidence_dir)
  table <- data.frame(
    scenario_key = "s1", row_key = "s1|neighbor|rmse", metric = "rmse",
    effect = 0.2, n_pairs = 100L, mcse = 0.03, lower = 0.14, upper = 0.26,
    stringsAsFactors = FALSE
  )
  utils::write.csv(table, file.path(evidence_dir, "nearest_neighbor_paired_contrasts.csv"), row.names = FALSE)
  claims <- data.frame(
    claim_id = "claim-1", scenario_key = "s1", row_key = "s1|neighbor|rmse",
    metric = "rmse", direction = "positive", denominator = "n_pairs",
    effect_artifact = "nearest_neighbor_paired_contrasts.csv", effect_column = "effect",
    mcse_column = "mcse", ci_lower_column = "lower", ci_upper_column = "upper",
    interval_support = "supported", wording_strength = "directional",
    stringsAsFactors = FALSE
  )
  expect_true(env$mvt_validate_phase2_claim_rows(claims, evidence_dir))
  bad <- claims; bad$denominator <- ""
  expect_error(env$mvt_validate_phase2_claim_rows(bad, evidence_dir), "blank")
  bad <- claims; bad$row_key <- "forged"
  expect_error(env$mvt_validate_phase2_claim_rows(bad, evidence_dir), "exactly one")
  bad <- claims; bad$effect_artifact <- "benchmark_summary.csv"
  expect_error(env$mvt_validate_phase2_claim_rows(bad, evidence_dir), "not allowlisted")
  bad_table <- table; bad_table$effect <- -0.2; bad_table$lower <- -0.3; bad_table$upper <- -0.1
  utils::write.csv(bad_table, file.path(evidence_dir, "nearest_neighbor_paired_contrasts.csv"), row.names = FALSE)
  expect_error(env$mvt_validate_phase2_claim_rows(claims, evidence_dir), "effect sign")
  bad_table <- table; bad_table$mcse <- -0.01
  utils::write.csv(bad_table, file.path(evidence_dir, "nearest_neighbor_paired_contrasts.csv"), row.names = FALSE)
  expect_error(env$mvt_validate_phase2_claim_rows(claims, evidence_dir), "MCSE")
  bad_table <- table; bad_table$lower <- 0.3; bad_table$upper <- 0.1
  utils::write.csv(bad_table, file.path(evidence_dir, "nearest_neighbor_paired_contrasts.csv"), row.names = FALSE)
  expect_error(env$mvt_validate_phase2_claim_rows(claims, evidence_dir), "inverted")
  bad_table <- table; bad_table$lower <- 0.21; bad_table$upper <- 0.3
  utils::write.csv(bad_table, file.path(evidence_dir, "nearest_neighbor_paired_contrasts.csv"), row.names = FALSE)
  expect_error(env$mvt_validate_phase2_claim_rows(claims, evidence_dir), "does not contain")
  bad_table <- table; bad_table$lower <- -0.1
  utils::write.csv(bad_table, file.path(evidence_dir, "nearest_neighbor_paired_contrasts.csv"), row.names = FALSE)
  expect_error(env$mvt_validate_phase2_claim_rows(claims, evidence_dir), "interval-support")
  for (wording in c("exact", "directional", "cautious", "descriptive")) {
    crossing <- claims
    crossing$interval_support <- "not_supported"
    crossing$wording_strength <- wording
    expect_error(
      env$mvt_validate_phase2_claim_rows(crossing, evidence_dir),
      "every directional expected direction requires", info = wording
    )
  }
  expect_error(env$mvt_validate_phase2_claim_outputs(claims, evidence_dir), "trust registry|snapshot")
})

test_that("mutated checkpoints are quarantined with a rejection ledger and rerun once", {
  env <- local_mvt_env()
  withr::local_envvar(c(
    GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm",
    GAMLSS_LONGITUDINAL_MVT_RESUME = "true"
  ))
  grid <- env$mvt_expand_grid(
    time_names = "t5", family_names = "gaussian",
    dependence_names = "external_exchangeable", reps = 1L
  )
  run_dir <- tempfile("mvt-quarantine-")
  marker <- tempfile("mvt-quarantine-runs-")
  fingerprints <- list(
    producer_fingerprint = "producer", code_fingerprint = "code",
    package_fingerprint = "packages"
  )
  runner <- function(row, seed_base, require_gamcopula, case_seed) {
    cat(case_seed, "\n", file = marker, append = TRUE)
    attempt <- data.frame(
      case_id = row$case_id, scenario = row$scenario, family = row$family_name,
      rep = row$rep, method = "glm", attempted = TRUE, success = TRUE,
      converged = TRUE, retained = TRUE, stop_reason = "converged",
      status_class = "success", elapsed_sec = 0,
      failure_reason_short = "none", error = NA_character_, stringsAsFactors = FALSE
    )
    list(
      fit_status = attempt, benchmark_results = attempt,
      coefficient_results = data.frame(), dependence_recovery = data.frame(),
      variogram_scores = data.frame(), runtime = attempt
    )
  }
  env$mvt_run_grid(
    grid, run_dir, seed_base = 815L, require_gamcopula = FALSE,
    workers = 1L, run_case_fun = runner, fingerprints = fingerprints
  )
  checkpoint_path <- env$mvt_case_checkpoint_path(run_dir, grid$case_id)
  mutated <- readRDS(checkpoint_path)
  mutated$result$benchmark_results$attempted <- NULL
  saveRDS(mutated, checkpoint_path)

  env$mvt_run_grid(
    grid, run_dir, seed_base = 815L, require_gamcopula = FALSE,
    workers = 1L, run_case_fun = runner, fingerprints = fingerprints
  )
  expect_length(readLines(marker), 2L)
  ledger <- utils::read.csv(file.path(run_dir, "checkpoint_rejections.csv"), stringsAsFactors = FALSE)
  expect_equal(nrow(ledger), 1L)
  expect_match(ledger$rejection_reason, "semantic_result_schema")
  expect_true(file.exists(ledger$quarantine_path))
  expect_true(env$mvt_case_checkpoint_valid(
    readRDS(checkpoint_path),
    env$mvt_prepare_tasks(grid, 815L, FALSE, fingerprints)[[1L]]
  ))
  env$mvt_run_grid(
    grid, run_dir, seed_base = 815L, require_gamcopula = FALSE,
    workers = 1L, run_case_fun = function(...) stop("valid resume executed twice"),
    fingerprints = fingerprints
  )
  expect_length(readLines(marker), 2L)
  persisted_ledger <- utils::read.csv(
    file.path(run_dir, "checkpoint_rejections.csv"), stringsAsFactors = FALSE
  )
  expect_equal(nrow(persisted_ledger), 1L)
  expect_match(persisted_ledger$rejection_reason, "semantic_result_schema")

  writeBin(as.raw(c(1, 2, 3, 4)), checkpoint_path)
  env$mvt_run_grid(
    grid, run_dir, seed_base = 815L, require_gamcopula = FALSE,
    workers = 1L, run_case_fun = runner, fingerprints = fingerprints
  )
  expect_length(readLines(marker), 3L)
  corrupt_ledger <- utils::read.csv(
    file.path(run_dir, "checkpoint_rejections.csv"), stringsAsFactors = FALSE
  )
  expect_equal(nrow(corrupt_ledger), 2L)
  expect_true(any(grepl("read_error", corrupt_ledger$rejection_reason, fixed = TRUE)))
})

test_that("checkpoint provenance values bind task order, worker ownership, and strict time bounds", {
  env <- local_mvt_env()
  withr::local_envvar(c(GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm"))
  grid <- env$mvt_expand_grid(
    time_names = "t5", family_names = "gaussian",
    dependence_names = "external_exchangeable", reps = 1L
  )
  fingerprints <- list(producer_fingerprint = "producer", code_fingerprint = "code", package_fingerprint = "packages")
  task <- env$mvt_prepare_tasks(grid, 711L, FALSE, fingerprints, workers_requested = 1L)[[1L]]
  attempt <- data.frame(
    case_id = task$case_id, scenario = task$row$scenario, family = task$row$family_name,
    rep = task$row$rep, method = "glm", attempted = TRUE, success = TRUE,
    converged = TRUE, retained = TRUE, stop_reason = "converged",
    status_class = "success", elapsed_sec = 0, failure_reason_short = "none",
    error = NA_character_, stringsAsFactors = FALSE
  )
  result <- list(
    fit_status = attempt, benchmark_results = attempt, coefficient_results = data.frame(),
    dependence_recovery = data.frame(), variogram_scores = data.frame(), runtime = attempt
  )
  path <- tempfile("strict-checkpoint-", fileext = ".rds")
  env$mvt_write_case_checkpoint_atomic(result, task, path)
  checkpoint <- readRDS(path)
  expect_true(env$mvt_case_checkpoint_valid(checkpoint, task))

  mutate_and_reject <- function(mutator) {
    forged <- checkpoint
    forged <- mutator(forged)
    expect_false(env$mvt_case_checkpoint_valid(forged, task))
  }
  mutate_and_reject(function(x) { x$metadata$task_index <- 2L; x })
  mutate_and_reject(function(x) { x$metadata$task_index <- 0L; x })
  mutate_and_reject(function(x) { x$metadata$worker_pid <- x$metadata$worker_pid + 1L; x })
  mutate_and_reject(function(x) { x$metadata$worker_pid <- -1L; x })
  mutate_and_reject(function(x) { x$metadata$created_at <- "2026-09-01 12:00:00 AEST"; x })
  mutate_and_reject(function(x) { x$metadata$task_started_at <- "2099-01-01T00:00:00.000Z"; x })
  mutate_and_reject(function(x) {
    x$metadata$worker_attestation$writer_role <- "psock_worker"
    x$metadata$worker_attestation$attestation_sha256 <- env$mvt_hash_object(
      x$metadata$worker_attestation[setdiff(names(x$metadata$worker_attestation), "attestation_sha256")]
    )
    x
  })
})

test_that("candidate and public integration APIs reject nonproduction and arbitrary row bundles", {
  env <- local_mvt_env()
  withr::local_envvar(c(
    GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm",
    GAMLSS_LONGITUDINAL_MVT_RESUME = "true"
  ))
  grid <- env$mvt_expand_grid(
    time_names = "t5", family_names = "gaussian",
    dependence_names = "external_exchangeable", reps = 1L
  )
  runner <- function(row, seed_base, require_gamcopula, case_seed) {
    attempt <- data.frame(
      case_id = row$case_id, scenario = row$scenario, family = row$family_name,
      rep = row$rep, method = "glm", attempted = TRUE, success = TRUE,
      converged = TRUE, retained = TRUE, stop_reason = "converged",
      status_class = "success", elapsed_sec = 0, failure_reason_short = "none",
      error = NA_character_, stringsAsFactors = FALSE
    )
    list(fit_status = attempt, benchmark_results = attempt,
      coefficient_results = data.frame(), dependence_recovery = data.frame(),
      variogram_scores = data.frame(), runtime = attempt)
  }
  run_dir <- tempfile("mvt-nonproduction-candidate-")
  fingerprints <- list(producer_fingerprint = "producer", code_fingerprint = "code", package_fingerprint = "packages")
  env$mvt_run_grid(grid, run_dir, 909L, FALSE, workers = 1L, run_case_fun = runner, fingerprints = fingerprints)
  snapshot <- readRDS(file.path(run_dir, "aggregate_snapshot.rds"))
  expect_error(env$mvt_write_phase2_snapshot_candidate(run_dir, snapshot), "exact production|not registered")

  forged <- tempfile("mvt-forged-public-"); dir.create(forged)
  utils::write.csv(data.frame(
    case_id = sprintf("case-%04d", rep(seq_len(1600L), each = 9L)),
    method = rep(sprintf("method-%02d", seq_len(9L)), 1600L),
    attempted = TRUE, success = TRUE, converged = TRUE, retained = TRUE,
    stringsAsFactors = FALSE
  ), file.path(forged, "benchmark_results_by_rep.csv"), row.names = FALSE)
  expect_equal(nrow(utils::read.csv(file.path(forged, "benchmark_results_by_rep.csv"))), 14400L)
  expect_error(
    env$mvt_integrate_approved_phase2_snapshot(forged, tempfile("mvt-public-destination-")),
    "valid immutable aggregate snapshot|snapshot"
  )
  expect_error(env$mvt_validate_phase2_claim_evidence(forged), "snapshot")

  source(file.path(local_mvt_repo_root(), "paper", "R", "public-paper-producers.R"), local = env)
  withr::local_envvar(GAMLSS_LONGITUDINAL_MVT_APPROVED_RUN_DIR = forged)
  settings <- list(
    public_data_dir = tempfile("unused-public-data-"),
    data_dir = tempfile("mvt-central-public-output-"),
    multivariate_benchmark_attestation = tempfile("unapproved-attestation-", fileext = ".bin"),
    multivariate_benchmark_signature = tempfile("unapproved-signature-", fileext = ".sig")
  )
  expect_error(
    env$jss_run_phase2_multivariate_benchmark(settings),
    "valid immutable aggregate snapshot|snapshot"
  )
})

test_that("committed aggregate snapshots reject mutable status artifacts", {
  env <- local_mvt_env()
  withr::local_envvar(c(
    GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm",
    GAMLSS_LONGITUDINAL_MVT_RESUME = "true"
  ))
  grid <- env$mvt_expand_grid(
    time_names = "t5", family_names = "gaussian",
    dependence_names = "external_exchangeable", reps = 1L
  )
  runner <- function(row, seed_base, require_gamcopula, case_seed) {
    attempt <- data.frame(
      case_id = row$case_id, scenario = row$scenario, family = row$family_name,
      rep = row$rep, method = "glm", attempted = TRUE, success = TRUE,
      converged = TRUE, retained = TRUE, stop_reason = "converged",
      status_class = "success", elapsed_sec = 0,
      failure_reason_short = "none", error = NA_character_, stringsAsFactors = FALSE
    )
    list(
      fit_status = attempt, benchmark_results = attempt,
      coefficient_results = data.frame(), dependence_recovery = data.frame(),
      variogram_scores = data.frame(), runtime = attempt
    )
  }
  run_dir <- tempfile("mvt-immutable-snapshot-")
  fingerprints <- list(producer_fingerprint = "producer", code_fingerprint = "code", package_fingerprint = "packages")
  env$mvt_run_grid(grid, run_dir, 909L, FALSE, workers = 1L, run_case_fun = runner, fingerprints = fingerprints)
  committed <- env$mvt_validate_committed_snapshot(run_dir)
  expect_identical(committed$status, "committed_immutable")
  trust_registry <- tempfile("mvt-trust-registry-", fileext = ".csv")
  utils::write.csv(data.frame(
    registry_version = 1L, snapshot_schema_version = committed$schema_version,
    profile = "paper", status = "approved",
    snapshot_sha256 = env$mvt_phase2_snapshot_trust_sha256(run_dir),
    producer_id = committed$producer_id, producer_version = committed$producer_version,
    configuration_fingerprint = committed$configuration_fingerprint,
    stringsAsFactors = FALSE
  ), trust_registry, row.names = FALSE)
  fake_signature <- tempfile("mvt-self-approval-", fileext = ".sig")
  writeBin(as.raw(rep(0L, 64L)), fake_signature)
  expect_error(env$mvt_read_signed_snapshot_approval(trust_registry, fake_signature),
    "valid detached production signature")

  worker_path <- file.path(run_dir, "worker_attestations.csv")
  worker_commit_path <- paste0(worker_path, ".commit.rds")
  worker_ownership_path <- paste0(worker_path, ".ownership.rds")
  original_commit <- readRDS(worker_commit_path)
  original_ownership <- readRDS(worker_ownership_path)
  for (field in c("writer_pid", "owner_role", "rows", "schema_sha256", "bytes")) {
    forged <- original_commit
    forged[[field]] <- switch(field,
      writer_pid = forged[[field]] + 1L,
      owner_role = "psock_worker",
      rows = forged[[field]] + 1L,
      schema_sha256 = paste(rep("0", 64), collapse = ""),
      bytes = forged[[field]] + 1
    )
    saveRDS(forged, worker_commit_path, version = 3)
    saveRDS(forged, worker_ownership_path, version = 3)
    expect_error(
      env$mvt_validate_aggregate_file_commit(
        worker_path, expected_nonce = original_commit$lease_nonce,
        expected_pid = original_commit$writer_pid
      ),
      "do not reconcile", info = field
    )
  }
  saveRDS(original_commit, worker_commit_path, version = 3)
  saveRDS(original_ownership, worker_ownership_path, version = 3)

  status_path <- file.path(run_dir, "fit_status_by_rep.csv")
  status_raw <- readBin(status_path, "raw", n = file.info(status_path)$size)
  status <- utils::read.csv(status_path, stringsAsFactors = FALSE)
  status$converged <- FALSE
  utils::write.csv(status, status_path, row.names = FALSE)
  expect_error(env$mvt_validate_committed_snapshot(run_dir), "missing or mutable")
  writeBin(status_raw, status_path)
  expect_identical(env$mvt_validate_committed_snapshot(run_dir)$status, "committed_immutable")

  worker <- utils::read.csv(worker_path, stringsAsFactors = FALSE, check.names = FALSE)
  worker$hostname <- paste0(worker$hostname, "-repinned")
  utils::write.csv(worker, worker_path, row.names = FALSE, na = "")
  repinned_commit <- original_commit
  repinned_commit$sha256 <- env$mvt_sha256_file(worker_path)
  repinned_commit$rows <- nrow(worker)
  repinned_commit$columns <- names(worker)
  repinned_commit$schema_sha256 <- env$mvt_hash_object(names(worker))
  repinned_commit$bytes <- as.numeric(file.info(worker_path)$size)
  saveRDS(repinned_commit, worker_commit_path, version = 3)
  saveRDS(repinned_commit, worker_ownership_path, version = 3)
  snapshot_path <- file.path(run_dir, "aggregate_snapshot.rds")
  snapshot <- readRDS(snapshot_path)
  idx <- match("worker_attestations.csv", snapshot$artifacts$file)
  snapshot$artifacts$sha256[[idx]] <- repinned_commit$sha256
  snapshot$artifacts$rows[[idx]] <- repinned_commit$rows
  snapshot$artifacts$schema_sha256[[idx]] <- repinned_commit$schema_sha256
  snapshot$artifacts$bytes[[idx]] <- repinned_commit$bytes
  snapshot$artifacts$commit_sha256[[idx]] <- env$mvt_sha256_file(worker_commit_path)
  snapshot$artifacts$ownership_sha256[[idx]] <- env$mvt_sha256_file(worker_ownership_path)
  saveRDS(snapshot, snapshot_path, version = 3)
  root_commit_path <- paste0(snapshot_path, ".commit.rds")
  root_ownership_path <- paste0(snapshot_path, ".ownership.rds")
  root_commit <- readRDS(root_commit_path)
  root_commit$sha256 <- env$mvt_sha256_file(snapshot_path)
  root_commit$bytes <- as.numeric(file.info(snapshot_path)$size)
  root_commit$snapshot_schema_sha256 <- env$mvt_hash_object(names(snapshot))
  root_commit$checkpoint_manifest_sha256 <- env$mvt_hash_object(snapshot$checkpoint_manifest)
  root_commit$artifacts_sha256 <- env$mvt_hash_object(snapshot$artifacts)
  saveRDS(root_commit, root_commit_path, version = 3)
  saveRDS(root_commit, root_ownership_path, version = 3)
  expect_identical(env$mvt_validate_committed_snapshot(run_dir)$status, "committed_immutable")
  expect_error(
    env$mvt_validate_approved_snapshot(run_dir, trust_registry, require_production = FALSE,
      signature_path = fake_signature),
    "production audit|valid detached production signature"
  )
})

test_that("SHA-256 identity and local child-process attestation defeat a fake v99 library", {
  env <- local_mvt_env()
  fingerprints <- env$mvt_checkpoint_fingerprints()
  expect_identical(fingerprints$algorithm, "SHA-256")
  expect_match(fingerprints$producer_fingerprint, "^[0-9a-f]{64}$")
  expect_match(fingerprints$code_fingerprint, "^[0-9a-f]{64}$")
  expect_match(fingerprints$package_fingerprint, "^[0-9a-f]{64}$")
  expect_identical(
    fingerprints$package_versions$namespace_path[
      fingerprints$package_versions$package == "gamlss.longitudinal"
    ],
    fingerprints$package_identity$checkout_path
  )

  fake_lib <- tempfile("fake-v99-lib-")
  fake_package <- file.path(fake_lib, "gamlss.longitudinal")
  dir.create(fake_package, recursive = TRUE)
  writeLines(c(
    "Package: gamlss.longitudinal", "Version: 99.0.0", "Title: Fake",
    "Description: Adversarial stale package.", "License: MIT"
  ), file.path(fake_package, "DESCRIPTION"))
  writeLines("exportPattern(\"^[[:alpha:]]+\")", file.path(fake_package, "NAMESPACE"))
  old_libpaths <- .libPaths()
  .libPaths(c(fake_lib, old_libpaths))
  withr::defer(.libPaths(old_libpaths), teardown_env())
  expect_identical(normalizePath(.libPaths()[[1L]], winslash = "/"), normalizePath(fake_lib, winslash = "/"))

  probe <- env$mvt_run_fit_with_timeout("mvt_subprocess_identity_probe", timeout = 60)
  expect_false(inherits(probe$value, "error"))
  expect_identical(probe$value$version, fingerprints$package_identity$version)
  expect_identical(probe$value$namespace_path, fingerprints$package_identity$checkout_path)
  expect_true(probe$subprocess_attestation$verified)
  expect_identical(
    probe$subprocess_attestation$verified_source_sha256,
    fingerprints$package_identity$source_sha256
  )

  withr::local_envvar(c(
    GAMLSS_LONGITUDINAL_MVT_SOURCE = "local",
    GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = "8",
    GAMLSS_LONGITUDINAL_MVT_GEE_TIMEOUT_SEC = "60"
  ))
  gee_grid <- env$mvt_expand_grid(
    time_names = "t5", family_names = "gaussian",
    dependence_names = "external_exchangeable", reps = 1L
  )
  gee_row <- gee_grid[1L, , drop = FALSE]
  gee_data <- env$mvt_simulate_case(
    gee_row, seed = env$mvt_case_seed(gee_row, seed_base = 817L)
  )
  gee <- env$mvt_run_one_gee(
    gee_data, env$mvt_family_specs()[["gaussian"]], gee_row, "independence"
  )
  expect_equal(nrow(gee$results), 1L)
  expect_true(gee$results$subprocess_package_verified)
  expect_identical(gee$results$subprocess_source_sha256, fingerprints$package_identity$source_sha256)
  expect_identical(gee$results$subprocess_namespace_path, fingerprints$package_identity$checkout_path)
})

test_that("actual seeded package fits are numerically equivalent with one and two workers", {
  env <- local_mvt_env()
  withr::local_envvar(c(
    GAMLSS_LONGITUDINAL_MVT_SOURCE = "local",
    GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = "10",
    GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "gamlss.longitudinal,gee_independence,gamCopula_markov",
    GAMLSS_LONGITUDINAL_MVT_RESUME = "true",
    GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM = "1",
    GAMLSS_LONGITUDINAL_MVT_COMPUTE_VCOV = "false",
    GAMLSS_LONGITUDINAL_MVT_SUMMARY_VCOV = "false",
    GAMLSS_LONGITUDINAL_MVT_MAX_INNER_ITER = "20",
    GAMLSS_LONGITUDINAL_MVT_MAX_OUTER_ITER = "20",
    GAMLSS_LONGITUDINAL_MVT_PRIMARY_MAX_ELAPSED_SEC = "60",
    GAMLSS_LONGITUDINAL_MVT_GEE_TIMEOUT_SEC = "60",
    GAMLSS_LONGITUDINAL_MVT_GAMCOPULA_MARKOV_TIMEOUT_SEC = "60"
  ))
  grid <- env$mvt_expand_grid(
    time_names = "t5", family_names = "gaussian",
    dependence_names = "external_exchangeable", reps = 1:2
  )
  serial_dir <- tempfile("mvt-real-serial-")
  parallel_dir <- tempfile("mvt-real-parallel-")
  env$mvt_run_grid(
    grid, serial_dir, seed_base = 816L, require_gamcopula = TRUE,
    workers = 1L
  )
  env$mvt_run_grid(
    grid, parallel_dir, seed_base = 816L, require_gamcopula = TRUE,
    workers = 2L
  )

  serial <- utils::read.csv(file.path(serial_dir, "benchmark_results_by_rep.csv"), stringsAsFactors = FALSE)
  parallel <- utils::read.csv(file.path(parallel_dir, "benchmark_results_by_rep.csv"), stringsAsFactors = FALSE)
  expect_setequal(unique(serial$method), c("gamlss.longitudinal", "gee_independence", "gamCopula_markov"))
  expect_true(all(serial$retained[serial$method == "gamCopula_markov"]))
  expect_true(all(serial$subprocess_full_verified))
  expect_identical(parallel[c("case_id", "scenario", "family", "rep", "method")], serial[c("case_id", "scenario", "family", "rep", "method")])
  deterministic_numeric <- setdiff(
    intersect(names(serial)[vapply(serial, is.numeric, logical(1L))], names(parallel)),
    "elapsed_sec"
  )
  expect_true(any(grepl("rmse|score|coverage", deterministic_numeric)))
  tolerance <- max(env$mvt_parallel_equivalence_tolerances())
  for (column in deterministic_numeric) {
    expect_equal(parallel[[column]], serial[[column]], tolerance = tolerance, info = column)
  }

  checkpoints <- lapply(
    list.files(env$mvt_case_checkpoint_dir(parallel_dir), pattern = "[.]rds$", full.names = TRUE),
    readRDS
  )
  expect_equal(length(checkpoints), nrow(grid))
  expect_true(all(vapply(checkpoints, function(x) x$metadata$writer_role == "psock_worker", logical(1L))))
  expect_equal(anyDuplicated(vapply(checkpoints, function(x) x$metadata$case_id, character(1L))), 0L)
  aggregate_owner <- readRDS(paste0(file.path(parallel_dir, "benchmark_results_by_rep.csv"), ".ownership.rds"))
  expect_identical(aggregate_owner$owner_role, "lease_parent")
  expect_identical(aggregate_owner$writer_pid, Sys.getpid())
  attest <- utils::read.csv(file.path(parallel_dir, "worker_attestations.csv"), stringsAsFactors = FALSE)
  expect_equal(sum(attest$role == "worker"), 2L)
  expect_true(all(attest$verified[attest$role == "worker"]))
  expect_true(all(attest$full_verified[attest$role == "worker"]))
  expect_true(all(c(
    "pid", "setup_path", "source_sha256", "loaded_namespace_path",
    "gc_peak_ncells", "gc_peak_vcells", "blas", "lapack"
  ) %in% names(attest)))
  metadata <- utils::read.csv(file.path(parallel_dir, "run_metadata.csv"), stringsAsFactors = FALSE)
  expect_true(all(c(
    "run_completed_at", "run_elapsed_sec", "hostname", "cpu_model",
    "logical_cores", "physical_cores", "ram_bytes", "blas", "lapack",
    "producer_id", "producer_version", "producer_fingerprint",
    "code_fingerprint", "package_fingerprint"
  ) %in% metadata$name))
  expect_true(nzchar(metadata$value[metadata$name == "run_completed_at"]))
  expect_true(as.numeric(metadata$value[metadata$name == "run_elapsed_sec"]) >= 0)
})

test_that("preflight records grid, package, resume, and timeout checks", {
  env <- local_mvt_env()
  run_dir <- file.path(tempdir(), paste0("mvt-preflight-", sample.int(1e6, 1)))
  dir.create(run_dir, recursive = TRUE)

  withr::local_envvar(GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm,gee_unstructured")
  grid <- env$mvt_expand_grid(
    time_names = "t50",
    family_names = "gaussian",
    dependence_names = "external_exchangeable",
    reps = 1L
  )
  preflight <- env$mvt_write_preflight(grid, run_dir, require_gamcopula = FALSE, resume = TRUE)

  expect_true(file.exists(file.path(run_dir, "preflight_checks.csv")))
  expect_true(file.exists(file.path(run_dir, "preflight_checks.md")))
  expect_true(all(c("grid_nonempty", "total_rows_within_cap", "gee_unstructured_timeout") %in% preflight$check))
  expect_false(any(preflight$status == "fail"))
})

test_that("review audit writes machine and markdown reports", {
  env <- local_mvt_env()
  run_dir <- file.path(tempdir(), paste0("mvt-audit-", sample.int(1e6, 1)))
  dir.create(run_dir, recursive = TRUE)

  grid <- data.frame(case_id = "case-a", stringsAsFactors = FALSE)
  audit_methods <- c("glm", "gamCopula_markov", "gamCopula_vine_simplified", "gamlss.longitudinal")
  status <- data.frame(
    case_id = "case-a",
    method = audit_methods,
    success = TRUE,
    status_class = "ok",
    stringsAsFactors = FALSE
  )
  bench <- status
  coefs <- data.frame(
    case_id = rep("case-a", length(audit_methods) * 4L),
    method = rep(audit_methods, each = 4),
    parameter = "mu",
    term = rep(c("intercept", "time", "x", "z"), length(audit_methods)),
    estimate = 0,
    std_error = 1,
    truth = rep(c(0, 1, 2, 3), length(audit_methods)),
    stringsAsFactors = FALSE
  )
  dep <- data.frame(
    case_id = "case-a",
    method = c("gamCopula_markov", "gamCopula_vine_simplified", "gamlss.longitudinal"),
    theta_mae = c(0.1, 0.1, 0.1),
    tau_mae = c(0.1, 0.1, 0.1),
    stringsAsFactors = FALSE
  )
  vario <- data.frame(
    case_id = "case-a",
    method = c("glm", "gamCopula_markov", "gamCopula_vine_simplified", "gamlss.longitudinal"),
    variogram_score_p05 = c(NA, 1, 1, 1),
    stringsAsFactors = FALSE
  )

  env$mvt_write_csv(grid, file.path(run_dir, "scenario_grid.csv"))
  env$mvt_write_csv(status, file.path(run_dir, "fit_status_by_rep.csv"))
  env$mvt_write_csv(bench, file.path(run_dir, "benchmark_results_by_rep.csv"))
  env$mvt_write_csv(coefs, file.path(run_dir, "coefficient_results_by_rep.csv"))
  env$mvt_write_csv(dep, file.path(run_dir, "dependence_recovery_by_rep.csv"))
  env$mvt_write_csv(vario, file.path(run_dir, "variogram_scores_by_rep.csv"))
  env$mvt_write_csv(status, file.path(run_dir, "runtime_by_rep.csv"))
  for (path in c(
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
    "case_method_completion_summary.csv",
    "pilot_feasibility_by_method.csv",
    "pilot_feasibility_by_scenario.csv",
    "pilot_feasibility_overall_method.csv",
    "preflight_checks.csv",
    "artifact_manifest.csv",
    "run_metadata.csv",
    "package_versions.csv"
  )) {
    env$mvt_write_csv(data.frame(x = "ok"), file.path(run_dir, path))
  }
  env$mvt_write_csv(data.frame(check = "fixture", status = "pass"), file.path(run_dir, "phase2_benchmark_audit.csv"))
  writeLines("pilot", file.path(run_dir, "pilot_feasibility.md"))
  writeLines("preflight", file.path(run_dir, "preflight_checks.md"))
  writeLines("manifest", file.path(run_dir, "artifact_manifest.md"))
  writeLines("phase2", file.path(run_dir, "phase2_benchmark_audit.md"))
  writeLines("session", file.path(run_dir, "session_info.txt"))

  audit <- env$mvt_audit_run_dir(run_dir)

  expect_true(file.exists(file.path(run_dir, "review_audit.csv")))
  expect_true(file.exists(file.path(run_dir, "review_audit.md")))
  expect_true(file.exists(file.path(run_dir, "artifact_manifest.csv")))
  expect_true(file.exists(file.path(run_dir, "artifact_manifest.md")))
  expect_false(any(audit$status == "fail"))
})

test_that("review bundle writes a human-facing index", {
  env <- local_mvt_env()
  run_dir <- file.path(tempdir(), paste0("mvt-bundle-", sample.int(1e6, 1)))
  dir.create(file.path(run_dir, "paper_tables"), recursive = TRUE)
  dir.create(file.path(run_dir, "figures"), recursive = TRUE)
  writeLines("audit", file.path(run_dir, "review_audit.md"))
  writeLines("preflight", file.path(run_dir, "preflight_checks.md"))
  env$mvt_write_csv(data.frame(x = 1), file.path(run_dir, "benchmark_summary.csv"))
  env$mvt_write_csv(data.frame(x = 1), file.path(run_dir, "paper_tables", "benchmark_summary.csv"))
  writeLines("png", file.path(run_dir, "figures", "runtime.png"))

  index <- env$mvt_write_review_bundle(run_dir)

  expect_true(file.exists(index))
  expect_true(file.exists(file.path(run_dir, "artifact_manifest.csv")))
  lines <- readLines(index, warn = FALSE)
  manifest <- utils::read.csv(file.path(run_dir, "artifact_manifest.csv"), stringsAsFactors = FALSE)
  expect_true(any(grepl("paper_tables/benchmark_summary.csv", lines, fixed = TRUE)))
  expect_true(any(grepl("figures/runtime.png", lines, fixed = TRUE)))
  expect_true("review_bundle/README.md" %in% manifest$path)
})

test_that("case method completion summary reports success and classified failures", {
  env <- local_mvt_env()
  status <- data.frame(
    case_id = "case-a",
    scenario = "scenario-a",
    generator = "external",
    dependence = "exchangeable",
    n_time = 5L,
    family = "gaussian",
    rep = 1L,
    method = c("glm", "gamlss.longitudinal", "gee_unstructured"),
    success = c(TRUE, TRUE, FALSE),
    status_class = c("ok", "warning", "timeout"),
    warning = c(NA_character_, "Hessian was not positive definite", NA_character_),
    error = c(NA_character_, NA_character_, "timed out after 30 seconds"),
    stringsAsFactors = FALSE
  )

  out <- env$mvt_case_method_completion_summary(status)

  expect_equal(nrow(out), 1L)
  expect_equal(out$n_methods_success, 2L)
  expect_equal(out$n_methods_warning, 1L)
  expect_equal(out$n_methods_timeout, 1L)
  expect_match(out$effective_methods_completed, "glm")
  expect_match(out$methods_with_warnings, "gamlss.longitudinal", fixed = TRUE)
  expect_match(out$methods_failed, "gee_unstructured", fixed = TRUE)
  expect_match(out$failure_reasons_short, "vcov_or_hessian_warning")
  expect_match(out$failure_reasons_short, "timeout")
})

test_that("optional empty result tables do not break summaries", {
  env <- local_mvt_env()
  run_dir <- file.path(tempdir(), paste0("mvt-empty-optional-", sample.int(1e6, 1)))
  dir.create(run_dir, recursive = TRUE)

  grid <- data.frame(
    case_id = "case-a",
    scenario = "scenario-a",
    generator = "external",
    dependence = "exchangeable",
    dependence_name = "external_exchangeable",
    correlation_level = "moderate",
    n_time = 5L,
    n_subject = 10L,
    total_rows = 50L,
    family = "gaussian",
    family_name = "gaussian",
    gamlss_family = "NO",
    rep = 1L,
    stringsAsFactors = FALSE
  )
  status <- data.frame(
    case_id = "case-a",
    scenario = "scenario-a",
    generator = "external",
    dependence = "exchangeable",
    correlation_level = "moderate",
    n_time = 5L,
    n_subject = 10L,
    total_rows = 50L,
    family = "gaussian",
    gamlss_family = "NO",
    rep = 1L,
    method = "glmm_slope",
    available = TRUE,
    success = TRUE,
    status_class = "ok",
    elapsed_sec = 1,
    nobs = 50L,
    warning = NA_character_,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
  bench <- cbind(status, mae = 1, rmse = 1, benchmark_mean_rmse = 1, benchmark_neg_log_score = 1)
  coef <- data.frame(
    case_id = "case-a",
    scenario = "scenario-a",
    n_time = 5L,
    family = "gaussian",
    method = "glmm_slope",
    parameter = "mu",
    term = "time",
    estimate = 0,
    std_error = 1,
    truth = 0,
    bias = 0,
    ci_width = 1,
    ci_covers_truth = TRUE,
    false_positive = FALSE,
    stringsAsFactors = FALSE
  )
  vario <- data.frame(
    case_id = "case-a",
    scenario = "scenario-a",
    n_time = 5L,
    family = "gaussian",
    method = "glmm_slope",
    variogram_score_p05 = NA_real_,
    variogram_score_p2 = NA_real_,
    stringsAsFactors = FALSE
  )

  env$mvt_write_csv(grid, file.path(run_dir, "scenario_grid.csv"))
  env$mvt_write_csv(status, file.path(run_dir, "fit_status_by_rep.csv"))
  env$mvt_write_csv(bench, file.path(run_dir, "benchmark_results_by_rep.csv"))
  env$mvt_write_csv(coef, file.path(run_dir, "coefficient_results_by_rep.csv"))
  writeLines('""', file.path(run_dir, "dependence_recovery_by_rep.csv"))
  env$mvt_write_csv(vario, file.path(run_dir, "variogram_scores_by_rep.csv"))
  env$mvt_write_csv(status[c("case_id", "scenario", "family", "method", "elapsed_sec", "success", "error")], file.path(run_dir, "runtime_by_rep.csv"))

  expect_equal(nrow(env$mvt_read_optional_csv(file.path(run_dir, "dependence_recovery_by_rep.csv"))), 0L)
  expect_error(env$mvt_summarise_results(run_dir), "valid immutable aggregate snapshot")
  expect_false(file.exists(file.path(run_dir, "case_method_completion_summary.csv")))
})

test_that("run shard merger combines scenario and by-replicate outputs", {
  env <- local_mvt_env()
  root <- file.path(tempdir(), paste0("mvt-merge-", sample.int(1e6, 1)))
  shard_dirs <- file.path(root, paste0("shard", 1:2))
  target <- file.path(root, "merged")
  dir.create(shard_dirs[[1L]], recursive = TRUE)
  dir.create(shard_dirs[[2L]], recursive = TRUE)

  write_shard <- function(dir, rep_id) {
    grid <- data.frame(
      case_id = sprintf("case-rep%03d", rep_id),
      scenario = "external_exchangeable_t20",
      generator = "external",
      dependence = "exchangeable",
      dependence_name = "external_exchangeable",
      correlation_level = "moderate",
      n_time = 20L,
      n_subject = 2L,
      total_rows = 40L,
      family = "gaussian",
      family_name = "gaussian",
      gamlss_family = "NO",
      rep = rep_id,
      stringsAsFactors = FALSE
    )
    status <- data.frame(
      case_id = grid$case_id,
      scenario = grid$scenario,
      generator = grid$generator,
      dependence = grid$dependence,
      correlation_level = grid$correlation_level,
      n_time = grid$n_time,
      n_subject = grid$n_subject,
      total_rows = grid$total_rows,
      family = grid$family,
      gamlss_family = grid$gamlss_family,
      rep = rep_id,
      method = "glm",
      available = TRUE,
      success = TRUE,
      status_class = "ok",
      elapsed_sec = 1,
      nobs = 40L,
      warning = NA_character_,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
    bench <- cbind(status, mae = 1, rmse = 1, benchmark_mean_rmse = 1, benchmark_neg_log_score = 1)
    coef <- data.frame(
      case_id = grid$case_id,
      scenario = grid$scenario,
      n_time = grid$n_time,
      family = grid$family,
      method = "glm",
      parameter = "mu",
      term = "time",
      estimate = 0,
      std_error = 1,
      truth = 0,
      bias = 0,
      ci_width = 1,
      ci_covers_truth = TRUE,
      false_positive = FALSE,
      stringsAsFactors = FALSE
    )
    dep <- data.frame(case_id = grid$case_id, scenario = grid$scenario, n_time = grid$n_time, family = grid$family, method = "gamCopula", theta_mae = 0, theta_rmse = 0, tau_mae = 0, tau_rmse = 0)
    vario <- data.frame(case_id = grid$case_id, scenario = grid$scenario, n_time = grid$n_time, family = grid$family, method = "gamCopula", variogram_score_p05 = 1, variogram_score_p2 = 1)
    runtime <- status[c("case_id", "scenario", "family", "method", "elapsed_sec", "success", "error")]
    env$mvt_write_csv(grid, file.path(dir, "scenario_grid.csv"))
    env$mvt_write_csv(status, file.path(dir, "fit_status_by_rep.csv"))
    env$mvt_write_csv(bench, file.path(dir, "benchmark_results_by_rep.csv"))
    env$mvt_write_csv(coef, file.path(dir, "coefficient_results_by_rep.csv"))
    env$mvt_write_csv(dep, file.path(dir, "dependence_recovery_by_rep.csv"))
    env$mvt_write_csv(vario, file.path(dir, "variogram_scores_by_rep.csv"))
    env$mvt_write_csv(runtime, file.path(dir, "runtime_by_rep.csv"))
  }
  write_shard(shard_dirs[[1L]], 1L)
  write_shard(shard_dirs[[2L]], 2L)
  for (name in env$mvt_result_names()) {
    file.rename(
      file.path(shard_dirs[[2L]], paste0(name, "_by_rep.csv")),
      file.path(shard_dirs[[2L]], paste0(name, "_checkpoint.csv"))
    )
  }

  env$mvt_merge_run_shards(shard_dirs, target)

  grid <- utils::read.csv(file.path(target, "scenario_grid.csv"), stringsAsFactors = FALSE)
  status <- utils::read.csv(file.path(target, "fit_status_by_rep.csv"), stringsAsFactors = FALSE)
  expect_equal(sort(grid$rep), 1:2)
  expect_equal(nrow(status), 2L)
  expect_true(file.exists(file.path(target, "NONPUBLICATION.txt")))
  expect_false(file.exists(file.path(target, "aggregate_snapshot.rds")))
  owner <- readRDS(paste0(file.path(target, "fit_status_by_rep.csv"), ".ownership.rds"))
  expect_identical(owner$owner_role, "lease_parent")
})

test_that("publication suite dry run writes a plan and restores environment variables", {
  repo_root <- local_mvt_repo_root()
  script <- file.path(repo_root, "paper", "R", "09-simulation-multivariate-longitudinal", "10-run-publication-suite.R")
  suite_dir <- file.path(tempdir(), paste0("mvt-suite-", sample.int(1e6, 1)))
  suite_env <- new.env(parent = globalenv())

  withr::local_dir(tempdir())
  withr::local_envvar(
    GAMLSS_LONGITUDINAL_MVT_SUITE_DIR = suite_dir,
    GAMLSS_LONGITUDINAL_MVT_RUN_PUBLICATION_SUITE = "false",
    GAMLSS_LONGITUDINAL_MVT_SUITE_INCLUDE_SPECIAL = "true",
    GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = "10"
  )
  source(script, local = suite_env)

  plan_path <- file.path(suite_dir, "publication_suite_plan.csv")
  expect_true(file.exists(plan_path))
  expect_true(file.exists(file.path(suite_dir, "README.md")))
  expect_true(file.exists(file.path(suite_dir, "publication_suite_artifacts.csv")))
  plan <- utils::read.csv(plan_path, stringsAsFactors = FALSE)
  expect_setequal(plan$role, c("pilot", "main_core", "appendix", "special_gamlss"))
  expect_true(all(c("estimated_cases", "max_total_rows") %in% names(plan)))
  expect_equal(plan$estimated_cases[plan$role == "main_core"], 1600L)
  expect_equal(plan$max_total_rows[plan$role == "appendix"], 5000L)

  readiness_env <- suite_env$mvt_suite_readiness_env(
    suite_env$mvt_suite_role_specs(suite_dir),
    c("pilot", "main_core")
  )
  expect_equal(readiness_env[["GAMLSS_LONGITUDINAL_MVT_APPENDIX_RUN_DIR"]], "")
  expect_match(readiness_env[["GAMLSS_LONGITUDINAL_MVT_MAIN_RUN_DIR"]], "main_core")

  Sys.setenv(MVT_TMP_CHECK = "old")
  suite_env$mvt_suite_with_env(c(MVT_TMP_CHECK = "new", MVT_TMP_NEW = "value"), {
    expect_equal(Sys.getenv("MVT_TMP_CHECK"), "new")
    expect_equal(Sys.getenv("MVT_TMP_NEW"), "value")
  })
  expect_equal(Sys.getenv("MVT_TMP_CHECK"), "old")
  expect_equal(Sys.getenv("MVT_TMP_NEW"), "")
})

test_that("publication suite strict mode fails early for partial publication runs", {
  repo_root <- local_mvt_repo_root()
  script <- file.path(repo_root, "paper", "R", "09-simulation-multivariate-longitudinal", "10-run-publication-suite.R")
  suite_dir <- file.path(tempdir(), paste0("mvt-suite-strict-", sample.int(1e6, 1)))

  withr::local_dir(tempdir())
  withr::local_envvar(
    GAMLSS_LONGITUDINAL_MVT_SUITE_DIR = suite_dir,
    GAMLSS_LONGITUDINAL_MVT_RUN_PUBLICATION_SUITE = "true",
    GAMLSS_LONGITUDINAL_MVT_SUITE_ROLES = "pilot",
    GAMLSS_LONGITUDINAL_MVT_SUITE_STRICT_READINESS = "true"
  )

  expect_error(
    source(script, local = new.env(parent = globalenv())),
    "Strict publication readiness requires suite role"
  )
  expect_true(file.exists(file.path(suite_dir, "publication_suite_plan.csv")))
  expect_false(file.exists(file.path(suite_dir, "pilot", "fit_status_by_rep.csv")))
})

test_that("publication suite preflight-only mode writes role preflight artifacts without fits", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  skip_if_not_installed("gamCopula")
  skip_if_not_installed("VineCopula")
  skip_if_not_installed("mvtnorm")
  repo_root <- local_mvt_repo_root()
  script <- file.path(repo_root, "paper", "R", "09-simulation-multivariate-longitudinal", "10-run-publication-suite.R")
  suite_dir <- file.path(tempdir(), paste0("mvt-suite-preflight-", sample.int(1e6, 1)))

  withr::local_dir(tempdir())
  withr::local_envvar(
    GAMLSS_LONGITUDINAL_MVT_SUITE_DIR = suite_dir,
    GAMLSS_LONGITUDINAL_MVT_SUITE_PREFLIGHT_ONLY = "true",
    GAMLSS_LONGITUDINAL_MVT_SUITE_ROLES = "special_gamlss",
    GAMLSS_LONGITUDINAL_MVT_SUITE_STRICT_READINESS = "false"
  )
  source(script, local = new.env(parent = globalenv()))

  role_dir <- file.path(suite_dir, "special_gamlss")
  expect_true(file.exists(file.path(suite_dir, "publication_suite_preflight.csv")))
  expect_true(file.exists(file.path(suite_dir, "publication_suite_preflight.md")))
  expect_true(file.exists(file.path(role_dir, "scenario_grid.csv")))
  expect_true(file.exists(file.path(role_dir, "preflight_checks.csv")))
  expect_true(file.exists(file.path(role_dir, "run_metadata.csv")))
  expect_false(file.exists(file.path(role_dir, "fit_status_by_rep.csv")))
  preflight <- utils::read.csv(file.path(suite_dir, "publication_suite_preflight.csv"), stringsAsFactors = FALSE)
  expect_equal(preflight$role, "special_gamlss")
})

test_that("study protocol generator writes design and evidence tables", {
  env <- local_mvt_env()
  output_dir <- file.path(tempdir(), paste0("mvt-protocol-", sample.int(1e6, 1)))

  protocol <- env$mvt_write_study_protocol(output_dir)

  expect_true(file.exists(protocol$md))
  expect_true(dir.exists(file.path(output_dir, "study_protocol")))
  time_design <- utils::read.csv(file.path(output_dir, "study_protocol", "time_design.csv"), stringsAsFactors = FALSE)
  family_design <- utils::read.csv(file.path(output_dir, "study_protocol", "family_design.csv"), stringsAsFactors = FALSE)
  comparators <- utils::read.csv(file.path(output_dir, "study_protocol", "comparators.csv"), stringsAsFactors = FALSE)
  readiness <- utils::read.csv(file.path(output_dir, "study_protocol", "publication_readiness.csv"), stringsAsFactors = FALSE)
  checklist <- utils::read.csv(file.path(output_dir, "study_protocol", "review_checklist.csv"), stringsAsFactors = FALSE)

  expect_true(50L %in% time_design$n_time)
  expect_true("gg_continuous" %in% family_design$family_name)
  expect_true(all(c("gamCopula_markov", "gamCopula_vine_simplified", "gamCopula_vine", "gamCopula") %in% comparators$method))
  expect_true(all(c("pilot", "main_core", "appendix") %in% readiness$role))
  expect_true(any(readiness$required_method[readiness$role == "main_core"] == "glm,glmm,gee_independence,gee_exchangeable,gee_ar1,gee_unstructured,gamCopula_markov,gamCopula_vine_simplified,gamlss.longitudinal"))
  expect_true(any(grepl("GJRM excluded", checklist$requirement, fixed = TRUE)))
  expect_true(any(grepl("No missingness/dropout", checklist$requirement, fixed = TRUE)))
  expect_true(all(c("implementation_evidence", "proof_after_run") %in% names(checklist)))
})

test_that("implementation status distinguishes preflight readiness from full fit evidence", {
  env <- local_mvt_env()
  output_dir <- file.path(tempdir(), paste0("mvt-status-", sample.int(1e6, 1)))
  suite_dir <- file.path(output_dir, "publication_suite_preflight_default")
  dir.create(suite_dir, recursive = TRUE)

  env$mvt_write_study_protocol(output_dir)
  env$mvt_write_csv(
    data.frame(
      role = c("pilot", "main_core", "appendix"),
      run_dir = file.path(suite_dir, c("pilot", "main_core", "appendix")),
      script = "script.R",
      reps = c(5L, 100L, 100L),
      estimated_cases = c(90L, 1600L, 4800L),
      max_total_rows = c(4000L, 4000L, 5000L),
      stringsAsFactors = FALSE
    ),
    file.path(suite_dir, "publication_suite_plan.csv")
  )
  env$mvt_write_csv(
    data.frame(
      role = c("pilot", "main_core", "appendix"),
      checks = 6L,
      failures = 0L,
      warnings = 0L,
      stringsAsFactors = FALSE
    ),
    file.path(suite_dir, "publication_suite_preflight.csv")
  )
  for (role in c("pilot", "main_core", "appendix")) {
    role_dir <- file.path(suite_dir, role)
    dir.create(role_dir, recursive = TRUE)
    for (path in c(
      "scenario_grid.csv",
      "preflight_checks.csv",
      "preflight_checks.md",
      "run_metadata.csv",
      "package_versions.csv",
      "session_info.txt",
      "artifact_manifest.csv"
    )) {
      writeLines("x", file.path(role_dir, path))
    }
  }

  status <- env$mvt_write_implementation_status(output_dir, suite_dir)
  checks <- status$checks

  expect_true(file.exists(status$csv))
  expect_true(file.exists(status$md))
  expect_false(status$ready)
  expect_true(all(checks$status[grepl("^preflight_|^default_suite|^study_protocol|^workflow", checks$check)] == "pass"))
  expect_equal(checks$status[checks$check == "full_method_stack_smoke"], "incomplete")
  expect_equal(checks$status[checks$check == "all_family_review_smoke"], "incomplete")
  expect_equal(checks$status[checks$check == "t20_full_method_smoke"], "incomplete")
  expect_equal(checks$status[checks$check == "pilot_real_size_rep1_shard"], "incomplete")
  expect_equal(checks$status[checks$check == "main_core_rep1_shard"], "incomplete")
  expect_equal(checks$status[checks$check == "full_fit_evidence"], "incomplete")
  expect_equal(checks$status[checks$check == "publication_review_ready"], "incomplete")
})

test_that("publication readiness audit fails clearly without full pilot and main evidence", {
  env <- local_mvt_env()
  output_dir <- file.path(tempdir(), paste0("mvt-readiness-", sample.int(1e6, 1)))
  dir.create(output_dir, recursive = TRUE)

  withr::local_envvar(
    GAMLSS_LONGITUDINAL_MVT_PILOT_RUN_DIR = "",
    GAMLSS_LONGITUDINAL_MVT_MAIN_RUN_DIR = "",
    GAMLSS_LONGITUDINAL_MVT_APPENDIX_RUN_DIR = "",
    GAMLSS_LONGITUDINAL_MVT_SPECIAL_RUN_DIR = "",
    GAMLSS_LONGITUDINAL_MVT_GLMM_SENSITIVITY_RUN_DIR = ""
  )
  readiness <- env$mvt_publication_readiness_audit(output_dir)

  expect_false(readiness$ready)
  expect_true(file.exists(file.path(output_dir, "publication_readiness_audit.csv")))
  expect_true(file.exists(file.path(output_dir, "publication_readiness_audit.md")))
  expect_true(any(readiness$audit$role == "pilot" & readiness$audit$status == "fail"))
  expect_true(any(readiness$audit$role == "main_core" & readiness$audit$status == "fail"))
})

test_that("coefficient terms normalize and map to true mean effects", {
  env <- local_mvt_env()
  spec <- env$mvt_family_specs()[["gaussian"]]

  terms <- env$mvt_normalize_term(c("(Intercept)", "mu.time", "time_covariate", "x", "z"))
  truth <- env$mvt_truth_for_term(terms, spec)

  expect_equal(terms, c("intercept", "time", "time", "x", "z"))
  expect_equal(unname(truth), unlist(spec$eta[c("intercept", "time", "time", "x", "z")], use.names = FALSE))
})

test_that("gamCopula margin coefficients are extracted from the fitted GAMLSS margin", {
  env <- local_mvt_env()
  skip_if_not_installed("gamlss")

  set.seed(1)
  dat <- data.frame(
    response = rnorm(30),
    time = rep(seq(0, 1, length.out = 5), 6),
    x = rep(rnorm(6), each = 5),
    z = rep(c(0, 1), each = 15)
  )
  fit <- gamlss::gamlss(
    response ~ time + x + z,
    sigma.formula = ~1,
    family = gamlss.dist::NO(),
    data = dat,
    trace = FALSE
  )

  out <- env$mvt_coef_table_gamlss_margin(fit, "gamCopula")

  expect_true(all(c("method", "parameter", "term", "estimate", "std_error") %in% names(out)))
  expect_equal(out$method, rep("gamCopula", nrow(out)))
  expect_true(all(c("(Intercept)", "time", "x", "z") %in% out$term))
  expect_true(any(is.finite(out$std_error)))
})

test_that("gamCopula pair tau helper handles simplified selected copula objects", {
  env <- local_mvt_env()
  skip_if_not_installed("VineCopula")

  tau <- env$mvt_gamcopula_pair_tau(list(family = 1, par = 0.4, par2 = 0))

  expect_true(is.finite(tau))
  expect_gt(tau, 0)
})

test_that("gamCopula adjacent Markov simulator returns balanced finite uniforms", {
  env <- local_mvt_env()
  skip_if_not_installed("VineCopula")

  dat <- data.frame(
    subject = rep(1:2, each = 3),
    time_index = rep(1:3, 2),
    time = rep(c(0, 0.5, 1), 2),
    response = 0,
    x = rep(c(-1, 1), each = 3),
    z = rep(c(0, 1), each = 3)
  )
  fit <- structure(
    list(
      pair_fits = list(
        list(family = 1, par = 0.3, par2 = 0),
        list(family = 1, par = 0.3, par2 = 0)
      ),
      pair_covariates = list(
        data.frame(subject = 1:2, time = c(0, 0), x = c(-1, 1), z_binary = c(0, 1)),
        data.frame(subject = 1:2, time = c(0.5, 0.5), x = c(-1, 1), z_binary = c(0, 1))
      )
    ),
    class = "mvt_gamcopula_fit"
  )

  u <- env$mvt_simulate_gamcopula_uniforms(fit, dat, nsim = 4L, seed = 1L)

  expect_equal(dim(u), c(6L, 4L))
  expect_true(all(is.finite(u)))
  expect_true(all(u > 0 & u < 1))
})

test_that("copula variogram coverage is required for successful fits only", {
  env <- local_mvt_env()
  status <- data.frame(
    case_id = c("ok-gc", "ok-gl", "timeout-gl"),
    method = c("gamCopula", "gamlss.longitudinal", "gamlss.longitudinal"),
    success = c(TRUE, TRUE, FALSE),
    status_class = c("ok", "ok", "timeout"),
    stringsAsFactors = FALSE
  )
  vario <- data.frame(
    case_id = c("ok-gc", "ok-gl", "timeout-gl"),
    method = c("gamCopula", "gamlss.longitudinal", "gamlss.longitudinal"),
    variogram_score_p05 = c(1, 2, NA),
    stringsAsFactors = FALSE
  )

  coverage <- env$mvt_copula_variogram_coverage(
    status,
    vario,
    c("gamCopula", "gamlss.longitudinal")
  )

  expect_equal(coverage$status, "pass")
  expect_equal(coverage$n, 2L)
})

test_that("joint-model fitted parameter extraction uses the package parameter argument", {
  env <- local_mvt_env()
  assign(
    "fitted.fake_mvt_gamlss_longitudinal",
    function(object, parameter = "mu", ...) object[[parameter]],
    envir = globalenv()
  )
  withr::defer(rm("fitted.fake_mvt_gamlss_longitudinal", envir = globalenv()), teardown_env())
  fit <- structure(
    list(
      mu = c(1, 2, 3),
      sigma = c(0.5, 0.6, 0.7),
      nu = c(0.1, 0.2, 0.3),
      tau = c(-0.1, -0.2, -0.3)
    ),
    class = "fake_mvt_gamlss_longitudinal"
  )
  dat <- data.frame(response = 1:3)

  params <- env$mvt_predict_gamlss_params(fit, dat, list(sigma = NA_real_))

  expect_equal(params$mu, c(1, 2, 3))
  expect_equal(params$sigma, c(0.5, 0.6, 0.7))
  expect_equal(params$nu, c(0.1, 0.2, 0.3))
  expect_equal(params$tau, c(-0.1, -0.2, -0.3))
})

test_that("tiny multivariate pilot writes expected output files when optional runtime packages are available", {
  env <- local_mvt_env()
  skip_on_cran()
  skip_if_not_installed("gamCopula")
  skip_if_not_installed("VineCopula")
  skip_if_not_installed("mvtnorm")
  skip_if_not_installed("gamlss")
  skip_if_not_installed("geepack")
  skip_if_not_installed("lme4")

  withr::local_envvar(
    GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = "10",
    GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm",
    GAMLSS_LONGITUDINAL_MVT_GEE_UNSTRUCTURED_TIMEOUT_SEC = "5",
    GAMLSS_LONGITUDINAL_MVT_PRIMARY_MAX_ELAPSED_SEC = "20",
    GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM = "1"
  )
  grid <- env$mvt_expand_grid(
    time_names = "t5",
    family_names = "gaussian",
    dependence_names = "external_exchangeable",
    reps = 1L
  )
  run_dir <- file.path(tempdir(), paste0("mvt-smoke-", sample.int(1e6, 1)))

  env$mvt_run_grid(grid, run_dir = run_dir, checkpoint_every = 1L, require_gamcopula = TRUE)

  expected <- file.path(
    run_dir,
    c(
      "scenario_grid.csv",
      "fit_status_by_rep.csv",
      "benchmark_results_by_rep.csv",
      "coefficient_results_by_rep.csv",
      "dependence_recovery_by_rep.csv",
      "variogram_scores_by_rep.csv",
      "runtime_by_rep.csv",
      "preflight_checks.csv",
      "preflight_checks.md",
      "run_metadata.csv",
      "package_versions.csv"
    )
  )
  expect_true(all(file.exists(expected)))
})
