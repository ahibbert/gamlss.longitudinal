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
  expect_true(all(c("intercept", "time", "x", "z") %in% out$coefficient_results$term))
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
  audit_methods <- c("glm", "gamCopula_markov", "gamCopula_vine", "gamlss.longitudinal")
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
    method = c("gamCopula_markov", "gamCopula_vine", "gamlss.longitudinal"),
    theta_mae = c(0.1, 0.1, 0.1),
    tau_mae = c(0.1, 0.1, 0.1),
    stringsAsFactors = FALSE
  )
  vario <- data.frame(
    case_id = "case-a",
    method = c("glm", "gamCopula_markov", "gamCopula_vine", "gamlss.longitudinal"),
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
  writeLines("pilot", file.path(run_dir, "pilot_feasibility.md"))
  writeLines("preflight", file.path(run_dir, "preflight_checks.md"))
  writeLines("manifest", file.path(run_dir, "artifact_manifest.md"))
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
  expect_silent(env$mvt_summarise_results(run_dir))
  expect_true(file.exists(file.path(run_dir, "case_method_completion_summary.csv")))
  expect_true(file.exists(file.path(run_dir, "variogram_summary.csv")))
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
  expect_true(file.exists(file.path(target, "benchmark_summary.csv")))
  expect_true(file.exists(file.path(target, "artifact_manifest.csv")))
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
  expect_true(all(c("gamCopula_markov", "gamCopula_vine", "gamCopula") %in% comparators$method))
  expect_true(all(c("pilot", "main_core", "appendix") %in% readiness$role))
  expect_true(any(readiness$required_method[readiness$role == "main_core"] == "glm,glmm,gee_independence,gee_exchangeable,gee_ar1,gee_unstructured,gamCopula_markov,gamCopula_vine,gamlss.longitudinal"))
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
