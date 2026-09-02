local_module07_env <- function() {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "07-gamma-copula-misspecification.R"), local = env)
  list(root = root, env = env)
}

test_that("Module 07 uses supported controls and an isolated v2 checkpoint namespace", {
  module <- local_module07_env()
  code <- paste(readLines(file.path(module$root, "paper", "R", "07-gamma-copula-misspecification.R"), warn = FALSE), collapse = "\n")
  expect_match(code, "optimizer_control = gamlss.longitudinal::gamlss_longitudinal_control", fixed = TRUE)
  expect_false(grepl("        max_outer_iter = config$max_outer_iter,\n        max_inner_iter", code, fixed = TRUE))
  paths <- module$env$jss_misspec_paths(list(data_dir = "data", tables_dir = "tables", figures_dir = "figures"))
  expect_identical(basename(paths$checkpoints), "07-gamma-copula-misspecification-checkpoints-ga-nc-v2")
})

test_that("Module 07 paired effects reconstruct analytic same-dataset uncertainty", {
  env <- local_module07_env()$env
  base <- data.frame(
    generating_copula = "N", tau_label = "high", target_tau = 0.55,
    n_subject = 500L, rep = 1:4, dataset_seed = 101:104,
    stringsAsFactors = FALSE
  )
  correct <- transform(base, fitted_copula = "N", retained = TRUE,
    joint_loglik = c(10, 11, 12, 13), margin_param_rmse = 0.1, tau_abs_error = 0.02)
  wrong <- transform(base, fitted_copula = "C", retained = TRUE,
    joint_loglik = c(8, 10, 9, 11), margin_param_rmse = c(.2, .3, .2, .4),
    tau_abs_error = c(.04, .05, .06, .03))
  effects <- env$jss_misspec_paired_effects(rbind(correct, wrong))
  row <- effects[effects$metric == "joint_loglik", , drop = FALSE]
  difference <- c(-2, -1, -3, -2)
  expected_mcse <- stats::sd(difference) / sqrt(length(difference))
  critical <- stats::qt(.975, df = 3)
  expect_equal(row$n_attempted, 4L)
  expect_equal(row$n_retained, 4L)
  expect_equal(row$estimate, mean(difference))
  expect_equal(row$mcse, expected_mcse)
  expect_equal(row$ci_lower, mean(difference) - critical * expected_mcse)
  expect_equal(row$ci_upper, mean(difference) + critical * expected_mcse)
  expect_identical(row$seed_pairing, "registered_dataset_seed_exact")
})

test_that("Module 07 warning audit distinguishes disclosed fallback from fatal unexpected events", {
  env <- local_module07_env()$env
  clean <- data.frame(
    fit_id = c("a", "b", "c"), retained = c(TRUE, TRUE, FALSE),
    warnings = c(NA, "Analytical Hessian unavailable; falling back to numerical Hessian.",
      "Model returned without satisfying the optimizer convergence contract (stop reason: max_iterations)."),
    stringsAsFactors = FALSE
  )
  audit <- env$jss_misspec_warning_audit(clean)
  expect_true(env$jss_misspec_validate_warning_audit(audit))
  expected <- audit[audit$warning_classification == "expected_numerical_fallback", ]
  expect_equal(expected$denominator_attempts, 3L)
  expect_equal(expected$n_affected_attempts, 1L)
  expect_equal(expected$n_retained_affected, 1L)
  expect_true(expected$warning_ci_lower <= expected$warning_rate)
  expect_true(expected$warning_ci_upper >= expected$warning_rate)
  nonconverged <- audit[audit$warning_classification == "expected_nonconvergence_signal", ]
  expect_equal(nonconverged$n_affected_attempts, 1L)
  expect_equal(nonconverged$n_retained_affected, 0L)

  bad <- clean
  bad$warnings[[1L]] <- "An unregistered warning"
  expect_error(
    env$jss_misspec_validate_warning_audit(env$jss_misspec_warning_audit(bad)),
    "unexpected warnings"
  )
})

test_that("Module 07 execution manifest binds generation-time source and artifact hashes", {
  module <- local_module07_env(); env <- module$env
  config <- env$jss_misspec_config(list(root = module$root, profile = "full", seed = 20260528L), stage = "smoke")
  grid <- env$jss_misspec_grid(config)
  out <- tempfile("m07-manifest-"); dir.create(out)
  paths <- c(results = file.path(out, "results.csv"), paired_effects = file.path(out, "paired-effects.csv"),
    warning_audit = file.path(out, "warning-audit.csv"), selection_attempts = file.path(out, "selection-attempts.csv"),
    selection = file.path(out, "selection.csv"))
  for (path in paths) env$jss_misspec_write_csv_atomic(data.frame(value = 1), path)
  git <- env$jss_misspec_git_identity(module$root)
  start <- list(
    started_at_utc = "2026-09-02T00:00:00Z", git_sha = git$sha, git_state = git$state,
    r_version = R.version.string, platform = R.version$platform,
    os = paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "),
    dependency_versions = env$jss_misspec_dependency_versions(), workers_requested = 2L
  )
  warning_audit <- env$jss_misspec_warning_audit(data.frame(
    fit_id = grid$fit_id, warnings = NA_character_, retained = TRUE, stringsAsFactors = FALSE
  ))
  env$jss_misspec_write_csv_atomic(warning_audit, paths[["warning_audit"]])
  manifest <- env$jss_misspec_execution_manifest(start, config, grid, grid, paths)
  expect_equal(manifest$workers_requested, 2L)
  expect_equal(manifest$dataset_seed_count, length(unique(grid$dataset_seed)))
  expect_equal(manifest$results_sha256, env$jss_misspec_sha256_file(paths[["results"]]))
  expect_equal(manifest$warning_audit_sha256, env$jss_misspec_sha256_file(paths[["warning_audit"]]))
  expect_match(manifest$dependency_versions_sha256, "^[0-9a-f]{64}$")
})

test_that("Module 07 claim mappings and visible figure metadata are explicit", {
  module <- local_module07_env(); env <- module$env
  claims <- env$jss_misspec_validate_claim_registry(module$root)
  expect_gte(sum(claims$claim_type == "paired_effect"), 2L)
  expect_gte(sum(claims$claim_type == "selection_rate"), 2L)
  expect_true(all(claims$effect_artifact_id[claims$claim_type == "paired_effect"] == "p2_copula_paired_effects"))
  metadata <- env$jss_misspec_figure_metadata(data.frame(
    n_subject = 500L, tau_label = "high", target_tau = 0.55
  ))
  expect_match(metadata$visible_label, "n = 500", fixed = TRUE)
  expect_match(metadata$visible_label, "tau = 0.55", fixed = TRUE)
  expect_match(metadata$visible_label, "R = 100", fixed = TRUE)
  expect_identical(metadata$tau_panel_title, "Mean absolute Kendall tau error")
  figure_code <- paste(deparse(body(env$jss_misspec_write_paper_summary_heatmap)), collapse = "\n")
  expect_match(figure_code, "grid::grid.text(title", fixed = TRUE)
})
