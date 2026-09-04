local_optimizer_lean_root <- function() {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  skip_if_not(file.exists(file.path(root, "paper", "R", "03-optimizer-publication-lean.R")))
  root
}

local_optimizer_lean_env <- function() {
  root <- local_optimizer_lean_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-optimizer-publication-lean.R"), local = env)
  env
}

local_optimizer_lean_attempts <- function(env, reps = 97L) {
  design <- env$jss_optimizer_lean_cases()
  grid <- expand.grid(
    case_id = design$case_id,
    joint_review_rep = seq_len(reps),
    method = c("rs_separate", "rs_joint"),
    stringsAsFactors = FALSE
  )
  grid <- grid[order(match(grid$case_id, design$case_id), grid$joint_review_rep,
    match(grid$method, c("rs_separate", "rs_joint"))), , drop = FALSE]
  case_index <- match(grid$case_id, design$case_id)
  method_shift <- ifelse(grid$method == "rs_joint", 0.1, 0)
  data.frame(
    schema_version = env$jss_optimizer_lean_version(),
    case_id = grid$case_id,
    case_role = design$case_role[case_index],
    joint_review_rep = as.integer(grid$joint_review_rep),
    paired_seed = as.integer(17000L + case_index * 1000L + grid$joint_review_rep),
    method = grid$method,
    success = TRUE,
    converged = TRUE,
    retained = TRUE,
    failure_reason = "none",
    stop_reason = "converged",
    outer_iterations = 5L,
    elapsed_sec = 2 + method_shift + grid$joint_review_rep / 1000,
    train_joint_loglik = -100 + method_shift,
    test_log_score_per_obs = -1 + method_shift / 10,
    heldout_variogram_score_p05 = 0.5 - method_shift / 10,
    train_rmse_mu = 0.2 - method_shift / 10,
    train_rmse_sigma = 0.3 - method_shift / 10,
    train_rmse_theta = 0.4 - method_shift / 10,
    train_rmse_tau = 0.5 - method_shift / 10,
    publication_scope = "supplementary",
    stringsAsFactors = FALSE
  )
}

local_optimizer_lean_bundle <- function(env, attempts) {
  root <- tempfile("optimizer-lean-bundle-")
  dirs <- file.path(root, c("data", "tables", "figures"))
  invisible(lapply(dirs, dir.create, recursive = TRUE))
  files <- env$jss_optimizer_lean_files()
  pairs <- env$jss_optimizer_lean_pairs(attempts)
  utils::write.csv(env$jss_optimizer_lean_cases(), file.path(dirs[[1L]], files$data[[1L]]), row.names = FALSE)
  utils::write.csv(env$jss_optimizer_lean_precision(), file.path(dirs[[1L]], files$data[[2L]]), row.names = FALSE)
  utils::write.csv(attempts, file.path(dirs[[1L]], files$data[[3L]]), row.names = FALSE, na = "")
  utils::write.csv(pairs, file.path(dirs[[1L]], files$data[[4L]]), row.names = FALSE, na = "")
  utils::write.csv(env$jss_optimizer_lean_summary(attempts),
    file.path(dirs[[2L]], files$tables[[1L]]), row.names = FALSE)
  utils::write.csv(env$jss_optimizer_lean_failures(attempts),
    file.path(dirs[[2L]], files$tables[[2L]]), row.names = FALSE)
  uncertainty <- env$jss_optimizer_lean_uncertainty(pairs)
  utils::write.csv(uncertainty, file.path(dirs[[2L]], files$tables[[3L]]), row.names = FALSE, na = "")
  env$jss_optimizer_lean_write_figure(uncertainty, file.path(dirs[[3L]], files$figures[[1L]]))
  list(root = root, dirs = dirs, files = files)
}

test_that("lean optimizer publication contract is three-scenario and supplementary", {
  env <- local_optimizer_lean_env()
  expect_identical(env$jss_optimizer_lean_cases()$case_id, c("JVS01", "JVS02", "JVS03"))
  expect_identical(env$jss_optimizer_lean_cases()$case_role,
    c("base", "dependence_strength", "scale_complexity"))
  expect_false(any(env$jss_optimizer_lean_cases()$case_id == "JVS04"))
  expect_true(all(env$jss_optimizer_lean_cases()$publication_scope == "supplementary"))
  expect_identical(env$jss_optimizer_lean_precision()$target_usable_pairs, 97L)
  forbidden <- c("checkpoint", "signature", "sha256", "host", "worker_pid", "process")
  expect_false(any(grepl(paste(forbidden, collapse = "|"),
    unlist(env$jss_optimizer_lean_files()), ignore.case = TRUE)))
})

test_that("lean optimizer bundle round trips and recomputes paired uncertainty", {
  env <- local_optimizer_lean_env()
  attempts <- local_optimizer_lean_attempts(env)
  bundle <- local_optimizer_lean_bundle(env, attempts)
  expect_silent(env$jss_optimizer_lean_validate(bundle$dirs[[1L]], bundle$dirs[[2L]], bundle$dirs[[3L]]))
  pairs <- utils::read.csv(file.path(bundle$dirs[[1L]], bundle$files$data[[4L]]), stringsAsFactors = FALSE)
  uncertainty <- utils::read.csv(file.path(bundle$dirs[[2L]], bundle$files$tables[[3L]]), stringsAsFactors = FALSE)
  expect_equal(nrow(attempts), 582L)
  expect_equal(nrow(pairs), 291L)
  expect_true(all(pairs$pair_retained))
  expect_setequal(unique(uncertainty$comparison_id),
    c("JVS01", "JVS02", "JVS03", "JVS02-JVS01", "JVS03-JVS01"))
  expect_equal(nrow(uncertainty), 45L)
  expect_true(all(uncertainty$usable_pairs == 97L))
  expect_true(all(uncertainty$publication_scope == "supplementary"))
  expect_true(all(uncertainty$interpretation[uncertainty$metric == "elapsed_sec"] ==
    "descriptive_order_confounded"))
  expect_true(all(uncertainty$interpretation[uncertainty$metric == "convergence_rate"] ==
    "paired_attempt_convergence"))
})

test_that("lean optimizer validator rejects hidden failures, nonfinite values, and stale tables", {
  env <- local_optimizer_lean_env()
  attempts <- local_optimizer_lean_attempts(env)
  hidden <- attempts
  hidden$retained[[1L]] <- FALSE
  expect_error(env$jss_optimizer_lean_validate_attempts(
    hidden, env$jss_optimizer_lean_cases(), env$jss_optimizer_lean_precision()),
    "lifecycle|failure")
  nonfinite <- attempts
  nonfinite$elapsed_sec[[1L]] <- Inf
  expect_error(env$jss_optimizer_lean_validate_attempts(
    nonfinite, env$jss_optimizer_lean_cases(), env$jss_optimizer_lean_precision()),
    "impossible numeric")
  visible <- attempts
  visible$success[[1L]] <- FALSE
  visible$converged[[1L]] <- FALSE
  visible$retained[[1L]] <- FALSE
  visible$failure_reason[[1L]] <- "fit_error"
  expect_equal(env$jss_optimizer_lean_failures(visible)$failure_count[
    env$jss_optimizer_lean_failures(visible)$failure_reason == "fit_error"], 1L)
  expect_error(env$jss_optimizer_lean_validate_attempts(
    visible, env$jss_optimizer_lean_cases(), env$jss_optimizer_lean_precision()),
    "precision target")
  bundle <- local_optimizer_lean_bundle(env, attempts)
  summary_path <- file.path(bundle$dirs[[2L]], bundle$files$tables[[1L]])
  summary <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
  summary$retained[[1L]] <- summary$retained[[1L]] - 1L
  utils::write.csv(summary, summary_path, row.names = FALSE)
  expect_error(env$jss_optimizer_lean_validate(
    bundle$dirs[[1L]], bundle$dirs[[2L]], bundle$dirs[[3L]]),
    "lifecycle summary")
})

test_that("lean optimizer exporter excludes retired JVS04 checkpoints", {
  env <- local_optimizer_lean_env()
  make_raw <- function(case_id, rep) data.frame(
    case_id = case_id, joint_review_rep = rep, paired_seed = 1000L + rep,
    method = c("rs_separate", "rs_joint"), success = TRUE, converged = TRUE,
    stop_reason = "converged", outer_iterations = 2L, elapsed_sec = c(1, 2),
    train_joint_loglik = c(-2, -1), test_log_score_per_obs = c(-1, -.9),
    heldout_variogram_score_p05 = c(.5, .4), train_rmse_mu = c(.2, .1),
    train_rmse_sigma = c(.2, .1), train_rmse_theta = c(.2, .1),
    train_rmse_tau = c(.2, .1), failure_type = "none", error = NA_character_,
    stringsAsFactors = FALSE
  )
  checkpoint_dir <- tempfile("optimizer-checkpoints-")
  dir.create(checkpoint_dir)
  saveRDS(list(result = make_raw("JVS01", 1L)), file.path(checkpoint_dir, "JVS01.rds"))
  saveRDS(list(result = make_raw("JVS04", 1L)), file.path(checkpoint_dir, "JVS04.rds"))
  exported <- env$jss_optimizer_lean_read_checkpoints(checkpoint_dir)
  expect_identical(unique(exported$case_id), "JVS01")
  expect_false(any(exported$case_id == "JVS04"))
})

test_that("paper target graph sources the lean optimizer layer after legacy helpers", {
  root <- local_optimizer_lean_root()
  code <- readLines(file.path(root, "paper", "_targets.R"), warn = FALSE)
  legacy <- grep("phase2-paper-evidence[.]R", code)
  lean <- grep("03-optimizer-publication-lean[.]R", code)
  expect_length(legacy, 1L)
  expect_length(lean, 1L)
  expect_gt(lean, legacy)
  exporter <- readLines(file.path(root, "paper", "R", "run-03-optimizer-publication-lean.R"),
    warn = FALSE)
  expect_true(any(grepl("jss_optimizer_lean_publish", exporter, fixed = TRUE)))
})

test_that("full-profile tolerance snapshots consume the lean uncertainty schema", {
  env <- local_optimizer_lean_env()
  root <- local_optimizer_lean_root()
  helpers <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "replication-helpers.R"), local = helpers)
  expect_false(exists("jss_phase2_public_approval_spec", envir = helpers, inherits = FALSE))
  signing_helper <- paste(readLines(file.path(root, "paper", "scripts", "phase2-evidence-approval.R"),
    warn = FALSE), collapse = "\n")
  expect_false(grepl('identical\\(study, "optimizer-benchmark"\\)', signing_helper))
  attempts <- local_optimizer_lean_attempts(env)
  bundle <- local_optimizer_lean_bundle(env, attempts)
  uncertainty_path <- file.path(bundle$dirs[[2L]], bundle$files$tables[[3L]])
  specs <- rbind(
    helpers$jss_metric_spec("estimate", tolerance_group = "paired_effect"),
    helpers$jss_metric_spec("mcse", tolerance_group = "paired_effect")
  )
  snapshots <- helpers$jss_metric_snapshots_from_csv(
    uncertainty_path, "optimizer_benchmark",
    c("comparison_id", "comparison_type", "metric"), specs
  )
  expect_equal(nrow(snapshots), 90L)
  expect_setequal(unique(snapshots$metric), c("estimate", "mcse"))
  expect_true(all(snapshots$tolerance_group == "paired_effect"))
})
