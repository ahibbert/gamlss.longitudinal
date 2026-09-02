local_phase2_repo_root <- function() {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  skip_if_not(file.exists(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R")))
  root
}

local_optimizer_manuscript_hits <- function(text) {
  normalized <- gsub("[[:space:]]+", " ", text)
  fixed_phrases <- c(
    cg_rs_steps_faster = "steps are substantially faster to run",
    cg_more_efficient_slower = "generally more efficient but slower to run each iteration",
    cg_speed_workflow = "use separate RS for exploratory screens when speed matters",
    cg_joint_outperform = "CG and RS joint (with cross terms), either match or outperform",
    joint_substantially_slower = "joint models is that they are substantially slower",
    separate_then_refine = "RS separate can be used as a",
    stale_six_examples = "six example cases",
    stale_broad_grid = "selected from a broad grid",
    stale_fixed_100 = "100 replications for each example",
    stale_fixed_final_attempts = "the number of attempts per cell is determined before fitting",
    unregistered_reasonable_fit = "will generally fit reasonably",
    unregistered_accurate_capture = "can be captured very accurately",
    retired_normal_table = "tables/normal-joint-vs-separate-six-case-median-iqr-table",
    retired_gamma_table = "tables/gamma-joint-vs-separate-six-case-median-iqr-table",
    retired_nbi_table = "tables/negative-binomial-joint-vs-separate-six-case-median-iqr-table"
  )
  exact_hits <- names(fixed_phrases)[vapply(
    fixed_phrases, grepl, logical(1), x = normalized, fixed = TRUE
  )]
  directional <- paste(
    "outperform(?:s|ed)?|underperform(?:s|ed)?|better|worse|superior|inferior|",
    "advantage|improv(?:e|es|ed)|increas(?:e|es|ed)|decreas(?:e|es|ed)|",
    "higher|lower|faster|slower|match(?:es|ed)?",
    sep = ""
  )
  semantic_patterns <- c(
    cg_followed_by_performance = paste0(
      "(?i)\\bCG\\b.{0,300}\\b(faster|slower|efficient|outperform|robust|",
      "difficult(?:y|ies)?|stuck|local opt(?:imum|ima)?)\\b"
    ),
    performance_followed_by_cg = paste0(
      "(?i)\\b(faster|slower|efficient|outperform|robust|",
      "difficult(?:y|ies)?|stuck|local opt(?:imum|ima)?)\\b.{0,300}\\bCG\\b"
    ),
    stale_fixed_100_optimizer = paste0(
      "(?i)\\b(optimizer|joint.{0,30}separate|separate.{0,30}joint)\\b.{0,300}",
      "\\b100[ -](replications|replicates|attempts)\\b"
    ),
    stale_fixed_100_optimizer_reverse = paste0(
      "(?i)\\b100[ -](replications|replicates|attempts)\\b.{0,300}",
      "\\b(optimizer|joint.{0,30}separate|separate.{0,30}joint)\\b"
    ),
    joint_directional_vs_separate = paste0(
      "(?i)\\b(?:joint(?:ly)?[ -](?:optim(?:i[sz](?:e|ed)|isation|ization)|RS)|joint model)\\b",
      ".{0,180}\\b(", directional, ")\\b.{0,180}",
      "\\b(?:separat(?:e|ely|ed)[ -](?:optim(?:i[sz](?:e|ed)|isation|ization)|RS)|separate model)\\b"
    ),
    separate_directional_vs_joint = paste0(
      "(?i)\\b(?:separat(?:e|ely|ed)[ -](?:optim(?:i[sz](?:e|ed)|isation|ization)|RS)|separate model)\\b",
      ".{0,180}\\b(", directional, ")\\b.{0,180}",
      "\\b(?:joint(?:ly)?[ -](?:optim(?:i[sz](?:e|ed)|isation|ization)|RS)|joint model)\\b"
    )
  )
  semantic_hits <- names(semantic_patterns)[vapply(
    semantic_patterns, grepl, logical(1), x = normalized, perl = TRUE
  )]
  unique(c(exact_hits, semantic_hits))
}

local_missingness_manuscript_hits <- function(text) {
  normalized <- gsub("[[:space:]]+", " ", text)
  patterns <- c(
    stale_appendix_wording = "(?i)missingness at random and by time",
    obsolete_time_mar = "(?i)\\btime_mar\\b",
    unfrozen_missingness_figure = "(?i)(fixed_margin_rmse_by_missingness|smooth_selected_recovery_curves)",
    directional_missingness_claim = paste0(
      "(?i)\\b(missingness|dropout)\\b.{0,220}\\b(outperform|better|worse|robust|",
      "improv(?:e|es|ed)|degrad(?:e|es|ed)|superior|inferior)\\b"
    )
  )
  names(patterns)[vapply(patterns, grepl, logical(1), x = normalized, perl = TRUE)]
}

test_that("optimizer evidence uses registered paired one-factor contrasts", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)

  cases <- env$jss_joint_case_definitions()
  expect_identical(cases$case_id, sprintf("JVS%02d", 1:4))
  expect_identical(cases$base_case_id, rep("JVS01", 4L))
  expect_identical(cases$family, c("NO", "NO", "NO", "ZIP"))
  expect_identical(cases$copula, rep("N", 4L))
  jvs04 <- cases[cases$case_id == "JVS04", , drop = FALSE]
  expect_identical(jvs04$contrast_factor, "family")
  expect_identical(jvs04$hypothesis_role, "count_margin")
  expect_equal(jvs04$mu_natural_reference, 4)
  expect_equal(jvs04$sigma_natural_reference, 0.20)
  expect_identical(jvs04$sigma_interpretation, "zero_inflation_probability")
  expect_equal(jvs04$mu_intercept_eta, log(4))
  expect_equal(jvs04$sigma_intercept_eta, stats::qlogis(0.20))
  zip <- env$jss_joint_margin_dist("ZIP")
  expect_equal(zip$mu.linkinv(jvs04$mu_intercept_eta), 4)
  expect_equal(zip$sigma.linkinv(jvs04$sigma_intercept_eta), 0.20)
  expect_true(env$jss_joint_validate_case_definitions(cases))
  expect_true(env$jss_joint_validate_registered_routes(cases))

  registry <- env$jss_joint_mc_precision_registry("full")
  expect_equal(registry$confidence_level, 0.95)
  expect_equal(registry$target_half_width, 0.10)
  expect_identical(registry$required_retained_pairs, 97L)
  expect_identical(registry$max_attempts, 194L)
  expect_identical(env$jss_joint_simulation_settings(list(profile = "smoke"))$method, "RS")

  invalid <- cases
  invalid$n[invalid$case_id == "JVS02"] <- invalid$n[1] + 1L
  expect_error(env$jss_joint_validate_case_definitions(invalid), "exactly its declared factor")

  invalid_truth <- cases
  invalid_truth$sigma_natural_reference[invalid_truth$case_id == "JVS04"] <- 0.30
  expect_error(env$jss_joint_validate_case_definitions(invalid_truth), "registered family truth")
})

test_that("optimizer route preflight rejects unsupported registered routes synchronously", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)

  cases <- env$jss_joint_case_definitions()
  invalid_route <- cases
  invalid_route$copula <- "C"
  expect_error(
    env$jss_joint_validate_registered_routes(invalid_route),
    "tested allowlist",
    class = "gamlss_longitudinal_unsupported_route_error"
  )

  out <- tempfile("jvs-route-preflight-")
  calls <- 0L
  mock_run <- function(...) {
    calls <<- calls + 1L
    stop("worker must not run", call. = FALSE)
  }
  settings <- list(profile = "smoke", seed = 101L, out_dir = out, root = root)
  expect_error(
    env$jss_joint_run_simulation(
      settings, invalid_route, list(reps = 1L, resume = TRUE), mock_run
    ),
    "tested allowlist",
    class = "gamlss_longitudinal_unsupported_route_error"
  )
  expect_identical(calls, 0L)
  expect_false(dir.exists(paste0(env$jss_joint_checkpoint_dir(settings), ".run-lock")))
  expect_false(dir.exists(env$jss_joint_checkpoint_dir(settings)))
  expect_error(
    env$jss_joint_run_simulation_fixed(
      settings, invalid_route, list(reps = 1L, resume = TRUE), mock_run,
      acquire_lock = FALSE, validate_design = FALSE
    ),
    "tested allowlist",
    class = "gamlss_longitudinal_unsupported_route_error"
  )
  expect_identical(calls, 0L)
})

test_that("optimizer ZIP/Gaussian case runs both RS modes with analytical discrete scoring", {
  skip_if_not_installed("gamlss.dist")
  skip_if_not_installed("scoringRules")
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)

  case <- env$jss_joint_case_definitions()[4L, , drop = FALSE]
  case$n <- 20L
  case$time_points <- 3L
  case$total_observations <- 60L
  cfg <- list(
    method = "RS", outer_stop_crit = 1, inner_stop_crit = 1,
    max_outer_iter = 2L, max_inner_iter = 2L, max_elapsed_sec = 60,
    warm_start_joint = TRUE, warm_start_joint_iter = 1L,
    discrete_score_method = "analytical", variogram_nsim = 1L
  )
  result <- env$jss_joint_run_case_rep(case, 1L, list(seed = 101L), cfg)

  expect_identical(result$method, c("rs_separate", "rs_joint"))
  expect_true(all(result$family == "ZIP" & result$copula == "N"))
  expect_true(all(result$success & result$converged & result$retained))
  expect_true(all(is.na(result$error)))
  expect_true(all(is.finite(result$test_log_score_per_obs)))
})

test_that("optimizer checkpoints are validated and resumed without rerunning cells", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)

  out <- tempfile("jvs-checkpoints-")
  dir.create(out)
  settings <- list(profile = "smoke", seed = 101L, out_dir = out, root = root)
  cases <- env$jss_joint_case_definitions()[1:2, , drop = FALSE]
  cfg <- list(reps = 2L, resume = TRUE)
  calls <- 0L
  contract_fixture <- env$jss_joint_contract_fixture
  mock_run <- function(case, rep_idx, settings, cfg) {
    calls <<- calls + 1L
    contract_fixture(case, rep_idx, settings)
  }

  first <- env$jss_joint_run_simulation(settings, cases, cfg, mock_run)
  expect_equal(calls, 4L)
  expect_false(any(attr(first, "checkpoint_status")$resumed))
  calls <- 0L
  second <- env$jss_joint_run_simulation(settings, cases, cfg, mock_run)
  expect_equal(calls, 0L)
  expect_true(all(attr(second, "checkpoint_status")$resumed))
  expect_equal(second, first, ignore_attr = TRUE)

  calls <- 0L
  changed_seed <- settings
  changed_seed$seed <- settings$seed + 1L
  third <- env$jss_joint_run_simulation(changed_seed, cases, cfg, mock_run)
  expect_equal(calls, 4L)
  expect_false(any(attr(third, "checkpoint_status")$resumed))

  calls <- 0L
  changed_cfg <- cfg
  changed_cfg$max_outer_iter <- 99L
  fourth <- env$jss_joint_run_simulation(changed_seed, cases, changed_cfg, mock_run)
  expect_equal(calls, 4L)
  expect_false(any(attr(fourth, "checkpoint_status")$resumed))

  forged_path <- env$jss_joint_checkpoint_path(changed_seed, cases[1, , drop = FALSE], 1L, changed_cfg)
  forged <- readRDS(forged_path)
  forged$result$converged[[1L]] <- FALSE
  forged$result$retained[[1L]] <- TRUE
  saveRDS(forged, forged_path)
  expect_false(env$jss_joint_checkpoint_valid(
    forged_path, changed_seed, cases[1, , drop = FALSE], 1L, changed_cfg
  ))
})

test_that("optimizer precision is conditional on usable pairs and hard-gated", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)
  registry <- env$jss_joint_mc_precision_registry("full")
  make_delta <- function(n, retained) data.frame(
    case_id = "JVS01", rs_joint_success = TRUE, rs_joint_converged = retained,
    rs_separate_success = TRUE, rs_separate_converged = retained,
    joint_review_rep = seq_len(n), stringsAsFactors = FALSE
  )
  failed <- env$jss_joint_precision_achievement(make_delta(97L, c(rep(TRUE, 96L), FALSE)), registry)
  achieved <- env$jss_joint_precision_achievement(make_delta(97L, rep(TRUE, 97L)), registry)
  expect_false(failed$precision_met)
  expect_true(achieved$precision_met)
  expect_equal(registry$max_attempts, 194L)
  withr::local_envvar(GAMLSS_LONGITUDINAL_JSS_JVS_MC_HALF_WIDTH = "Inf")
  expect_error(env$jss_joint_mc_precision_registry("full"), "finite and strictly positive")
})

test_that("optimizer interval decisions use literal zero with exact boundary behavior", {
  root <- local_phase2_repo_root(); env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)
  epsilon <- .Machine$double.eps
  expect_identical(env$jss_joint_interval_direction(0, 1), 0L)
  expect_identical(env$jss_joint_interval_direction(-1, 0), 0L)
  expect_identical(env$jss_joint_interval_direction(epsilon, 1), 1L)
  expect_identical(env$jss_joint_interval_direction(-1, -epsilon), -1L)
})

test_that("optimizer adaptively tops up only deficient paired cells", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)
  cases <- env$jss_joint_case_definitions()[1:2, , drop = FALSE]
  out <- tempfile("jvs-adaptive-"); dir.create(out)
  settings <- list(profile = "smoke", seed = 909L, out_dir = out, root = root,
    logs_dir = file.path(out, "logs"))
  precision <- env$jss_joint_mc_precision_registry("smoke")
  precision$required_retained_pairs <- 2L; precision$initial_attempts <- 2L
  precision$top_up_batch <- 1L; precision$max_attempts <- 4L
  precision$target_half_width <- .8
  calls <- character(); fixture <- env$jss_joint_contract_fixture
  run <- function(case, rep_idx, settings, cfg) {
    calls <<- c(calls, paste(case$case_id[[1L]], rep_idx))
    x <- fixture(case, rep_idx, settings)
    if (case$case_id[[1L]] == "JVS02" && rep_idx == 1L) {
      x$converged <- FALSE; x$retained <- FALSE
      x$stop_reason <- "max_iterations"
      x$failure_type <- "optimizer_nonconvergence:max_iterations"
    }
    x
  }
  result <- env$jss_joint_run_simulation(settings, cases,
    list(reps = 4L, resume = TRUE, workers = 1L, precision = precision), run)
  expect_equal(table(result$case_id, result$method)["JVS01", ], c(rs_joint = 2L, rs_separate = 2L))
  expect_equal(table(result$case_id, result$method)["JVS02", ], c(rs_joint = 3L, rs_separate = 3L))
  expect_equal(length(calls), 5L)
  expect_equal(sum(grepl("JVS01", calls)), 2L)
  expect_true(all(attr(result, "precision_achievement")$precision_met))
})

test_that("optimizer claim evidence separates the fixed initial plan from deficient top-ups", {
  root <- local_phase2_repo_root(); env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "phase2-paper-evidence.R"), local = env)
  cases <- sprintf("JVS%02d", 1:4)
  attempts <- do.call(rbind, lapply(seq_along(cases), function(i) {
    reps <- if (i == 4L) 107L else 97L
    expand.grid(case_id = cases[[i]], joint_review_rep = seq_len(reps),
      method = c("rs_separate", "rs_joint"), stringsAsFactors = FALSE)
  }))
  attempts$success <- TRUE; attempts$converged <- TRUE
  estimates <- env$jss_phase2_optimizer_claim_estimates(attempts,
    data.frame(initial_attempts = 97L))
  expect_identical(unname(estimates[["registered_scenario_cells"]]), 4L)
  expect_identical(unname(estimates[["registered_method_cells"]]), 8L)
  expect_identical(unname(estimates[["planned_initial_pairs"]]), 388L)
  expect_identical(unname(estimates[["actual_attempt_rows"]]), 796L)
  expect_identical(unname(estimates[["achieved_min_retained_pairs"]]), 97L)
  expect_false(estimates[["actual_attempt_rows"]] == 2L * estimates[["planned_initial_pairs"]])
})

test_that("optimizer output lock admits exactly one concurrent launcher", {
  root <- local_phase2_repo_root()
  out <- tempfile("optimizer-concurrent-output-")
  dir.create(out)
  settings <- list(out_dir = out)
  cl <- parallel::makePSOCKcluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  result <- parallel::clusterCall(cl, function(settings, source_file) {
    source(source_file, local = .GlobalEnv)
    tryCatch({
      lock <- jss_joint_acquire_run_lock(settings)
      on.exit(jss_joint_release_run_lock(lock), add = TRUE)
      Sys.sleep(0.5)
      TRUE
    }, error = function(e) FALSE)
  }, settings = settings,
  source_file = file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"))
  expect_identical(sort(unlist(result)), c(FALSE, TRUE))
  expect_false(dir.exists(paste0(file.path(out, "checkpoints", "03-joint-vs-separate-optimization",
    "paired-one-factor-v7"), ".run-lock")))
})

test_that("optimizer checkpoint contract rejects forged metrics and labels", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)
  settings <- list(profile = "smoke", seed = 1L, out_dir = tempfile("jvs-forgery-"), root = root)
  case <- env$jss_joint_case_definitions()[1, , drop = FALSE]
  cfg <- list(reps = 1L, resume = TRUE, workers = 1L)
  valid <- env$jss_joint_contract_fixture(case, 1L, settings)
  path <- env$jss_joint_checkpoint_path(settings, case, 1L, cfg)
  forged <- valid
  forged$train_joint_loglik[[1L]] <- 1e300
  expect_error(env$jss_joint_write_checkpoint(forged, path, settings, case, 1L, cfg), "out of contract")
  forged <- valid
  forged$converged[[1L]] <- FALSE
  forged$retained[[1L]] <- FALSE
  forged$stop_reason[[1L]] <- "made_up_stop"
  forged$failure_type[[1L]] <- "made_up_failure"
  expect_error(env$jss_joint_write_checkpoint(forged, path, settings, case, 1L, cfg), "unregistered")
  forged <- valid
  forged$joint_review_rep <- as.character(forged$joint_review_rep)
  expect_error(env$jss_joint_write_checkpoint(forged, path, settings, case, 1L, cfg), "integer typed")
  forged <- valid
  forged$train_joint_loglik <- as.character(forged$train_joint_loglik)
  expect_error(env$jss_joint_write_checkpoint(forged, path, settings, case, 1L, cfg),
    "must be numeric")
  forged <- valid
  forged$fabricated_note <- "not in checkpoint schema"
  expect_error(env$jss_joint_write_checkpoint(forged, path, settings, case, 1L, cfg),
    "extra, missing, or reordered")
  forged <- valid
  forged$success[[1L]] <- FALSE
  expect_error(env$jss_joint_write_checkpoint(forged, path, settings, case, 1L, cfg),
    "cannot be marked unsuccessful|unregistered")
})

test_that("optimizer public figure guard rejects dummy PNG payloads", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "phase2-paper-evidence.R"), local = env)
  good <- tempfile(fileext = ".png")
  grDevices::png(good, width = 400, height = 300)
  graphics::plot(1:10, (1:10)^2)
  grDevices::dev.off()
  expect_true(env$jss_phase2_validate_png(good))
  dummy <- tempfile(fileext = ".png")
  writeBin(as.raw(c(137, 80, 78, 71, 13, 10, 26, 10, rep(0, 2000))), dummy)
  expect_false(env$jss_phase2_validate_png(dummy))
})

test_that("production attestation rejects self-forged signatures and invalid UTC order", {
  root <- local_phase2_repo_root(); env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)
  source(file.path(root, "paper", "R", "phase2-paper-evidence.R"), local = env)
  checkpoints <- data.frame(task_id = 1:2, case_id = c("JVS01", "JVS02"),
    checkpoint = c("checkpoints/a.rds", "checkpoints/b.rds"),
    result_content_sha256 = c(strrep("a", 64), strrep("b", 64)),
    checkpoint_timestamp_utc = c("2026-09-01T00:00:00Z", "2026-09-01T00:00:01Z"),
    stringsAsFactors = FALSE)
  make_attestation <- function(manifest = checkpoints, approved = "2026-09-01T00:00:02Z") list(
    schema_version = 1L, study = "optimizer-benchmark", bundle_sha256 = strrep("c", 64),
    package_source_sha256 = strrep("d", 64), producer_sha256 = strrep("e", 64),
    approved_at_utc = approved, approver = "external-reviewer", checkpoint_manifest = manifest)
  path <- tempfile(fileext = ".rds")
  signature <- tempfile(fileext = ".sig")
  check <- function(x) {
    writeBin(serialize(x, NULL, version = 3L), path)
    key <- sodium::sig_keygen()
    message <- readBin(path, "raw", n = file.info(path)$size)
    writeBin(sodium::sig_sign(message, key), signature)
    env$jss_phase2_require_external_attestation("optimizer-benchmark", strrep("c", 64),
      strrep("d", 64), strrep("e", 64), checkpoints, root, path,
      signature)
  }
  expect_error(check(make_attestation()), "valid detached production signature")
  expect_error(env$jss_phase2_validate_approval_order("2026-08-31T23:59:59Z",
    checkpoints$checkpoint_timestamp_utc), "after every checkpoint")
  expect_error(env$jss_phase2_parse_rfc3339_utc("2026-02-30T00:00:00Z"),
    "not a real RFC3339")
  expect_false(identical(checkpoints, checkpoints[2:1, ]))

  missing_env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "04-missingness-dropout-sensitivity.R"), local = missing_env)
  missing_checkpoints <- data.frame(timestamp_utc = "2026-09-01T00:00:01Z")
  missing_attestation <- make_attestation(missing_checkpoints)
  missing_attestation$study <- "missingness"
  names(missing_attestation) <- c("schema_version", "study", "bundle_sha256",
    "package_source_sha256", "producer_sha256", "approved_at_utc", "approver",
    "checkpoint_manifest")
  writeBin(serialize(missing_attestation, NULL, version = 3L), path)
  test_key <- sodium::sig_keygen()
  message <- readBin(path, "raw", n = file.info(path)$size)
  writeBin(sodium::sig_sign(message, test_key), signature)
  expect_error(missing_env$jss_missing_require_external_attestation(strrep("c", 64),
    strrep("d", 64), strrep("e", 64), missing_checkpoints, root, path, signature),
    "valid detached production signature")
})

test_that("approval helper exposes only the pinned public key and documents independent rotation", {
  root <- local_phase2_repo_root(); env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "scripts", "phase2-evidence-approval.R"), local = env)
  expect_identical(digest::digest(env$phase2_pinned_public_key(), "sha256", serialize = FALSE),
    "cb73c05cede55bfd56357b1780b90c1bf254413d82765c7a682fc9db3a0d8587")
  expect_false(identical(sodium::sig_pubkey(sodium::sig_keygen()), env$phase2_pinned_public_key()))
  ceremony <- paste(readLines(file.path(root, "paper", "notes", "phase2-evidence-signing.md"),
    warn = FALSE), collapse = "\n")
  expect_match(ceremony, "access-restricted, checkout-external", fixed = TRUE)
  expect_match(ceremony, "all five validators", fixed = TRUE)
  expect_match(ceremony, "code-only change", fixed = TRUE)
  expect_match(ceremony, "separate data-only review/commit", fixed = TRUE)
})

test_that("public validators detect source mutation while copying immutable snapshots", {
  root <- local_phase2_repo_root()
  optimizer <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = optimizer)
  source(file.path(root, "paper", "R", "phase2-paper-evidence.R"), local = optimizer)
  base <- tempfile("optimizer-toctou-"); dirs <- file.path(base, c("data", "tables", "figures"))
  invisible(lapply(dirs, dir.create, recursive = TRUE))
  files <- optimizer$jss_phase2_optimizer_files()
  paths <- c(file.path(dirs[[1L]], files$data), file.path(dirs[[2L]], files$tables),
    file.path(dirs[[3L]], files$figures))
  invisible(vapply(paths, function(path) file.create(path), logical(1)))
  real_sha <- optimizer$jss_joint_sha256_file; calls <- 0L
  optimizer$jss_joint_sha256_file <- function(path) {
    calls <<- calls + 1L
    if (calls == length(paths) + 1L) write("mutation", paths[[1L]], append = TRUE)
    real_sha(path)
  }
  expect_error(optimizer$jss_phase2_optimizer_validate(dirs[[1L]], dirs[[2L]], dirs[[3L]], root),
    "changed while creating")

  missing <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "04-missingness-dropout-sensitivity.R"), local = missing)
  input <- tempfile("missingness-toctou-"); dir.create(input)
  missing_paths <- file.path(input, missing$jss_missingness_public_artifacts())
  invisible(vapply(missing_paths, function(path) file.create(path), logical(1)))
  old_options <- options(jss_missing_toctou_calls = 0L,
    jss_missing_toctou_paths = missing_paths)
  on.exit(options(old_options), add = TRUE)
  mutating_missing_sha <- function(path) {
    calls <- getOption("jss_missing_toctou_calls") + 1L
    options(jss_missing_toctou_calls = calls)
    paths <- getOption("jss_missing_toctou_paths")
    if (calls == length(paths) + 1L) write("mutation", paths[[1L]], append = TRUE)
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }
  missing$jss_missing_sha256_file <- mutating_missing_sha
  expect_error(missing$jss_missingness_validate_public_bundle(input, root), "changed while creating")
})

test_that("full 4x97 optimizer public bundle validates across named-key CSV boundaries", {
  root <- local_phase2_repo_root(); env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)
  source(file.path(root, "paper", "R", "phase2-paper-evidence.R"), local = env)
  base <- tempfile("optimizer-public-"); data_dir <- file.path(base, "data")
  table_dir <- file.path(base, "tables"); figure_dir <- file.path(base, "figures")
  dir.create(data_dir, recursive = TRUE); dir.create(table_dir); dir.create(figure_dir)
  settings <- list(seed = 17L); cases <- env$jss_joint_case_definitions()
  precision <- env$jss_joint_mc_precision_registry("full")
  design <- cases; for (nm in names(precision)) design[[nm]] <- precision[[nm]][[1L]]
  rows <- list(); archive <- list(); k <- 0L
  for (i in seq_len(nrow(cases))) for (rep in seq_len(97L)) {
    k <- k + 1L; x <- env$jss_joint_contract_fixture(cases[i, , drop = FALSE], rep, settings,
      deterministic_offset = i + rep / 1000)
    x$elapsed_sec <- 1 + i / 10 + rep / 10000 + c(0, .001)
    archive[[k]] <- x; portable <- env$jss_joint_portable_result_sha256(x)
    x$result_portable_sha256 <- portable; x$profile <- "full"; rows[[k]] <- x
  }
  attempts <- do.call(rbind, rows); rownames(attempts) <- NULL
  deltas <- env$jss_joint_delta_table(attempts)
  raw_deltas <- deltas
  achieved <- env$jss_joint_precision_achievement(deltas, precision)
  precision$achieved_min_retained_pairs <- min(achieved$retained_pairs)
  precision$achieved_max_half_width <- max(achieved$achieved_worst_case_half_width)
  precision$all_cells_precision_met <- all(achieved$precision_met)
  paths <- env$jss_joint_output_paths(list(data_dir = data_dir, tables_dir = table_dir,
    figures_dir = figure_dir))
  write.csv(design, paths$candidate_selection, row.names = FALSE)
  write.csv(precision, paths$precision_registry, row.names = FALSE)
  write.csv(attempts, paths$results, row.names = FALSE)
  saveRDS(archive, paths$checkpoint_payload_archive, version = 3L)
  manifest <- do.call(rbind, lapply(archive, function(x) data.frame(case_id = x$case_id[[1L]],
    joint_review_rep = x$joint_review_rep[[1L]], result_content_sha256 = env$jss_joint_content_sha256(x),
    result_portable_sha256 = env$jss_joint_portable_result_sha256(x))))
  write.csv(manifest, paths$checkpoint_content_manifest, row.names = FALSE)
  identity <- env$jss_joint_checkout_package_identity(list(root = root))
  producer <- env$jss_joint_producer_fingerprint(list(root = root))
  runtime <- env$jss_joint_runtime_identity()
  rlibs_sha <- digest::digest(gsub("\\\\", "/", Sys.getenv("R_LIBS_USER", unset = "")),
    "sha256", serialize = FALSE)
  libpaths_sha <- digest::digest(paste(gsub("\\\\", "/", .libPaths()), collapse = ";"),
    "sha256", serialize = FALSE)
  status <- data.frame(task_id = seq_len(388L), case_id = manifest$case_id,
    joint_review_rep = manifest$joint_review_rep,
    paired_seed = 17L + 3000L + manifest$joint_review_rep * 100L,
    checkpoint = sprintf("checkpoints/03-joint-vs-separate-optimization/paired-one-factor-v7/spec/%s-rep-%d.rds",
      manifest$case_id, manifest$joint_review_rep), checkpoint_schema_version = 7L,
    result_content_sha256 = manifest$result_content_sha256,
    result_portable_sha256 = manifest$result_portable_sha256,
    producer_fingerprint = producer, producer_fingerprint_algorithm = "SHA-256",
    package_version = identity$version, package_fingerprint_scope = identity$fingerprint_scope,
    package_source_file_count = identity$source_file_count, package_checkout_path = "checkout",
    package_source_sha256 = identity$source_sha256, verified_package_path = "checkout",
    verified_package_version = identity$version, verified_source_sha256 = identity$source_sha256,
    verified_producer_path = "checkout/paper/R/03-joint-vs-separate-optimization.R",
    verified_producer_sha256 = producer, package_identity_verified = TRUE,
    verified_rlibs_user = rlibs_sha,
    verified_libpaths = libpaths_sha, package_load_strategy = "pkgload_checkout",
    resumed = FALSE, executed = TRUE, workers_requested = 2L, workers_used = 2L,
    execution_mode = "psock", worker_pid = 1000L + seq_len(388L),
    checkpoint_timestamp_utc = "2026-09-01T00:00:00Z", execution_host = runtime$host,
    execution_os = runtime$os, execution_r_version = runtime$r_version,
    execution_platform = runtime$platform,
    execution_rng_kind = runtime$rng_kind, execution_blas = runtime$blas,
    execution_lapack = runtime$lapack, adaptive_round = 1L,
    registered_initial_attempts = 97L, registered_top_up_batch = 10L,
    registered_hard_cap = 194L, actual_attempts_in_cell = 97L, stringsAsFactors = FALSE)
  write.csv(status, paths$checkpoint_status, row.names = FALSE)
  canonical <- env$jss_joint_read_canonical_results_csv(paths$results)
  expect_true(is.character(canonical$warnings))
  expect_true(is.character(canonical$error))
  expect_true(all(canonical$warnings == ""))
  expect_true(all(is.na(canonical$error)))
  deltas <- env$jss_joint_delta_table(canonical)
  expect_equal(deltas, raw_deltas, tolerance = 1e-12)
  summary <- env$jss_joint_summary_table(canonical, deltas, design)
  wins <- env$jss_joint_metric_wins(canonical); failures <- env$jss_joint_failure_summary(canonical)
  uncertainty <- env$jss_joint_difference_uncertainty(deltas)
  hypotheses <- env$jss_joint_hypothesis_summary(deltas)
  evidence <- env$jss_joint_hypothesis_evidence(hypotheses, uncertainty)
  write.csv(summary, paths$summary, row.names = FALSE); write.csv(wins, paths$metric_wins, row.names = FALSE)
  write.csv(hypotheses, paths$hypothesis_summary, row.names = FALSE)
  write.csv(failures, paths$failure_summary, row.names = FALSE)
  write.csv(uncertainty, paths$difference_uncertainty, row.names = FALSE)
  write.csv(evidence, paths$hypothesis_evidence, row.names = FALSE)
  env$jss_joint_write_delta_figure(deltas, paths$deltas_figure)
  env$jss_joint_write_metric_dashboard(wins, paths$metric_dashboard)
  registry <- data.frame(figure = basename(c(paths$deltas_figure, paths$metric_dashboard)),
    png_sha256 = vapply(c(paths$deltas_figure, paths$metric_dashboard), env$jss_joint_sha256_file, character(1)),
    plotted_data_sha256 = c(env$jss_joint_portable_frame_sha256(deltas), env$jss_joint_portable_frame_sha256(wins)),
    plot_spec_sha256 = c(digest::digest("delta density+zero line+case facet; width=10 height=7 dpi=320", "sha256", serialize = FALSE),
      digest::digest("metric win-or-tie bars; width=9 height=6 dpi=320", "sha256", serialize = FALSE)))
  write.csv(registry, paths$figure_registry, row.names = FALSE)
  files <- env$jss_phase2_optimizer_files()
  required <- c(file.path(data_dir, files$data), file.path(table_dir, files$tables), file.path(figure_dir, files$figures))
  validate_candidate <- function() env$jss_phase2_optimizer_validate_candidate(
    data_dir, table_dir, figure_dir, root, require_promotion = FALSE)
  expect_silent(validate_candidate())
  forged_runtime <- status
  forged_runtime$execution_host[[1L]] <- paste0(runtime$host, "-forged")
  write.csv(forged_runtime, paths$checkpoint_status, row.names = FALSE)
  expect_error(validate_candidate(), "runtime provenance")
  write.csv(status, paths$checkpoint_status, row.names = FALSE)
  expect_error(env$jss_phase2_optimizer_validate(data_dir, table_dir, figure_dir, root),
    "detached production promotion signature")
  named <- canonical; names(named$method) <- paste0("row", seq_len(nrow(named)))
  expect_silent(env$jss_phase2_compare_frame(named[1:4, ], canonical[1:4, ],
    c("case_id", "method"), "named key"))
  forged <- canonical; forged$fabricated_note <- "self-approved"
  write.csv(forged, paths$results, row.names = FALSE)
  expect_error(validate_candidate(),
    "extra or missing columns")
  write.csv(canonical, paths$results, row.names = FALSE)

  forged_status <- status
  forged_status$checkpoint[[1L]] <- normalizePath(paths$results, winslash = "/", mustWork = TRUE)
  write.csv(forged_status, paths$checkpoint_status, row.names = FALSE)
  expect_error(validate_candidate(),
    "portable paths|canonical checkpoint paths")
  write.csv(status, paths$checkpoint_status, row.names = FALSE)

  forged_status <- status
  forged_status$checkpoint[1:2] <- rev(forged_status$checkpoint[1:2])
  write.csv(forged_status, paths$checkpoint_status, row.names = FALSE)
  expect_error(validate_candidate(), "canonical checkpoint paths")
  write.csv(status, paths$checkpoint_status, row.names = FALSE)

  forged_status <- status
  forged_status$adaptive_round[[98L]] <- 3L
  write.csv(forged_status, paths$checkpoint_status, row.names = FALSE)
  expect_error(validate_candidate(), "precision eligibility")
  write.csv(status, paths$checkpoint_status, row.names = FALSE)

  forged_manifest <- manifest
  forged_manifest$result_content_sha256[[1L]] <- strrep("d", 64L)
  write.csv(forged_manifest, paths$checkpoint_content_manifest, row.names = FALSE)
  expect_error(validate_candidate(),
    "content manifest")
  write.csv(manifest, paths$checkpoint_content_manifest, row.names = FALSE)

  original_figure <- tempfile(fileext = ".png")
  file.copy(paths$deltas_figure, original_figure, overwrite = TRUE)
  grDevices::png(paths$deltas_figure, width = 3200, height = 2240, res = 320)
  graphics::plot.new(); graphics::title("placeholder")
  grDevices::dev.off()
  placeholder_registry <- registry
  placeholder_registry$png_sha256[placeholder_registry$figure == basename(paths$deltas_figure)] <-
    env$jss_joint_sha256_file(paths$deltas_figure)
  write.csv(placeholder_registry, paths$figure_registry, row.names = FALSE)
  expect_error(validate_candidate(),
    "deterministic rerendering")
  file.copy(original_figure, paths$deltas_figure, overwrite = TRUE)
  write.csv(registry, paths$figure_registry, row.names = FALSE)
})

test_that("optimizer producer fingerprint invalidates checkpoints after code mutation", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)
  producer_copy <- tempfile("optimizer-producer-", fileext = ".R")
  file.copy(
    file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"),
    producer_copy
  )
  package_copy <- tempfile("optimizer-package-")
  dir.create(file.path(package_copy, "R"), recursive = TRUE)
  file.copy(file.path(root, c("DESCRIPTION", "NAMESPACE")), package_copy)
  file.copy(
    file.path(root, "R", "benchmark-adoption-scenarios.R"),
    file.path(package_copy, "R", "benchmark-adoption-scenarios.R")
  )
  settings <- list(
    profile = "smoke", seed = 101L, out_dir = tempfile("fingerprint-checkpoint-"),
    root = root, producer_source_path = producer_copy, package_source_root = package_copy
  )
  case <- env$jss_joint_case_definitions()[1, , drop = FALSE]
  cfg <- list(reps = 1L, resume = TRUE, workers = 1L)
  before <- env$jss_joint_checkpoint_spec(settings, case, 1L, cfg)
  checkpoint <- env$jss_joint_checkpoint_path(settings, case, 1L, cfg)
  result <- env$jss_joint_contract_fixture(case, 1L, settings)
  env$jss_joint_write_checkpoint(result, checkpoint, settings, case, 1L, cfg)
  expect_true(env$jss_joint_checkpoint_valid(checkpoint, settings, case, 1L, cfg))
  write(
    "# package-code mutation",
    file = file.path(package_copy, "R", "benchmark-adoption-scenarios.R"),
    append = TRUE
  )
  after_package <- env$jss_joint_checkpoint_spec(settings, case, 1L, cfg)
  expect_identical(before$schema_version, 7L)
  expect_false(identical(
    before$package_identity$source_sha256,
    after_package$package_identity$source_sha256
  ))
  expect_false(env$jss_joint_checkpoint_valid(checkpoint, settings, case, 1L, cfg))
  env$jss_joint_write_checkpoint(result, checkpoint, settings, case, 1L, cfg)
  expect_true(env$jss_joint_checkpoint_valid(checkpoint, settings, case, 1L, cfg))
  attested <- settings
  attested$checkpoint_package_identity <- env$jss_joint_checkout_package_identity(settings)
  attested$checkpoint_producer_fingerprint <- env$jss_joint_producer_fingerprint(settings)
  write("# producer mutation", file = producer_copy, append = TRUE)
  after <- env$jss_joint_checkpoint_spec(settings, case, 1L, cfg)
  expect_false(identical(after_package$producer_fingerprint, after$producer_fingerprint))
  expect_false(identical(after_package, after))
  expect_false(env$jss_joint_checkpoint_valid(checkpoint, settings, case, 1L, cfg))
  expect_error(env$jss_joint_reverify_checkpoint_or_quarantine(attested, checkpoint),
    "changed after attestation")
  expect_false(file.exists(checkpoint))
  expect_true(length(Sys.glob(paste0(checkpoint, ".source-identity-changed-*"))) == 1L)
})

test_that("optimizer worker refuses a parent-to-worker producer mutation", {
  skip_if_not_installed("pkgload")
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)
  identity <- env$jss_joint_checkout_package_identity(list(root = root))
  expect_error(env$jss_joint_initialize_worker(
    root, .libPaths(), Sys.getenv("R_LIBS_USER", unset = ""), identity,
    expected_producer_fingerprint = strrep("0", 64L)
  ), "changed between parent fingerprinting")
})

test_that("optimizer PSOCK execution preserves order and resumes without duplicates", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)
  cases <- env$jss_joint_case_definitions()[1:2, , drop = FALSE]
  serial_out <- tempfile("jvs-serial-")
  parallel_out <- tempfile("jvs-parallel-")
  marker_dir <- tempfile("jvs-markers-")
  dir.create(serial_out)
  dir.create(parallel_out)
  dir.create(marker_dir)
  stale_library <- tempfile("stale-library-")
  stale_package <- file.path(stale_library, "gamlss.longitudinal")
  dir.create(file.path(stale_package, "R"), recursive = TRUE)
  writeLines(
    c(
      "Package: gamlss.longitudinal", "Version: 99.0.0",
      "Title: Deliberately stale test package", "Description: Adversarial stale package.",
      "License: GPL-3"
    ),
    file.path(stale_package, "DESCRIPTION")
  )
  writeLines("exportPattern(\"^[[:alpha:]]+\")", file.path(stale_package, "NAMESPACE"))
  withr::local_libpaths(c(stale_library, .libPaths()))
  expect_identical(normalizePath(.libPaths()[[1L]], winslash = "/"), normalizePath(stale_library, winslash = "/"))
  expect_identical(
    normalizePath(find.package("gamlss.longitudinal", lib.loc = .libPaths(), quiet = TRUE), winslash = "/"),
    normalizePath(stale_package, winslash = "/")
  )
  contract_fixture <- env$jss_joint_contract_fixture
  mock_run <- function(case, rep_idx, settings, cfg) {
    marker <- file.path(
      settings$marker_dir,
      sprintf("%s-rep-%d.txt", case$case_id[[1L]], rep_idx)
    )
    write(as.character(Sys.getpid()), file = marker, append = TRUE)
    contract_fixture(case, rep_idx, settings,
      deterministic_offset = case$n[[1L]] + rep_idx)
  }
  serial_settings <- list(
    profile = "smoke", seed = 202L, out_dir = serial_out, root = root,
    logs_dir = file.path(serial_out, "logs"), marker_dir = tempfile("unused-serial-markers-")
  )
  dir.create(serial_settings$marker_dir)
  parallel_settings <- list(
    profile = "smoke", seed = 202L, out_dir = parallel_out, root = root,
    logs_dir = file.path(parallel_out, "logs"), marker_dir = marker_dir
  )
  serial <- env$jss_joint_run_simulation(
    serial_settings, cases, list(reps = 2L, resume = TRUE, workers = 1L), mock_run
  )
  parallel_result <- env$jss_joint_run_simulation(
    parallel_settings, cases, list(reps = 2L, resume = TRUE, workers = 2L), mock_run
  )
  expect_equal(parallel_result, serial, ignore_attr = TRUE)
  parallel_status <- attr(parallel_result, "checkpoint_status")
  expect_identical(parallel_status$task_id, seq_len(4L))
  expect_true(all(parallel_status$checkpoint_schema_version == 7L))
  expect_true(all(nzchar(parallel_status$producer_fingerprint)))
  expect_true(all(parallel_status$package_source_sha256 == parallel_status$verified_source_sha256))
  expect_true(all(parallel_status$package_identity_verified))
  expect_true(all(parallel_status$producer_fingerprint_algorithm == "SHA-256"))
  expect_true(all(normalizePath(parallel_status$verified_package_path, winslash = "/") == normalizePath(root, winslash = "/")))
  expect_true(all(parallel_status$verified_package_version != "99.0.0"))
  expect_true(all(parallel_status$execution_mode == "psock"))
  expect_true(all(parallel_status$workers_used == 2L))
  markers <- list.files(marker_dir, pattern = "[.]txt$", full.names = TRUE)
  expect_length(markers, 4L)
  expect_true(all(vapply(markers, function(path) length(readLines(path, warn = FALSE)) == 1L, logical(1))))

  resumed <- env$jss_joint_run_simulation(
    parallel_settings, cases, list(reps = 2L, resume = TRUE, workers = 2L), mock_run
  )
  expect_equal(resumed, parallel_result, ignore_attr = TRUE)
  resumed_status <- attr(resumed, "checkpoint_status")
  expect_true(all(resumed_status$resumed))
  expect_true(all(resumed_status$execution_mode == "resume_only"))
  expect_true(all(resumed_status$workers_used == 0L))
  expect_true(all(vapply(markers, function(path) length(readLines(path, warn = FALSE)) == 1L, logical(1))))
})

test_that("optimizer workers quarantine forged checkpoints before returning", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)
  out <- tempfile("jvs-forged-worker-")
  dir.create(out)
  settings <- list(
    profile = "smoke", seed = 303L, out_dir = out, root = root,
    logs_dir = file.path(out, "logs")
  )
  cases <- env$jss_joint_case_definitions()[1, , drop = FALSE]
  contract_fixture <- env$jss_joint_contract_fixture
  forged_run <- function(case, rep_idx, settings, cfg) {
    out <- contract_fixture(case, rep_idx, settings)
    out$converged <- c(TRUE, FALSE)
    out$retained <- TRUE
    out
  }
  expect_error(
    env$jss_joint_run_simulation(
      settings, cases, list(reps = 2L, resume = TRUE, workers = 2L), forged_run
    ),
    "invalid checkpoint|invalid payload"
  )
  checkpoint_dir <- env$jss_joint_checkpoint_dir(settings)
  expect_length(list.files(checkpoint_dir, pattern = "[.]rds$", recursive = TRUE), 0L)
  expect_equal(length(list.files(checkpoint_dir, pattern = "[.]rds$", recursive = TRUE)), 0L)
})

test_that("actual optimizer runner has serial and PSOCK warning/content parity", {
  skip_if_not_installed("pkgload")
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)
  cases <- env$jss_joint_case_definitions()[1:2, , drop = FALSE]
  cases$n <- 10L
  cases$time_points <- 3L
  cases$total_observations <- 30L
  cfg <- env$jss_joint_simulation_settings(list(profile = "smoke"))
  cfg$reps <- 1L
  cfg$resume <- FALSE
  cfg$max_outer_iter <- 1L
  cfg$max_inner_iter <- 1L
  cfg$max_elapsed_sec <- 3
  cfg$variogram_nsim <- 1L
  serial_settings <- list(profile = "smoke", seed = 717L, out_dir = tempfile("jvs-real-serial-"),
    logs_dir = tempfile("jvs-real-serial-log-"), root = root)
  psock_settings <- list(profile = "smoke", seed = 717L, out_dir = tempfile("jvs-real-psock-"),
    logs_dir = tempfile("jvs-real-psock-log-"), root = root)
  dir.create(serial_settings$out_dir)
  dir.create(psock_settings$out_dir)
  cfg$workers <- 1L
  serial <- env$jss_joint_run_simulation(serial_settings, cases, cfg, env$jss_joint_run_case_rep)
  cfg$workers <- 2L
  psock <- env$jss_joint_run_simulation(psock_settings, cases, cfg, env$jss_joint_run_case_rep)
  normalize <- function(x) {
    x$elapsed_sec <- 0
    x$result_portable_sha256 <- NULL
    attr(x, "checkpoint_status") <- NULL
    rownames(x) <- NULL
    x
  }
  expect_identical(as.character(psock$warnings), as.character(serial$warnings))
  expect_equal(normalize(psock), normalize(serial), tolerance = 1e-10)
  expect_false(any(grepl("Flat optimizer arguments are deprecated", serial$warnings, fixed = TRUE)))
})

test_that("monotone dropout has no observations after subject dropout", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "missingness-study-helpers.R"), local = env)

  set.seed(42)
  subjects <- data.frame(
    id = 1:300,
    x1 = stats::rnorm(300),
    x2 = stats::rbinom(300, 1, 0.5),
    s1 = stats::runif(300)
  )
  dat <- merge(subjects, data.frame(time = 1:4), by = NULL)
  dat <- dat[order(dat$id, dat$time), , drop = FALSE]
  dat$t <- (dat$time - 1) / 3
  dat$response <- stats::rnorm(nrow(dat))

  dropout <- env$jss_missing_apply(dat, 0.30, "monotone_dropout", seed = 123)
  summary <- attr(dropout, "missingness_summary")
  expect_true(summary$no_observations_after_dropout)
  expect_identical(summary$n_subjects_with_interior_gaps, 0L)
  expect_identical(summary$analysis_role, "headline")
  expect_lt(abs(summary$observed_missing_rate - 0.30), 0.08)

  intermittent <- env$jss_missing_apply(
    dat,
    0.30,
    "time_dependent_intermittent_mar",
    seed = 123
  )
  intermittent_summary <- attr(intermittent, "missingness_summary")
  expect_identical(intermittent_summary$missingness_label, "time-dependent intermittent MAR")
  expect_gt(intermittent_summary$n_subjects_with_interior_gaps, 0L)
  expect_error(env$jss_missing_validate_mechanisms("time_mar"), "renamed")
})

test_that("missingness optional diagnostics normalize non-finite scalars to typed missing values", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "missingness-study-helpers.R"), local = env)

  expect_identical(env$jss_missing_finite_scalar_or_na(1.25), 1.25)
  expect_identical(env$jss_missing_finite_scalar_or_na(-Inf), NA_real_)
  expect_identical(env$jss_missing_finite_scalar_or_na(Inf), NA_real_)
  expect_identical(env$jss_missing_finite_scalar_or_na(NaN), NA_real_)
  expect_identical(env$jss_missing_finite_scalar_or_na(NA_real_), NA_real_)
  expect_identical(env$jss_missing_finite_scalar_or_na(3, integer = TRUE), 3L)
  expect_identical(env$jss_missing_finite_scalar_or_na(Inf, integer = TRUE), NA_integer_)
})

test_that("missingness producer encodes an unavailable fit as an explicit fit error", {
  root <- local_phase2_repo_root()
  producer <- parse(file.path(
    root, "paper", "scripts", "final-simulations", "missingness",
    "run_missingness_study.R"
  ))
  is_extractor <- vapply(producer, function(expr) {
    is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      identical(expr[[2L]], as.name("extract_convergence_info"))
  }, logical(1))
  expect_equal(sum(is_extractor), 1L)
  env <- new.env(parent = globalenv())
  eval(producer[[which(is_extractor)]], envir = env)

  unavailable <- env$extract_convergence_info(NULL)
  expect_identical(unavailable$converged, FALSE)
  expect_identical(unavailable$stop_reason, "fit_error")
  expect_identical(unavailable$outer_iterations, NA_integer_)
})

test_that("missingness registered Cartesian grid rejects deletion substitution and seed mutation", {
  root <- local_phase2_repo_root(); env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "missingness-study-helpers.R"), local = env)
  source(file.path(root, "paper", "R", "04-missingness-dropout-sensitivity.R"), local = env)
  grid <- env$jss_missing_registered_task_grid()
  runs <- grid[rep(seq_len(nrow(grid)), each = 2L), , drop = FALSE]
  runs$model <- rep(c("gamlss.longitudinal", "gamlss2"), nrow(grid))
  checkpoints <- grid[c("scenario", "rep")]
  expect_silent(env$jss_missing_validate_registered_design(runs, checkpoints))
  csv_roundtrip <- runs
  csv_roundtrip$target_missing_rate[
    abs(csv_roundtrip$target_missing_rate - 0.3) < 1e-12
  ] <- as.numeric("0.3")
  expect_silent(env$jss_missing_validate_registered_design(csv_roundtrip, checkpoints))
  expect_error(env$jss_missing_validate_registered_design(runs[-1L, ], checkpoints),
    "240-task/480-model")
  substituted <- runs
  substituted$target_missing_rate[1:2] <- .25
  expect_error(env$jss_missing_validate_registered_design(substituted, checkpoints),
    "Cartesian task grid")
  meaningful_drift <- csv_roundtrip
  affected <- which(abs(meaningful_drift$target_missing_rate - 0.3) < 1e-12)[1:2]
  meaningful_drift$target_missing_rate[affected] <-
    meaningful_drift$target_missing_rate[affected] + 1e-8
  expect_error(env$jss_missing_validate_registered_design(meaningful_drift, checkpoints),
    "Cartesian task grid")
  mutated_seed <- runs
  mutated_seed$missingness_seed[[1L]] <- mutated_seed$missingness_seed[[1L]] + 1L
  expect_error(env$jss_missing_validate_registered_design(mutated_seed, checkpoints),
    "Cartesian task grid")
})

test_that("missingness checkpoints reject stale task metadata", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "missingness-study-helpers.R"), local = env)
  task <- data.frame(
    scenario_id = 1L,
    n = 500L,
    d = 4L,
    rep = 2L,
    missing_mechanism = "monotone_dropout",
    missing_rate = 0.3
  )
  configuration <- list(method = "RS", max_outer = 100L)
  spec <- env$jss_missing_checkpoint_spec(task, configuration)
  intermittent_task <- task
  intermittent_task$missing_mechanism <- "time_dependent_intermittent_mar"
  intermittent_spec <- env$jss_missing_checkpoint_spec(intermittent_task, configuration)
  expect_false(identical(spec$missingness_seed, intermittent_spec$missingness_seed))
  result <- list(
    checkpoint_schema_version = 5L,
    checkpoint_configuration_key = "test-v5",
    checkpoint_spec = spec,
    checkpoint_provenance = env$jss_missing_runtime_identity(),
    runs = data.frame(
      scenario = spec$scenario,
      n = 500L,
      d = 4L,
      rep = 2L,
      missing_mechanism = "monotone_dropout",
      target_missing_rate = 0.3,
      simulation_seed = spec$simulation_seed,
      missingness_seed = spec$missingness_seed,
      model = c("gamlss2", "gamlss.longitudinal"),
      success = TRUE,
      converged = TRUE,
      retained = TRUE,
      stop_reason = "converged",
      failure_type = "none",
      logLik = -10,
      df = 5,
      elapsed_sec = 1,
      error = NA_character_
    ),
    fixed = data.frame(
      scenario = spec$scenario, model = "gamlss.longitudinal", rep = 2L,
      n = 500L, d = 4L, missing_mechanism = "monotone_dropout",
      target_missing_rate = 0.3, parameter = "mu", estimate = 999, true_value = 1
    ),
    smooth = NULL,
    joint = NULL,
    predictive = NULL,
    missingness = data.frame(
      scenario = spec$scenario, n = 500L, d = 4L, rep = 2L,
      missing_mechanism = "monotone_dropout", target_missing_rate = 0.3,
      missingness_label = "subject-level monotone dropout MAR",
      missingness_pattern = "monotone_dropout", analysis_role = "headline",
      observed_missing_rate = 0.5, n_rows = 2000L, n_observed_rows = 1000L,
      n_subjects = 500L, n_complete_subjects = 0L, n_dropout_subjects = 500L,
      n_monotone_incomplete_subjects = 500L, n_subjects_with_interior_gaps = 0L,
      no_observations_after_dropout = TRUE, total_adjacent_pairs = 1500L,
      complete_adjacent_pairs = 500L, complete_adjacent_pair_rate = 1 / 3
    ),
    missingness_pattern = data.frame(
      scenario = spec$scenario, n = 500L, d = 4L, rep = 2L,
      missing_mechanism = "monotone_dropout", id = rep(seq_len(500L), each = 4L),
      time = rep(as.numeric(1:4), 500L),
      response_observed = rep(c(TRUE, TRUE, FALSE, FALSE), 500L)
    ),
    warning_events = env$jss_missing_empty_warning_events()
  )
  result$checkpoint_content_sha256 <- env$jss_missing_content_sha256(
    env$jss_missing_checkpoint_content(result)
  )
  expect_true(env$jss_missing_checkpoint_valid(result, task, configuration))
  failed_fit <- result
  failed_row <- failed_fit$runs$model == "gamlss.longitudinal"
  failed_fit$runs$success[failed_row] <- FALSE
  failed_fit$runs$converged[failed_row] <- FALSE
  failed_fit$runs$retained[failed_row] <- FALSE
  failed_fit$runs$stop_reason[failed_row] <- "time_limit"
  failed_fit$runs$failure_type[failed_row] <- "fit_error:time_limit"
  failed_fit$runs$logLik[failed_row] <- NA_real_
  failed_fit$runs$df[failed_row] <- NA_real_
  failed_fit$runs$error[failed_row] <- paste(
    "Model exceeded max_elapsed_sec during RS outer iteration",
    "(elapsed 183.5 sec > 180.0 sec)."
  )
  failed_fit["fixed"] <- list(NULL)
  failed_fit$checkpoint_content_sha256 <- env$jss_missing_content_sha256(
    env$jss_missing_checkpoint_content(failed_fit)
  )
  expect_true(env$jss_missing_checkpoint_valid(failed_fit, task, configuration))
  expect_match(failed_fit$runs$error[failed_row], "max_elapsed_sec", fixed = TRUE)
  expect_identical(env$jss_missing_failure_stop_reason(failed_fit$runs$error[failed_row]),
    "time_limit")
  rs_sentinel <- result
  rs_sentinel$runs$best_raw_loglik <- c(-Inf, NA_real_)
  rs_sentinel$checkpoint_content_sha256 <- env$jss_missing_content_sha256(
    env$jss_missing_checkpoint_content(rs_sentinel)
  )
  expect_false(env$jss_missing_checkpoint_valid(rs_sentinel, task, configuration))
  rs_sentinel$runs$best_raw_loglik <- vapply(
    rs_sentinel$runs$best_raw_loglik,
    env$jss_missing_finite_scalar_or_na,
    numeric(1L)
  )
  rs_sentinel$checkpoint_content_sha256 <- env$jss_missing_content_sha256(
    env$jss_missing_checkpoint_content(rs_sentinel)
  )
  expect_true(env$jss_missing_checkpoint_valid(rs_sentinel, task, configuration))
  durable_record <- env$jss_missing_checkpoint_archive_record(result)
  payload_mutation <- result
  payload_mutation$fixed$estimate <- payload_mutation$fixed$estimate + 0.25
  mutated_record <- env$jss_missing_checkpoint_archive_record(payload_mutation)
  expect_false(identical(durable_record$checkpoint_content_sha256,
    mutated_record$checkpoint_content_sha256))
  expect_false(identical(durable_record$public_payload_sha256,
    mutated_record$public_payload_sha256))
  warning_bound <- result
  captured_warning <- list(value = NULL, messages = "NaNs produced",
    condition_classes = "simpleWarning/warning/condition")
  warning_bound$warning_events <- env$jss_missing_warning_events(captured_warning,
    spec$scenario, spec$n, spec$d, spec$replicate, spec$missing_mechanism,
    spec$missing_rate, "gamlss2")
  warning_bound$checkpoint_content_sha256 <- env$jss_missing_content_sha256(
    env$jss_missing_checkpoint_content(warning_bound))
  expect_true(env$jss_missing_checkpoint_valid(warning_bound, task, configuration))
  warning_tamper <- warning_bound
  warning_tamper$warning_events$warning_message <- "Algorithm RS has not yet converged"
  expect_false(env$jss_missing_checkpoint_valid(warning_tamper, task, configuration))
  warning_tamper$checkpoint_content_sha256 <- env$jss_missing_content_sha256(
    env$jss_missing_checkpoint_content(warning_tamper))
  expect_false(env$jss_missing_checkpoint_valid(warning_tamper, task, configuration))
  forged_metric <- result
  forged_metric$fixed$estimate <- 1e300
  forged_metric$checkpoint_content_sha256 <- env$jss_missing_content_sha256(
    env$jss_missing_checkpoint_content(forged_metric)
  )
  expect_false(env$jss_missing_checkpoint_valid(forged_metric, task, configuration))
  character_key <- result
  character_key$runs$rep <- as.character(character_key$runs$rep)
  character_key$checkpoint_content_sha256 <- env$jss_missing_content_sha256(
    env$jss_missing_checkpoint_content(character_key))
  expect_false(env$jss_missing_checkpoint_valid(character_key, task, configuration))
  character_metric <- result
  character_metric$fixed$estimate <- as.character(character_metric$fixed$estimate)
  character_metric$checkpoint_content_sha256 <- env$jss_missing_content_sha256(
    env$jss_missing_checkpoint_content(character_metric))
  expect_false(env$jss_missing_checkpoint_valid(character_metric, task, configuration))
  contradictory <- result
  contradictory$runs$success[[1L]] <- FALSE
  contradictory$runs$converged[[1L]] <- TRUE
  contradictory$checkpoint_content_sha256 <- env$jss_missing_content_sha256(
    env$jss_missing_checkpoint_content(contradictory))
  expect_false(env$jss_missing_checkpoint_valid(contradictory, task, configuration))
  forged_provenance <- result
  forged_provenance$checkpoint_provenance$host <- paste0(
    forged_provenance$checkpoint_provenance$host, "-forged")
  forged_provenance$checkpoint_content_sha256 <- env$jss_missing_content_sha256(
    env$jss_missing_checkpoint_content(forged_provenance))
  expect_false(env$jss_missing_checkpoint_valid(forged_provenance, task, configuration))
  failed_with_metrics <- result
  failed_with_metrics$runs$success[[1L]] <- FALSE
  failed_with_metrics$runs$converged[[1L]] <- FALSE
  failed_with_metrics$runs$retained[[1L]] <- FALSE
  failed_with_metrics$runs$stop_reason[[1L]] <- "fit_error"
  failed_with_metrics$runs$failure_type[[1L]] <- "fit_error"
  failed_with_metrics$runs$error[[1L]] <- "synthetic fit error"
  failed_with_metrics$checkpoint_content_sha256 <- env$jss_missing_content_sha256(
    env$jss_missing_checkpoint_content(failed_with_metrics))
  expect_false(env$jss_missing_checkpoint_valid(failed_with_metrics, task, configuration))
  forged_summary <- result
  forged_summary$missingness$complete_adjacent_pairs <-
    forged_summary$missingness$complete_adjacent_pairs + 1L
  forged_summary$checkpoint_content_sha256 <- env$jss_missing_content_sha256(
    env$jss_missing_checkpoint_content(forged_summary))
  expect_false(env$jss_missing_checkpoint_valid(forged_summary, task, configuration))
  forged_pattern <- result
  forged_pattern$missingness_pattern$response_observed[1:4] <- c(TRUE, FALSE, TRUE, FALSE)
  forged_pattern$checkpoint_content_sha256 <- env$jss_missing_content_sha256(
    env$jss_missing_checkpoint_content(forged_pattern))
  expect_false(env$jss_missing_checkpoint_valid(forged_pattern, task, configuration))
  forged <- result
  forged$runs$converged[forged$runs$model == "gamlss.longitudinal"] <- FALSE
  forged$runs$retained[forged$runs$model == "gamlss.longitudinal"] <- TRUE
  expect_false(env$jss_missing_checkpoint_valid(forged, task, configuration))
  expect_equal(nrow(env$jss_missing_filter_payload(forged$fixed, forged$runs)), 0L)
  stale_payload <- result
  stale_payload$fixed$rep <- 999L
  expect_false(env$jss_missing_checkpoint_valid(stale_payload, task, configuration))
  result$runs$missing_mechanism <- "time_dependent_intermittent_mar"
  expect_false(env$jss_missing_checkpoint_valid(result, task, configuration))
  result$runs$missing_mechanism <- "monotone_dropout"
  result$runs$n <- 499L
  expect_false(env$jss_missing_checkpoint_valid(result, task, configuration))
  result$runs$n <- 500L
  expect_false(env$jss_missing_checkpoint_valid(result, task, list(method = "RS", max_outer = 101L)))
})

test_that("missingness fit warnings are classified, captured, and fail closed", {
  root <- local_phase2_repo_root(); env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "missingness-study-helpers.R"), local = env)
  captured <- env$jss_missing_capture_warnings({
    warning(paste(
      "Intermittent observation gaps were detected. Using the segmented likelihood",
      "requested by missingness = \"segment\": observations in different contiguous",
      "segments are treated as independent."
    ))
    warning("Algorithm RS has not yet converged")
    warning("NaNs produced")
    42L
  })
  expect_identical(captured$value, 42L)
  events <- env$jss_missing_warning_events(captured, "scenario", 500L, 4L, 1L,
    "monotone_dropout", 0.3, "gamlss.longitudinal")
  expect_identical(events$classification,
    c("segmented_likelihood_assumption", "optimizer_nonconvergence", "numerical_domain"))
  expect_true(all(events$expected))
  expect_silent(env$jss_missing_assert_warning_policy(events))

  unknown <- env$jss_missing_capture_warnings({ warning("unregistered warning text"); NULL })
  unknown_events <- env$jss_missing_warning_events(unknown, "scenario", 500L, 4L, 1L,
    "monotone_dropout", 0.3, "gamlss2")
  expect_identical(unknown_events$classification, "unexpected")
  expect_error(env$jss_missing_assert_warning_policy(unknown_events),
    "1 unexpected warning event")
})

test_that("missingness warning aggregation de-duplicates resume rows and rejects conflicts", {
  root <- local_phase2_repo_root(); env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "missingness-study-helpers.R"), local = env)
  captured <- list(value = NULL, messages = "NaNs produced", condition_classes = "simpleWarning/warning/condition")
  events <- env$jss_missing_warning_events(captured, "scenario", 500L, 4L, 1L,
    "monotone_dropout", 0.3, "gamlss2")
  result <- list(warning_events = events)
  expect_identical(nrow(env$jss_missing_collect_warning_events(list(result, result))), 1L)
  conflict <- result
  conflict$warning_events$warning_message <- "Algorithm RS has not yet converged"
  expect_error(env$jss_missing_collect_warning_events(list(result, conflict)),
    "conflicting duplicate event keys")
})

test_that("missingness warning evidence is content-bound and semantically validated", {
  root <- local_phase2_repo_root(); env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "missingness-study-helpers.R"), local = env)
  captured <- list(value = NULL, messages = "NaNs produced", condition_classes = "simpleWarning/warning/condition")
  event <- env$jss_missing_warning_events(captured, "scenario", 500L, 4L, 1L,
    "monotone_dropout", 0.3, "gamlss2")
  base_hash <- env$jss_missing_portable_task_sha256(data.frame(a = 1), NULL, NULL,
    data.frame(a = 1), warning_events = event)
  changed <- event; changed$warning_message <- "Algorithm RS has not yet converged"
  changed_hash <- env$jss_missing_portable_task_sha256(data.frame(a = 1), NULL, NULL,
    data.frame(a = 1), warning_events = changed)
  expect_false(identical(base_hash, changed_hash))

  classified <- env$jss_missing_classify_warning(changed$warning_message)
  expect_false(identical(changed$classification, classified$classification))
})

test_that("missingness coverage plot filters and audits every omitted row", {
  root <- local_phase2_repo_root(); env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "missingness-study-helpers.R"), local = env)
  fixed <- data.frame(
    target_missing_pct = c(0, 10, 20, 30),
    mean_coverage_95 = c(.95, NA_real_, .92, Inf),
    model = c("gamlss2", "gamlss.longitudinal", "gamlss2", "gamlss.longitudinal"),
    parameter = c("mu", "mu", "theta", "sigma"),
    missing_mechanism = "monotone_dropout", stringsAsFactors = FALSE)
  payload <- env$jss_missing_coverage_plot_payload(fixed, c("mu", "sigma"))
  expect_identical(payload$candidate_rows, 3L)
  expect_identical(payload$omitted_rows, 2L)
  expect_true(all(is.finite(payload$data$mean_coverage_95)))
  audit <- env$jss_missing_warning_audit(env$jss_missing_empty_warning_events(),
    data.frame(model = c("gamlss2", "gamlss.longitudinal")), payload)
  plot_audit <- audit[audit$audit_type == "plot_omission", , drop = FALSE]
  expect_identical(plot_audit$omitted_rows, 2L)
  expect_match(plot_audit$uncertainty, "candidates=3; plotted=1")
})

test_that("missingness PSOCK workers close cleanly after successful scheduling", {
  log <- tempfile("missingness-psock-", fileext = ".log")
  cl <- parallel::makePSOCKcluster(1L, outfile = log)
  closed <- FALSE
  on.exit(if (!closed) parallel::stopCluster(cl), add = TRUE)
  result <- parallel::parLapplyLB(cl, 1:2, function(x) x + 1L)
  expect_identical(unlist(result), 2:3)
  expect_silent(parallel::stopCluster(cl))
  closed <- TRUE
  worker_log <- if (file.exists(log)) readLines(log, warn = FALSE) else character()
  expect_false(any(grepl("unserialize|error reading from connection", worker_log,
    ignore.case = TRUE)))
})

test_that("missingness estimands use a fixed population centering target", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "missingness-study-helpers.R"), local = env)
  registry <- env$jss_missing_estimand_registry()
  expect_false(any(registry$observed_data_projection))
  truth <- function(s1) 0.3 * cos(2 * pi * s1) - 0.1 * (s1 - 0.5)^2
  population_mean <- env$jss_missing_population_smooth_mean(truth)
  grid <- seq(0, 1, length.out = 10001L)
  expect_equal(mean(truth(grid) - population_mean), 0, tolerance = 1e-12)
  set.seed(99)
  observed_projection <- mean(truth(stats::rbeta(500, 2, 8)))
  expect_gt(abs(observed_projection - population_mean), 0.01)
  producer <- paste(readLines(file.path(root, "paper", "scripts", "final-simulations",
    "missingness", "run_missingness_study.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("fitted_population_smooth_mean_longitudinal", producer, fixed = TRUE))
  expect_true(grepl("fitted_population_level_gamlss2", producer, fixed = TRUE))
  expect_true(grepl("smooth_truth[[p]](s1_grid) - jss_missing_population_smooth_mean", producer, fixed = TRUE))
  expect_true(grepl("not_available_without_joint_fixed_smooth_covariance", producer, fixed = TRUE))
  expect_true(grepl("if (population_intercept) NA_real_", producer, fixed = TRUE))
})

test_that("missingness longitudinal intercept metadata identifies only smoothed parameters", {
  root <- local_phase2_repo_root()
  producer <- parse(file.path(
    root, "paper", "scripts", "final-simulations", "missingness",
    "run_missingness_study.R"
  ))
  is_extractor <- vapply(producer, function(expr) {
    is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      identical(expr[[2L]], as.name("extract_fixed_estimates_longitudinal"))
  }, logical(1))
  expect_equal(sum(is_extractor), 1L)

  env <- new.env(parent = globalenv())
  env$compute_se <- FALSE
  env$params_all <- c("mu", "sigma", "nu", "tau", "theta", "zeta")
  env$fixed_terms <- "intercept"
  env$smooth_truth <- list(
    mu = identity,
    sigma = identity,
    theta = identity
  )
  env$true_beta <- setNames(
    lapply(env$params_all, function(parameter) c(intercept = 0)),
    env$params_all
  )
  env$calc_smooth_mean <- function(data_used, parameter) 0
  env$extract_one_longitudinal_term <- function(par_vec, parameter, term) 0
  eval(producer[[which(is_extractor)]], envir = env)

  extracted <- env$extract_fixed_estimates_longitudinal(
    list(par = numeric(), model_matrix = list(s = list()), par_s = list()),
    data.frame(s1 = 0.5)
  )
  smooth_intercepts <- extracted$parameter %in% c("mu", "sigma", "theta")
  expect_identical(extracted$intercept_includes_fitted_smooth_mean, smooth_intercepts)
  expect_identical(
    extracted$inference_status[smooth_intercepts],
    rep("not_available_without_joint_fixed_smooth_covariance", 3L)
  )
  expect_identical(
    extracted$inference_status[!smooth_intercepts],
    rep("coefficient_covariance_available_when_finite", 3L)
  )
})

test_that("missingness package source mutation changes the checkpoint contract", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "missingness-study-helpers.R"), local = env)
  copy <- tempfile("missingness-checkout-")
  dir.create(file.path(copy, "R"), recursive = TRUE)
  file.copy(list.files(file.path(root, "R"), full.names = TRUE), file.path(copy, "R"))
  file.copy(file.path(root, c("DESCRIPTION", "NAMESPACE")), copy)
  before <- env$jss_missing_checkout_identity(copy)
  target <- list.files(file.path(copy, "R"), full.names = TRUE)[[1L]]
  write("# adversarial mutation", file = target, append = TRUE)
  after <- env$jss_missing_checkout_identity(copy)
  expect_false(identical(before$source_sha256, after$source_sha256))
  producer_copy <- tempfile(fileext = ".R")
  file.copy(file.path(root, "paper", "R", "missingness-study-helpers.R"), producer_copy)
  producer_before <- env$jss_missing_producer_sha256(producer_copy)
  expect_silent(env$jss_missing_reverify_sources(copy, after, producer_copy, producer_before))
  write("# post-attestation producer mutation", producer_copy, append = TRUE)
  expect_error(env$jss_missing_reverify_sources(copy, after, producer_copy, producer_before),
    "changed after attestation")
})

test_that("public missingness loader rejects dummy summaries and forged identity", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "missingness-study-helpers.R"), local = env)
  source(file.path(root, "paper", "R", "04-missingness-dropout-sensitivity.R"), local = env)
  input <- tempfile("missingness-public-")
  dir.create(input)
  mechanisms <- c("monotone_dropout", "time_dependent_intermittent_mar")
  scenarios <- paste0("scenario_", seq_along(mechanisms))
  registry <- data.frame(missing_mechanism = mechanisms,
    display_name = c("subject-level monotone dropout MAR", "time-dependent intermittent MAR"),
    pattern = c("monotone_dropout", "intermittent"), analysis_role = c("headline", "sensitivity"))
  estimands <- data.frame(component = c("fixed_intercept", "smooth_curve"),
    estimand = c("population coefficient plus smooth mean", "population centered smooth"),
    observed_data_projection = FALSE, reference_distribution = c("Uniform", "fixed grid"))
  penalties <- data.frame(penalty = 1, role = "primary", scale = "link", application = "failed fit",
    justification = "registered sensitivity")
  runs <- do.call(rbind, lapply(seq_along(mechanisms), function(i) data.frame(
    scenario = scenarios[[i]], model = c("gamlss.longitudinal", "gamlss2"), rep = 1L,
    missing_mechanism = mechanisms[[i]], analysis_role = registry$analysis_role[[i]],
    missingness_pattern = registry$pattern[[i]], success = TRUE, converged = TRUE, retained = TRUE,
    failure_type = "none", no_observations_after_dropout = mechanisms[[i]] == "monotone_dropout"
  )))
  attempts <- data.frame(scenario = runs$scenario, model = runs$model,
    attempted = 1L, fit_successful = 1L, converged = 1L, retained = 1L, failed = 0L)
  failures <- data.frame(scenario = runs$scenario, model = runs$model,
    failure_reason = "retained", attempts = 1L)
  missing_by_rep <- data.frame(scenario = scenarios, rep = 1L, missing_mechanism = mechanisms,
    no_observations_after_dropout = c(TRUE, FALSE))
  missingness_pattern <- do.call(rbind, lapply(seq_along(scenarios), function(i) data.frame(
    scenario = scenarios[[i]], n = 10L, d = 2L, rep = 1L,
    missing_mechanism = mechanisms[[i]], id = 1L, time = c(1, 2),
    response_observed = if (i == 1L) c(TRUE, FALSE) else c(FALSE, TRUE)
  )))
  payload <- runs[c("scenario", "model", "rep")]
  payload$parameter <- "mu"
  fixed_raw <- transform(payload, term = "x1", estimate = 0.1, true_value = 0.1)
  smooth_raw <- transform(payload, s1 = 0.5, smooth_hat = 0, smooth_true = 0)
  headline <- data.frame(scenario = runs$scenario, model = runs$model, parameter = "mu",
    missing_mechanism = runs$missing_mechanism, missingness_pattern = runs$missingness_pattern,
    attempted = 1L, retained = 1L, failed = 0L, conditional_rmse_mcse = 0,
    conditional_rmse_conf_low = 0, conditional_rmse_conf_high = 0,
    failure_penalty = 1, failure_inclusive_rmse_mcse = 0,
    failure_inclusive_rmse_conf_low = 0, failure_inclusive_rmse_conf_high = 0,
    inference_status = ifelse(runs$model == "gamlss.longitudinal" & runs$missingness_pattern == "intermittent",
      "not_applicable_segmented_model_hessian", "conditional_on_retained_fits"))
  public_hash <- vapply(scenarios, function(scenario) env$jss_missing_portable_task_sha256(
    runs[runs$scenario == scenario, , drop = FALSE],
    fixed_raw[fixed_raw$scenario == scenario, , drop = FALSE],
    smooth_raw[smooth_raw$scenario == scenario, , drop = FALSE],
    missing_by_rep[missing_by_rep$scenario == scenario, , drop = FALSE],
    missingness_pattern[missingness_pattern$scenario == scenario, , drop = FALSE]
  ), character(1))
  hash_map <- setNames(public_hash, scenarios)
  runs$public_payload_sha256 <- unname(hash_map[runs$scenario])
  fixed_raw$public_payload_sha256 <- unname(hash_map[fixed_raw$scenario])
  smooth_raw$public_payload_sha256 <- unname(hash_map[smooth_raw$scenario])
  missing_by_rep$public_payload_sha256 <- unname(hash_map[missing_by_rep$scenario])
  missingness_pattern$public_payload_sha256 <- unname(hash_map[missingness_pattern$scenario])
  checkpoints <- data.frame(checkpoint_schema_version = 5L, scenario = scenarios, rep = 1L,
    checkpoint = paste0("rep_results/", scenarios, ".rds"),
    checkpoint_content_sha256 = strrep("a", 64), package_source_sha256 = strrep("b", 64),
    public_payload_sha256 = public_hash, producer_sha256 = strrep("c", 64), package_identity_verified = TRUE,
    timestamp_utc = "2026-01-01T00:00:00Z", worker_pid = 1234L,
    host = "host", os = "Windows", platform = R.version$platform, r_version = R.version.string,
    rng_kind = paste(RNGkind(), collapse = "/"), blas = "blas", lapack = "lapack",
    rlibs_user_sha256 = strrep("d", 64), libpaths_sha256 = strrep("e", 64))
  objects <- list(
    missingness_design_registry.csv = registry, missingness_estimand_registry.csv = estimands,
    missingness_sensitivity_registry.csv = penalties, missingness_checkpoint_status.csv = checkpoints,
    missingness_checkpoint_content_manifest.csv = data.frame(scenario = scenarios, rep = 1L,
      checkpoint = paste0("rep_results/", scenarios, ".rds"),
      checkpoint_content_sha256 = strrep("a", 64), public_payload_sha256 = public_hash),
    missingness_warning_events.csv = transform(env$jss_missing_empty_warning_events(),
      public_payload_sha256 = character()),
    missingness_warning_audit.csv = data.frame(
      audit_type = character(), classification = character(), policy_action = character(),
      expected = logical(), events = integer(), affected_fits = integer(),
      attempted_fits = integer(), omitted_rows = integer(),
      reconciliation_status = character(), uncertainty = character()),
    fit_run_log.csv = runs, attempt_failure_summary.csv = attempts,
    failure_reason_summary.csv = failures, missingness_by_rep.csv = missing_by_rep,
    missingness_pattern_by_subject_visit.csv = missingness_pattern,
    fixed_effects_by_rep.csv = fixed_raw, smooth_estimates_by_rep.csv = smooth_raw,
    missingness_headline_summary.csv = headline,
    fixed_term_summary_by_missingness.csv = data.frame(dummy = 1),
    smooth_irmse_summary.csv = data.frame(dummy = 1),
    smooth_selected_plot_data.csv = data.frame(dummy = 1)
  )
  for (nm in names(objects)) utils::write.csv(objects[[nm]], file.path(input, nm), row.names = FALSE)
  saveRDS(list(), file.path(input, "missingness_checkpoint_payloads.rds"))
  expect_error(env$jss_missingness_validate_public_bundle(input),
    "extra or missing columns|estimand registry|canonical reconstruction|package/source")
  attempts$attempted[[1L]] <- 2L
  utils::write.csv(attempts, file.path(input, "attempt_failure_summary.csv"), row.names = FALSE)
  expect_error(env$jss_missingness_validate_public_bundle(input))
})

test_that("public Monte Carlo interval guards reject impossible bounds and negative MCSE", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "04-missingness-dropout-sensitivity.R"), local = env)
  valid <- data.frame(effect = 1, effect_mcse = .1, effect_conf_low = .8, effect_conf_high = 1.2)
  expect_silent(env$jss_missing_validate_intervals(valid, "effect", "test"))
  negative <- valid; negative$effect_mcse <- -.1
  expect_error(env$jss_missing_validate_intervals(negative, "effect", "test"), "impossible")
  reversed <- valid; reversed$effect_conf_low <- 1.1
  expect_error(env$jss_missing_validate_intervals(reversed, "effect", "test"), "impossible")
})

test_that("nonconverged optimizer pairs are failures and have no conditional deltas", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)
  cases <- env$jss_joint_case_definitions()[1, , drop = FALSE]
  make_row <- function(method, converged, value) data.frame(
    case_id = cases$case_id, base_case_id = cases$base_case_id,
    contrast_factor = cases$contrast_factor, contrast_label = cases$contrast_label,
    contrast_level = cases$contrast_level, hypothesis_role = cases$hypothesis_role,
    joint_review_rep = 1L, family = cases$family, copula = cases$copula,
    design = "paired_one_factor", n = cases$n, time_points = cases$time_points,
    n_subject = cases$n, n_time = cases$time_points, total_observations = cases$total_observations,
    mu_strength = cases$mu_strength, sigma_strength = cases$sigma_strength,
    theta_strength = cases$theta_strength, time_shape = cases$time_shape,
    dependence = "case_specific", missingness = "none", start_mode = "default",
    paired_seed = 123L, method = method, success = TRUE, converged = converged,
    retained = converged, train_joint_loglik = value, elapsed_sec = 1,
    failure_type = if (converged) "none" else "nonconverged:outer_limit",
    stop_reason = if (converged) "converged" else "outer_limit",
    error = NA_character_, stringsAsFactors = FALSE
  )
  results <- rbind(make_row("rs_separate", TRUE, 1), make_row("rs_joint", FALSE, 2))
  results$retained[results$method == "rs_joint"] <- TRUE
  delta <- env$jss_joint_delta_table(results)
  expect_false(delta$rs_joint_retained)
  expect_true(is.na(delta$delta_train_joint_loglik))
  failures <- env$jss_joint_failure_summary(results)
  expect_true(any(grepl("optimizer_nonconvergence", failures$failure_reason)))
  expect_equal(nrow(env$jss_joint_metric_wins(results[results$method == "rs_joint", , drop = FALSE])), 0L)
  forged_delta <- delta
  forged_delta$rs_joint_retained <- TRUE
  forged_delta$delta_train_joint_loglik <- 999
  forged_evidence <- env$jss_joint_difference_uncertainty(forged_delta)
  expect_identical(
    forged_evidence$retained_pairs[forged_evidence$metric == "train_joint_loglik"],
    0L
  )
  expect_true(is.na(
    forged_evidence$conditional_difference[forged_evidence$metric == "train_joint_loglik"]
  ))
})

test_that("every registered optimizer difference has MC uncertainty and sign uncertainty", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)
  registry <- env$jss_joint_difference_metric_registry()
  deltas <- data.frame(
    case_id = rep("JVS01", 4), base_case_id = "JVS01", contrast_factor = "base",
    rs_joint_success = TRUE, rs_joint_converged = TRUE,
    rs_separate_success = TRUE, rs_separate_converged = TRUE,
    rs_joint_retained = TRUE, rs_separate_retained = TRUE,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(registry))) deltas[[paste0("delta_", registry$metric[[i]])]] <- c(-2, -1, 1, 2) * i
  evidence <- env$jss_joint_difference_uncertainty(deltas)
  expect_setequal(evidence$metric, registry$metric)
  expect_true(all(evidence$retained_pairs == 4L))
  expect_true(all(is.finite(evidence$difference_mcse)))
  expect_true(all(is.finite(evidence$sign_probability_mcse)))
  expect_true(all(is.finite(evidence$sign_probability_conf_low)))
  expect_true(all(is.finite(evidence$sign_probability_conf_high)))
  expect_true(all(evidence$metric_finite_pairs == 4L))
  expect_true(all(evidence$metric_failed_pairs == 0L))
  expect_true(all(is.finite(evidence$retention_difference_mcse)))
  expect_true(all(is.finite(evidence$convergence_difference_mcse)))
  expect_true(all(is.finite(evidence$failure_inclusive_mcse)))
  hypotheses <- data.frame(
    hypothesis = "registered contrast", focal_case = "JVS01", comparator_case = "JVS01", metric = registry$metric,
    decision = "pending production", rationale = "test", stringsAsFactors = FALSE
  )
  hypothesis_evidence <- env$jss_joint_hypothesis_evidence(hypotheses, evidence)
  expect_true(all(c(
    "focal_difference_mcse", "focal_difference_conf_low", "focal_difference_conf_high",
    "comparator_difference_mcse", "comparator_difference_conf_low", "comparator_difference_conf_high"
  ) %in% names(hypothesis_evidence)))
})

test_that("public missingness loader rejects the tracked legacy bundle", {
  root <- local_phase2_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "04-missingness-dropout-sensitivity.R"), local = env)
  settings <- list(
    public_data_dir = file.path(root, "paper", "data", "public-derived"),
    figures_dir = tempfile("figures-")
  )
  dir.create(settings$figures_dir)
  expect_error(env$jss_run_04_missingness_from_public_inputs(settings), "legacy")
  expect_true(file.exists(file.path(root, "paper", "data", "public-derived", "missingness", "README.md")))
})

test_that("missingness producer writes public plot data and registered sensitivity evidence", {
  root <- local_phase2_repo_root()
  code <- paste(readLines(
    file.path(root, "paper", "scripts", "final-simulations", "missingness", "run_missingness_study.R"),
    warn = FALSE
  ), collapse = "\n")
  expect_true(grepl('"smooth_selected_plot_data.csv"', code, fixed = TRUE))
  expect_true(grepl('"smooth_irmse_summary.csv"', code, fixed = TRUE))
  expect_true(grepl('"missingness_headline_summary.csv"', code, fixed = TRUE))
  expect_true(grepl('"missingness_sensitivity_registry.csv"', code, fixed = TRUE))
  expect_true(grepl('"missingness_pattern_by_subject_visit.csv"', code, fixed = TRUE))
  expect_true(grepl('"missingness_checkpoint_payloads.rds"', code, fixed = TRUE))
  expect_true(grepl('"missingness_checkpoint_content_manifest.csv"', code, fixed = TRUE))
  expect_true(grepl("not_applicable_segmented_model_hessian", code, fixed = TRUE))
  expect_true(grepl('Sys.getenv("RESUME_CHECKPOINTS", unset = "1")', code, fixed = TRUE))
  expect_true(grepl('design_version = "missingness-v5"', code, fixed = TRUE))
  expect_true(grepl("jss_missing_verify_checkout", code, fixed = TRUE))
  expect_true(grepl("jss_missing_portable_task_sha256", code, fixed = TRUE))
  expect_true(grepl("checkpoint_content_sha256", code, fixed = TRUE))
  expect_true(grepl("jss_missing_reverify_sources", code, fixed = TRUE))
  expect_true(grepl(
    "runs_all$retained <- runs_all$success %in% TRUE & runs_all$converged %in% TRUE",
    code,
    fixed = TRUE
  ))
  expect_true(grepl("retained_df <- df[df$retained %in% TRUE", code, fixed = TRUE))
  expect_true(all(vapply(
    c("fixed", "smooth", "joint", "predictive"),
    function(payload) grepl(
      paste0(payload, "_all <- jss_missing_filter_payload"), code, fixed = TRUE
    ),
    logical(1)
  )))
})

test_that("CG performance guidance is explicitly withdrawn", {
  root <- local_phase2_repo_root()
  checklist <- paste(readLines(
    file.path(root, "paper", "notes", "03-optimizer-guidance-checklist.md"),
    warn = FALSE
  ), collapse = "\n")
  expect_match(checklist, "empirical CG performance guidance is withdrawn", fixed = TRUE)
  expect_match(checklist, "Do not claim that CG and joint RS are interchangeable", fixed = TRUE)
  expect_match(checklist, "test-p0-optimizer-cg-helpers.R", fixed = TRUE)
  expect_match(checklist, "The RS method uses first derivatives", fixed = TRUE)
  expect_match(checklist, "CG is not included in this empirical comparison", fixed = TRUE)
  expect_match(checklist, "The CG implementation constrains each simultaneous Hessian-based update", fixed = TRUE)
})

test_that("protected manuscript contains no unsupported optimizer performance guidance", {
  root <- local_phase2_repo_root()
  manuscript <- paste(readLines(
    file.path(root, "paper", "manuscript", "main.tex"),
    warn = FALSE
  ), collapse = " ")
  hits <- local_optimizer_manuscript_hits(manuscript)
  expect_true(
    length(hits) == 0L,
    info = paste(
      "Protected manuscript still contains unsupported CG/RS performance guidance:",
      paste(hits, collapse = " | "),
      "See paper/notes/03-optimizer-guidance-checklist.md for exact removals."
    )
  )
})

test_that("protected manuscript contains no unfrozen missingness claims or figures", {
  root <- local_phase2_repo_root()
  manuscript <- paste(readLines(file.path(root, "paper", "manuscript", "main.tex"), warn = FALSE), collapse = " ")
  hits <- local_missingness_manuscript_hits(manuscript)
  expect_identical(hits, character(), info = paste(
    "Protected manuscript must stay neutral until the schema-v4 missingness bundle validates:",
    paste(hits, collapse = " | ")
  ))
  normalized <- gsub("[[:space:]]+", " ", manuscript)
  expect_match(normalized, "full-data population", ignore.case = TRUE)
  expect_match(normalized, "not observed-data projections", ignore.case = TRUE)
})

test_that("missingness manuscript guard catches stale wording and unfrozen figures", {
  stale <- paste(
    "Sensitivity to missingness at random and by time.",
    "The dropout analysis improves robustness.",
    "\\includegraphics{smooth_selected_recovery_curves.png}"
  )
  hits <- local_missingness_manuscript_hits(stale)
  expect_true(all(c("stale_appendix_wording", "unfrozen_missingness_figure",
    "directional_missingness_claim") %in% hits))
})

test_that("optimizer manuscript guard permits registered pending wording and catches stale claims", {
  neutral <- paste(
    "The registered comparison evaluates RS with separately optimized margin and copula blocks",
    "against RS optimization of the joint likelihood. CG is not included in this empirical",
    "comparison. Production is in progress; no comparative optimizer-performance recommendation",
    "is made until outputs are frozen. The registered design uses one reference scenario and",
    "three paired, one-factor contrasts, with attempts determined from the Monte Carlo precision target."
  )
  expect_identical(local_optimizer_manuscript_hits(neutral), character())

  synthetic_stale <- c(
    "We developed six example cases selected from a broad grid.",
    "We performed 100 replications for each example in the joint versus separate optimizer study.",
    "\\input{tables/gamma-joint-vs-separate-six-case-median-iqr-table}",
    "The jointly optimized model generally outperforms the separately optimized model.",
    "Joint optimization substantially increases runtime over separate optimization.",
    "Auto-regressive structures will generally fit reasonably.",
    "Covariate-dependent adjacent dependence can be captured very accurately."
  )
  stale_hits <- local_optimizer_manuscript_hits(paste(synthetic_stale, collapse = " "))
  expect_true(all(c(
    "stale_six_examples", "stale_broad_grid", "stale_fixed_100",
    "retired_gamma_table", "joint_directional_vs_separate",
    "unregistered_reasonable_fit", "unregistered_accurate_capture"
  ) %in% stale_hits))
})
