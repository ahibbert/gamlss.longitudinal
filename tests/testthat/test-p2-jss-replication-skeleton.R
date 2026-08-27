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
  expect_false(grepl("05-application-lipid", graph, fixed = TRUE))
  expect_false(grepl("06-application-rand", graph, fixed = TRUE))
  expect_false(grepl("GAMLSS_LONGITUDINAL_LIPID_DATA", graph, fixed = TRUE))
})

test_that("LIPID recipe works only against a generic contract fixture", {
  root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "05-application-lipid.R"), local = env)
  fixture <- data.frame(subject = rep(1:3, each = 2), time = rep(0:1, 3), response = c(1, 2, 2, 3, 3, 4), treatment = rep(c(0, 1, 0), each = 2))
  expect_true(env$jss_validate_lipid_input(fixture))
  recipe <- env$jss_lipid_analysis_recipe(fixture)
  expect_identical(recipe$fit, "gamlss_longitudinal")
  expect_error(env$jss_validate_lipid_input(fixture[-1]), "missing")
})

test_that("checkpoint resumption excludes completed fits", {
  root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(root, "paper", "R", "07-gamma-copula-misspecification.R"), local = env)
  settings <- list(profile = "full", seed = 20260528L)
  config <- env$jss_misspec_config(settings, stage = "smoke")
  grid <- env$jss_misspec_grid(config)
  checkpoint_dir <- tempfile("checkpoints-"); dir.create(checkpoint_dir)
  utils::write.csv(data.frame(fit_id = grid$fit_id[[1]], success = TRUE), env$jss_misspec_checkpoint_path(checkpoint_dir, grid$fit_id[[1]]), row.names = FALSE)
  expect_equal(nrow(env$jss_misspec_pending_grid(grid, checkpoint_dir)), nrow(grid) - 1L)
})

test_that("copula misspecification checkpoints are atomic and worker-aware", {
  root <- local_jss_repo_root()
  code <- paste(readLines(file.path(root, "paper", "R", "07-gamma-copula-misspecification.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("file.rename(temporary, path)", code, fixed = TRUE))
  expect_true(grepl("parallel::parLapplyLB", code, fixed = TRUE))
  expect_true(grepl("workers = settings$workers", code, fixed = TRUE))
})

test_that("promoted Monte Carlo runners validate and atomically replace checkpoints", {
  root <- local_jss_repo_root()
  scripts <- file.path(root, "paper", "scripts", "final-simulations", c(
    "missingness/run_missingness_study.R",
    "bcpe-t/simulation_bcpe_t_gamlss_comparison.R"
  ))
  for (script in scripts) {
    code <- paste(readLines(script, warn = FALSE), collapse = "\n")
    expect_true(grepl("tryCatch(readRDS(path)", code, fixed = TRUE))
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
  expect_false(any(grepl("lipid|rand", manifest$artifact_id, ignore.case = TRUE)))
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
