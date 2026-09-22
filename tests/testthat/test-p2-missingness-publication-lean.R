local_lean_missingness_env <- function() {
  root <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE)
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "paper", "R", "paper-plot-style.R"), envir = env)
  sys.source(file.path(root, "paper", "R", "missingness-study-helpers.R"), envir = env)
  sys.source(file.path(root, "paper", "R", "04-missingness-dropout-sensitivity.R"), envir = env)
  sys.source(file.path(root, "paper", "R", "04-missingness-publication-lean.R"), envir = env)
  env
}

test_that("lean metric summaries distinguish zero, one, and estimable states", {
  env <- local_lean_missingness_env()
  zero <- env$jss_missingness_lean_metric(numeric())
  one <- env$jss_missingness_lean_metric(2)
  many <- env$jss_missingness_lean_metric(c(1, 2, 3))
  expect_identical(unname(zero[["status"]]), "no_success")
  expect_true(all(is.na(as.numeric(zero[c("estimate", "mcse", "conf_low", "conf_high")]))))
  expect_identical(unname(one[["status"]]), "one_success_no_mcse")
  expect_equal(as.numeric(one[["estimate"]]), 2)
  expect_true(all(is.na(as.numeric(one[c("mcse", "conf_low", "conf_high")]))))
  expect_identical(unname(many[["status"]]), "estimated")
  expect_true(all(is.finite(as.numeric(many[c("estimate", "mcse", "conf_low", "conf_high")]))))
})

test_that("interior gaps and delayed entry are separate registered concepts", {
  env <- local_lean_missingness_env()
  pattern <- data.frame(
    id = rep(1:4, each = 4), time = rep(1:4, 4),
    response_observed = c(TRUE, FALSE, TRUE, TRUE,
      FALSE, TRUE, TRUE, TRUE,
      TRUE, TRUE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE)
  )
  audit <- env$jss_missingness_lean_subject_audit(pattern)
  expect_equal(unname(audit[["n_subjects_with_interior_gaps"]]), 1)
  expect_equal(unname(audit[["n_subjects_with_delayed_entry"]]), 1)
  expect_equal(unname(audit[["n_subjects_requiring_segmented_likelihood"]]), 2)
})

test_that("metric-state validation rejects impossible zero and one-success rows", {
  env <- local_lean_missingness_env()
  valid <- data.frame(
    n_metric_successful = c(0L, 1L, 2L),
    metric_status = c("no_success", "one_success_no_mcse", "estimated"),
    estimate = c(NA, 2, 2), mcse = c(NA, NA, .1),
    conf_low = c(NA, NA, 1.8), conf_high = c(NA, NA, 2.2))
  expect_invisible(env$jss_missingness_lean_validate_metric_state(valid,
    "n_metric_successful", "estimate", c("mcse", "conf_low", "conf_high")))
  invalid_zero <- valid; invalid_zero$estimate[[1L]] <- 0
  expect_error(env$jss_missingness_lean_validate_metric_state(invalid_zero,
    "n_metric_successful", "estimate", c("mcse", "conf_low", "conf_high")),
    "inconsistent")
  invalid_one <- valid; invalid_one$mcse[[2L]] <- 0
  expect_error(env$jss_missingness_lean_validate_metric_state(invalid_one,
    "n_metric_successful", "estimate", c("mcse", "conf_low", "conf_high")),
    "uncertainty")
})

test_that("registered design comparison rejects a meaningful grid mutation", {
  env <- local_lean_missingness_env()
  expected <- env$jss_missing_registered_task_grid()
  mutated <- expected
  mutated$target_missing_rate[[1L]] <- mutated$target_missing_rate[[1L]] + 1e-6
  expect_error(env$jss_missing_compare_frame(mutated, expected,
    c("scenario", "rep"), "registered design", tolerance = 1e-15),
    "target_missing_rate")
})

test_that("warning policy is narrow and unknown text is fatal", {
  env <- local_lean_missingness_env()
  known <- data.frame(
    scenario = "s", rep = 1L, missing_mechanism = "monotone_dropout",
    target_missing_rate = 0, model = "gamlss.longitudinal", warning_index = 1L,
    normalized_message = paste0(
      "Hessian-based inference is unavailable (fitted_gradient_too_large). ",
      "The fitted model is unchanged."), stringsAsFactors = FALSE)
  classified <- env$jss_missingness_lean_warning_events(known)
  expect_identical(classified$classification, "inference_unavailable")
  expect_false(classified$coverage_eligible)
  unknown <- known
  unknown$normalized_message <- "a new numerical warning"
  expect_error(env$jss_missingness_lean_warning_events(unknown), "fatal")
})

test_that("missingness presentation retains failure denominators", {
  env <- local_lean_missingness_env()
  attempts <- data.frame(
    missing_mechanism = c("monotone_dropout", "time_dependent_intermittent_mar"),
    target_missing_rate = c(.3, .4),
    model = c("gamlss.longitudinal", "gamlss2"),
    attempted = c(20L, 20L), retained = c(13L, 20L),
    retention_rate_conf_low = c(.4329, .8389),
    retention_rate_conf_high = c(.8188, 1),
    stringsAsFactors = FALSE
  )
  labels <- c(
    monotone_dropout = "Monotone dropout (headline)",
    time_dependent_intermittent_mar = "Intermittent MAR (sensitivity)"
  )
  plot_data <- env$jss_missingness_retention_plot_data(attempts, labels)
  expect_equal(plot_data$retention_pct, c(65, 100))
  expect_identical(plot_data$retention_label, c("13/20", "20/20"))
  expect_identical(as.character(plot_data$mechanism_panel), unname(labels))

  invalid <- attempts
  invalid$retained[[1L]] <- 21L
  expect_error(
    env$jss_missingness_retention_plot_data(invalid, labels),
    "invalid mechanism or denominator"
  )
})

test_that("reviewer artifacts are CSV-only and omit release forensics", {
  env <- local_lean_missingness_env()
  files <- env$jss_missingness_lean_artifacts()
  expect_true(all(grepl("[.]csv$", files)))
  expect_false(any(grepl("checkpoint|attestation|signature|operational|[.]rds$", files,
    ignore.case = TRUE)))
  expect_false(env$jss_missingness_lean_no_private_location(
    data.frame(value = "C:/Users/reviewer/private/results.csv")))
  expect_false(env$jss_missingness_lean_no_private_location(
    data.frame(value = "/home/reviewer/results.csv")))
  expect_true(env$jss_missingness_lean_no_private_location(
    data.frame(value = "monotone_dropout")))

  root <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE)
  public <- file.path(root, "paper", "data", "public-derived", "missingness")
  tables <- lapply(file.path(public, files), utils::read.csv,
    stringsAsFactors = FALSE, check.names = FALSE)
  expect_true(all(vapply(tables, env$jss_missingness_lean_no_private_location,
    logical(1))))
  expect_false(any(grepl("host|pid|path|checkpoint|signature|attestation|quarantine",
    unlist(lapply(tables, names), use.names = FALSE), ignore.case = TRUE)))
})

test_that("manifest and artifact contracts expose the same lean missingness files", {
  env <- local_lean_missingness_env()
  root <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE)
  manifest <- utils::read.csv(file.path(root, "paper", "manifest.csv"), stringsAsFactors = FALSE)
  contract <- utils::read.csv(file.path(root, "paper", "phase2-artifact-contract.csv"),
    stringsAsFactors = FALSE)
  manifest <- manifest[manifest$producer == "module_04_missingness" &
    manifest$artifact_type != "figure", , drop = FALSE]
  contract <- contract[contract$producer == "module_04_missingness" &
    contract$artifact_type != "figure", , drop = FALSE]
  expect_setequal(basename(manifest$generated_path), env$jss_missingness_lean_artifacts())
  expect_setequal(manifest$artifact_id, contract$artifact_id)
  expect_false(any(manifest$artifact_type == "checkpoint_archive"))
})
