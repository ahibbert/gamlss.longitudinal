local_lean_env <- function() {
  env <- new.env(parent = globalenv())
  root <- normalizePath(test_path("..", ".."), winslash = "/", mustWork = TRUE)
  sys.source(file.path(root, "paper", "R", "replication-helpers.R"), envir = env)
  sys.source(file.path(root, "paper", "R", "submission-replication.R"), envir = env)
  list(env = env, root = root)
}

test_that("submission registry names exactly the five active studies", {
  setup <- local_lean_env()
  registry <- setup$env$jss_submission_input_registry(setup$root)
  expect_identical(registry$schema_version, rep("jss-submission-inputs-1", 5L))
  expect_identical(registry$study, c(
    "main-recovery", "optimizer-benchmark", "missingness",
    "copula-misspecification", "fit-scaling"
  ))
})

test_that("registered inputs reject duplicate keys and invalid rates", {
  setup <- local_lean_env()
  registry <- setup$env$jss_submission_input_registry(setup$root)
  scratch <- tempfile("lean-schema-"); dir.create(file.path(scratch, "paper"), recursive = TRUE)
  file.copy(file.path(setup$root, "paper", "submission-input-schemas.csv"),
    file.path(scratch, "paper", "submission-input-schemas.csv"))
  valid <- data.frame(
    scenario_id = "a", attempt_id = 1L, seed = 10L, attempted = TRUE,
    fit_returned = TRUE, converged = TRUE, retained = TRUE,
    failure_reason = "", elapsed_sec = 1
  )
  input <- file.path(scratch, "attempts.csv")
  utils::write.csv(valid, input, row.names = FALSE, na = "")
  expect_s3_class(setup$env$jss_submission_read_registered_input(
    "fit-scaling", input, scratch), "data.frame")
  utils::write.csv(rbind(valid, valid), input, row.names = FALSE, na = "")
  expect_error(setup$env$jss_submission_read_registered_input(
    "fit-scaling", input, scratch), "duplicate keys")

  invalid_lifecycle <- valid
  invalid_lifecycle$converged <- FALSE
  utils::write.csv(invalid_lifecycle, input, row.names = FALSE, na = "")
  expect_error(setup$env$jss_submission_read_registered_input(
    "fit-scaling", input, scratch), "lifecycle")

  rate_registry <- registry[registry$study == "copula-misspecification", , drop = FALSE]
  rate_registry$min_rows <- 1L; rate_registry$exact_rows <- NA
  utils::write.csv(rate_registry, file.path(scratch, "paper", "submission-input-schemas.csv"),
    row.names = FALSE, na = "")
  bad <- data.frame(
    fit_id = "x", rep = 1L, generating_copula = "N", fitted_copula = "N",
    tau_label = "high", n_subject = 500L, n_time = 4L,
    attempted = TRUE, success = TRUE, converged = TRUE,
    retained = TRUE, failure_type = "", target_tau = 1.1
  )
  utils::write.csv(bad, input, row.names = FALSE, na = "")
  expect_error(setup$env$jss_submission_read_registered_input(
    "copula-misspecification", input, scratch), "outside \\[0, 1\\]")
})

test_that("Module 09 is deferred from every submission interface", {
  setup <- local_lean_env()
  manifest <- utils::read.csv(file.path(setup$root, "paper", "manifest.csv"), stringsAsFactors = FALSE)
  contract <- utils::read.csv(file.path(setup$root, "paper", "phase2-artifact-contract.csv"), stringsAsFactors = FALSE)
  target_code <- paste(readLines(file.path(setup$root, "paper", "_targets.R"), warn = FALSE), collapse = "\n")
  rows <- grepl("^p2_mvt_", manifest$artifact_id)
  expect_true(any(rows))
  expect_true(all(manifest$publication_status[rows] == "deferred"))
  expect_true(all(manifest$profiles[rows] == ""))
  expect_false(any(grepl("^p2_mvt_", contract$artifact_id)))
  expect_false(grepl("module_09_multivariate", target_code, fixed = TRUE))
  expect_true(file.exists(file.path(setup$root, "paper", "DEFERRED-MULTIVARIATE-STUDY.md")))
})
