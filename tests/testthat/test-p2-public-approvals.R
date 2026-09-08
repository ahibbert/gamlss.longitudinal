test_that("reviewer entry points do not require approval credentials", {
  replicate_code <- paste(readLines(test_path("..", "..", "paper", "replicate.R"), warn = FALSE), collapse = "\n")
  target_code <- paste(readLines(test_path("..", "..", "paper", "_targets.R"), warn = FALSE), collapse = "\n")

  expect_false(grepl("materialize_public_approvals", replicate_code, fixed = TRUE))
  expect_false(grepl("verify_materialized_public_approvals", target_code, fixed = TRUE))
  expect_false(grepl("APPROVAL", target_code, fixed = TRUE))
})

test_that("ordinary settings contain no approval identity", {
  env <- new.env(parent = globalenv())
  sys.source(test_path("..", "..", "paper", "R", "replication-helpers.R"), envir = env)
  root <- normalizePath(test_path("..", ".."), winslash = "/", mustWork = TRUE)
  old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(root)
  settings <- env$jss_settings(create = FALSE)
  expect_false("approval_bundle_sha256" %in% names(settings))
})
