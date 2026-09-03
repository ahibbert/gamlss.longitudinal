local_public_approval_env <- function() {
  env <- new.env(parent = globalenv())
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  source(file.path(root, "paper", "R", "replication-helpers.R"), local = env)
  list(env = env, root = root)
}

local_clear_public_approval_env <- function(spec) {
  vars <- c(spec$attestation_env, spec$signature_env,
    "GAMLSS_LONGITUDINAL_JSS_APPROVAL_BUNDLE_SHA256")
  old <- Sys.getenv(vars, unset = NA_character_)
  Sys.unsetenv(vars)
  function() {
    Sys.unsetenv(vars)
    present <- !is.na(old)
    if (any(present)) do.call(Sys.setenv, as.list(stats::setNames(old[present], vars[present])))
  }
}

local_write_public_approval_fixture <- function(env, root) {
  approval_dir <- file.path(root, "paper", "approvals")
  dir.create(approval_dir, recursive = TRUE)
  spec <- env$jss_phase2_public_approval_spec()
  writeLines("fixture", file.path(approval_dir, "README.md"))
  for (i in seq_len(nrow(spec))) {
    writeBin(as.raw(c(i, 1:16)), file.path(approval_dir, spec$attestation_file[[i]]))
    writeBin(as.raw(c(i, 17:32)), file.path(approval_dir, spec$signature_file[[i]]))
  }
  manifest <- data.frame(
    study = spec$study,
    attestation_file = spec$attestation_file,
    attestation_sha256 = vapply(file.path(approval_dir, spec$attestation_file),
      digest::digest, character(1), file = TRUE, algo = "sha256", serialize = FALSE),
    signature_file = spec$signature_file,
    signature_sha256 = vapply(file.path(approval_dir, spec$signature_file),
      digest::digest, character(1), file = TRUE, algo = "sha256", serialize = FALSE),
    public_key_sha256 = rep("cb73c05cede55bfd56357b1780b90c1bf254413d82765c7a682fc9db3a0d8587", nrow(spec)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, file.path(approval_dir, "manifest.csv"), row.names = FALSE, na = "")
  approval_dir
}

test_that("smoke replication does not require Phase 2 approvals", {
  setup <- local_public_approval_env()
  expect_identical(setup$env$jss_phase2_materialize_public_approvals(setup$root, "smoke"), character())
})

test_that("reviewer entry point materializes approvals before targets run", {
  setup <- local_public_approval_env()
  code <- readLines(file.path(setup$root, "paper", "replicate.R"), warn = FALSE)
  materialize <- grep("jss_phase2_materialize_public_approvals", code, fixed = TRUE)
  make <- grep("targets::tar_make", code, fixed = TRUE)
  expect_length(materialize, 1L)
  expect_true(length(make) >= 1L)
  expect_lt(materialize, min(make))
})

test_that("public approvals are hash checked and materialized outside checkout", {
  skip_if_not_installed("digest")
  setup <- local_public_approval_env()
  spec <- setup$env$jss_phase2_public_approval_spec()
  restore <- local_clear_public_approval_env(spec); on.exit(restore(), add = TRUE)
  fixture_root <- tempfile("approval-checkout-"); dir.create(fixture_root)
  approval_dir <- local_write_public_approval_fixture(setup$env, fixture_root)
  external <- tempfile("approval-external-")
  paths <- setup$env$jss_phase2_materialize_public_approvals(
    fixture_root, "paper", approval_dir = approval_dir, external_dir = external
  )
  expect_length(paths, 2L * nrow(spec))
  expect_true(all(file.exists(paths)))
  expect_true(all(!startsWith(
    tolower(normalizePath(paths, winslash = "/")),
    paste0(tolower(normalizePath(fixture_root, winslash = "/")), "/")
  )))
  expect_identical(unname(normalizePath(Sys.getenv(spec$attestation_env), winslash = "/")),
    unname(normalizePath(paths[seq_len(nrow(spec))], winslash = "/")))
  expect_identical(unname(normalizePath(Sys.getenv(spec$signature_env), winslash = "/")),
    unname(normalizePath(paths[nrow(spec) + seq_len(nrow(spec))], winslash = "/")))
  first_identity <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_APPROVAL_BUNDLE_SHA256")
  Sys.unsetenv(c(spec$attestation_env, spec$signature_env,
    "GAMLSS_LONGITUDINAL_JSS_APPROVAL_BUNDLE_SHA256"))
  setup$env$jss_phase2_materialize_public_approvals(
    fixture_root, "paper", approval_dir = approval_dir,
    external_dir = tempfile("approval-external-")
  )
  expect_identical(Sys.getenv("GAMLSS_LONGITUDINAL_JSS_APPROVAL_BUNDLE_SHA256"),
    first_identity)

  writeBin(as.raw(0:10), file.path(approval_dir, spec$signature_file[[1L]]))
  Sys.unsetenv(c(spec$attestation_env, spec$signature_env))
  expect_error(setup$env$jss_phase2_materialize_public_approvals(
    fixture_root, "full", approval_dir = approval_dir,
    external_dir = tempfile("approval-external-")
  ), "does not match")
})

test_that("approval bundle rejects partial overrides and private-key-like files", {
  skip_if_not_installed("digest")
  setup <- local_public_approval_env()
  spec <- setup$env$jss_phase2_public_approval_spec()
  restore <- local_clear_public_approval_env(spec); on.exit(restore(), add = TRUE)
  Sys.setenv(GAMLSS_LONGITUDINAL_JSS_MAIN_RECOVERY_ATTESTATION = "only-one-half")
  expect_error(setup$env$jss_phase2_materialize_public_approvals(setup$root, "paper"),
    "requires both")
  Sys.unsetenv(c(spec$attestation_env, spec$signature_env))

  fixture_root <- tempfile("approval-checkout-"); dir.create(fixture_root)
  approval_dir <- local_write_public_approval_fixture(setup$env, fixture_root)
  writeBin(as.raw(1:8), file.path(approval_dir, "production-private.key"))
  expect_error(setup$env$jss_phase2_materialize_public_approvals(
    fixture_root, "paper", approval_dir = approval_dir,
    external_dir = tempfile("approval-external-")
  ), "private-key-like")
})

test_that("complete per-study external overrides remain supported", {
  skip_if_not_installed("digest")
  setup <- local_public_approval_env()
  spec <- setup$env$jss_phase2_public_approval_spec()
  restore <- local_clear_public_approval_env(spec); on.exit(restore(), add = TRUE)
  fixture_root <- tempfile("approval-checkout-"); dir.create(fixture_root)
  approval_dir <- local_write_public_approval_fixture(setup$env, fixture_root)
  override <- c(tempfile("external-attestation-"), tempfile("external-signature-"))
  writeBin(as.raw(1:8), override[[1L]]); writeBin(as.raw(9:16), override[[2L]])
  do.call(Sys.setenv, stats::setNames(as.list(override),
    c(spec$attestation_env[[1L]], spec$signature_env[[1L]])))
  paths <- setup$env$jss_phase2_materialize_public_approvals(
    fixture_root, "paper", approval_dir = approval_dir,
    external_dir = tempfile("approval-external-")
  )
  expect_identical(normalizePath(paths[[1L]], winslash = "/"),
    normalizePath(override[[1L]], winslash = "/"))
  expect_identical(normalizePath(paths[[nrow(spec) + 1L]], winslash = "/"),
    normalizePath(override[[2L]], winslash = "/"))
  expect_match(Sys.getenv("GAMLSS_LONGITUDINAL_JSS_APPROVAL_BUNDLE_SHA256"),
    "^[0-9a-f]{64}$")
})

test_that("external overrides are preflighted before targets run", {
  skip_if_not_installed("digest")
  setup <- local_public_approval_env()
  spec <- setup$env$jss_phase2_public_approval_spec()
  restore <- local_clear_public_approval_env(spec); on.exit(restore(), add = TRUE)
  fixture_root <- tempfile("approval-checkout-"); dir.create(fixture_root)
  approval_dir <- local_write_public_approval_fixture(setup$env, fixture_root)

  inside <- file.path(fixture_root, "inside.attestation")
  writeBin(as.raw(1:8), inside)
  outside <- tempfile("outside-signature-"); writeBin(as.raw(9:16), outside)
  do.call(Sys.setenv, stats::setNames(as.list(c(inside, outside)),
    c(spec$attestation_env[[1L]], spec$signature_env[[1L]])))
  expect_error(setup$env$jss_phase2_materialize_public_approvals(
    fixture_root, "paper", approval_dir = approval_dir,
    external_dir = tempfile("approval-external-")), "external to the editable checkout")

  Sys.unsetenv(c(spec$attestation_env, spec$signature_env))
  complete <- vapply(seq_len(2L * nrow(spec)), function(i) {
    path <- tempfile("complete-override-"); writeBin(as.raw(i), path); path
  }, character(1L))
  unlink(complete[[length(complete)]])
  do.call(Sys.setenv, stats::setNames(as.list(complete[seq_len(nrow(spec))]), spec$attestation_env))
  do.call(Sys.setenv, stats::setNames(as.list(complete[nrow(spec) + seq_len(nrow(spec))]), spec$signature_env))
  expect_error(setup$env$jss_phase2_materialize_public_approvals(
    fixture_root, "full"), "does not exist or is not a regular file")
})

test_that("a complete external override set needs no bundled manifest", {
  skip_if_not_installed("digest")
  setup <- local_public_approval_env()
  spec <- setup$env$jss_phase2_public_approval_spec()
  restore <- local_clear_public_approval_env(spec); on.exit(restore(), add = TRUE)
  fixture_root <- tempfile("approval-checkout-"); dir.create(fixture_root)
  complete <- vapply(seq_len(2L * nrow(spec)), function(i) {
    path <- tempfile("complete-override-"); writeBin(as.raw(i), path); path
  }, character(1L))
  do.call(Sys.setenv, stats::setNames(as.list(complete[seq_len(nrow(spec))]), spec$attestation_env))
  do.call(Sys.setenv, stats::setNames(as.list(complete[nrow(spec) + seq_len(nrow(spec))]), spec$signature_env))
  paths <- setup$env$jss_phase2_materialize_public_approvals(fixture_root, "paper")
  expect_identical(unname(normalizePath(paths, winslash = "/")),
    unname(normalizePath(complete, winslash = "/")))
  expect_match(Sys.getenv("GAMLSS_LONGITUDINAL_JSS_APPROVAL_BUNDLE_SHA256"),
    "^[0-9a-f]{64}$")
})

test_that("effective approval-byte changes invalidate cached targets", {
  skip_if_not_installed("digest")
  skip_if_not_installed("targets")
  skip_if_not_installed("callr")
  setup <- local_public_approval_env()
  spec <- setup$env$jss_phase2_public_approval_spec()
  restore <- local_clear_public_approval_env(spec); on.exit(restore(), add = TRUE)
  fixture_root <- tempfile("approval-checkout-"); dir.create(fixture_root)
  approval_dir <- local_write_public_approval_fixture(setup$env, fixture_root)
  work <- tempfile("approval-target-cache-"); dir.create(work)
  script <- file.path(work, "_targets.R")
  store <- file.path(work, "_targets")
  helper <- normalizePath(file.path(setup$root, "paper", "R", "replication-helpers.R"),
    winslash = "/", mustWork = TRUE)
  writeLines(c(
    "library(targets)",
    paste0("source(", dQuote(helper), ")"),
    "list(",
    "  tar_target(settings, jss_settings(create = FALSE), cue = tar_cue(mode = 'always')),",
    "  tar_target(approval_identity, settings$approval_bundle_sha256)",
    ")"
  ), script)
  old <- Sys.getenv(c("GAMLSS_LONGITUDINAL_JSS_PROFILE", "GAMLSS_LONGITUDINAL_JSS_ROOT"),
    unset = NA_character_)
  on.exit({
    Sys.unsetenv(c("GAMLSS_LONGITUDINAL_JSS_PROFILE", "GAMLSS_LONGITUDINAL_JSS_ROOT"))
    present <- !is.na(old)
    if (any(present)) do.call(Sys.setenv, as.list(stats::setNames(old[present], names(old)[present])))
  }, add = TRUE)
  Sys.setenv(GAMLSS_LONGITUDINAL_JSS_PROFILE = "paper",
    GAMLSS_LONGITUDINAL_JSS_ROOT = fixture_root)

  setup$env$jss_phase2_materialize_public_approvals(fixture_root, "paper",
    approval_dir = approval_dir, external_dir = tempfile("approval-external-"))
  first_identity <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_APPROVAL_BUNDLE_SHA256")
  targets::tar_make(script = script, store = store, reporter = "silent")
  expect_identical(targets::tar_read(approval_identity, store = store), first_identity)

  writeBin(as.raw(c(255L, 1:16)), file.path(approval_dir, spec$attestation_file[[1L]]))
  manifest <- utils::read.csv(file.path(approval_dir, "manifest.csv"), stringsAsFactors = FALSE)
  manifest$attestation_sha256[[1L]] <- digest::digest(
    file.path(approval_dir, spec$attestation_file[[1L]]), file = TRUE,
    algo = "sha256", serialize = FALSE)
  utils::write.csv(manifest, file.path(approval_dir, "manifest.csv"), row.names = FALSE, na = "")
  Sys.unsetenv(c(spec$attestation_env, spec$signature_env,
    "GAMLSS_LONGITUDINAL_JSS_APPROVAL_BUNDLE_SHA256"))
  setup$env$jss_phase2_materialize_public_approvals(fixture_root, "paper",
    approval_dir = approval_dir, external_dir = tempfile("approval-external-"))
  second_identity <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_APPROVAL_BUNDLE_SHA256")
  expect_false(identical(second_identity, first_identity))
  targets::tar_make(script = script, store = store, reporter = "silent")
  expect_identical(targets::tar_read(approval_identity, store = store), second_identity)
})

test_that("targets preflight recomputes approval identity and rejects tampering", {
  skip_if_not_installed("digest")
  setup <- local_public_approval_env()
  spec <- setup$env$jss_phase2_public_approval_spec()
  restore <- local_clear_public_approval_env(spec); on.exit(restore(), add = TRUE)
  fixture_root <- tempfile("approval-checkout-"); dir.create(fixture_root)
  approval_dir <- local_write_public_approval_fixture(setup$env, fixture_root)

  Sys.setenv(GAMLSS_LONGITUDINAL_JSS_APPROVAL_BUNDLE_SHA256 = paste(rep("a", 64L), collapse = ""))
  expect_error(setup$env$jss_phase2_verify_materialized_public_approvals(
    fixture_root, "paper"), "require all Phase 2 approvals")

  external <- tempfile("approval-external-")
  paths <- setup$env$jss_phase2_materialize_public_approvals(
    fixture_root, "paper", approval_dir = approval_dir, external_dir = external
  )
  expect_match(setup$env$jss_phase2_verify_materialized_public_approvals(
    fixture_root, "paper"), "^[0-9a-f]{64}$")

  writeBin(as.raw(c(255L, 1:8)), paths[[1L]])
  expect_error(setup$env$jss_phase2_verify_materialized_public_approvals(
    fixture_root, "paper"), "do not match their declared bundle identity")
})
