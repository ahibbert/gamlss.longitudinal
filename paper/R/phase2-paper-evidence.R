# Central Phase 2 staging and gate helpers. These functions keep paper-profile
# execution read-only with respect to tracked public-derived inputs.

jss_phase2_optimizer_files <- function() {
  list(
    data = c(
      "03-joint-vs-separate-optimization-candidate-selection.csv",
      "03-joint-vs-separate-optimization-mc-precision.csv",
      "03-joint-vs-separate-optimization-checkpoint-status.csv",
      "03-joint-vs-separate-optimization-results.csv",
      "03-joint-vs-separate-optimization-figure-registry.csv",
      "03-joint-vs-separate-optimization-checkpoint-payloads.rds",
      "03-joint-vs-separate-optimization-checkpoint-content-manifest.csv"
    ),
    tables = c(
      "03-joint-vs-separate-optimization-summary.csv",
      "03-joint-vs-separate-optimization-metric-wins.csv",
      "03-joint-vs-separate-optimization-hypothesis-summary.csv",
      "03-joint-vs-separate-optimization-failures.csv",
      "03-joint-vs-separate-optimization-difference-uncertainty.csv",
      "03-joint-vs-separate-optimization-hypothesis-evidence.csv"
    ),
    figures = c(
      "03-joint-vs-separate-optimization-deltas.png",
      "03-joint-vs-separate-optimization-metric-dashboard.png"
    )
  )
}

jss_phase2_copy_exact <- function(source, destination) {
  if (length(source) != length(destination) || any(!file.exists(source))) {
    missing <- source[!file.exists(source)]
    stop("Phase 2 source bundle is incomplete: ", paste(basename(missing), collapse = ", "), call. = FALSE)
  }
  invisible(lapply(unique(dirname(destination)), dir.create, recursive = TRUE, showWarnings = FALSE))
  copied <- file.copy(source, destination, overwrite = TRUE)
  if (any(!copied)) stop("Could not stage Phase 2 artifact(s): ", paste(basename(source[!copied]), collapse = ", "), call. = FALSE)
  destination
}

jss_phase2_compare_frame <- function(actual, expected, keys, label, tolerance = 1e-10) {
  if (!setequal(names(actual), names(expected))) stop(label, " schema is not canonical.", call. = FALSE)
  sort_frame <- function(x) {
    ord <- if (!nrow(x)) integer() else if (!length(keys)) seq_len(nrow(x)) else
      do.call(order, c(unname(lapply(x[keys], function(value) unname(as.character(value)))),
        list(na.last = TRUE)))
    x[ord, names(expected), drop = FALSE]
  }
  actual <- sort_frame(actual); expected <- sort_frame(expected)
  rownames(actual) <- NULL; rownames(expected) <- NULL
  if (nrow(actual) != nrow(expected)) stop(label, " row count is not canonical.", call. = FALSE)
  for (nm in names(expected)) {
    ok <- if (is.numeric(expected[[nm]])) {
      a <- as.numeric(actual[[nm]]); e <- as.numeric(expected[[nm]])
      a[is.nan(a)] <- NA_real_; e[is.nan(e)] <- NA_real_
      isTRUE(all.equal(a, e,
        tolerance = tolerance, check.attributes = FALSE))
    } else identical(as.character(actual[[nm]]), as.character(expected[[nm]]))
    if (!ok) stop(label, " disagrees with attempt rows in column ", nm, ".", call. = FALSE)
  }
  invisible(TRUE)
}

jss_phase2_parse_rfc3339_utc <- function(x, label = "timestamp") {
  if (!is.character(x) || anyNA(x) ||
      any(!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$", x))) {
    stop(label, " is not an RFC3339 UTC instant.", call. = FALSE)
  }
  parsed <- as.POSIXct(strptime(x, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  roundtrip <- format(parsed, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  if (anyNA(parsed) || !identical(roundtrip, x)) {
    stop(label, " is not a real RFC3339 UTC instant.", call. = FALSE)
  }
  parsed
}

jss_phase2_validate_approval_order <- function(approved_at_utc, checkpoint_timestamps) {
  approved <- jss_phase2_parse_rfc3339_utc(approved_at_utc, "Approval timestamp")
  executed <- jss_phase2_parse_rfc3339_utc(checkpoint_timestamps, "Checkpoint timestamp")
  if (approved <= max(executed)) stop("Promotion approval must occur after every checkpoint execution.", call. = FALSE)
  invisible(TRUE)
}

jss_phase2_require_external_attestation <- function(study, bundle_sha256, package_sha256,
    producer_sha256, checkpoints, root, attestation_path, signature_path) {
  if (!requireNamespace("sodium", quietly = TRUE)) {
    stop("sodium is required for production promotion-signature verification.", call. = FALSE)
  }
  if (!nzchar(attestation_path) || !nzchar(signature_path)) {
    stop(study, " public bundle lacks a detached production promotion signature.", call. = FALSE)
  }
  root <- paste0(normalizePath(root, winslash = "/", mustWork = TRUE), "/")
  path <- normalizePath(attestation_path, winslash = "/", mustWork = TRUE)
  signature <- normalizePath(signature_path, winslash = "/", mustWork = TRUE)
  if (startsWith(tolower(paste0(path, "/")), tolower(root)) ||
      startsWith(tolower(paste0(signature, "/")), tolower(root))) {
    stop(study, " promotion attestation and signature must be external to the editable checkout.", call. = FALSE)
  }
  public_key <- as.raw(c(66, 91, 71, 233, 21, 247, 172, 45, 215, 202, 170, 0, 64,
    43, 83, 206, 23, 50, 48, 154, 25, 217, 178, 37, 252, 59, 158, 195, 237, 0, 31, 216))
  message_raw <- readBin(path, "raw", n = file.info(path)$size)
  signature_raw <- readBin(signature, "raw", n = file.info(signature)$size)
  verified <- tryCatch(sodium::sig_verify(message_raw, signature_raw, public_key),
    error = function(e) FALSE)
  if (!isTRUE(verified)) {
    stop(study, " promotion attestation lacks a valid detached production signature.", call. = FALSE)
  }
  attestation <- tryCatch(unserialize(message_raw), error = function(e) NULL)
  expected_names <- c("schema_version", "study", "bundle_sha256", "package_source_sha256",
    "producer_sha256", "approved_at_utc", "approver", "checkpoint_manifest")
  if (!is.list(attestation) || !identical(names(attestation), expected_names) ||
      !identical(attestation$schema_version, 1L) || !identical(attestation$study, study) ||
      !identical(attestation$bundle_sha256, bundle_sha256) ||
      !identical(attestation$package_source_sha256, package_sha256) ||
      !identical(attestation$producer_sha256, producer_sha256) ||
      !is.character(attestation$approver) || length(attestation$approver) != 1L ||
      !nzchar(attestation$approver) || !identical(attestation$checkpoint_manifest, checkpoints)) {
    stop(study, " external promotion attestation does not bind the canonical checkpoint manifest.", call. = FALSE)
  }
  jss_phase2_validate_approval_order(attestation$approved_at_utc,
    checkpoints$checkpoint_timestamp_utc)
  invisible(TRUE)
}

jss_phase2_validate_png <- function(path) {
  bytes <- readBin(path, "raw", n = file.info(path)$size)
  signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  if (length(bytes) < 1000L || !identical(bytes[seq_along(signature)], signature)) return(FALSE)
  u32 <- function(pos) sum(as.integer(bytes[pos + 0:3]) * 256^(3:0))
  pos <- 9L; idat <- raw(); width <- height <- bit_depth <- colour_type <- NA_real_; ended <- FALSE
  while (pos + 11L <= length(bytes)) {
    n <- u32(pos); type <- rawToChar(bytes[(pos + 4L):(pos + 7L)])
    data_start <- pos + 8L; data_end <- data_start + n - 1L
    if (n < 0 || data_end + 4L > length(bytes)) return(FALSE)
    payload <- if (n) bytes[data_start:data_end] else raw()
    if (type == "IHDR" && n == 13L) {
      width <- sum(as.integer(payload[1:4]) * 256^(3:0)); height <- sum(as.integer(payload[5:8]) * 256^(3:0))
      bit_depth <- as.integer(payload[[9L]]); colour_type <- as.integer(payload[[10L]])
    } else if (type == "IDAT") idat <- c(idat, payload) else if (type == "IEND") { ended <- TRUE; break }
    pos <- data_end + 5L
  }
  if (!ended || !is.finite(width) || width < 2 || height < 2 || !length(idat)) return(FALSE)
  decoded <- tryCatch(memDecompress(idat, type = "gzip"), error = function(e) raw())
  channels <- c(`0` = 1, `2` = 3, `3` = 1, `4` = 2, `6` = 4)[as.character(colour_type)]
  expected_min <- height * (1 + ceiling(width * bit_depth * channels / 8))
  length(decoded) >= expected_min && length(unique(as.integer(decoded))) > 2L
}

jss_phase2_optimizer_normalize_attempt_csv <- function(attempts) {
  if (exists("jss_joint_normalize_results_csv", mode = "function")) {
    return(jss_joint_normalize_results_csv(attempts))
  }
  normalize <- function(x, empty) {
    if (is.logical(x) && all(is.na(x))) return(if (empty) rep("", length(x)) else rep(NA_character_, length(x)))
    if (!is.character(x)) stop("Optimizer CSV diagnostic columns have invalid declared types.", call. = FALSE)
    if (empty) x[is.na(x)] <- ""
    x
  }
  attempts$warnings <- normalize(attempts$warnings, TRUE)
  attempts$error <- normalize(attempts$error, FALSE)
  attempts
}

jss_phase2_optimizer_validate_candidate <- function(data_dir, tables_dir, figures_dir, root = getwd(),
    attestation_path = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_OPTIMIZER_ATTESTATION", unset = ""),
    signature_path = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_OPTIMIZER_ATTESTATION_SIGNATURE", unset = ""),
    require_promotion = TRUE) {
  files <- jss_phase2_optimizer_files()
  required <- c(file.path(data_dir, files$data), file.path(tables_dir, files$tables), file.path(figures_dir, files$figures))
  if (anyDuplicated(normalizePath(required, winslash = "/", mustWork = FALSE)) || anyDuplicated(basename(required))) {
    stop("Optimizer bundle artifact paths are not uniquely canonical.", call. = FALSE)
  }
  if (any(!file.exists(required))) {
    stop("Optimizer benchmark bundle is incomplete: ", paste(basename(required[!file.exists(required)]), collapse = ", "), call. = FALSE)
  }
  bundle_manifest <- paste(sort(paste(basename(required),
    vapply(required, jss_joint_sha256_file, character(1)), sep = "\t")), collapse = "\n")
  bundle_sha256 <- digest::digest(bundle_manifest, algo = "sha256", serialize = FALSE)
  design <- utils::read.csv(file.path(data_dir, files$data[[1L]]), stringsAsFactors = FALSE)
  precision <- utils::read.csv(file.path(data_dir, files$data[[2L]]), stringsAsFactors = FALSE)
  checkpoints <- utils::read.csv(file.path(data_dir, files$data[[3L]]), stringsAsFactors = FALSE)
  attempts <- utils::read.csv(file.path(data_dir, files$data[[4L]]), stringsAsFactors = FALSE)
  uncertainty <- utils::read.csv(file.path(tables_dir, files$tables[[5L]]), stringsAsFactors = FALSE)
  failures <- utils::read.csv(file.path(tables_dir, files$tables[[4L]]), stringsAsFactors = FALSE)
  summary_reported <- utils::read.csv(file.path(tables_dir, files$tables[[1L]]), stringsAsFactors = FALSE)
  wins_reported <- utils::read.csv(file.path(tables_dir, files$tables[[2L]]), stringsAsFactors = FALSE)
  hypothesis_reported <- utils::read.csv(file.path(tables_dir, files$tables[[3L]]), stringsAsFactors = FALSE)
  evidence_reported <- utils::read.csv(file.path(tables_dir, files$tables[[6L]]), stringsAsFactors = FALSE)
  figure_registry <- utils::read.csv(file.path(data_dir, files$data[[5L]]), stringsAsFactors = FALSE)
  payload_archive <- readRDS(file.path(data_dir, files$data[[6L]]))
  content_manifest <- utils::read.csv(file.path(data_dir, files$data[[7L]]), stringsAsFactors = FALSE)
  expected_attempt_names <- c(names(jss_joint_contract_fixture(jss_joint_case_definitions()[1L, , drop = FALSE],
    1L, list(seed = 1L))), "result_portable_sha256", "profile")
  if (!setequal(names(attempts), expected_attempt_names)) {
    stop("Optimizer public attempt CSV has extra or missing columns.", call. = FALSE)
  }
  attempts <- jss_phase2_optimizer_normalize_attempt_csv(attempts)
  logical_attempt_fields <- c("success", "converged", "retained", "hit_outer_limit")
  integer_attempt_fields <- c("joint_review_rep", "paired_seed", "outer_iterations", "heldout_variogram_nsim")
  if (!all(vapply(attempts[logical_attempt_fields], is.logical, logical(1))) ||
      !all(vapply(attempts[integer_attempt_fields], is.integer, logical(1)))) {
    stop("Optimizer public attempt CSV violates exact logical/integer types.", call. = FALSE)
  }
  expected_status_names <- c(
    "task_id", "case_id", "joint_review_rep", "paired_seed", "checkpoint",
    "checkpoint_schema_version", "result_content_sha256", "result_portable_sha256",
    "producer_fingerprint", "producer_fingerprint_algorithm", "package_version",
    "package_fingerprint_scope", "package_source_file_count", "package_checkout_path",
    "package_source_sha256", "verified_package_path", "verified_package_version",
    "verified_source_sha256", "verified_producer_path", "verified_producer_sha256",
    "package_identity_verified", "verified_rlibs_user", "verified_libpaths",
    "package_load_strategy", "resumed", "executed", "workers_requested", "workers_used",
    "execution_mode", "worker_pid", "checkpoint_timestamp_utc", "execution_host",
    "execution_os", "execution_platform", "execution_r_version", "execution_rng_kind", "execution_blas",
    "execution_lapack", "adaptive_round", "registered_initial_attempts",
    "registered_top_up_batch", "registered_hard_cap", "actual_attempts_in_cell")
  if (!setequal(names(checkpoints), expected_status_names)) {
    stop("Optimizer checkpoint-status CSV has extra or missing columns.", call. = FALSE)
  }
  checkpoint_integer_fields <- c("task_id", "joint_review_rep", "paired_seed",
    "checkpoint_schema_version", "package_source_file_count", "workers_requested",
    "workers_used", "worker_pid", "adaptive_round", "registered_initial_attempts",
    "registered_top_up_batch", "registered_hard_cap", "actual_attempts_in_cell")
  checkpoint_logical_fields <- c("package_identity_verified", "resumed", "executed")
  checkpoint_character_fields <- setdiff(expected_status_names,
    c(checkpoint_integer_fields, checkpoint_logical_fields))
  if (!all(vapply(checkpoints[checkpoint_integer_fields], is.integer, logical(1))) ||
      !all(vapply(checkpoints[checkpoint_logical_fields], is.logical, logical(1))) ||
      !all(vapply(checkpoints[checkpoint_character_fields], is.character, logical(1)))) {
    stop("Optimizer checkpoint-status CSV violates exact declared types.", call. = FALSE)
  }

  expected_precision <- jss_joint_mc_precision_registry("full")
  expected_design <- jss_joint_case_definitions()
  for (nm in names(expected_precision)) expected_design[[nm]] <- expected_precision[[nm]][[1L]]
  jss_phase2_compare_frame(design, expected_design, "case_id", "Optimizer candidate design")
  jss_phase2_compare_frame(precision[names(expected_precision)], expected_precision,
    character(), "Optimizer precision registry")
  if (!identical(as.character(design$case_id), sprintf("JVS%02d", 1:4)) ||
      nrow(design) != 4L ||
      !isTRUE(as.logical(precision$all_cells_precision_met[[1L]]))) {
    stop("Optimizer benchmark does not match the registered four-scenario retained-pair precision design.", call. = FALSE)
  }
  required_attempt_fields <- c("case_id", "joint_review_rep", "method", "success", "converged", "retained", "paired_seed")
  attempts_per_cell <- table(attempts$case_id, attempts$method)
  initial_attempts <- as.integer(precision$initial_attempts[[1L]])
  top_up_batch <- as.integer(precision$top_up_batch[[1L]])
  max_attempts <- as.integer(precision$max_attempts[[1L]])
  if (!all(required_attempt_fields %in% names(attempts)) ||
      !setequal(unique(attempts$case_id), design$case_id) ||
      !setequal(unique(attempts$method), c("rs_separate", "rs_joint")) ||
      any(attempts_per_cell[, "rs_joint"] != attempts_per_cell[, "rs_separate"]) ||
      any(attempts_per_cell < initial_attempts | attempts_per_cell > max_attempts) ||
      any(!(attempts_per_cell == max_attempts |
        (attempts_per_cell - initial_attempts) %% top_up_batch == 0L))) {
    stop("Optimizer attempt rows do not reconcile to registered adaptive batches.", call. = FALSE)
  }
  key <- paste(attempts$case_id, attempts$joint_review_rep, attempts$method, sep = "\r")
  if (anyDuplicated(key)) stop("Optimizer attempt keys are duplicated.", call. = FALSE)
  if (any(as.logical(attempts$retained) != (as.logical(attempts$success) & as.logical(attempts$converged)))) {
    stop("Optimizer retained flags do not equal success and convergence.", call. = FALSE)
  }
  provenance_fields <- c("checkpoint_schema_version", "result_content_sha256", "result_portable_sha256",
    "producer_fingerprint", "package_source_sha256", "verified_source_sha256",
    "verified_producer_path", "verified_producer_sha256", "package_identity_verified",
    "package_version", "package_fingerprint_scope", "package_source_file_count",
    "checkpoint_timestamp_utc", "execution_host", "execution_os", "execution_r_version",
    "execution_platform", "execution_rng_kind", "execution_blas", "execution_lapack")
  expected_checkpoint_rows <- sum(attempts_per_cell[, "rs_joint"])
  expected_package <- jss_joint_checkout_package_identity(list(root = root))
  expected_producer <- jss_joint_producer_fingerprint(list(root = root))
  observed_runtime <- jss_joint_runtime_identity()
  observed_rlibs_sha256 <- digest::digest(
    gsub("\\\\", "/", Sys.getenv("R_LIBS_USER", unset = "")), "sha256", serialize = FALSE)
  observed_libpaths_sha256 <- digest::digest(
    paste(gsub("\\\\", "/", .libPaths()), collapse = ";"), "sha256", serialize = FALSE)
  hex <- "^[0-9a-f]{64}$"
  portable_path_fields <- intersect(c("checkpoint", "package_checkout_path", "verified_package_path",
    "verified_producer_path"), names(checkpoints))
  canonical_status_order <- order(match(checkpoints$case_id, design$case_id), checkpoints$joint_review_rep)
  status_seed_check <- merge(
    unique(attempts[c("case_id", "joint_review_rep", "paired_seed")]),
    checkpoints[c("case_id", "joint_review_rep", "paired_seed")],
    by = c("case_id", "joint_review_rep"), suffixes = c("_attempt", "_status"), all = TRUE)
  if (nrow(checkpoints) != expected_checkpoint_rows ||
      anyDuplicated(checkpoints$checkpoint) || !identical(as.integer(checkpoints$task_id), seq_len(nrow(checkpoints))) ||
      !identical(canonical_status_order, seq_len(nrow(checkpoints))) ||
      nrow(status_seed_check) != nrow(checkpoints) ||
      any(status_seed_check$paired_seed_attempt != status_seed_check$paired_seed_status) ||
      any(checkpoints$resumed == checkpoints$executed) ||
      any(!checkpoints$resumed & checkpoints$workers_used < 1L) ||
      any(checkpoints$resumed & checkpoints$workers_used != 0L) ||
      any(checkpoints$workers_used > checkpoints$workers_requested) ||
      any(!checkpoints$execution_mode %in% c("serial", "psock", "resume_only")) ||
      any(checkpoints$execution_mode == "resume_only" & !checkpoints$resumed) ||
      any(checkpoints$registered_initial_attempts != initial_attempts) ||
      any(checkpoints$registered_top_up_batch != top_up_batch) ||
      any(checkpoints$registered_hard_cap != max_attempts) ||
      any(checkpoints$actual_attempts_in_cell != attempts_per_cell[cbind(checkpoints$case_id, rep("rs_joint", nrow(checkpoints)))]) ||
      any(vapply(checkpoints[portable_path_fields], function(x)
        any(grepl("^([A-Za-z]:[/\\\\]|/)", as.character(x))), logical(1))) ||
      !all(provenance_fields %in% names(checkpoints)) || any(checkpoints$checkpoint_schema_version != 6L) ||
      any(!as.logical(checkpoints$package_identity_verified)) ||
      any(!grepl(hex, checkpoints$result_content_sha256)) ||
      any(!grepl(hex, checkpoints$result_portable_sha256)) ||
      any(!grepl(hex, checkpoints$package_source_sha256)) ||
      any(!grepl(hex, checkpoints$producer_fingerprint)) ||
      any(checkpoints$package_source_sha256 != expected_package$source_sha256) ||
      any(checkpoints$package_version != expected_package$version) ||
      any(checkpoints$package_fingerprint_scope != expected_package$fingerprint_scope) ||
      any(checkpoints$package_source_file_count != expected_package$source_file_count) ||
      any(checkpoints$producer_fingerprint != expected_producer) ||
      any(checkpoints$package_checkout_path != "checkout") ||
      any(checkpoints$verified_package_path != "checkout") ||
      any(checkpoints$verified_producer_path != "checkout/paper/R/03-joint-vs-separate-optimization.R") ||
      any(checkpoints$verified_rlibs_user != observed_rlibs_sha256) ||
      any(checkpoints$verified_libpaths != observed_libpaths_sha256) ||
      any(checkpoints$package_load_strategy != "pkgload_checkout") ||
      any(checkpoints$package_source_sha256 != checkpoints$verified_source_sha256) ||
      any(checkpoints$producer_fingerprint != checkpoints$verified_producer_sha256) ||
      any(vapply(checkpoints[c("checkpoint_timestamp_utc", "execution_host", "execution_os",
        "execution_platform", "execution_r_version", "execution_rng_kind", "execution_blas", "execution_lapack")],
        function(x) any(is.na(x) | !nzchar(as.character(x))), logical(1)))) {
    stop("Optimizer checkpoint audit lacks portable paths, hashes, or verified runtime provenance.", call. = FALSE)
  }
  timestamp_ok <- grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$",
    checkpoints$checkpoint_timestamp_utc)
  checkpoint_name_ok <- grepl("checkpoints/03-joint-vs-separate-optimization/paired-one-factor-v6/.+/JVS[0-9]{2}-rep-[0-9]+[.]rds$",
    gsub("\\\\", "/", checkpoints$checkpoint))
  checkpoint_embedded_key_ok <- endsWith(gsub("\\\\", "/", checkpoints$checkpoint),
    sprintf("/%s-rep-%d.rds", checkpoints$case_id, checkpoints$joint_review_rep))
  if (!all(timestamp_ok) || !all(checkpoint_name_ok) ||
      !all(checkpoint_embedded_key_ok) ||
      any(checkpoints$execution_host != observed_runtime$host) ||
      any(checkpoints$execution_os != observed_runtime$os) ||
      any(checkpoints$execution_platform != observed_runtime$platform) ||
      any(checkpoints$execution_r_version != observed_runtime$r_version) ||
      any(checkpoints$execution_rng_kind != observed_runtime$rng_kind) ||
      any(checkpoints$execution_blas != observed_runtime$blas) ||
      any(checkpoints$execution_lapack != observed_runtime$lapack) ||
      any(grepl("^([A-Za-z]:[/\\\\]|/)", checkpoints$execution_blas)) ||
      any(grepl("^([A-Za-z]:[/\\\\]|/)", checkpoints$execution_lapack)) ||
      any(!is.finite(checkpoints$worker_pid) | checkpoints$worker_pid < 1 |
        checkpoints$worker_pid != floor(checkpoints$worker_pid))) {
    stop("Optimizer runtime provenance formats or canonical checkpoint paths are invalid.", call. = FALSE)
  }
  jss_phase2_parse_rfc3339_utc(checkpoints$checkpoint_timestamp_utc,
    "Optimizer checkpoint timestamp")
  result_hash <- unique(attempts[c("case_id", "joint_review_rep", "result_portable_sha256")])
  status_hash <- checkpoints[c("case_id", "joint_review_rep", "result_portable_sha256")]
  hashes <- merge(result_hash, status_hash, by = c("case_id", "joint_review_rep"), suffixes = c("_results", "_status"), all = TRUE)
  if (nrow(hashes) != nrow(checkpoints) || any(hashes$result_portable_sha256_results != hashes$result_portable_sha256_status)) {
    stop("Optimizer results do not reconcile to checkpoint content hashes.", call. = FALSE)
  }
  if (!exists("jss_joint_portable_result_sha256", mode = "function")) {
    stop("Optimizer portable hash validator is unavailable.", call. = FALSE)
  }
  grouped_attempts <- split(attempts, interaction(attempts$case_id, attempts$joint_review_rep, drop = TRUE))
  recalculated <- vapply(grouped_attempts, jss_joint_portable_result_sha256, character(1))
  recorded <- vapply(grouped_attempts, function(x) unique(as.character(x$result_portable_sha256))[[1L]], character(1))
  if (!identical(unname(recalculated), unname(recorded))) {
    stop("Optimizer result content hashes do not validate against the public attempt rows.", call. = FALSE)
  }
  if (!is.list(payload_archive) || length(payload_archive) != nrow(checkpoints) ||
      !setequal(names(content_manifest), c("case_id", "joint_review_rep", "result_content_sha256",
        "result_portable_sha256")) || nrow(content_manifest) != nrow(checkpoints)) {
    stop("Optimizer durable checkpoint payload archive is incomplete.", call. = FALSE)
  }
  if (!is.character(content_manifest$case_id) || !is.integer(content_manifest$joint_review_rep) ||
      !is.character(content_manifest$result_content_sha256) ||
      !is.character(content_manifest$result_portable_sha256) ||
      any(!grepl("^[0-9a-f]{64}$", content_manifest$result_content_sha256)) ||
      any(!grepl("^[0-9a-f]{64}$", content_manifest$result_portable_sha256))) {
    stop("Optimizer checkpoint content manifest violates exact declared types or hashes.", call. = FALSE)
  }
  archive_check <- do.call(rbind, lapply(payload_archive, function(x) data.frame(
    case_id = x$case_id[[1L]], joint_review_rep = x$joint_review_rep[[1L]],
    result_content_sha256 = jss_joint_content_sha256(x),
    result_portable_sha256 = jss_joint_portable_result_sha256(x), stringsAsFactors = FALSE)))
  jss_phase2_compare_frame(content_manifest, archive_check,
    c("case_id", "joint_review_rep"), "Optimizer checkpoint content manifest")
  status_manifest <- checkpoints[c("case_id", "joint_review_rep", "result_content_sha256",
    "result_portable_sha256")]
  jss_phase2_compare_frame(status_manifest, archive_check,
    c("case_id", "joint_review_rep"), "Optimizer checkpoint attestations")
  payload_checks <- lapply(grouped_attempts, function(x) {
    case <- expected_design[expected_design$case_id == x$case_id[[1L]], , drop = FALSE]
    expected_spec <- list(paired_seed = as.integer(unique(x$paired_seed)[[1L]]))
    checkpoint_x <- x[, setdiff(names(x), c("result_portable_sha256", "profile")), drop = FALSE]
    jss_joint_validate_result_payload(checkpoint_x, expected_spec, case,
      as.integer(x$joint_review_rep[[1L]]))
  })
  payload_valid <- vapply(payload_checks, function(x) isTRUE(x$valid), logical(1))
  inferred_seed <- unique(as.integer(attempts$paired_seed) - 3000L -
    as.integer(attempts$joint_review_rep) * 100L)
  sequential_reps <- vapply(split(attempts$joint_review_rep, attempts$case_id), function(x)
    identical(sort(unique(as.integer(x))), seq_len(length(unique(x)))), logical(1))
  design_fields <- intersect(c("base_case_id", "contrast_factor", "contrast_label", "contrast_level",
    "hypothesis_role", "family", "copula", "n", "time_points", "total_observations",
    "mu_strength", "sigma_strength", "theta_strength", "time_shape"), names(attempts))
  metadata_ok <- all(vapply(seq_len(nrow(design)), function(i) {
    x <- attempts[attempts$case_id == design$case_id[[i]], , drop = FALSE]
    all(vapply(design_fields, function(nm)
      all(as.character(x[[nm]]) == as.character(design[[nm]][[i]])), logical(1)))
  }, logical(1)))
  variation_metrics <- intersect(c("train_joint_loglik", "test_log_score_per_obs",
    "heldout_variogram_score_p05", "train_rmse_theta"), names(attempts))
  informative_metrics <- all(vapply(variation_metrics, function(nm) {
    x <- attempts[[nm]][as.logical(attempts$retained) & is.finite(attempts[[nm]])]
    length(unique(signif(x, 12L))) > 1L
  }, logical(1)))
  if (!all(payload_valid) || length(inferred_seed) != 1L || !all(sequential_reps) ||
      !metadata_ok || length(variation_metrics) < 4L || !informative_metrics) {
    invalid_reasons <- unique(unlist(lapply(payload_checks[!payload_valid], `[[`, "reasons")))
    stop("Optimizer attempt payloads violate the exact typed/status/metric contract",
      if (length(invalid_reasons)) paste0(": ", paste(invalid_reasons, collapse = "; ")) else ".",
      call. = FALSE)
  }
  deltas <- jss_joint_delta_table(attempts)
  achieved <- jss_joint_precision_achievement(deltas, precision)
  adaptive_contract <- vapply(achieved$case_id, function(case_id) {
    rows <- checkpoints[checkpoints$case_id == case_id, , drop = FALSE]
    rows <- rows[order(rows$joint_review_rep), , drop = FALSE]
    rounds <- sort(unique(rows$adaptive_round))
    if (!identical(rounds, seq_len(max(rounds))) || rounds[[1L]] != 1L) return(FALSE)
    for (round in rounds) {
      members <- rows$joint_review_rep[rows$adaptive_round == round]
      expected_start <- if (round == 1L) 1L else max(rows$joint_review_rep[rows$adaptive_round < round]) + 1L
      expected_size <- if (round == 1L) initial_attempts else
        min(top_up_batch, max_attempts - expected_start + 1L)
      if (!identical(members, seq.int(expected_start, length.out = expected_size))) return(FALSE)
      if (round > 1L) {
        prior <- deltas[deltas$case_id == case_id & deltas$joint_review_rep < expected_start, , drop = FALSE]
        if (isTRUE(jss_joint_precision_achievement(prior, precision)$precision_met[[1L]])) return(FALSE)
      }
    }
    TRUE
  }, logical(1))
  if (nrow(achieved) != nrow(design) || any(!achieved$precision_met) ||
      !all(adaptive_contract) ||
      min(achieved$retained_pairs) != precision$achieved_min_retained_pairs[[1L]] ||
      !isTRUE(all.equal(max(achieved$achieved_worst_case_half_width),
        precision$achieved_max_half_width[[1L]], tolerance = 1e-12))) {
    stop("Optimizer production precision eligibility does not reconstruct from attempt rows.", call. = FALSE)
  }
  canonical_failure <- jss_joint_failure_summary(attempts)
  canonical_uncertainty <- jss_joint_difference_uncertainty(deltas)
  canonical_summary <- jss_joint_summary_table(attempts, deltas, design)
  canonical_wins <- jss_joint_metric_wins(attempts)
  canonical_hypothesis <- jss_joint_hypothesis_summary(deltas)
  canonical_evidence <- jss_joint_hypothesis_evidence(canonical_hypothesis, canonical_uncertainty)
  jss_phase2_compare_frame(failures, canonical_failure,
    c("case_id", "method", "failure_reason"), "Optimizer failure summary")
  jss_phase2_compare_frame(uncertainty, canonical_uncertainty,
    c("case_id", "metric"), "Optimizer uncertainty summary")
  jss_phase2_compare_frame(summary_reported, canonical_summary,
    intersect(c("case_id", "method"), names(canonical_summary)), "Optimizer summary")
  jss_phase2_compare_frame(wins_reported, canonical_wins,
    intersect(c("case_id", "metric", "method"), names(canonical_wins)), "Optimizer metric-win table")
  jss_phase2_compare_frame(hypothesis_reported, canonical_hypothesis,
    intersect(c("hypothesis", "focal_case", "metric"), names(canonical_hypothesis)), "Optimizer hypothesis summary")
  jss_phase2_compare_frame(evidence_reported, canonical_evidence,
    intersect(c("hypothesis", "focal_case", "metric"), names(canonical_evidence)), "Optimizer hypothesis evidence")
  uncertainty_fields <- c(
    "conditional_difference", "difference_mcse", "difference_conf_low",
    "difference_conf_high", "retained_pairs", "metric_finite_pairs", "metric_failed_pairs",
    "retention_difference_mcse", "retention_difference_conf_low", "retention_difference_conf_high",
    "convergence_difference_mcse", "convergence_difference_conf_low", "convergence_difference_conf_high",
    "failure_inclusive_difference", "failure_inclusive_mcse", "failure_inclusive_conf_low",
    "failure_inclusive_conf_high"
  )
  if (!nrow(uncertainty) || !all(uncertainty_fields %in% names(uncertainty)) ||
      any(!is.finite(uncertainty$difference_mcse[uncertainty$retained_pairs > 1L])) ||
      any(uncertainty$metric_finite_pairs < 0 | uncertainty$metric_failed_pairs < 0 |
        uncertainty$metric_finite_pairs + uncertainty$metric_failed_pairs != uncertainty$attempted_pairs)) {
    stop("Optimizer differences lack complete Monte Carlo uncertainty.", call. = FALSE)
  }
  interval_sets <- list(
    c("conditional_difference", "difference_mcse", "difference_conf_low", "difference_conf_high"),
    c("failure_inclusive_difference", "failure_inclusive_mcse", "failure_inclusive_conf_low", "failure_inclusive_conf_high"),
    c("retention_difference", "retention_difference_mcse", "retention_difference_conf_low", "retention_difference_conf_high"),
    c("convergence_difference", "convergence_difference_mcse", "convergence_difference_conf_low", "convergence_difference_conf_high")
  )
  if (any(vapply(interval_sets, function(fields) {
    est <- uncertainty[[fields[[1L]]]]; se <- uncertainty[[fields[[2L]]]]
    lo <- uncertainty[[fields[[3L]]]]; hi <- uncertainty[[fields[[4L]]]]
    complete <- is.finite(est) | is.finite(se) | is.finite(lo) | is.finite(hi)
    any(complete & (!is.finite(est) | !is.finite(se) | se < 0 |
      !is.finite(lo) | !is.finite(hi) | lo > est | est > hi))
  }, logical(1)))) stop("Optimizer uncertainty contains an impossible interval.", call. = FALSE)
  failure_totals <- stats::aggregate(attempts ~ case_id + method, failures, sum)
  if (!nrow(failures) || !all(c("case_id", "method", "failure_reason", "attempts", "cell_attempts") %in% names(failures)) ||
      any(failures$cell_attempts != attempts_per_cell[cbind(failures$case_id, failures$method)]) ||
      any(failure_totals$attempts != attempts_per_cell[cbind(failure_totals$case_id, failure_totals$method)])) {
    stop("Optimizer failure denominators are absent.", call. = FALSE)
  }
  figure_paths <- file.path(figures_dir, files$figures)
  expected_figure_registry <- data.frame(
    figure = basename(figure_paths),
    png_sha256 = vapply(figure_paths, jss_joint_sha256_file, character(1)),
    plotted_data_sha256 = c(jss_joint_portable_frame_sha256(deltas),
      jss_joint_portable_frame_sha256(canonical_wins)),
    plot_spec_sha256 = c(
      digest::digest("delta density+zero line+case facet; width=10 height=7 dpi=320", "sha256", serialize = FALSE),
      digest::digest("metric win-or-tie bars; width=9 height=6 dpi=320", "sha256", serialize = FALSE)
    ), stringsAsFactors = FALSE)
  jss_phase2_compare_frame(figure_registry, expected_figure_registry, "figure",
    "Optimizer figure registry")
  if (!all(vapply(figure_paths, jss_phase2_validate_png, logical(1)))) {
    stop("Optimizer figures are not decodable nonblank PNG files.", call. = FALSE)
  }
  rerender_dir <- tempfile("optimizer-rerender-"); dir.create(rerender_dir)
  on.exit(unlink(rerender_dir, recursive = TRUE, force = TRUE), add = TRUE)
  rerendered <- file.path(rerender_dir, basename(figure_paths))
  jss_joint_write_delta_figure(deltas, rerendered[[1L]])
  jss_joint_write_metric_dashboard(canonical_wins, rerendered[[2L]])
  if (!identical(unname(vapply(rerendered, jss_joint_sha256_file, character(1))),
      unname(vapply(figure_paths, jss_joint_sha256_file, character(1))))) {
    stop("Optimizer figures do not match deterministic rerendering from registered data/specification.", call. = FALSE)
  }
  if (isTRUE(require_promotion)) {
    jss_phase2_require_external_attestation("optimizer-benchmark", bundle_sha256,
      expected_package$source_sha256, expected_producer, checkpoints, root,
      attestation_path, signature_path)
  }
  final_manifest <- paste(sort(paste(basename(required),
    vapply(required, jss_joint_sha256_file, character(1)), sep = "\t")), collapse = "\n")
  if (!identical(digest::digest(final_manifest, "sha256", serialize = FALSE), bundle_sha256)) {
    stop("Optimizer staged snapshot changed during validation.", call. = FALSE)
  }
  final_package <- jss_joint_checkout_package_identity(list(root = root))
  final_producer <- jss_joint_producer_fingerprint(list(root = root))
  if (!identical(final_package$source_sha256, expected_package$source_sha256) ||
      !identical(final_producer, expected_producer)) {
    stop("Optimizer checkout/producer identity changed during validation.", call. = FALSE)
  }
  invisible(required)
}

jss_phase2_optimizer_validate <- function(data_dir, tables_dir, figures_dir, root = getwd(),
    attestation_path = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_OPTIMIZER_ATTESTATION", unset = ""),
    signature_path = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_OPTIMIZER_ATTESTATION_SIGNATURE", unset = "")) {
  files <- jss_phase2_optimizer_files()
  original <- c(file.path(data_dir, files$data), file.path(tables_dir, files$tables),
    file.path(figures_dir, files$figures))
  before <- vapply(original, jss_joint_sha256_file, character(1))
  stage <- tempfile("optimizer-immutable-snapshot-"); dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
  staged_dirs <- file.path(stage, c("data", "tables", "figures"))
  invisible(lapply(staged_dirs, dir.create, recursive = TRUE))
  staged <- c(file.path(staged_dirs[[1L]], files$data), file.path(staged_dirs[[2L]], files$tables),
    file.path(staged_dirs[[3L]], files$figures))
  if (!all(file.copy(original, staged, overwrite = FALSE)) ||
      !identical(before, vapply(original, jss_joint_sha256_file, character(1))) ||
      !identical(unname(before), unname(vapply(staged, jss_joint_sha256_file, character(1))))) {
    stop("Optimizer source bundle changed while creating its immutable validation snapshot.", call. = FALSE)
  }
  jss_phase2_optimizer_validate_candidate(staged_dirs[[1L]], staged_dirs[[2L]], staged_dirs[[3L]], root,
    attestation_path, signature_path, require_promotion = TRUE)
}

jss_phase2_optimizer_relativize_checkpoints <- function(settings) {
  path <- jss_joint_output_paths(settings)$checkpoint_status
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  root <- paste0(normalizePath(settings$out_dir, winslash = "/", mustWork = TRUE), "/")
  path_fields <- intersect(c("checkpoint", "package_checkout_path", "verified_package_path",
    "verified_producer_path"), names(x))
  for (field in path_fields) {
    normalized <- gsub("\\\\", "/", as.character(x[[field]]))
    x[[field]] <- ifelse(startsWith(tolower(normalized), tolower(root)),
      substring(normalized, nchar(root) + 1L), normalized)
    root_checkout <- paste0(normalizePath(settings$root, winslash = "/", mustWork = TRUE), "/")
    x[[field]] <- ifelse(startsWith(tolower(x[[field]]), tolower(root_checkout)),
      paste0("checkout/", substring(x[[field]], nchar(root_checkout) + 1L)), x[[field]])
    x[[field]] <- ifelse(tolower(x[[field]]) == tolower(sub("/$", "", root_checkout)),
      "checkout", x[[field]])
  }
  for (field in intersect(c("verified_rlibs_user", "verified_libpaths"), names(x))) {
    x[[field]] <- vapply(as.character(x[[field]]), function(value)
      digest::digest(gsub("\\\\", "/", value), "sha256", serialize = FALSE), character(1))
  }
  if (any(vapply(x[path_fields], function(value)
      any(grepl("^([A-Za-z]:[/\\\\]|/)", value)), logical(1)))) {
    stop("Optimizer checkpoint status contains paths outside the replication output root.", call. = FALSE)
  }
  utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}

jss_run_phase2_optimizer <- function(settings) {
  files <- jss_phase2_optimizer_files()
  if (identical(settings$profile, "paper")) {
    tracked <- file.path(settings$public_data_dir, "optimizer-benchmark")
    source_data <- file.path(tracked, "data", files$data)
    source_tables <- file.path(tracked, "tables", files$tables)
    source_figures <- file.path(tracked, "figures", files$figures)
    data <- jss_phase2_copy_exact(source_data, file.path(settings$data_dir, files$data))
    tables <- jss_phase2_copy_exact(source_tables, file.path(settings$tables_dir, files$tables))
    figures <- jss_phase2_copy_exact(source_figures, file.path(settings$figures_dir, files$figures))
  } else if (identical(settings$profile, "full")) {
    result <- jss_run_03_joint_vs_separate(settings)
    jss_phase2_optimizer_relativize_checkpoints(settings)
    data <- result$data
    tables <- result$tables
    figures <- result$figures
  } else {
    return(list(module_id = "03-joint-vs-separate-optimization", status = "not_run_in_smoke", data = character(), tables = character(), figures = character()))
  }
  jss_phase2_optimizer_validate(settings$data_dir, settings$tables_dir, settings$figures_dir)
  list(
    module_id = "03-joint-vs-separate-optimization", status = "current",
    data = data, tables = tables, figures = figures,
    notes = "Registered paired RS-separate versus RS-joint benchmark; CG excluded."
  )
}

jss_phase2_gate_row <- function(check_id, pass, detail, source = NA_character_) {
  data.frame(check_id = check_id, pass = isTRUE(pass), detail = as.character(detail), source = as.character(source), stringsAsFactors = FALSE)
}

jss_phase2_claim_source_paths <- function(settings) {
  c(
    `main-recovery` = file.path(settings$data_dir, "main-recovery", "attempt_metadata.csv"),
    `optimizer-benchmark` = file.path(settings$data_dir, "03-joint-vs-separate-optimization-results.csv"),
    missingness = file.path(settings$data_dir, "missingness", "fit_run_log.csv"),
    `copula-misspecification` = file.path(settings$data_dir, "07-gamma-copula-misspecification-results.csv"),
    `multivariate-benchmark` = file.path(settings$data_dir, "multivariate-benchmark", "fit_status_by_rep.csv"),
    `fit-scaling` = file.path(settings$data_dir, "10-fit-scaling-attempts.csv")
  )
}

jss_phase2_claim_lifecycle <- function(x, study_id) {
  choose <- function(candidates, default = NULL) {
    hit <- candidates[candidates %in% names(x)]
    if (length(hit)) x[[hit[[1L]]]] else default
  }
  attempted <- choose(c("attempted", "execution_attempted"), rep(TRUE, nrow(x))) %in% TRUE
  converged <- choose(c("converged", "convergence_eligible", "fit_converged"), rep(FALSE, nrow(x))) %in% TRUE
  retained <- choose(c("retained", "publication_retained", "success"), converged) %in% TRUE
  retained <- retained & converged
  if (!all(attempted) || any(retained & !converged)) stop("Claim source lifecycle is contradictory: ", study_id, call. = FALSE)
  failed <- !retained
  reason <- as.character(choose(
    c("failure_reason", "failure_reason_short", "failure_type", "stop_reason", "error"),
    ifelse(failed, "unclassified_failure", "none")
  ))
  reason[is.na(reason) | !nzchar(trimws(reason))] <- ifelse(failed[is.na(reason) | !nzchar(trimws(reason))], "unclassified_failure", "none")
  if (any(failed & reason == "none")) stop("Claim source hides one or more failures: ", study_id, call. = FALSE)
  reasons <- sort(table(reason[failed]), decreasing = TRUE)
  list(
    attempted = sum(attempted), converged = sum(converged), retained = sum(retained),
    failed = sum(failed),
    failure_reasons = if (length(reasons)) paste(names(reasons), as.integer(reasons), sep = ":", collapse = ";") else "none"
  )
}

jss_phase2_optimizer_claim_estimates <- function(attempts, precision) {
  required <- c("case_id", "joint_review_rep", "method", "success", "converged")
  if (!all(required %in% names(attempts)) || !nrow(attempts) ||
      is.null(precision$initial_attempts) || length(precision$initial_attempts) != 1L) {
    stop("Optimizer claim evidence lacks its registered attempt/precision fields.", call. = FALSE)
  }
  pair_key <- unique(attempts[c("case_id", "joint_review_rep")])
  pair_eligible <- stats::aggregate(
    as.logical(attempts$success) & as.logical(attempts$converged),
    attempts[c("case_id", "joint_review_rep")], all
  )
  retained_by_case <- stats::aggregate(pair_eligible$x, pair_eligible["case_id"], sum)
  c(
    registered_scenario_cells = nrow(unique(attempts["case_id"])),
    registered_method_cells = nrow(unique(attempts[c("case_id", "method")])),
    planned_initial_pairs = nrow(unique(attempts["case_id"])) * as.integer(precision$initial_attempts[[1L]]),
    actual_attempt_rows = nrow(attempts),
    achieved_min_retained_pairs = min(retained_by_case$x)
  )
}

jss_build_phase2_claim_evidence <- function(settings) {
  source_identities <- jss_phase2_validate_claim_sources(settings)
  claims_path <- file.path(settings$root, "paper", "phase2-claims.csv")
  claims <- utils::read.csv(claims_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(
    "claim_id", "tex_anchor", "study_id", "scenario_id", "claim_type", "metric",
    "expected_direction", "attempt_artifact_id", "effect_artifact_id",
    "source_row_key", "wording_strength", "expected_value"
  )
  if (!all(required %in% names(claims)) || anyDuplicated(claims$claim_id)) {
    stop("Phase 2 claim register is malformed or has duplicate IDs.", call. = FALSE)
  }
  manuscript <- paste(readLines(file.path(settings$root, "paper", "manuscript", "main.tex"), warn = FALSE), collapse = "\n")
  missing_anchor <- claims$tex_anchor[!vapply(
    claims$tex_anchor,
    function(anchor) grepl(paste0("\\label{", anchor, "}"), manuscript, fixed = TRUE),
    logical(1)
  )]
  if (length(missing_anchor)) stop("Claim register names missing TeX anchors: ", paste(unique(missing_anchor), collapse = ", "), call. = FALSE)

  manifest <- utils::read.csv(file.path(settings$root, "paper", "manifest.csv"), stringsAsFactors = FALSE)
  if (any(!claims$attempt_artifact_id %in% manifest$artifact_id)) {
    stop("Claim register references unregistered attempt artifacts.", call. = FALSE)
  }
  sources <- jss_phase2_claim_source_paths(settings)
  rows <- lapply(seq_len(nrow(claims)), function(i) {
    claim <- claims[i, , drop = FALSE]
    source <- unname(sources[[claim$study_id]])
    if (is.null(source) || !file.exists(source)) {
      stop("Claim source is missing for ", claim$claim_id, ": ", source, call. = FALSE)
    }
    x <- utils::read.csv(source, stringsAsFactors = FALSE)
    scenario <- as.character(claim$scenario_id)
    if (!identical(scenario, "all")) {
      if (!"scenario_id" %in% names(x)) stop("Claim source lacks scenario_id: ", claim$claim_id, call. = FALSE)
      x <- x[as.character(x$scenario_id) == scenario, , drop = FALSE]
    }
    row_key <- as.character(claim$source_row_key)
    if (!identical(row_key, "all")) {
      pieces <- strsplit(row_key, ";", fixed = TRUE)[[1L]]
      for (piece in pieces) {
        pair <- strsplit(piece, "=", fixed = TRUE)[[1L]]
        if (length(pair) != 2L || !pair[[1L]] %in% names(x)) {
          stop("Claim source_row_key is not an exact registered field=value selector: ", claim$claim_id, call. = FALSE)
        }
        x <- x[as.character(x[[pair[[1L]]]]) == pair[[2L]], , drop = FALSE]
      }
    }
    if (!nrow(x)) stop("Claim scenario/source_row_key selects no approved attempt rows: ", claim$claim_id, call. = FALSE)
    lifecycle <- jss_phase2_claim_lifecycle(x, as.character(claim$study_id))
    estimate <- if (identical(claim$study_id, "main-recovery")) {
      groups <- unique(x[c("study_id", "scenario_id", "method", "target_replicates")])
      sum(as.integer(groups$target_replicates))
    } else if (identical(claim$study_id, "optimizer-benchmark")) {
      precision <- utils::read.csv(file.path(settings$data_dir,
        "03-joint-vs-separate-optimization-mc-precision.csv"), stringsAsFactors = FALSE)
      optimizer_estimates <- jss_phase2_optimizer_claim_estimates(x, precision)
      metric <- as.character(claim$metric)
      if (!metric %in% names(optimizer_estimates)) {
        stop("Unregistered optimizer claim metric: ", metric, call. = FALSE)
      }
      unname(optimizer_estimates[[metric]])
    } else if (identical(claim$study_id, "fit-scaling") && identical(as.character(claim$metric), "scenario_cells")) {
      length(unique(as.character(x$scenario_id)))
    } else nrow(x)
    dynamic <- identical(as.character(claim$wording_strength), "dynamic")
    if (!dynamic && !isTRUE(all.equal(as.numeric(estimate), as.numeric(claim$expected_value)))) {
      stop("Claim count does not match its registered value: ", claim$claim_id, call. = FALSE)
    }
    data.frame(
      claim_id = claim$claim_id,
      tex_anchor = claim$tex_anchor,
      study_id = claim$study_id,
      scenario_id = claim$scenario_id,
      claim_type = claim$claim_type,
      metric = claim$metric,
      expected_direction = claim$expected_direction,
      wording_strength = claim$wording_strength,
      attempt_artifact_id = claim$attempt_artifact_id,
      effect_artifact_id = claim$effect_artifact_id,
      source_row_key = claim$source_row_key,
      source_file = substring(normalizePath(source, winslash = "/", mustWork = TRUE), nchar(normalizePath(settings$out_dir, winslash = "/", mustWork = TRUE)) + 2L),
      source_sha256 = jss_file_sha256(source),
      effect_sha256 = jss_file_sha256(source),
      estimate = as.numeric(estimate),
      mcse = 0,
      lower = as.numeric(estimate),
      upper = as.numeric(estimate),
      denominator = lifecycle$attempted,
      attempted = lifecycle$attempted,
      converged = lifecycle$converged,
      retained = lifecycle$retained,
      failed = lifecycle$failed,
      failure_reasons = lifecycle$failure_reasons,
      source_identity_sha256 = digest::digest(
        source_identities[[as.character(claim$study_id)]], algo = "sha256", serialize = TRUE
      ),
      status = if (dynamic) "verified_exact_evidence_count" else "verified_exact_design_count",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

jss_write_phase2_claim_evidence <- function(settings) {
  evidence <- jss_build_phase2_claim_evidence(settings)
  path <- file.path(settings$tables_dir, "phase2-claim-evidence.csv")
  utils::write.csv(evidence, path, row.names = FALSE)
  installed <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  registered <- utils::read.csv(file.path(settings$root, "paper", "phase2-claims.csv"), stringsAsFactors = FALSE)
  jss_phase2_validate_claim_evidence_rows(installed, registered, settings = settings)
  path
}

jss_write_phase2_gate_audit <- function(settings, modules, claim_evidence) {
  status <- vapply(modules, function(x) identical(x$status, "current") || identical(x$status, "regenerated"), logical(1))
  rows <- list(jss_phase2_gate_row(
    "all_phase2_modules_current", all(status),
    paste(names(status), status, sep = "=", collapse = ";"), "target module return values"
  ))
  production_error <- NULL
  production_ok <- tryCatch({
    jss_phase2_validate_production_modules(settings, modules)
    TRUE
  }, error = function(e) { production_error <<- conditionMessage(e); FALSE })
  rows[[2L]] <- jss_phase2_gate_row(
    "all_module_production_gates_pass", production_ok,
    if (production_ok) "optimizer, missingness, recovery, Module07, Module09, and scaling gates passed" else production_error,
    "registered module production validators"
  )
  manifest_error <- NULL
  manifest_ok <- tryCatch({
    jss_phase2_validate_manifest_allowlist(settings, modules, claim_evidence)
    TRUE
  }, error = function(e) { manifest_error <<- conditionMessage(e); FALSE })
  rows[[3L]] <- jss_phase2_gate_row(
    "phase2_manifest_exact_bidirectional", manifest_ok,
    if (manifest_ok) paste("exact_paths", length(jss_phase2_module_output_allowlist(settings, modules, claim_evidence)), sep = "=") else manifest_error,
    "paper/manifest.csv"
  )
  manuscript <- paste(readLines(file.path(settings$root, "paper", "manuscript", "main.tex"), warn = FALSE), collapse = " ")
  stale <- c(
    "run 1,000 simulations", "approximately 20x slower", "58/60 Clayton", "six example cases",
    "will generally fit reasonably", "can be captured very accurately",
    "tables/08-simulation-sensitivity-correlation-misspecification"
  )
  hits <- stale[vapply(stale, grepl, logical(1), x = manuscript, fixed = TRUE)]
  rows[[4L]] <- jss_phase2_gate_row(
    "no_registered_stale_claims", !length(hits),
    if (length(hits)) paste(hits, collapse = ";") else "no stale registered phrases", "paper/manuscript/main.tex"
  )
  claims <- utils::read.csv(claim_evidence, stringsAsFactors = FALSE)
  registered_claims <- utils::read.csv(file.path(settings$root, "paper", "phase2-claims.csv"), stringsAsFactors = FALSE)
  claim_error <- NULL
  claim_ok <- tryCatch({
    jss_phase2_validate_claim_evidence_rows(claims, registered_claims, settings = settings)
    TRUE
  }, error = function(e) { claim_error <<- conditionMessage(e); FALSE })
  rows[[5L]] <- jss_phase2_gate_row(
    "all_registered_claims_verified", nrow(claims) == nrow(registered_claims) && setequal(claims$claim_id, registered_claims$claim_id) &&
      claim_ok,
    if (claim_ok) paste("verified_claims", nrow(claims), sep = "=") else claim_error,
    "paper/phase2-claims.csv"
  )
  main_bundle <- file.path(settings$data_dir, "main-recovery")
  validation_path <- file.path(main_bundle, "evidence_validation.csv")
  validation <- if (file.exists(validation_path)) utils::read.csv(validation_path, stringsAsFactors = FALSE) else data.frame(pass = FALSE)
  rows[[6L]] <- jss_phase2_gate_row(
    "main_recovery_code_gate", nrow(validation) > 0L && all(as.logical(validation$pass)),
    "structural and recomputation contract checks; distinct from production completion", validation_path
  )
  recovery_attestation <- if (!is.null(settings$main_recovery_attestation)) settings$main_recovery_attestation else Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MAIN_RECOVERY_ATTESTATION", unset = "")
  recovery_signature <- if (!is.null(settings$main_recovery_attestation_signature)) settings$main_recovery_attestation_signature else Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MAIN_RECOVERY_ATTESTATION_SIGNATURE", unset = "")
  tracked_bundle <- file.path(settings$public_data_dir, "main-recovery")
  trust_error <- NULL
  trust_ok <- tryCatch({
    jss_main_recovery_validate_public_bundle(tracked_bundle, repo_root = settings$root,
      attestation_path = recovery_attestation, signature_path = recovery_signature, require_attestation = TRUE)
    TRUE
  }, error = function(e) { trust_error <<- conditionMessage(e); FALSE })
  rows[[7L]] <- jss_phase2_gate_row(
    "main_recovery_production_evidence_gate", trust_ok,
    if (trust_ok) "checkout-external pinned-key detached attestation validates the immutable bundle" else paste("external attestation gate failed", trust_error, sep = ": "), recovery_attestation
  )
  audit <- do.call(rbind, rows)
  if (!identical(as.character(audit$check_id), jss_phase2_central_gate_registry()) ||
      anyDuplicated(audit$check_id)) {
    stop("Phase 2 central audit registry is missing, reordered, duplicated, or extended.", call. = FALSE)
  }
  path <- file.path(settings$tables_dir, "phase2-gate-audit.csv")
  utils::write.csv(audit, path, row.names = FALSE)
  if (any(!audit$pass)) stop("Phase 2 central gate failed: ", paste(audit$check_id[!audit$pass], collapse = ", "), call. = FALSE)
  path
}
