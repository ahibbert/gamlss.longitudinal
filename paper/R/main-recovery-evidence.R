# Evidence helpers for the JSS-004 main recovery study.
#
# This file is sourced by the production evidence runner and focused tests; it
# is not package API. SHA-256 provenance requires the pinned digest package.

jss_recovery_schema_version <- function() "2026-09-02.10"

jss_recovery_phase1_contract_version <- function() {
  "likelihood-jss001-2026-09-01|inference-2026.1|capability-2026.2"
}

jss_recovery_gaussian_limit_zeta <- function() 7.5

jss_recovery_sha256 <- function(paths) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for main-recovery SHA-256 provenance.", call. = FALSE)
  }
  vapply(paths, digest::digest, character(1), file = TRUE, algo = "sha256", serialize = FALSE)
}

jss_recovery_package_source_files <- function(repo_root) {
  paths <- c(
    file.path(repo_root, "DESCRIPTION"),
    file.path(repo_root, "NAMESPACE"),
    sort(list.files(file.path(repo_root, "R"), pattern = "[.]R$", full.names = TRUE))
  )
  paths[file.exists(paths)]
}

jss_recovery_package_source_sha256 <- function(repo_root) {
  paths <- jss_recovery_package_source_files(repo_root)
  root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  normalized <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  relative <- substring(normalized, nchar(root) + 2L)
  individual <- unname(jss_recovery_sha256(normalized))
  digest::digest(paste(relative, individual, sep = "\t", collapse = "\n"), algo = "sha256", serialize = FALSE)
}

jss_recovery_git_identity <- function(repo_root) {
  sha <- tryCatch(system2("git", c("-C", shQuote(repo_root), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[[1]], error = function(e) NA_character_)
  state <- tryCatch(if (length(system2("git", c("-C", shQuote(repo_root), "status", "--porcelain"), stdout = TRUE, stderr = FALSE))) "dirty" else "clean", error = function(e) "unknown")
  if (length(sha) != 1L || !grepl("^[0-9a-f]{40}$", sha) || !state %in% c("clean", "dirty")) stop("Could not determine canonical Git identity.", call. = FALSE)
  tracked <- tryCatch(system2("git", c("-C", shQuote(repo_root), "ls-files"), stdout = TRUE, stderr = FALSE), error = function(e) character())
  if (!length(tracked)) stop("Could not enumerate the tracked worktree.", call. = FALSE)
  tracked <- sort(gsub("\\\\", "/", tracked))
  tracked_paths <- file.path(repo_root, tracked)
  hashes <- vapply(tracked_paths, function(path) if (file.exists(path)) unname(jss_recovery_sha256(path)) else "<deleted>", character(1))
  worktree_sha <- digest::digest(paste(tracked, hashes, sep = "\t", collapse = "\n"), algo = "sha256", serialize = FALSE)
  c(git_sha = sha, git_state = state, tracked_worktree_sha256 = worktree_sha)
}

jss_recovery_control_source_paths <- function(repo_root) {
  relative <- c(
    "paper/R/main-recovery-evidence.R",
    "paper/scripts/final-simulations/bcpe-t/simulation_bcpe_t_gamlss_comparison.R",
    "paper/scripts/final-simulations/nbi-clayton/compare_gamlss_ours_nbi_sigma_smooth.R",
    "paper/scripts/final-simulations/main-recovery/run_main_recovery_evidence.R",
    "paper/scripts/final-simulations/main-recovery/README.md",
    "paper/R/phase2-paper-evidence.R", "paper/R/phase2-central-integration.R",
    "paper/R/replication-helpers.R", "paper/_targets.R", "paper/manifest.csv",
    "paper/phase2-claims.csv", "paper/scripts/phase2-evidence-approval.R",
    "paper/notes/phase2-evidence-signing.md"
  )
  paths <- file.path(repo_root, relative)
  if (any(!file.exists(paths))) stop("A registered main-recovery control source is missing: ", paste(relative[!file.exists(paths)], collapse = ", "), call. = FALSE)
  stats::setNames(paths, relative)
}

jss_recovery_control_source_identity <- function(repo_root) {
  paths <- jss_recovery_control_source_paths(repo_root)
  data.frame(
    schema_version = jss_recovery_schema_version(), source_path = names(paths),
    sha256 = unname(jss_recovery_sha256(paths)), bytes = unname(file.info(paths)$size),
    stringsAsFactors = FALSE
  )
}

jss_recovery_producer_sha256 <- function(repo_root = normalizePath(".", winslash = "/", mustWork = TRUE)) {
  identity <- jss_recovery_control_source_identity(repo_root)
  digest::digest(paste(identity$source_path, identity$sha256, identity$bytes, sep = "\t", collapse = "\n"), algo = "sha256", serialize = FALSE)
}

jss_recovery_portable_source_identity <- function(x) {
  value <- gsub("\\\\", "/", trimws(as.character(x)))
  value <- sub("/+$", "", value)
  value[value %in% c("", ".")] <- "."
  absolute <- grepl("^/|^[A-Za-z]:/", value)
  outside <- grepl("(^|/)\\.\\.(/|$)", value)
  value[absolute | outside | value != "."] <- NA_character_
  value
}

jss_recovery_raw_convergence_fields <- function() c(
  "raw_convergence_schema", "raw_convergence_api", "raw_convergence_basis", "raw_convergence_status",
  "raw_convergence_indicator_name", "raw_convergence_indicator_value", "raw_convergence_loglik",
  "raw_convergence_deviance", "raw_convergence_coefficient_count", "raw_convergence_coefficient_nonfinite_count",
  "raw_convergence_fitted_count", "raw_convergence_fitted_nonfinite_count", "raw_convergence_iteration_count",
  "raw_convergence_iteration_cap", "raw_convergence_hit_outer_limit", "raw_convergence_hit_max_stall",
  "raw_convergence_hit_raw_loglik_deterioration"
)

jss_recovery_convergence_method_contract <- function(margin_family, method) {
  key <- paste(margin_family, method, sep = "/")
  switch(key,
    "BCPE/gamlss2" = list(api = "gamlss2::RS_result", basis = "registered_gamlss2_rs_v1:finite_logLik_and_deviance+finite_coefficients+finite_fitted_values+iterations_below_control_maxit", indicator = "not_exposed_by_gamlss2"),
    "BCPE/gamlss.longitudinal" = list(api = "gamlss.longitudinal::fit_convergence", basis = "gamlss.longitudinal_explicit_convergence_v1", indicator = "fit$convergence$converged"),
    "NBI/gamlss" = list(api = "gamlss::gamlss_result", basis = "gamlss_explicit_converged_v1", indicator = "fit$converged"),
    "NBI/ours_rs_joint" = list(api = "gamlss.longitudinal::fit_convergence", basis = "gamlss.longitudinal_explicit_convergence_v1", indicator = "fit$convergence$converged"),
    NULL
  )
}

jss_recovery_recompute_convergence <- function(ledger, reject_extra = TRUE) {
  fields <- jss_recovery_raw_convergence_fields()
  present <- names(ledger)[startsWith(names(ledger), "raw_convergence_")]
  if (!identical(sort(present), sort(fields))) {
    stop("Raw convergence evidence does not exactly match the registered schema; missing/extra fields: ",
      paste(c(setdiff(fields, present), setdiff(present, fields)), collapse = ", "), call. = FALSE)
  }
  if (any(as.character(ledger$raw_convergence_schema) != "raw-convergence-2026-09-01.1")) stop("Raw convergence evidence schema version mismatch.", call. = FALSE)
  logical_value <- function(x) jss_recovery_logical(x)
  result <- lapply(seq_len(nrow(ledger)), function(i) {
    row <- ledger[i, , drop = FALSE]
    contract <- jss_recovery_convergence_method_contract(as.character(row$margin_family), as.character(row$method))
    if (is.null(contract)) stop("No registered raw convergence contract for ", row$margin_family, "/", row$method, ".", call. = FALSE)
    if (!identical(as.character(row$raw_convergence_api), contract$api) ||
        !identical(as.character(row$raw_convergence_basis), contract$basis) ||
        !identical(as.character(row$raw_convergence_indicator_name), contract$indicator)) {
      stop("Raw convergence API/basis/indicator contradiction for attempt ", row$attempt_id, ".", call. = FALSE)
    }
    execution <- isTRUE(logical_value(row$execution_success))
    loglik <- suppressWarnings(as.numeric(row$raw_convergence_loglik)); deviance <- suppressWarnings(as.numeric(row$raw_convergence_deviance))
    coefficient_count <- suppressWarnings(as.integer(row$raw_convergence_coefficient_count)); coefficient_bad <- suppressWarnings(as.integer(row$raw_convergence_coefficient_nonfinite_count))
    fitted_count <- suppressWarnings(as.integer(row$raw_convergence_fitted_count)); fitted_bad <- suppressWarnings(as.integer(row$raw_convergence_fitted_nonfinite_count))
    iterations <- suppressWarnings(as.integer(row$raw_convergence_iteration_count)); cap <- suppressWarnings(as.integer(row$raw_convergence_iteration_cap))
    indicator <- logical_value(row$raw_convergence_indicator_value)
    hit_outer <- logical_value(row$raw_convergence_hit_outer_limit); hit_stall <- logical_value(row$raw_convergence_hit_max_stall)
    hit_drop <- logical_value(row$raw_convergence_hit_raw_loglik_deterioration)
    if (!execution) {
      status <- "fit_execution_failed"; converged <- FALSE
    } else {
      finite_payload <- is.finite(loglik) && is.finite(deviance) && is.finite(coefficient_count) && coefficient_count > 0L && coefficient_bad == 0L &&
        is.finite(fitted_count) && fitted_count > 0L && fitted_bad == 0L
      iteration_valid <- is.finite(iterations) && iterations >= 1L && is.finite(cap) && cap >= 1L
      if (!is.finite(loglik) || !is.finite(deviance)) status <- "nonfinite_objective"
      else if (!is.finite(coefficient_count) || coefficient_count < 1L || !is.finite(coefficient_bad) || coefficient_bad != 0L) status <- "nonfinite_or_missing_coefficients"
      else if (!is.finite(fitted_count) || fitted_count < 1L || !is.finite(fitted_bad) || fitted_bad != 0L) status <- "nonfinite_or_missing_fitted_values"
      else if (!iteration_valid) status <- "missing_or_invalid_iterations"
      else if (as.character(row$method) == "gamlss2" && iterations >= cap) status <- "outer_iteration_cap_reached_or_unverified"
      else if (as.character(row$method) != "gamlss2" && (isTRUE(hit_outer) || iterations >= cap)) status <- "outer_iteration_cap_reached"
      else if (as.character(row$method) != "gamlss2" && isTRUE(hit_stall)) status <- "maximum_stall_reached"
      else if (as.character(row$method) != "gamlss2" && isTRUE(hit_drop)) status <- "raw_loglik_deterioration"
      else if (as.character(row$method) != "gamlss2" && is.na(indicator)) status <- "missing_explicit_convergence_indicator"
      else if (as.character(row$method) != "gamlss2" && !isTRUE(indicator)) status <- "explicit_optimizer_nonconvergence"
      else status <- if (as.character(row$method) == "gamlss2") "finite_fit_below_outer_iteration_cap" else "explicit_optimizer_convergence"
      converged <- finite_payload && iteration_valid && if (as.character(row$method) == "gamlss2") iterations < cap else isTRUE(indicator) && !isTRUE(hit_outer) && !isTRUE(hit_stall) && !isTRUE(hit_drop) && iterations < cap
    }
    if (!identical(as.character(row$raw_convergence_status), status)) stop("Raw convergence named status contradicts its numeric/indicator evidence for attempt ", row$attempt_id, ".", call. = FALSE)
    data.frame(converged = converged, convergence_eligible = execution && converged,
      convergence_status = if (!execution) "not_run_to_completion" else if (converged) "converged" else "nonconverged", stringsAsFactors = FALSE)
  })
  do.call(rbind, result)
}

jss_recovery_require_clean_checkout <- function(git_identity) {
  if (!identical(unname(git_identity[["git_state"]]), "clean")) stop("Publication production requires a clean checkout with no modified or untracked files.", call. = FALSE)
  invisible(TRUE)
}

jss_recovery_require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop(label, " is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(x)
}

jss_recovery_logical <- function(x) {
  if (is.logical(x)) return(x)
  value <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(value))
  out[value %in% c("true", "t", "1", "yes")] <- TRUE
  out[value %in% c("false", "f", "0", "no")] <- FALSE
  out
}

jss_recovery_nonempty <- function(x, fallback = NA_character_) {
  x <- as.character(x)
  bad <- is.na(x) | !nzchar(trimws(x))
  x[bad] <- fallback
  x
}

jss_recovery_atomic_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(x, temporary, row.names = FALSE, na = "")
  if (file.exists(path)) {
    unlink(path)
    if (file.exists(path)) stop("Could not replace output: ", path, call. = FALSE)
  }
  if (!file.rename(temporary, path)) stop("Could not atomically install output: ", path, call. = FALSE)
  invisible(path)
}

jss_recovery_frame_identity <- function(frames) {
  rows <- lapply(names(frames), function(name) {
    temporary <- tempfile("normalized-input-", fileext = ".csv")
    on.exit(unlink(temporary), add = TRUE)
    utils::write.csv(frames[[name]], temporary, row.names = FALSE, na = "")
    data.frame(schema_version = jss_recovery_schema_version(), input_id = name,
      path = paste0(name, ".csv"), sha256 = unname(jss_recovery_sha256(temporary)),
      bytes = unname(file.info(temporary)$size), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

jss_recovery_public_input_map <- function() c(
  attempt_metadata = "attempt_metadata.csv", fixed_by_attempt = "fixed_by_attempt.csv",
  smooth_by_attempt = "smooth_by_attempt.csv", predictive_by_attempt = "predictive_by_attempt.csv",
  diagnostic_by_attempt = "diagnostic_by_attempt.csv", runner_settings_identity = "runner_settings_identity.csv",
  runner_settings_bcpe = "runner_settings_bcpe.csv", runner_settings_nbi = "runner_settings_nbi.csv"
)

jss_recovery_bundle_checkpoint_path <- function(output_dir) file.path(output_dir, ".main-recovery-checkpoint.rds")

jss_recovery_bundle_is_current <- function(output_dir, expected_outputs, provenance) {
  checkpoint_path <- jss_recovery_bundle_checkpoint_path(output_dir)
  if (!file.exists(checkpoint_path) || !all(file.exists(expected_outputs))) return(FALSE)
  checkpoint <- tryCatch(readRDS(checkpoint_path), error = function(e) NULL)
  if (!is.list(checkpoint) || !identical(checkpoint$schema_version, jss_recovery_schema_version())) return(FALSE)
  comparable <- c("schema_version", "producer_sha256", "analysis_config", "input_id", "sha256", "bytes")
  if (!all(comparable %in% names(checkpoint$input_provenance)) ||
      !isTRUE(all.equal(checkpoint$input_provenance[comparable], provenance[comparable], check.attributes = FALSE))) return(FALSE)
  current <- data.frame(
    file = basename(expected_outputs), sha256 = unname(jss_recovery_sha256(expected_outputs)),
    bytes = unname(file.info(expected_outputs)$size), stringsAsFactors = FALSE
  )
  isTRUE(all.equal(checkpoint$output_files, current, check.attributes = FALSE))
}

jss_recovery_install_bundle <- function(outputs, output_dir, provenance, repo_root) {
  parent <- dirname(output_dir)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  staging <- file.path(parent, paste0(basename(output_dir), ".staging-", Sys.getpid()))
  backup <- file.path(parent, paste0(basename(output_dir), ".backup-", Sys.getpid()))
  if (file.exists(staging) || file.exists(backup)) stop("Staging or backup path already exists.", call. = FALSE)
  dir.create(staging, recursive = FALSE)
  installed <- FALSE
  on.exit({
    if (!installed && dir.exists(staging)) unlink(staging, recursive = TRUE, force = TRUE)
    if (!installed && dir.exists(backup) && !dir.exists(output_dir)) file.rename(backup, output_dir)
  }, add = TRUE)
  if (dir.exists(output_dir)) {
    extras <- list.files(output_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    extras <- extras[!basename(extras) %in% c(paste0(names(outputs), ".csv"), ".main-recovery-checkpoint.rds")]
    if (length(extras) && !all(file.copy(extras, staging, recursive = TRUE, copy.mode = TRUE))) stop("Could not stage existing non-generated bundle files.", call. = FALSE)
  }
  staged_outputs <- vapply(names(outputs), function(name) {
    path <- file.path(staging, paste0(name, ".csv"))
    jss_recovery_atomic_csv(outputs[[name]], path)
    path
  }, character(1))
  output_manifest_path <- file.path(staging, "output_manifest.csv")
  output_manifest <- data.frame(
    schema_version = jss_recovery_schema_version(),
    artifact = basename(staged_outputs),
    sha256 = unname(jss_recovery_sha256(staged_outputs)),
    bytes = unname(file.info(staged_outputs)$size),
    stringsAsFactors = FALSE
  )
  jss_recovery_atomic_csv(output_manifest, output_manifest_path)
  portable_checkpoint_path <- file.path(staging, "bundle_checkpoint.csv")
  public_provenance <- outputs$input_provenance
  input_digest <- digest::digest(paste(public_provenance$input_id, public_provenance$sha256, sep = "=", collapse = "\n"), algo = "sha256", serialize = FALSE)
  portable_checkpoint <- data.frame(
    schema_version = jss_recovery_schema_version(),
    producer_sha256 = jss_recovery_producer_sha256(repo_root),
    input_provenance_sha256 = input_digest,
    output_manifest_sha256 = unname(jss_recovery_sha256(output_manifest_path)),
    artifact_count = nrow(output_manifest), stringsAsFactors = FALSE
  )
  jss_recovery_atomic_csv(portable_checkpoint, portable_checkpoint_path)
  staged_all <- c(staged_outputs, output_manifest = output_manifest_path, bundle_checkpoint = portable_checkpoint_path)
  checkpoint <- list(
    schema_version = jss_recovery_schema_version(), input_provenance = provenance,
    output_files = data.frame(file = basename(staged_all), sha256 = unname(jss_recovery_sha256(staged_all)),
      bytes = unname(file.info(staged_all)$size), stringsAsFactors = FALSE)
  )
  saveRDS(checkpoint, jss_recovery_bundle_checkpoint_path(staging))
  if (dir.exists(output_dir) && !file.rename(output_dir, backup)) stop("Could not move the previous bundle aside.", call. = FALSE)
  if (!file.rename(staging, output_dir)) {
    if (dir.exists(backup)) file.rename(backup, output_dir)
    stop("Could not atomically promote the staged evidence bundle.", call. = FALSE)
  }
  installed <- TRUE
  if (dir.exists(backup)) unlink(backup, recursive = TRUE, force = TRUE)
  file.path(output_dir, basename(staged_all))
}

jss_recovery_mean_stats <- function(x, conf.level = 0.95) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  estimate <- if (n) mean(x) else NA_real_
  empirical_sd <- if (n > 1L) stats::sd(x) else NA_real_
  mcse <- if (n > 1L) empirical_sd / sqrt(n) else NA_real_
  critical <- if (n > 1L) stats::qt(1 - (1 - conf.level) / 2, df = n - 1L) else NA_real_
  c(
    estimate = estimate,
    empirical_sd = empirical_sd,
    mcse = mcse,
    ci_low = if (is.finite(mcse)) estimate - critical * mcse else NA_real_,
    ci_high = if (is.finite(mcse)) estimate + critical * mcse else NA_real_,
    n_metric = n
  )
}

jss_recovery_rmse_stats <- function(error, conf.level = 0.95) {
  error <- as.numeric(error)
  error <- error[is.finite(error)]
  n <- length(error)
  if (!n) return(c(estimate = NA_real_, mcse = NA_real_, ci_low = NA_real_, ci_high = NA_real_, n_metric = 0))
  squared <- error^2
  mse <- mean(squared)
  estimate <- sqrt(mse)
  mse_mcse <- if (n > 1L) stats::sd(squared) / sqrt(n) else NA_real_
  mcse <- if (is.finite(mse_mcse) && estimate > 0) mse_mcse / (2 * estimate) else NA_real_
  critical <- if (n > 1L) stats::qt(1 - (1 - conf.level) / 2, df = n - 1L) else NA_real_
  mse_low <- if (is.finite(mse_mcse)) max(0, mse - critical * mse_mcse) else NA_real_
  mse_high <- if (is.finite(mse_mcse)) max(0, mse + critical * mse_mcse) else NA_real_
  c(
    estimate = estimate,
    mcse = mcse,
    ci_low = sqrt(mse_low),
    ci_high = sqrt(mse_high),
    n_metric = n
  )
}

jss_recovery_sd_stats <- function(x, conf.level = 0.95) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n <= 1L) return(c(estimate = NA_real_, mcse = NA_real_, ci_low = NA_real_, ci_high = NA_real_, n_metric = n))
  estimate <- stats::sd(x)
  alpha <- 1 - conf.level
  c(
    estimate = estimate,
    mcse = estimate / sqrt(2 * (n - 1L)),
    ci_low = sqrt((n - 1L) * estimate^2 / stats::qchisq(1 - alpha / 2, df = n - 1L)),
    ci_high = sqrt((n - 1L) * estimate^2 / stats::qchisq(alpha / 2, df = n - 1L)),
    n_metric = n
  )
}

jss_recovery_bootstrap_quantile <- function(x, probability, conf.level = 0.95, seed = 20260901L, B = 2000L) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (!n) return(c(estimate = NA_real_, mcse = NA_real_, ci_low = NA_real_, ci_high = NA_real_, n_metric = 0L))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))
  boot <- replicate(B, stats::quantile(sample(x, n, replace = TRUE), probability, names = FALSE, type = 8))
  alpha <- 1 - conf.level
  c(
    estimate = stats::quantile(x, probability, names = FALSE, type = 8),
    mcse = stats::sd(boot),
    ci_low = stats::quantile(boot, alpha / 2, names = FALSE, type = 8),
    ci_high = stats::quantile(boot, 1 - alpha / 2, names = FALSE, type = 8),
    n_metric = n
  )
}

jss_recovery_column <- function(x, name, fallback) {
  if (name %in% names(x)) return(x[[name]])
  if (length(fallback) == nrow(x)) fallback else rep(fallback, nrow(x))
}

jss_recovery_wilson <- function(success, n, conf.level = 0.95) {
  if (!is.finite(n) || n <= 0) return(c(estimate = NA_real_, mcse = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  p <- success / n
  z <- stats::qnorm(1 - (1 - conf.level) / 2)
  denominator <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denominator
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denominator
  c(
    estimate = p,
    mcse = sqrt(p * (1 - p) / n),
    ci_low = max(0, centre - half),
    ci_high = min(1, centre + half)
  )
}

jss_recovery_split_apply <- function(x, keys, fun) {
  if (!nrow(x)) return(data.frame())
  key <- interaction(x[keys], drop = TRUE, lex.order = TRUE, sep = "\r")
  rows <- lapply(split(x, key, drop = TRUE), fun)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

jss_recovery_cell_key <- function(x) paste(x$study_id, x$scenario_id, x$method, x$replicate, sep = "\r")

jss_recovery_output_attempt_key <- function(x) {
  if ("attempt_id" %in% names(x)) return(as.character(x$attempt_id))
  if ("retry_index" %in% names(x)) return(paste(jss_recovery_cell_key(x), x$retry_index, sep = "\r"))
  jss_recovery_cell_key(x)
}

jss_recovery_attach_attempt_identity <- function(x, ledger) {
  if (!nrow(x)) return(x)
  if (!"retry_index" %in% names(x)) x$retry_index <- NA_integer_
  if (anyNA(x$retry_index)) {
    candidates <- split(seq_len(nrow(ledger)), jss_recovery_cell_key(ledger))
    selected <- vapply(candidates, function(index) {
      eligible <- index[ledger$execution_success[index] %in% TRUE]
      if (!length(eligible)) eligible <- index
      eligible[which.max(ledger$retry_index[eligible])]
    }, integer(1))
    retry_map <- stats::setNames(ledger$retry_index[selected], names(selected))
    missing <- is.na(x$retry_index)
    x$retry_index[missing] <- unname(retry_map[jss_recovery_cell_key(x)[missing]])
  }
  attempt_key <- paste(jss_recovery_cell_key(ledger), ledger$retry_index, sep = "\r")
  attempt_map <- stats::setNames(ledger$attempt_id, attempt_key)
  x$attempt_id <- unname(attempt_map[paste(jss_recovery_cell_key(x), x$retry_index, sep = "\r")])
  x
}

jss_recovery_attempt_ledger <- function(bcpe_log, nbi_log) {
  jss_recovery_require_columns(
    bcpe_log, c("scenario", "n", "d", "rep", "model", "success", "converged", "elapsed_sec", "error"),
    "BCPE/t fit log"
  )
  raw_fields <- jss_recovery_raw_convergence_fields()
  for (input in list(bcpe_log = bcpe_log, nbi_log = nbi_log)) {
    production_input <- "evidence_status" %in% names(input) && any(as.character(input$evidence_status) == "post_phase1_production")
    present <- names(input)[startsWith(names(input), "raw_convergence_")]
    if (production_input && !identical(sort(present), sort(raw_fields))) stop("Production runner log raw convergence schema is missing fields or contains extras.", call. = FALSE)
  }
  jss_recovery_require_columns(
    nbi_log, c("rep", "engine", "success", "converged", "elapsed_sec", "error"),
    "NBI/Clayton fit log"
  )

  bcpe_seed <- if ("seed" %in% names(bcpe_log)) {
    as.integer(bcpe_log$seed)
  } else {
    110000L + as.integer(bcpe_log$rep)
  }
  bcpe_seed_source <- if ("seed" %in% names(bcpe_log)) "runner_metadata" else "reconstructed_runner_rule_110000_plus_rep"
  bcpe_stop <- if ("stop_reason" %in% names(bcpe_log)) bcpe_log$stop_reason else NA_character_
  bcpe_target <- if ("target_replicates" %in% names(bcpe_log)) as.integer(bcpe_log$target_replicates) else 100L
  bcpe_retry <- if ("retry_index" %in% names(bcpe_log)) {
    as.integer(bcpe_log$retry_index)
  } else {
    stats::ave(seq_len(nrow(bcpe_log)), bcpe_log$scenario, bcpe_log$model, bcpe_log$rep, FUN = seq_along)
  }
  bcpe_success <- jss_recovery_logical(bcpe_log$success)
  bcpe_converged <- jss_recovery_logical(bcpe_log$converged)

  bcpe <- data.frame(
    schema_version = jss_recovery_schema_version(),
    evidence_status = jss_recovery_column(bcpe_log, "evidence_status", "legacy_pre_phase1_reconciliation"),
    study_id = jss_recovery_column(bcpe_log, "study_id", "bcpe_t_legacy_reconciliation"),
    scenario_id = as.character(bcpe_log$scenario),
    margin_family = jss_recovery_column(bcpe_log, "margin_family", "BCPE"),
    copula_family = jss_recovery_column(bcpe_log, "copula_family", "t"),
    copula_code = jss_recovery_column(bcpe_log, "copula_code", "t"),
    runner_contract_version = jss_recovery_column(bcpe_log, "runner_contract_version", "legacy_unversioned"),
    phase1_contract_version = jss_recovery_column(bcpe_log, "phase1_contract_version", "legacy_pre_phase1"),
    runner_settings_signature = jss_recovery_column(bcpe_log, "runner_settings_signature", "legacy_unrecorded"),
    runner_settings_sha256 = jss_recovery_column(bcpe_log, "runner_settings_sha256", "legacy_unrecorded"),
    runner_sha256 = jss_recovery_column(bcpe_log, "runner_sha256", "legacy_unrecorded"),
    package_source_path = jss_recovery_column(bcpe_log, "package_source_path", "legacy_unrecorded"),
    package_version = jss_recovery_column(bcpe_log, "package_version", "legacy_unrecorded"),
    package_source_sha256 = jss_recovery_column(bcpe_log, "package_source_sha256", "legacy_unrecorded"),
    n_subjects = as.integer(bcpe_log$n),
    n_time = as.integer(bcpe_log$d),
    target_replicates = bcpe_target,
    replicate = as.integer(bcpe_log$rep),
    data_seed = bcpe_seed,
    training_covariate_seed = as.integer(jss_recovery_column(bcpe_log, "training_covariate_seed", NA_integer_)),
    training_response_seed = as.integer(jss_recovery_column(bcpe_log, "training_response_seed", bcpe_seed)),
    test_response_seed = as.integer(jss_recovery_column(bcpe_log, "test_response_seed", bcpe_seed + 500000L)),
    diagnostic_seed = as.integer(jss_recovery_column(bcpe_log, "diagnostic_seed", NA_integer_)),
    predictive_seed = as.integer(jss_recovery_column(bcpe_log, "predictive_seed", NA_integer_)),
    execution_completed_at_utc = as.character(jss_recovery_column(bcpe_log, "execution_completed_at_utc", NA_character_)),
    seed_source = bcpe_seed_source,
    method = as.character(bcpe_log$model),
    attempted = TRUE,
    success = bcpe_success,
    execution_success = jss_recovery_logical(jss_recovery_column(bcpe_log, "execution_success", bcpe_success)),
    converged = bcpe_converged,
    descriptive_outputs_complete = jss_recovery_logical(jss_recovery_column(bcpe_log, "descriptive_outputs_complete", NA)),
    output_failure_reason = jss_recovery_nonempty(jss_recovery_column(bcpe_log, "output_failure_reason", NA_character_)),
    runtime_n_cores = as.integer(jss_recovery_column(bcpe_log, "runtime_n_cores", NA_integer_)),
    runtime_backend = jss_recovery_nonempty(jss_recovery_column(bcpe_log, "runtime_backend", NA_character_)),
    runtime_rscript_sha256 = jss_recovery_nonempty(jss_recovery_column(bcpe_log, "runtime_rscript_sha256", NA_character_)),
    elapsed_sec = as.numeric(bcpe_log$elapsed_sec),
    stop_reason = jss_recovery_nonempty(bcpe_stop),
    error = jss_recovery_nonempty(bcpe_log$error),
    retry_index = as.integer(bcpe_retry),
    stringsAsFactors = FALSE
  )

  nbi_success <- jss_recovery_logical(nbi_log$success)
  nbi_converged <- jss_recovery_logical(nbi_log$converged)
  nbi_scenario <- if ("scenario" %in% names(nbi_log)) as.character(nbi_log$scenario) else "n500_d4_nbi_signal2"
  nbi_n <- if ("n" %in% names(nbi_log)) as.integer(nbi_log$n) else 500L
  nbi_d <- if ("d" %in% names(nbi_log)) as.integer(nbi_log$d) else 4L
  nbi_seed <- if ("seed" %in% names(nbi_log)) as.integer(nbi_log$seed) else 20260529L + as.integer(nbi_log$rep)
  nbi_seed_source <- if ("seed" %in% names(nbi_log)) "runner_metadata" else "reconstructed_runner_rule_20260529_plus_rep"
  nbi_stop <- if ("stop_reason" %in% names(nbi_log)) nbi_log$stop_reason else NA_character_
  nbi_target <- if ("target_replicates" %in% names(nbi_log)) as.integer(nbi_log$target_replicates) else 100L
  nbi_retry <- if ("retry_index" %in% names(nbi_log)) {
    as.integer(nbi_log$retry_index)
  } else {
    stats::ave(seq_len(nrow(nbi_log)), nbi_scenario, nbi_log$engine, nbi_log$rep, FUN = seq_along)
  }

  nbi <- data.frame(
    schema_version = jss_recovery_schema_version(),
    evidence_status = jss_recovery_column(nbi_log, "evidence_status", "legacy_pre_phase1_reconciliation"),
    study_id = jss_recovery_column(nbi_log, "study_id", "nbi_gaussian_legacy_reconciliation"),
    scenario_id = nbi_scenario,
    margin_family = jss_recovery_column(nbi_log, "margin_family", "NBI"),
    copula_family = jss_recovery_column(nbi_log, "copula_family", "Gaussian"),
    copula_code = jss_recovery_column(nbi_log, "copula_code", "N"),
    runner_contract_version = jss_recovery_column(nbi_log, "runner_contract_version", "legacy_unversioned_gaussian_runner"),
    phase1_contract_version = jss_recovery_column(nbi_log, "phase1_contract_version", "legacy_pre_phase1"),
    runner_settings_signature = jss_recovery_column(nbi_log, "runner_settings_signature", "legacy_unrecorded"),
    runner_settings_sha256 = jss_recovery_column(nbi_log, "runner_settings_sha256", "legacy_unrecorded"),
    runner_sha256 = jss_recovery_column(nbi_log, "runner_sha256", "legacy_unrecorded"),
    package_source_path = jss_recovery_column(nbi_log, "package_source_path", "legacy_unrecorded"),
    package_version = jss_recovery_column(nbi_log, "package_version", "legacy_unrecorded"),
    package_source_sha256 = jss_recovery_column(nbi_log, "package_source_sha256", "legacy_unrecorded"),
    n_subjects = nbi_n,
    n_time = nbi_d,
    target_replicates = nbi_target,
    replicate = as.integer(nbi_log$rep),
    data_seed = nbi_seed,
    training_covariate_seed = as.integer(jss_recovery_column(nbi_log, "training_covariate_seed", NA_integer_)),
    training_response_seed = as.integer(jss_recovery_column(nbi_log, "training_response_seed", nbi_seed)),
    test_response_seed = as.integer(jss_recovery_column(nbi_log, "test_response_seed", nbi_seed + 500000L)),
    diagnostic_seed = as.integer(jss_recovery_column(nbi_log, "diagnostic_seed", NA_integer_)),
    predictive_seed = as.integer(jss_recovery_column(nbi_log, "predictive_seed", NA_integer_)),
    execution_completed_at_utc = as.character(jss_recovery_column(nbi_log, "execution_completed_at_utc", NA_character_)),
    seed_source = nbi_seed_source,
    method = as.character(nbi_log$engine),
    attempted = TRUE,
    success = nbi_success,
    execution_success = jss_recovery_logical(jss_recovery_column(nbi_log, "execution_success", nbi_success)),
    converged = nbi_converged,
    descriptive_outputs_complete = jss_recovery_logical(jss_recovery_column(nbi_log, "descriptive_outputs_complete", NA)),
    output_failure_reason = jss_recovery_nonempty(jss_recovery_column(nbi_log, "output_failure_reason", NA_character_)),
    runtime_n_cores = as.integer(jss_recovery_column(nbi_log, "runtime_n_cores", NA_integer_)),
    runtime_backend = jss_recovery_nonempty(jss_recovery_column(nbi_log, "runtime_backend", NA_character_)),
    runtime_rscript_sha256 = jss_recovery_nonempty(jss_recovery_column(nbi_log, "runtime_rscript_sha256", NA_character_)),
    elapsed_sec = as.numeric(nbi_log$elapsed_sec),
    stop_reason = jss_recovery_nonempty(nbi_stop),
    error = jss_recovery_nonempty(nbi_log$error),
    retry_index = as.integer(nbi_retry),
    stringsAsFactors = FALSE
  )
  for (field in raw_fields) {
    bcpe[[field]] <- jss_recovery_column(bcpe_log, field, NA)
    nbi[[field]] <- jss_recovery_column(nbi_log, field, NA)
  }

  common <- union(names(bcpe), names(nbi))
  for (column in setdiff(common, names(bcpe))) bcpe[[column]] <- NA
  for (column in setdiff(common, names(nbi))) nbi[[column]] <- NA
  out <- rbind(bcpe[common], nbi[common])
  out$success <- out$execution_success
  out$convergence_eligible <- out$execution_success %in% TRUE & out$converged %in% TRUE
  out$convergence_status <- ifelse(
    !out$execution_success, "not_run_to_completion",
    ifelse(is.na(out$converged), "not_reported", ifelse(out$converged, "converged", "nonconverged"))
  )
  out$failure_reason <- ifelse(
    !out$execution_success,
    ifelse(is.na(out$error), "execution_error_unspecified", paste0("execution_error: ", out$error)),
    ifelse(out$converged %in% FALSE,
      ifelse(is.na(out$stop_reason), "optimizer_nonconvergence", paste0("optimizer_nonconvergence: ", out$stop_reason)),
      ifelse(is.na(out$converged), "convergence_not_reported",
        ifelse(!is.na(out$output_failure_reason), paste0("incomplete_postfit_outputs: ", out$output_failure_reason), "none"))
    )
  )
  # Publication retention is intentionally stricter than execution. It is
  # reconciled against exact descriptive-output cardinality below.
  out$retained <- out$convergence_eligible
  out$retention_reason <- ifelse(out$retained, "provisional_converged_pending_output_reconciliation", out$failure_reason)
  out$attempt_id <- paste(out$study_id, out$scenario_id, out$method, sprintf("rep%03d", out$replicate), sprintf("try%02d", out$retry_index), sep = "::")
  production <- out$evidence_status == "post_phase1_production"
  if (any(production)) {
    recomputed <- jss_recovery_recompute_convergence(out[production, , drop = FALSE])
    contradiction <- !identical(as.logical(out$converged[production]), as.logical(recomputed$converged)) ||
      !identical(as.logical(out$convergence_eligible[production]), as.logical(recomputed$convergence_eligible)) ||
      !identical(as.character(out$convergence_status[production]), as.character(recomputed$convergence_status))
    if (contradiction) stop("Runner-derived convergence fields contradict strict raw convergence evidence.", call. = FALSE)
    out$converged[production] <- recomputed$converged
    out$convergence_eligible[production] <- recomputed$convergence_eligible
    out$convergence_status[production] <- recomputed$convergence_status
  }
  out <- out[order(out$study_id, out$scenario_id, out$method, out$replicate, out$retry_index), ]
  rownames(out) <- NULL
  out
}

jss_recovery_counts <- function(ledger) {
  jss_recovery_split_apply(ledger, c("study_id", "scenario_id", "method"), function(df) {
    collapse_recorded <- function(x) {
      value <- sort(unique(as.character(x[!is.na(x)])))
      if (length(value)) paste(value, collapse = "|") else "not_recorded"
    }
    planned <- unique(as.integer(df$target_replicates))
    if (length(planned) != 1L || !is.finite(planned)) planned <- length(unique(df$replicate))
    attempted_cells <- length(unique(df$replicate[df$attempted %in% TRUE]))
    cell_any <- function(flag) length(unique(df$replicate[flag %in% TRUE]))
    execution_cells <- cell_any(df$execution_success)
    converged_cells <- cell_any(df$convergence_eligible)
    retained_cells <- cell_any(df$retained)
    failed_cells <- planned - retained_cells
    attempted_rate <- jss_recovery_wilson(attempted_cells, planned)
    execution_rate <- jss_recovery_wilson(execution_cells, planned)
    converged_rate <- jss_recovery_wilson(converged_cells, planned)
    retained_rate <- jss_recovery_wilson(retained_cells, planned)
    failure_rate <- jss_recovery_wilson(failed_cells, planned)
    data.frame(
      study_id = df$study_id[[1]], scenario_id = df$scenario_id[[1]], method = df$method[[1]],
      runtime_n_cores = collapse_recorded(df$runtime_n_cores[is.finite(df$runtime_n_cores)]),
      runtime_backend = collapse_recorded(df$runtime_backend),
      runtime_rscript_sha256 = collapse_recorded(df$runtime_rscript_sha256),
      planned_method_replicate_cells = planned,
      actual_attempt_rows = sum(df$attempted %in% TRUE),
      retry_rows = sum(df$attempted %in% TRUE) - attempted_cells,
      attempted = attempted_cells,
      successful = execution_cells,
      converged = converged_cells,
      convergence_not_reported = sum(is.na(df$converged) & df$success %in% TRUE),
      convergence_not_reported_attempt_rows = sum(is.na(df$converged) & df$success %in% TRUE),
      retained = retained_cells,
      failed = failed_cells,
      attempted_rate = attempted_rate[["estimate"]], attempted_rate_mcse = attempted_rate[["mcse"]], attempted_rate_ci_low = attempted_rate[["ci_low"]], attempted_rate_ci_high = attempted_rate[["ci_high"]],
      execution_success_rate = execution_rate[["estimate"]], execution_success_rate_mcse = execution_rate[["mcse"]], execution_success_rate_ci_low = execution_rate[["ci_low"]], execution_success_rate_ci_high = execution_rate[["ci_high"]],
      convergence_rate = converged_rate[["estimate"]], convergence_rate_mcse = converged_rate[["mcse"]], convergence_rate_ci_low = converged_rate[["ci_low"]], convergence_rate_ci_high = converged_rate[["ci_high"]],
      retention_rate = retained_rate[["estimate"]], retention_rate_mcse = retained_rate[["mcse"]], retention_rate_ci_low = retained_rate[["ci_low"]], retention_rate_ci_high = retained_rate[["ci_high"]],
      failure_rate = failure_rate[["estimate"]], failure_rate_mcse = failure_rate[["mcse"]], failure_rate_ci_low = failure_rate[["ci_low"]], failure_rate_ci_high = failure_rate[["ci_high"]],
      rate_denominator = planned,
      stringsAsFactors = FALSE
    )
  })
}

jss_recovery_reconcile_retention <- function(ledger, fixed, smooth, predictive, diagnostic) {
  fixed <- jss_recovery_attach_attempt_identity(fixed, ledger)
  smooth <- jss_recovery_attach_attempt_identity(smooth, ledger)
  predictive <- jss_recovery_attach_attempt_identity(predictive, ledger)
  diagnostic <- jss_recovery_attach_attempt_identity(diagnostic, ledger)
  contract_check <- lapply(seq_len(nrow(ledger)), function(i) {
    contract <- jss_recovery_metric_contract(ledger$margin_family[[i]], ledger$method[[i]])
    if (is.null(contract)) return(c(fixed = FALSE, smooth = FALSE, predictive = FALSE, diagnostic = FALSE))
    fixed_i <- fixed[fixed$attempt_id == ledger$attempt_id[[i]], , drop = FALSE]
    smooth_i <- smooth[smooth$attempt_id == ledger$attempt_id[[i]], , drop = FALSE]
    predictive_i <- predictive[predictive$attempt_id == ledger$attempt_id[[i]], , drop = FALSE]
    diagnostic_i <- diagnostic[diagnostic$attempt_id == ledger$attempt_id[[i]], , drop = FALSE]
    jss_recovery_attempt_output_contract(fixed_i, smooth_i, predictive_i, diagnostic_i, contract)
  })
  contract_matrix <- do.call(rbind, contract_check)
  ledger$has_fixed_output <- contract_matrix[, "fixed"]
  ledger$has_smooth_output <- contract_matrix[, "smooth"]
  ledger$has_predictive_output <- contract_matrix[, "predictive"]
  ledger$has_diagnostic_output <- contract_matrix[, "diagnostic"]
  complete <- apply(contract_matrix, 1L, all)
  declared_complete <- is.na(ledger$descriptive_outputs_complete) | ledger$descriptive_outputs_complete %in% TRUE
  complete <- complete & declared_complete
  ledger$descriptive_outputs_complete <- complete
  candidate <- ledger$execution_success %in% TRUE & ledger$convergence_eligible %in% TRUE & complete
  group <- interaction(ledger[c("study_id", "scenario_id", "method", "replicate")], drop = TRUE, lex.order = TRUE)
  ledger$retained <- FALSE
  for (indices in split(seq_len(nrow(ledger)), group)) {
    eligible <- indices[candidate[indices]]
    if (length(eligible)) ledger$retained[eligible[which.max(ledger$retry_index[eligible])]] <- TRUE
  }
  missing_outputs <- vapply(seq_len(nrow(ledger)), function(i) {
    missing <- c(
      if (!ledger$has_fixed_output[[i]]) "fixed" else NULL,
      if (!ledger$has_smooth_output[[i]]) "smooth" else NULL,
      if (!ledger$has_predictive_output[[i]]) "predictive" else NULL,
      if (!ledger$has_diagnostic_output[[i]]) "diagnostic" else NULL
    )
    paste(missing, collapse = "|")
  }, character(1))
  ledger$retention_reason <- ifelse(
    !(ledger$execution_success %in% TRUE), ledger$failure_reason,
    ifelse(!(ledger$convergence_eligible %in% TRUE),
      ifelse(ledger$converged %in% FALSE,
        ifelse(is.na(ledger$stop_reason), "optimizer_nonconvergence", paste0("optimizer_nonconvergence: ", ledger$stop_reason)),
        "convergence_not_reported"),
      ifelse(!complete, paste0("incomplete_postfit_outputs: ", ifelse(!is.na(ledger$output_failure_reason), ledger$output_failure_reason, missing_outputs)),
        ifelse(ledger$retained, "publication_retained_converged_complete", "superseded_eligible_retry")))
  )
  failed <- !ledger$retained & !grepl("^superseded_", ledger$retention_reason)
  ledger$failure_reason[failed] <- ledger$retention_reason[failed]
  ledger
}

jss_recovery_design_table <- function(ledger) {
  jss_recovery_split_apply(ledger, c("study_id", "scenario_id"), function(df) {
    replicas <- sort(unique(df$replicate[df$attempted %in% TRUE]))
    seeds <- sort(unique(df$data_seed[df$attempted %in% TRUE & is.finite(df$data_seed)]))
    methods <- sort(unique(df$method))
    data.frame(
      schema_version = unique(df$schema_version)[[1]],
      evidence_status = unique(df$evidence_status)[[1]],
      study_id = df$study_id[[1]], scenario_id = df$scenario_id[[1]],
      margin_family = unique(df$margin_family)[[1]], copula_family = unique(df$copula_family)[[1]], copula_code = unique(df$copula_code)[[1]],
      n_subjects = unique(df$n_subjects)[[1]], n_time = unique(df$n_time)[[1]],
      target_replicates = unique(df$target_replicates)[[1]],
      attempted_replicates = length(replicas),
      replicate_min = min(replicas), replicate_max = max(replicas),
      data_seed_min = min(seeds), data_seed_max = max(seeds),
      seed_source = paste(sort(unique(df$seed_source)), collapse = "|"),
      methods = paste(methods, collapse = "|"),
      attempted_method_fits = sum(df$attempted %in% TRUE),
      stringsAsFactors = FALSE
    )
  })
}

jss_recovery_failure_summary <- function(ledger, conf.level = 0.95) {
  jss_recovery_split_apply(ledger, c("study_id", "scenario_id", "method", "failure_reason"), function(df) {
    denominator <- sum(ledger$study_id == df$study_id[[1]] & ledger$scenario_id == df$scenario_id[[1]] & ledger$method == df$method[[1]])
    interval <- jss_recovery_wilson(nrow(df), denominator, conf.level)
    data.frame(
      study_id = df$study_id[[1]], scenario_id = df$scenario_id[[1]], method = df$method[[1]],
      failure_reason = df$failure_reason[[1]], count = nrow(df),
      attempt_row_denominator = denominator,
      proportion_of_attempts = interval[["estimate"]], rate_mcse = interval[["mcse"]],
      rate_ci_low = interval[["ci_low"]], rate_ci_high = interval[["ci_high"]],
      stringsAsFactors = FALSE
    )
  })
}

jss_recovery_study_identity <- function(ledger, margin_family) {
  x <- unique(ledger[ledger$margin_family == margin_family, c("study_id", "scenario_id")])
  if (nrow(x) != 1L) stop("Expected one study/scenario identity for margin ", margin_family, ".", call. = FALSE)
  x
}

jss_recovery_adapt_fixed <- function(bcpe_fixed, nbi_fixed, ledger) {
  required <- c("model", "parameter", "term", "rep", "estimate", "std_error", "true_value")
  jss_recovery_require_columns(bcpe_fixed, c("scenario", required), "BCPE/t fixed-effect data")
  jss_recovery_require_columns(nbi_fixed, required, "NBI/Clayton fixed-effect data")
  bcpe_id <- jss_recovery_study_identity(ledger, "BCPE")
  nbi_id <- jss_recovery_study_identity(ledger, "NBI")
  bcpe <- data.frame(
    study_id = bcpe_id$study_id, scenario_id = as.character(bcpe_fixed$scenario),
    method = as.character(bcpe_fixed$model), replicate = as.integer(bcpe_fixed$rep),
    parameter = as.character(bcpe_fixed$parameter), term = as.character(bcpe_fixed$term),
    estimate = as.numeric(bcpe_fixed$estimate), std_error = as.numeric(bcpe_fixed$std_error),
    true_value = as.numeric(bcpe_fixed$true_value),
    inference_status = as.character(jss_recovery_column(bcpe_fixed, "inference_status", ifelse(is.finite(as.numeric(bcpe_fixed$std_error)), "available", "unavailable_unspecified"))),
    inference_denominator = as.numeric(jss_recovery_column(bcpe_fixed, "inference_denominator", ifelse(is.finite(as.numeric(bcpe_fixed$std_error)), 1, 0))),
    retry_index = as.integer(jss_recovery_column(bcpe_fixed, "retry_index", NA_integer_)), stringsAsFactors = FALSE
  )
  nbi <- data.frame(
    study_id = nbi_id$study_id, scenario_id = nbi_id$scenario_id,
    method = as.character(nbi_fixed$model), replicate = as.integer(nbi_fixed$rep),
    parameter = as.character(nbi_fixed$parameter), term = as.character(nbi_fixed$term),
    estimate = as.numeric(nbi_fixed$estimate), std_error = as.numeric(nbi_fixed$std_error),
    true_value = as.numeric(nbi_fixed$true_value),
    inference_status = as.character(jss_recovery_column(nbi_fixed, "inference_status", ifelse(is.finite(as.numeric(nbi_fixed$std_error)), "available", "unavailable_unspecified"))),
    inference_denominator = as.numeric(jss_recovery_column(nbi_fixed, "inference_denominator", ifelse(is.finite(as.numeric(nbi_fixed$std_error)), 1, 0))),
    retry_index = as.integer(jss_recovery_column(nbi_fixed, "retry_index", NA_integer_)), stringsAsFactors = FALSE
  )
  out <- rbind(bcpe, nbi)
  if (all(is.na(out$retry_index))) out$retry_index <- NULL
  jss_recovery_attach_attempt_identity(out, ledger)
}

jss_recovery_fixed_summary <- function(fixed, ledger, conf.level = 0.95) {
  counts <- jss_recovery_counts(ledger)
  fixed <- jss_recovery_attach_attempt_identity(fixed, ledger)
  if (!"inference_status" %in% names(fixed)) fixed$inference_status <- ifelse(is.finite(fixed$std_error), "available", "unavailable_unspecified")
  if (!"inference_denominator" %in% names(fixed)) fixed$inference_denominator <- as.integer(is.finite(fixed$std_error))
  retained <- ledger[ledger$retained %in% TRUE, "attempt_id", drop = FALSE]
  x <- merge(fixed, retained, by = "attempt_id", all = FALSE)
  out <- jss_recovery_split_apply(x, c("study_id", "scenario_id", "method", "parameter", "term"), function(df) {
    error <- df$estimate - df$true_value
    bias <- jss_recovery_mean_stats(error, conf.level)
    rmse <- jss_recovery_rmse_stats(error, conf.level)
    estimate <- jss_recovery_mean_stats(df$estimate, conf.level)
    estimate_sd <- jss_recovery_sd_stats(df$estimate, conf.level)
    se <- jss_recovery_mean_stats(df$std_error, conf.level)
    valid_ci <- is.finite(df$estimate) & is.finite(df$std_error) & is.finite(df$true_value)
    covered <- valid_ci & df$true_value >= df$estimate - stats::qnorm(1 - (1 - conf.level) / 2) * df$std_error &
      df$true_value <= df$estimate + stats::qnorm(1 - (1 - conf.level) / 2) * df$std_error
    coverage <- jss_recovery_wilson(sum(covered[valid_ci]), sum(valid_ci), conf.level)
    data.frame(
      study_id = df$study_id[[1]], scenario_id = df$scenario_id[[1]], method = df$method[[1]],
      parameter = df$parameter[[1]], term = df$term[[1]],
      n_metric = as.integer(rmse[["n_metric"]]), n_se = sum(is.finite(df$std_error)),
      inference_denominator = sum(as.numeric(df$inference_denominator), na.rm = TRUE),
      n_inference_available = sum(df$inference_status == "available", na.rm = TRUE),
      n_inference_unavailable = sum(df$inference_status != "available" | is.na(df$inference_status)),
      inference_status = paste(sort(unique(as.character(df$inference_status))), collapse = "|"),
      estimand_status = if (rmse[["n_metric"]] > 0) "estimated" else "not_applicable_or_unavailable",
      true_value_mean = mean(df$true_value, na.rm = TRUE), mean_estimate = estimate[["estimate"]],
      bias = bias[["estimate"]], bias_mcse = bias[["mcse"]], bias_ci_low = bias[["ci_low"]], bias_ci_high = bias[["ci_high"]],
      rmse = rmse[["estimate"]], rmse_mcse = rmse[["mcse"]], rmse_ci_low = rmse[["ci_low"]], rmse_ci_high = rmse[["ci_high"]],
      empirical_sd = estimate_sd[["estimate"]], empirical_sd_mcse = estimate_sd[["mcse"]],
      empirical_sd_ci_low = estimate_sd[["ci_low"]], empirical_sd_ci_high = estimate_sd[["ci_high"]],
      mean_se = se[["estimate"]], mean_se_mcse = se[["mcse"]], mean_se_ci_low = se[["ci_low"]], mean_se_ci_high = se[["ci_high"]],
      se_to_empirical_sd = se[["estimate"]] / estimate[["empirical_sd"]],
      coverage = coverage[["estimate"]], coverage_mcse = coverage[["mcse"]], coverage_ci_low = coverage[["ci_low"]], coverage_ci_high = coverage[["ci_high"]],
      stringsAsFactors = FALSE
    )
  })
  merge(counts, out, by = c("study_id", "scenario_id", "method"), all.y = TRUE, sort = FALSE)
}

jss_recovery_adapt_smooth <- function(bcpe_smooth, nbi_smooth, ledger) {
  required <- c("scenario", "model", "parameter", "rep", "bias_abs_integrated", "irmse")
  jss_recovery_require_columns(bcpe_smooth, required, "BCPE/t smooth data")
  jss_recovery_require_columns(nbi_smooth, required, "NBI/Clayton smooth data")
  bcpe_smooth$study_id <- jss_recovery_study_identity(ledger, "BCPE")$study_id
  nbi_smooth$study_id <- jss_recovery_study_identity(ledger, "NBI")$study_id
  x <- rbind(
    bcpe_smooth[c("study_id", "scenario", "model", "parameter", "rep", "bias_abs_integrated", "irmse")],
    nbi_smooth[c("study_id", "scenario", "model", "parameter", "rep", "bias_abs_integrated", "irmse")]
  )
  names(x)[names(x) == "scenario"] <- "scenario_id"
  names(x)[names(x) == "model"] <- "method"
  names(x)[names(x) == "rep"] <- "replicate"
  if ("retry_index" %in% names(bcpe_smooth) || "retry_index" %in% names(nbi_smooth)) {
    # rbind above preserves retry_index only when it is part of the selected set.
    x$retry_index <- c(
      as.integer(jss_recovery_column(bcpe_smooth, "retry_index", NA_integer_)),
      as.integer(jss_recovery_column(nbi_smooth, "retry_index", NA_integer_))
    )
    if (all(is.na(x$retry_index))) x$retry_index <- NULL
  }
  jss_recovery_attach_attempt_identity(x, ledger)
}

jss_recovery_smooth_summary <- function(smooth, ledger, conf.level = 0.95) {
  counts <- jss_recovery_counts(ledger)
  smooth <- jss_recovery_attach_attempt_identity(smooth, ledger)
  retained <- ledger[ledger$retained %in% TRUE, "attempt_id", drop = FALSE]
  x <- merge(smooth, retained, by = "attempt_id", all = FALSE)
  out <- jss_recovery_split_apply(x, c("study_id", "scenario_id", "method", "parameter"), function(df) {
    irmse <- jss_recovery_rmse_stats(df$irmse, conf.level)
    irmse_sd <- jss_recovery_sd_stats(df$irmse, conf.level)
    abs_bias <- jss_recovery_mean_stats(df$bias_abs_integrated, conf.level)
    data.frame(
      study_id = df$study_id[[1]], scenario_id = df$scenario_id[[1]], method = df$method[[1]], parameter = df$parameter[[1]],
      n_metric = as.integer(irmse[["n_metric"]]),
      irmse = irmse[["estimate"]], irmse_mcse = irmse[["mcse"]], irmse_ci_low = irmse[["ci_low"]], irmse_ci_high = irmse[["ci_high"]],
      replicate_irmse_empirical_sd = irmse_sd[["estimate"]], replicate_irmse_sd_mcse = irmse_sd[["mcse"]],
      replicate_irmse_sd_ci_low = irmse_sd[["ci_low"]], replicate_irmse_sd_ci_high = irmse_sd[["ci_high"]],
      mean_integrated_abs_bias = abs_bias[["estimate"]], integrated_abs_bias_mcse = abs_bias[["mcse"]],
      integrated_abs_bias_ci_low = abs_bias[["ci_low"]], integrated_abs_bias_ci_high = abs_bias[["ci_high"]],
      stringsAsFactors = FALSE
    )
  })
  merge(counts, out, by = c("study_id", "scenario_id", "method"), all.y = TRUE, sort = FALSE)
}

jss_recovery_runtime_summary <- function(ledger, conf.level = 0.95) {
  counts <- jss_recovery_counts(ledger)
  x <- ledger[ledger$retained %in% TRUE, ]
  out <- jss_recovery_split_apply(x, c("study_id", "scenario_id", "method"), function(df) {
    stats <- jss_recovery_mean_stats(df$elapsed_sec, conf.level)
    runtime_sd <- jss_recovery_sd_stats(df$elapsed_sec, conf.level)
    finite <- df$elapsed_sec[is.finite(df$elapsed_sec)]
    seed_base <- 20260901L + sum(utf8ToInt(paste(df$study_id[[1]], df$method[[1]])))
    q25 <- jss_recovery_bootstrap_quantile(finite, 0.25, conf.level, seed_base + 1L)
    median <- jss_recovery_bootstrap_quantile(finite, 0.5, conf.level, seed_base + 2L)
    q75 <- jss_recovery_bootstrap_quantile(finite, 0.75, conf.level, seed_base + 3L)
    data.frame(
      study_id = df$study_id[[1]], scenario_id = df$scenario_id[[1]], method = df$method[[1]],
      n_runtime = length(finite), mean_runtime_sec = stats[["estimate"]], runtime_empirical_sd = runtime_sd[["estimate"]],
      runtime_empirical_sd_mcse = runtime_sd[["mcse"]], runtime_empirical_sd_ci_low = runtime_sd[["ci_low"]], runtime_empirical_sd_ci_high = runtime_sd[["ci_high"]],
      runtime_mcse = stats[["mcse"]], runtime_ci_low = stats[["ci_low"]], runtime_ci_high = stats[["ci_high"]],
      median_runtime_sec = median[["estimate"]], median_runtime_mcse = median[["mcse"]], median_runtime_ci_low = median[["ci_low"]], median_runtime_ci_high = median[["ci_high"]],
      runtime_q25 = q25[["estimate"]], runtime_q25_mcse = q25[["mcse"]], runtime_q25_ci_low = q25[["ci_low"]], runtime_q25_ci_high = q25[["ci_high"]],
      runtime_q75 = q75[["estimate"]], runtime_q75_mcse = q75[["mcse"]], runtime_q75_ci_low = q75[["ci_low"]], runtime_q75_ci_high = q75[["ci_high"]],
      stringsAsFactors = FALSE
    )
  })
  merge(counts, out, by = c("study_id", "scenario_id", "method"), all.y = TRUE, sort = FALSE)
}

jss_recovery_adapt_by_rep <- function(bcpe, nbi, study_ids, ledger = NULL) {
  required <- c("scenario", "model", "rep")
  jss_recovery_require_columns(bcpe, required, "BCPE/t per-replicate metrics")
  jss_recovery_require_columns(nbi, required, "NBI/Clayton per-replicate metrics")
  bcpe$study_id <- study_ids[[1]]
  nbi$study_id <- study_ids[[2]]
  common <- union(names(bcpe), names(nbi))
  for (column in setdiff(common, names(bcpe))) bcpe[[column]] <- NA
  for (column in setdiff(common, names(nbi))) nbi[[column]] <- NA
  x <- rbind(bcpe[common], nbi[common])
  names(x)[names(x) == "scenario"] <- "scenario_id"
  names(x)[names(x) == "model"] <- "method"
  names(x)[names(x) == "rep"] <- "replicate"
  if (!is.null(ledger)) jss_recovery_attach_attempt_identity(x, ledger) else x
}

jss_recovery_long_metric_summary <- function(x, ledger, metrics, conf.level = 0.95, extra_keys = character()) {
  counts <- jss_recovery_counts(ledger)
  x <- jss_recovery_attach_attempt_identity(x, ledger)
  retained <- ledger[ledger$retained %in% TRUE, "attempt_id", drop = FALSE]
  x <- merge(x, retained, by = "attempt_id", all = FALSE)
  rows <- list()
  for (metric in metrics) {
    if (!metric %in% names(x)) next
    keys <- c("study_id", "scenario_id", "method", extra_keys)
    tmp <- jss_recovery_split_apply(x, keys, function(df) {
      stats <- jss_recovery_mean_stats(df[[metric]], conf.level)
      metric_sd <- jss_recovery_sd_stats(df[[metric]], conf.level)
      key_values <- df[1, keys, drop = FALSE]
      data.frame(
        key_values, metric = metric, n_metric = as.integer(stats[["n_metric"]]),
        mean = stats[["estimate"]], empirical_sd = metric_sd[["estimate"]],
        empirical_sd_mcse = metric_sd[["mcse"]], empirical_sd_ci_low = metric_sd[["ci_low"]], empirical_sd_ci_high = metric_sd[["ci_high"]],
        mcse = stats[["mcse"]],
        ci_low = stats[["ci_low"]], ci_high = stats[["ci_high"]], stringsAsFactors = FALSE, check.names = FALSE
      )
    })
    rows[[length(rows) + 1L]] <- tmp
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  merge(counts, out, by = c("study_id", "scenario_id", "method"), all.y = TRUE, sort = FALSE)
}

jss_recovery_predictive_summary <- function(predictive, ledger, conf.level = 0.95) {
  non_variogram <- c("test_log_score_joint", "test_log_score_marginal", "test_log_score_copula", "test_log_score_per_obs")
  identity_cols <- c("study_id", "scenario_id", "method", "replicate", non_variogram)
  unique_scores <- unique(predictive[intersect(identity_cols, names(predictive))])
  score_summary <- jss_recovery_long_metric_summary(unique_scores, ledger, non_variogram, conf.level)
  variogram_summary <- if (all(c("variogram_score", "variogram_p") %in% names(predictive))) {
    jss_recovery_long_metric_summary(predictive, ledger, "variogram_score", conf.level, extra_keys = "variogram_p")
  } else data.frame()
  if (nrow(score_summary) && !"variogram_p" %in% names(score_summary)) score_summary$variogram_p <- NA_real_
  if (nrow(variogram_summary)) {
    missing <- setdiff(names(score_summary), names(variogram_summary))
    for (column in missing) variogram_summary[[column]] <- NA
    missing <- setdiff(names(variogram_summary), names(score_summary))
    for (column in missing) score_summary[[column]] <- NA
  }
  out <- rbind(score_summary, variogram_summary[names(score_summary)])
  rownames(out) <- NULL
  out
}

jss_recovery_paired_method_differences <- function(ledger, fixed, smooth, predictive, diagnostic, conf.level = 0.95) {
  fixed <- jss_recovery_attach_attempt_identity(fixed, ledger)
  smooth <- jss_recovery_attach_attempt_identity(smooth, ledger)
  predictive <- jss_recovery_attach_attempt_identity(predictive, ledger)
  diagnostic <- jss_recovery_attach_attempt_identity(diagnostic, ledger)
  retained_ids <- ledger$attempt_id[ledger$retained %in% TRUE]
  values <- list()
  add_values <- function(x) values[[length(values) + 1L]] <<- x
  retained_ledger <- ledger[ledger$retained %in% TRUE, ]
  add_values(data.frame(retained_ledger[c("study_id", "scenario_id", "method", "replicate")], metric = "runtime_sec", stratum = "overall", value = retained_ledger$elapsed_sec))
  if (nrow(fixed)) {
    x <- fixed[fixed$attempt_id %in% retained_ids, ]
    add_values(data.frame(x[c("study_id", "scenario_id", "method", "replicate")], metric = "fixed_estimation_error", stratum = paste(x$parameter, x$term, sep = ":"), value = x$estimate - x$true_value))
    add_values(data.frame(x[c("study_id", "scenario_id", "method", "replicate")], metric = "fixed_absolute_error", stratum = paste(x$parameter, x$term, sep = ":"), value = abs(x$estimate - x$true_value)))
    add_values(data.frame(x[c("study_id", "scenario_id", "method", "replicate")], metric = "fixed_squared_error", stratum = paste(x$parameter, x$term, sep = ":"), value = (x$estimate - x$true_value)^2))
    add_values(data.frame(x[c("study_id", "scenario_id", "method", "replicate")], metric = "fixed_standard_error", stratum = paste(x$parameter, x$term, sep = ":"), value = x$std_error))
    valid <- is.finite(x$estimate) & is.finite(x$std_error) & is.finite(x$true_value)
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    coverage <- rep(NA_real_, nrow(x))
    coverage[valid] <- as.numeric(x$true_value[valid] >= x$estimate[valid] - z * x$std_error[valid] &
      x$true_value[valid] <= x$estimate[valid] + z * x$std_error[valid])
    add_values(data.frame(x[c("study_id", "scenario_id", "method", "replicate")], metric = "fixed_coverage", stratum = paste(x$parameter, x$term, sep = ":"), value = coverage))
  }
  if (nrow(smooth)) {
    x <- smooth[smooth$attempt_id %in% retained_ids, ]
    add_values(data.frame(x[c("study_id", "scenario_id", "method", "replicate")], metric = "smooth_irmse", stratum = x$parameter, value = x$irmse))
  }
  sources <- list(predictive = predictive, diagnostic = diagnostic)
  for (source_name in names(sources)) {
    x <- sources[[source_name]][sources[[source_name]]$attempt_id %in% retained_ids, , drop = FALSE]
    identity <- c("study_id", "scenario_id", "method", "replicate", "attempt_id", "retry_index", "scenario", "model", "rep")
    numeric_metrics <- setdiff(names(x)[vapply(x, is.numeric, logical(1))], c(identity, "variogram_p", "predictive_nsim", "n", "d"))
    for (metric in numeric_metrics) {
      stratum <- if (identical(source_name, "predictive") && "variogram_p" %in% names(x) && metric == "variogram_score") paste0("p=", x$variogram_p) else "overall"
      add_values(data.frame(x[c("study_id", "scenario_id", "method", "replicate")], metric = paste0(source_name, ":", metric), stratum = stratum, value = x[[metric]]))
    }
  }
  numeric_values <- do.call(rbind, values)
  numeric_values <- numeric_values[is.finite(numeric_values$value), ]

  cell <- unique(ledger[c("study_id", "scenario_id", "method", "replicate")])
  cell$retained <- vapply(seq_len(nrow(cell)), function(i) any(ledger$retained & jss_recovery_cell_key(ledger) == jss_recovery_cell_key(cell[i, ])), logical(1))
  cell$converged <- vapply(seq_len(nrow(cell)), function(i) any(ledger$convergence_eligible & jss_recovery_cell_key(ledger) == jss_recovery_cell_key(cell[i, ])), logical(1))
  cell$execution_success <- vapply(seq_len(nrow(cell)), function(i) any(ledger$execution_success & jss_recovery_cell_key(ledger) == jss_recovery_cell_key(cell[i, ])), logical(1))
  for (metric in c("execution_success", "converged", "retained")) {
    add <- data.frame(cell[1:4], metric = paste0("rate:", metric), stratum = "overall", value = as.numeric(cell[[metric]]))
    numeric_values <- rbind(numeric_values, add)
  }

  cells <- split(numeric_values, interaction(numeric_values[c("study_id", "scenario_id", "metric", "stratum")], drop = TRUE, lex.order = TRUE))
  rows <- list()
  for (df in cells) {
    methods <- sort(unique(ledger$method[ledger$study_id == df$study_id[[1]] & ledger$scenario_id == df$scenario_id[[1]]]))
    if (length(methods) != 2L) next
    collapsed <- stats::aggregate(value ~ method + replicate, df, mean)
    wide <- stats::reshape(collapsed, idvar = "replicate", timevar = "method", direction = "wide")
    columns <- paste0("value.", methods)
    if (!all(columns %in% names(wide))) {
      available <- methods[columns %in% names(wide)]
      unavailable <- methods[!columns %in% names(wide)]
      rows[[length(rows) + 1L]] <- data.frame(
        study_id = df$study_id[[1]], scenario_id = df$scenario_id[[1]], metric = df$metric[[1]], stratum = df$stratum[[1]],
        method_a = methods[[1]], method_b = methods[[2]], contrast = "not_computed_noncommon_estimand",
        orientation = "not_applicable", analysis = "not_applicable_noncommon_estimand",
        planned_pair_denominator = unique(ledger$target_replicates[ledger$study_id == df$study_id[[1]] & ledger$scenario_id == df$scenario_id[[1]]])[[1]],
        paired_denominator = 0L, missing_pair_cells = NA_integer_, estimate = NA_real_, mcse = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        interval_method = paste0("available=", paste(available, collapse = "|"), ";not_applicable=", paste(unavailable, collapse = "|")), stringsAsFactors = FALSE)
      next
    }
    complete <- is.finite(wide[[columns[[1]]]]) & is.finite(wide[[columns[[2]]]])
    differences <- wide[[columns[[1]]]][complete] - wide[[columns[[2]]]][complete]
    conditional <- jss_recovery_mean_stats(differences, conf.level)
    planned <- unique(ledger$target_replicates[ledger$study_id == df$study_id[[1]] & ledger$scenario_id == df$scenario_id[[1]]])
    if (length(planned) != 1L) planned <- length(unique(ledger$replicate[ledger$study_id == df$study_id[[1]]]))
    is_rate <- startsWith(df$metric[[1]], "rate:")
    higher_is_better <- df$metric[[1]] %in% c("fixed_coverage") || grepl("predictive:test_log_score", df$metric[[1]], fixed = TRUE) || is_rate
    orientation <- if (df$metric[[1]] == "fixed_estimation_error") "signed_bias_difference" else if (higher_is_better) "higher_is_better" else "lower_is_better"
    rows[[length(rows) + 1L]] <- data.frame(
      study_id = df$study_id[[1]], scenario_id = df$scenario_id[[1]], metric = df$metric[[1]], stratum = df$stratum[[1]],
      method_a = methods[[1]], method_b = methods[[2]], contrast = "method_a_minus_method_b",
      orientation = orientation,
      analysis = "complete_pair", planned_pair_denominator = planned, paired_denominator = length(differences), missing_pair_cells = planned - length(differences),
      estimate = conditional[["estimate"]], mcse = conditional[["mcse"]], ci_low = conditional[["ci_low"]], ci_high = conditional[["ci_high"]],
      interval_method = if (is_rate) "paired_difference_t_interval" else "paired_seed_t_interval", stringsAsFactors = FALSE
    )
    if (!is_rate) {
      penalty <- if (orientation == "signed_bias_difference") max(abs(df$value), na.rm = TRUE) else if (higher_is_better) min(df$value, na.rm = TRUE) else max(df$value, na.rm = TRUE)
      reps <- seq_len(planned)
      filled <- merge(expand.grid(replicate = reps, method = methods, stringsAsFactors = FALSE), collapsed, by = c("replicate", "method"), all.x = TRUE)
      filled$value[!is.finite(filled$value)] <- penalty
      wide_filled <- stats::reshape(filled, idvar = "replicate", timevar = "method", direction = "wide")
      inclusive <- wide_filled[[columns[[1]]]] - wide_filled[[columns[[2]]]]
      inclusive_stats <- jss_recovery_mean_stats(inclusive, conf.level)
      rows[[length(rows) + 1L]] <- transform(rows[[length(rows)]],
        analysis = "failure_inclusive_penalized", paired_denominator = planned, missing_pair_cells = 0L,
        estimate = inclusive_stats[["estimate"]], mcse = inclusive_stats[["mcse"]], ci_low = inclusive_stats[["ci_low"]], ci_high = inclusive_stats[["ci_high"]],
        interval_method = "paired_seed_t_interval_failure_penalty")
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  rownames(out) <- NULL
  out
}

jss_recovery_failure_sensitivity <- function(ledger, fixed, smooth, predictive, conf.level = 0.95) {
  values <- list()
  fixed <- jss_recovery_attach_attempt_identity(fixed, ledger)
  smooth <- jss_recovery_attach_attempt_identity(smooth, ledger)
  predictive <- jss_recovery_attach_attempt_identity(predictive, ledger)
  fixed$error_sq <- (fixed$estimate - fixed$true_value)^2
  fixed_rep <- jss_recovery_split_apply(fixed[is.finite(fixed$error_sq), ], c("attempt_id"), function(df) {
    data.frame(df[1, c("attempt_id", "study_id", "scenario_id", "method", "replicate")], metric = "fixed_effect_mse", loss = mean(df$error_sq), stringsAsFactors = FALSE)
  })
  values[[length(values) + 1L]] <- fixed_rep
  smooth_rep <- jss_recovery_split_apply(smooth[is.finite(smooth$irmse), ], c("attempt_id"), function(df) {
    data.frame(df[1, c("attempt_id", "study_id", "scenario_id", "method", "replicate")], metric = "mean_smooth_irmse", loss = mean(df$irmse), stringsAsFactors = FALSE)
  })
  values[[length(values) + 1L]] <- smooth_rep
  runtime <- ledger[ledger$retained %in% TRUE & is.finite(ledger$elapsed_sec), c("attempt_id", "study_id", "scenario_id", "method", "replicate", "elapsed_sec")]
  if (nrow(runtime)) values[[length(values) + 1L]] <- data.frame(runtime[1:5], metric = "runtime_sec", loss = runtime$elapsed_sec, stringsAsFactors = FALSE)
  score_cols <- c("attempt_id", "study_id", "scenario_id", "method", "replicate", "test_log_score_per_obs")
  if (all(score_cols %in% names(predictive))) {
    score <- unique(predictive[score_cols])
    score <- score[is.finite(score$test_log_score_per_obs), ]
    values[[length(values) + 1L]] <- data.frame(score[1:5], metric = "negative_test_log_score_per_obs", loss = -score$test_log_score_per_obs, stringsAsFactors = FALSE)
  }
  if (all(c("variogram_score", "variogram_p") %in% names(predictive))) {
    variogram <- predictive[is.finite(predictive$variogram_score), ]
    variogram$metric <- paste0("variogram_score_p", format(variogram$variogram_p, trim = TRUE, scientific = FALSE))
    values[[length(values) + 1L]] <- variogram[c("attempt_id", "study_id", "scenario_id", "method", "replicate", "metric", "variogram_score")]
    names(values[[length(values)]])[7] <- "loss"
  }
  metric_values <- do.call(rbind, values)
  retained_keys <- unique(ledger[ledger$retained %in% TRUE, "attempt_id", drop = FALSE])
  metric_values <- merge(metric_values, retained_keys, by = "attempt_id", all = FALSE)
  method_cells <- unique(ledger[c("study_id", "scenario_id", "method")])
  metric_cells <- merge(method_cells, data.frame(metric = sort(unique(metric_values$metric)), stringsAsFactors = FALSE), by = NULL)
  observed_cells <- unique(metric_values[c("study_id", "scenario_id", "method", "metric")])
  observed_cells$observed <- TRUE
  missing_cells <- merge(metric_cells, observed_cells, by = c("study_id", "scenario_id", "method", "metric"), all.x = TRUE, all.y = FALSE)
  missing_cells <- missing_cells[is.na(missing_cells$observed), c("study_id", "scenario_id", "method", "metric"), drop = FALSE]
  if (nrow(missing_cells)) {
    missing_cells$replicate <- NA_integer_
    missing_cells$attempt_id <- NA_character_
    missing_cells$loss <- NA_real_
    metric_values <- rbind(metric_values, missing_cells[names(metric_values)])
  }
  out <- jss_recovery_split_apply(metric_values, c("study_id", "scenario_id", "method", "metric"), function(df) {
    attempts <- ledger[ledger$study_id == df$study_id[[1]] & ledger$scenario_id == df$scenario_id[[1]] & ledger$method == df$method[[1]], ]
    target <- unique(attempts$target_replicates)
    if (length(target) != 1L) target <- length(unique(attempts$replicate))
    observed <- df$loss[is.finite(df$loss)]
    n_observed <- length(unique(df$replicate[is.finite(df$loss)]))
    n_failed_or_missing <- max(0L, target - n_observed)
    conditional <- jss_recovery_mean_stats(observed, conf.level)
    fallback_losses <- metric_values$loss[metric_values$study_id == df$study_id[[1]] & metric_values$scenario_id == df$scenario_id[[1]] & metric_values$metric == df$metric[[1]] & is.finite(metric_values$loss)]
    penalty <- if (length(observed)) max(observed) else if (length(fallback_losses)) max(fallback_losses) else NA_real_
    penalty_source <- if (length(observed)) "same_method_maximum_observed" else if (length(fallback_losses)) "same_scenario_metric_maximum_observed" else "unavailable_no_observed_loss"
    failure_rate <- jss_recovery_wilson(n_failed_or_missing, target, conf.level)
    data.frame(
      study_id = df$study_id[[1]], scenario_id = df$scenario_id[[1]], method = df$method[[1]], metric = df$metric[[1]],
      orientation = "lower_is_better", attempted_replicates = target, observed_replicates = n_observed,
      failed_or_missing_replicates = n_failed_or_missing,
      conditional_mean = conditional[["estimate"]], conditional_mcse = conditional[["mcse"]],
      conditional_ci_low = conditional[["ci_low"]], conditional_ci_high = conditional[["ci_high"]],
      optimistic_failure_inclusive_mean = sum(observed) / target,
      failure_penalty = penalty,
      failure_penalty_source = penalty_source,
      failure_penalized_mean = if (is.finite(penalty)) (sum(observed) + n_failed_or_missing * penalty) / target else NA_real_,
      failure_rate = failure_rate[["estimate"]], failure_rate_mcse = failure_rate[["mcse"]],
      failure_rate_ci_low = failure_rate[["ci_low"]], failure_rate_ci_high = failure_rate[["ci_high"]],
      penalty_rule = "missing loss uses the same method/metric maximum; an all-failed method uses the same scenario/metric maximum",
      stringsAsFactors = FALSE
    )
  })
  out
}

jss_recovery_t_tail_dependence <- function(rho, degrees_freedom) {
  rho <- pmin(pmax(as.numeric(rho), -1 + 1e-12), 1 - 1e-12)
  degrees_freedom <- as.numeric(degrees_freedom)
  2 * stats::pt(-sqrt((degrees_freedom + 1) * (1 - rho) / (1 + rho)), df = degrees_freedom + 1)
}

jss_recovery_t_shape_summary <- function(fixed, ledger, conf.level = 0.95, gaussian_limit_zeta = jss_recovery_gaussian_limit_zeta()) {
  fixed <- jss_recovery_attach_attempt_identity(fixed, ledger)
  t_studies <- unique(ledger$study_id[ledger$margin_family == "BCPE" & ledger$copula_family %in% c("t", "Student t", "Student's t")])
  if (length(t_studies) != 1L) stop("Expected one BCPE/Student t study for shape recovery.", call. = FALSE)
  x <- fixed[
    fixed$study_id == t_studies[[1]] & fixed$method == "gamlss.longitudinal" &
      fixed$term == "intercept" & fixed$parameter %in% c("theta", "zeta"),
  ]
  theta <- x[x$parameter == "theta", c("attempt_id", "study_id", "scenario_id", "method", "replicate", "estimate", "std_error", "true_value")]
  zeta <- x[x$parameter == "zeta", c("attempt_id", "study_id", "scenario_id", "method", "replicate", "estimate", "std_error", "true_value")]
  names(theta)[6:8] <- c("theta_estimate", "theta_se", "theta_true")
  names(zeta)[6:8] <- c("zeta_estimate", "zeta_se", "zeta_true")
  paired <- merge(theta, zeta, by = c("attempt_id", "study_id", "scenario_id", "method", "replicate"), all = FALSE)
  retained <- ledger[ledger$retained %in% TRUE, "attempt_id", drop = FALSE]
  paired <- merge(paired, retained, by = "attempt_id", all = FALSE)
  paired$rho_estimate <- tanh(paired$theta_estimate)
  paired$rho_true <- tanh(paired$theta_true)
  paired$df_estimate_overflow <- is.finite(paired$zeta_estimate) & paired$zeta_estimate > log(.Machine$double.xmax)
  paired$df_estimate <- 2 + exp(paired$zeta_estimate)
  paired$df_true <- 2 + exp(paired$zeta_true)
  paired$tail_dependence_estimate <- jss_recovery_t_tail_dependence(paired$rho_estimate, paired$df_estimate)
  paired$tail_dependence_true <- jss_recovery_t_tail_dependence(paired$rho_true, paired$df_true)
  paired$zeta_interval_width <- 2 * stats::qnorm(1 - (1 - conf.level) / 2) * paired$zeta_se
  zcrit <- stats::qnorm(1 - (1 - conf.level) / 2)
  paired$df_interval_lower_link <- paired$zeta_estimate - zcrit * paired$zeta_se
  paired$df_interval_upper_link <- paired$zeta_estimate + zcrit * paired$zeta_se
  paired$df_interval_width_overflow <- is.finite(paired$df_interval_upper_link) &
    paired$df_interval_upper_link > log(.Machine$double.xmax)
  paired$df_interval_width <- exp(paired$df_interval_upper_link) - exp(paired$df_interval_lower_link)
  # Inf - Inf yields NaN; the pre-transform flag keeps this as an overflow
  # event rather than misclassifying it as an ordinary missing value.
  paired$df_interval_width[paired$df_interval_width_overflow] <- Inf
  paired$near_gaussian_limit <- paired$zeta_estimate >= gaussian_limit_zeta

  quantities <- list(
    zeta_link = c("zeta_estimate", "zeta_true"),
    degrees_of_freedom = c("df_estimate", "df_true"),
    rho = c("rho_estimate", "rho_true"),
    lower_tail_dependence = c("tail_dependence_estimate", "tail_dependence_true")
  )
  rows <- lapply(names(quantities), function(quantity) {
    columns <- quantities[[quantity]]
    error <- paired[[columns[[1]]]] - paired[[columns[[2]]]]
    bias <- jss_recovery_mean_stats(error, conf.level)
    rmse <- jss_recovery_rmse_stats(error, conf.level)
    estimate <- jss_recovery_mean_stats(paired[[columns[[1]]]], conf.level)
    estimate_sd <- jss_recovery_sd_stats(paired[[columns[[1]]]], conf.level)
    data.frame(
      study_id = t_studies[[1]], scenario_id = unique(paired$scenario_id)[[1]], method = "gamlss.longitudinal",
      reference_profile = "intercept predictor: t=x1=x2=0; fitted smooth centering retained",
      quantity = quantity, n_metric = as.integer(rmse[["n_metric"]]), true_value_mean = mean(paired[[columns[[2]]]], na.rm = TRUE),
      mean_estimate = estimate[["estimate"]], bias = bias[["estimate"]], bias_mcse = bias[["mcse"]], bias_ci_low = bias[["ci_low"]], bias_ci_high = bias[["ci_high"]],
      rmse = rmse[["estimate"]], rmse_mcse = rmse[["mcse"]], rmse_ci_low = rmse[["ci_low"]], rmse_ci_high = rmse[["ci_high"]],
      empirical_sd = estimate_sd[["estimate"]], empirical_sd_mcse = estimate_sd[["mcse"]],
      empirical_sd_ci_low = estimate_sd[["ci_low"]], empirical_sd_ci_high = estimate_sd[["ci_high"]],
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  zeta_width <- jss_recovery_mean_stats(paired$zeta_interval_width, conf.level)
  df_width <- jss_recovery_mean_stats(paired$df_interval_width, conf.level)
  zeta_width_median <- jss_recovery_bootstrap_quantile(paired$zeta_interval_width, 0.5, conf.level, 20260941L)
  df_width_median <- jss_recovery_bootstrap_quantile(paired$df_interval_width, 0.5, conf.level, 20260942L)
  for (prefix in c("zeta", "df")) {
    width <- if (prefix == "zeta") zeta_width else df_width
    median <- if (prefix == "zeta") zeta_width_median else df_width_median
    summary[[paste0("mean_", prefix, "_interval_width")]] <- width[["estimate"]]
    summary[[paste0("mean_", prefix, "_interval_width_mcse")]] <- width[["mcse"]]
    summary[[paste0("mean_", prefix, "_interval_width_ci_low")]] <- width[["ci_low"]]
    summary[[paste0("mean_", prefix, "_interval_width_ci_high")]] <- width[["ci_high"]]
    summary[[paste0("median_", prefix, "_interval_width")]] <- median[["estimate"]]
    summary[[paste0("median_", prefix, "_interval_width_mcse")]] <- median[["mcse"]]
    summary[[paste0("median_", prefix, "_interval_width_ci_low")]] <- median[["ci_low"]]
    summary[[paste0("median_", prefix, "_interval_width_ci_high")]] <- median[["ci_high"]]
  }
  boundary <- jss_recovery_wilson(sum(paired$near_gaussian_limit, na.rm = TRUE), sum(!is.na(paired$near_gaussian_limit)), conf.level)
  event_counts <- function(value, overflow = rep(FALSE, length(value))) c(
    total = length(value), finite = sum(is.finite(value) & !overflow), missing = sum(is.na(value) & !overflow),
    infinite = sum(is.infinite(value) & !overflow), overflow = sum(overflow, na.rm = TRUE)
  )
  zeta_events <- event_counts(paired$zeta_estimate)
  df_events <- event_counts(paired$df_estimate, paired$df_estimate_overflow)
  zeta_width_events <- event_counts(paired$zeta_interval_width)
  df_width_overflow <- paired$df_interval_width_overflow
  df_width_events <- event_counts(paired$df_interval_width, df_width_overflow)
  for (name in names(zeta_events)) summary[[paste0("zeta_", name, "_count")]] <- unname(zeta_events[[name]])
  for (name in names(df_events)) summary[[paste0("df_", name, "_count")]] <- unname(df_events[[name]])
  for (name in names(zeta_width_events)) summary[[paste0("zeta_interval_width_", name, "_count")]] <- unname(zeta_width_events[[name]])
  for (name in names(df_width_events)) summary[[paste0("df_interval_width_", name, "_count")]] <- unname(df_width_events[[name]])
  summary$near_gaussian_zeta_threshold <- gaussian_limit_zeta
  summary$near_gaussian_threshold_contract <- "registered descriptive threshold: fitted zeta >= 7.5; event retained, not discarded"
  summary$near_gaussian_limit_rate <- boundary[["estimate"]]
  summary$near_gaussian_limit_mcse <- boundary[["mcse"]]
  summary$near_gaussian_limit_ci_low <- boundary[["ci_low"]]
  summary$near_gaussian_limit_ci_high <- boundary[["ci_high"]]
  summary
}

jss_recovery_copula_family_consistent <- function(code, family) {
  code <- as.character(code)
  family <- tolower(trimws(as.character(family)))
  (code == "N" & family %in% c("gaussian", "normal")) |
    (code == "C" & family == "clayton") |
    (code == "t" & family %in% c("t", "student t", "student's t")) |
    (code == "F" & family == "frank") |
    (code == "G" & family == "gumbel") |
    (code == "J" & family == "joe")
}

jss_recovery_bundle_status <- function(ledger) {
  routes <- unique(ledger[c("margin_family", "copula_code")])
  correct_routes <- nrow(routes) == 2L && any(routes$margin_family == "BCPE" & routes$copula_code == "t") &&
    any(routes$margin_family == "NBI" & routes$copula_code == "C")
  methods <- split(ledger$method, ledger$margin_family)
  correct_methods <- setequal(unique(methods$BCPE), c("gamlss.longitudinal", "gamlss2")) &&
    setequal(unique(methods$NBI), c("gamlss", "ours_rs_joint"))
  production <- all(ledger$evidence_status == "post_phase1_production") &&
    all(ledger$phase1_contract_version == jss_recovery_phase1_contract_version()) &&
    all(jss_recovery_copula_family_consistent(ledger$copula_code, ledger$copula_family)) &&
    all((ledger$margin_family == "BCPE" & ledger$runner_contract_version == "bcpe-t-main-recovery-2026-09-02.7") |
      (ledger$margin_family == "NBI" & ledger$runner_contract_version == "nbi-clayton-main-recovery-2026-09-02.5")) &&
    all(!is.na(ledger$runner_settings_signature) & nzchar(ledger$runner_settings_signature) &
      !grepl("legacy|unrecorded", ledger$runner_settings_signature, ignore.case = TRUE)) &&
    all(vapply(c("runner_settings_sha256", "runner_sha256", "package_source_path", "package_version", "package_source_sha256"), function(name)
      all(!is.na(ledger[[name]]) & nzchar(ledger[[name]]) & !grepl("legacy|unrecorded", ledger[[name]], ignore.case = TRUE)), logical(1)))
  if (production && correct_routes && correct_methods) "authoritative_post_phase1" else "legacy_reconciliation_not_authoritative"
}

jss_recovery_metric_contract <- function(margin_family, method) {
  if (margin_family == "BCPE" && method %in% c("gamlss.longitudinal", "gamlss2")) {
    parameters <- if (method == "gamlss2") c("mu", "sigma", "nu", "tau") else c("mu", "sigma", "nu", "tau", "theta", "zeta")
    return(list(
      fixed = as.vector(outer(parameters, c("intercept", "x1", "x2", "t"), paste, sep = "\r")),
      smooth = if (method == "gamlss.longitudinal") c("mu", "sigma", "theta") else c("mu", "sigma"),
      variogram_p = c(0.5, 2),
      predictive = c("test_log_score_joint", "test_log_score_marginal", "test_log_score_copula", "test_log_score_per_obs", "variogram_score", "predictive_nsim"),
      diagnostic = c("logLik", "rosenblatt_ks", "rosenblatt_cvm", "abs_rosenblatt_lag1_cor", "abs_rosenblatt_normal_lag1_cor", "rosenblatt_mean_abs_time_mean", "rosenblatt_normal_mean_abs_time_mean")
    ))
  }
  if (margin_family == "NBI" && method == "gamlss") {
    return(list(fixed = c("mu\rintercept", "mu\rx1", "mu\rx2", "mu\rtime_scaled", "sigma\rintercept", "sigma\rx1", "sigma\rtime_scaled"), smooth = c("mu", "sigma"), variogram_p = c(0.5, 2), predictive = c("test_log_score_joint", "test_log_score_marginal", "test_log_score_copula", "test_log_score_per_obs", "variogram_score", "predictive_nsim"), diagnostic = c("logLik", "rosenblatt_ks", "rosenblatt_cvm", "abs_rosenblatt_lag1_cor", "abs_rosenblatt_normal_lag1_cor", "rosenblatt_mean_abs_time_mean", "rosenblatt_normal_mean_abs_time_mean")))
  }
  if (margin_family == "NBI" && startsWith(method, "ours_rs_")) {
    return(list(fixed = c("mu\rintercept", "mu\rx1", "mu\rx2", "mu\rtime_scaled", "sigma\rintercept", "sigma\rx1", "sigma\rtime_scaled", "theta\rintercept", "theta\rtime_scaled"), smooth = c("mu", "sigma", "theta"), variogram_p = c(0.5, 2), predictive = c("test_log_score_joint", "test_log_score_marginal", "test_log_score_copula", "test_log_score_per_obs", "variogram_score", "predictive_nsim"), diagnostic = c("logLik", "rosenblatt_ks", "rosenblatt_cvm", "abs_rosenblatt_lag1_cor", "abs_rosenblatt_normal_lag1_cor", "rosenblatt_mean_abs_time_mean", "rosenblatt_normal_mean_abs_time_mean")))
  }
  NULL
}

jss_recovery_bcpe_settings_fields <- function() c(
  "runner_contract_version", "phase1_contract_version", "n_fits", "rep_ids", "max_attempts_per_fit",
  "n_cores", "verbose_fits", "rscript_path", "rscript_sha256", "rscript_version", "parallel_backend",
  "parallel_scheduler", "psock_setup_strategy", "worker_count_rule", "rscript_args", "omp_num_threads",
  "openblas_num_threads", "mkl_num_threads", "veclib_maximum_threads", "blis_num_threads",
  "n_subject", "n_time", "smooth_k", "theta_effect_mode", "theta_binary_effect", "theta_binary_prob",
  "compute_se", "save_fits", "vcov_method_longitudinal", "include_dlcopdpar", "optim_method",
  "max_outer_iter", "max_inner_iter", "max_elapsed_sec", "outer_stop_crit", "inner_stop_crit",
  "use_backtracking", "backtracking_max_halves", "start_step_size", "step_adjustment", "max_steps",
  "lambda_start", "warm_start_joint", "cg_max_delta", "cg_armijo_c1", "cg_max_stall", "cg_update_lambda",
  "cg_line_search", "cg_max_line_search_evals", "cg_gradient_method", "cg_zeta_hessian", "cg_lambda_update_every",
  "cg_max_lambda_updates", "cg_raw_loglik_drop_tol", "compute_predictive_scores", "predictive_nsim", "variogram_p_values",
  "dgp_global_seed", "dgp_seed_registry", "dgp_covariates", "dgp_dependence", "dgp_uniform_clip", "smooth_grid_points",
  "dgp_true_beta_mu", "dgp_true_beta_sigma", "dgp_true_beta_nu", "dgp_true_beta_tau", "dgp_true_beta_theta", "dgp_true_beta_zeta",
  "dgp_smooth_mu", "dgp_smooth_sigma", "dgp_smooth_theta"
)

jss_recovery_settings_signature <- function(settings, margin_family) {
  one <- function(name) as.character(settings[[name]][[1]])
  if (margin_family == "BCPE") {
    fields <- jss_recovery_bcpe_settings_fields()
    if (!all(fields %in% names(settings))) stop("BCPE canonical settings row is incomplete.", call. = FALSE)
    return(paste(fields, vapply(fields, one, character(1)), sep = "=", collapse = "\n"))
  } else if (margin_family == "NBI") {
    values <- c(one("runner_contract_version"), one("phase1_contract_version"), one("n_subject"),
      gsub("[|]", ",", one("times")), one("reps"), one("base_seed"), one("sigma_signal_multiplier"),
      gsub("[|]", ",", one("engines")), one("max_elapsed_sec"), one("max_outer_iter"), one("max_inner_iter"),
      one("start_step_size"), one("step_adjustment_env"), one("lambda_start"), one("rs_update_lambda"), one("warm_start_joint_iter"),
      one("compute_se"), one("vcov_method_longitudinal"), one("compute_predictive_scores"), one("predictive_nsim"),
      gsub("[|]", ",", one("variogram_p_values")), one("theta_intercept"), one("theta_time_coef"), one("max_attempts_per_fit"),
      one("runtime_n_cores"), one("runtime_backend"), one("rscript_path"), one("rscript_sha256"), one("rscript_version"))
  } else stop("Unknown runner settings family.", call. = FALSE)
  paste(values, collapse = "|")
}

jss_recovery_attempt_output_contract <- function(fixed, smooth, predictive, diagnostic, contract) {
  fixed_keys <- if (all(c("parameter", "term") %in% names(fixed))) paste(fixed$parameter, fixed$term, sep = "\r") else character()
  fixed_schema <- nrow(fixed) == length(contract$fixed) && !anyDuplicated(fixed_keys) && identical(sort(fixed_keys), sort(contract$fixed))
  fixed_values <- fixed_schema && all(is.finite(fixed$estimate)) && all(is.finite(fixed$true_value)) &&
    all(is.na(fixed$std_error) | (is.finite(fixed$std_error) & fixed$std_error >= 0))
  # Missing standard errors are allowed only as an explicitly denominated,
  # named inferential status. They are never silently treated as zero rows.
  if (fixed_values && anyNA(fixed$std_error)) {
    fixed_values <- all(c("inference_status", "inference_denominator") %in% names(fixed)) &&
      all(!is.na(fixed$inference_status[is.na(fixed$std_error)]) & nzchar(fixed$inference_status[is.na(fixed$std_error)])) &&
      all(fixed$inference_status[is.na(fixed$std_error)] != "available") &&
      all(is.finite(fixed$inference_denominator[is.na(fixed$std_error)]) & fixed$inference_denominator[is.na(fixed$std_error)] == 0)
  }
  smooth_schema <- nrow(smooth) == length(contract$smooth) && !anyDuplicated(smooth$parameter) && identical(sort(as.character(smooth$parameter)), sort(contract$smooth))
  smooth_values <- smooth_schema && all(is.finite(smooth$bias_abs_integrated) & smooth$bias_abs_integrated >= 0) && all(is.finite(smooth$irmse) & smooth$irmse >= 0)
  predictive_schema <- nrow(predictive) == length(contract$variogram_p) && !anyDuplicated(as.numeric(predictive$variogram_p)) &&
    identical(sort(as.numeric(predictive$variogram_p)), sort(contract$variogram_p)) && all(contract$predictive %in% names(predictive))
  predictive_values <- predictive_schema && all(vapply(contract$predictive, function(name) all(is.finite(as.numeric(predictive[[name]]))), logical(1))) &&
    all(predictive$variogram_score >= 0) && all(predictive$predictive_nsim > 0)
  diagnostic_schema <- nrow(diagnostic) == 1L && all(contract$diagnostic %in% names(diagnostic))
  diagnostic_values <- diagnostic_schema && all(vapply(contract$diagnostic, function(name) all(is.finite(as.numeric(diagnostic[[name]]))), logical(1))) &&
    diagnostic$rosenblatt_ks >= 0 && diagnostic$rosenblatt_ks <= 1 && diagnostic$rosenblatt_cvm >= 0 &&
    diagnostic$abs_rosenblatt_lag1_cor >= 0 && diagnostic$abs_rosenblatt_lag1_cor <= 1 &&
    diagnostic$abs_rosenblatt_normal_lag1_cor >= 0 && diagnostic$abs_rosenblatt_normal_lag1_cor <= 1 &&
    diagnostic$rosenblatt_mean_abs_time_mean >= 0 && diagnostic$rosenblatt_mean_abs_time_mean <= 1 &&
    diagnostic$rosenblatt_normal_mean_abs_time_mean >= 0
  result <- c(fixed = fixed_schema && fixed_values, smooth = smooth_schema && smooth_values,
    predictive = predictive_schema && predictive_values, diagnostic = diagnostic_schema && diagnostic_values)
  attr(result, "reason") <- paste(names(result)[!result], collapse = "|")
  result
}

jss_recovery_validate <- function(ledger, fixed, smooth, predictive, diagnostic) {
  fixed <- jss_recovery_attach_attempt_identity(fixed, ledger)
  smooth <- jss_recovery_attach_attempt_identity(smooth, ledger)
  predictive <- jss_recovery_attach_attempt_identity(predictive, ledger)
  diagnostic <- jss_recovery_attach_attempt_identity(diagnostic, ledger)
  checks <- list()
  add <- function(check_id, pass, detail) {
    checks[[length(checks) + 1L]] <<- data.frame(check_id = check_id, pass = isTRUE(pass), detail = detail, stringsAsFactors = FALSE)
  }
  duplicate_attempt <- duplicated(ledger$attempt_id)
  add("unique_attempt_ids", !any(duplicate_attempt), paste(sum(duplicate_attempt), "duplicate attempt IDs"))
  metadata_columns <- c("n_subjects", "n_time", "margin_family", "copula_family", "copula_code", "target_replicates", "evidence_status", "phase1_contract_version", "runner_settings_signature", "runner_settings_sha256", "runner_sha256", "package_source_path", "package_version", "package_source_sha256")
  scenario_groups <- split(ledger, interaction(ledger[c("study_id", "scenario_id")], drop = TRUE, lex.order = TRUE))
  metadata_consistent <- all(vapply(scenario_groups, function(df) all(vapply(df[metadata_columns], function(x) length(unique(x)) == 1L, logical(1))), logical(1)))
  add("scenario_metadata_consistent", metadata_consistent, "n/T/family/target/status must be single-valued within scenario")
  retry_groups <- split(ledger, interaction(ledger[c("study_id", "scenario_id", "method", "replicate")], drop = TRUE, lex.order = TRUE))
  retry_valid <- all(vapply(retry_groups, function(df) identical(sort(df$retry_index), seq_len(nrow(df))), logical(1)))
  add("retry_indices_contiguous", retry_valid, "retry indices must be unique and contiguous from one")
  add("one_retained_attempt_per_method_replicate", all(vapply(retry_groups, function(df) sum(df$retained) <= 1L, logical(1))), "at most one retained retry")
  design <- jss_recovery_design_table(ledger)
  add("two_registered_scenarios", nrow(design) == 2L, paste(nrow(design), "scenario rows"))
  add("n_is_500", all(design$n_subjects == 500L), paste(unique(design$n_subjects), collapse = ","))
  add("T_is_4", all(design$n_time == 4L), paste(unique(design$n_time), collapse = ","))
  add("R_is_100", all(design$attempted_replicates == 100L & design$target_replicates == 100L), paste(design$attempted_replicates, collapse = ","))
  expected_reps <- seq_len(100L)
  by_method <- split(ledger, interaction(ledger[c("study_id", "scenario_id", "method")], drop = TRUE, lex.order = TRUE))
  add("method_replicate_cardinality", all(vapply(by_method, function(x) identical(sort(unique(x$replicate)), expected_reps), logical(1))), "every method must contain replicate IDs 1:100")
  method_sets <- split(ledger$method, ledger$margin_family)
  add("registered_method_sets", setequal(unique(method_sets$BCPE), c("gamlss.longitudinal", "gamlss2")) && setequal(unique(method_sets$NBI), c("gamlss", "ours_rs_joint")), "BCPE and NBI studies must contain their two registered comparators")
  planned_groups <- unique(ledger[c("study_id", "scenario_id", "method", "target_replicates")])
  add("exactly_400_planned_method_replicate_cells", sum(as.integer(planned_groups$target_replicates)) == 400L,
    paste("planned_cells=", sum(as.integer(planned_groups$target_replicates)), "; actual_attempt_rows=", nrow(ledger), sep = ""))
  route_rows <- unique(ledger[c("margin_family", "copula_code", "evidence_status")])
  route_ok <- (route_rows$margin_family == "BCPE" & route_rows$copula_code == "t") |
    (route_rows$margin_family == "NBI" & route_rows$copula_code == "C") |
    route_rows$evidence_status == "legacy_pre_phase1_reconciliation"
  add("registered_routes_or_explicit_legacy", all(route_ok), paste(route_rows$margin_family, route_rows$copula_code, route_rows$evidence_status, collapse = ";"))
  add("copula_code_family_consistent", all(jss_recovery_copula_family_consistent(ledger$copula_code, ledger$copula_family)), "copula code and declared family must agree")
  production_rows <- ledger$evidence_status == "post_phase1_production"
  add("production_phase1_contract_current", all(!production_rows | ledger$phase1_contract_version == jss_recovery_phase1_contract_version()), "post-Phase1 attempts must declare the current contract")
  expected_runner <- ifelse(ledger$margin_family == "BCPE", "bcpe-t-main-recovery-2026-09-02.7",
    ifelse(ledger$margin_family == "NBI", "nbi-clayton-main-recovery-2026-09-02.5", NA_character_))
  add("production_runner_contract_current", all(!production_rows | (!is.na(expected_runner) & ledger$runner_contract_version == expected_runner)), "post-Phase1 attempts must use the current route runner")
  signature_ok <- !is.na(ledger$runner_settings_signature) & nzchar(ledger$runner_settings_signature) & !grepl("legacy|unrecorded", ledger$runner_settings_signature, ignore.case = TRUE)
  add("production_settings_signature_recorded", all(!production_rows | signature_ok), "post-Phase1 attempts require a nonlegacy settings signature")
  identity_fields <- c("runner_settings_sha256", "runner_sha256", "package_source_path", "package_version", "package_source_sha256")
  identity_ok <- vapply(identity_fields, function(name) all(!production_rows | (!is.na(ledger[[name]]) & nzchar(ledger[[name]]) & !grepl("legacy|unrecorded", ledger[[name]], ignore.case = TRUE))), logical(1))
  add("production_source_identity_recorded", all(identity_ok), paste(identity_fields[!identity_ok], collapse = ","))
  key <- jss_recovery_cell_key
  ledger_keys <- unique(ledger$attempt_id)
  retained_keys <- unique(ledger$attempt_id[ledger$retained])
  add("fixed_rows_reference_attempts", all(fixed$attempt_id %in% ledger_keys), paste(nrow(fixed), "fixed rows"))
  add("smooth_rows_reference_attempts", all(smooth$attempt_id %in% ledger_keys), paste(nrow(smooth), "smooth rows"))
  add("predictive_rows_reference_attempts", all(predictive$attempt_id %in% ledger_keys), paste(nrow(predictive), "predictive rows"))
  add("diagnostic_rows_reference_attempts", all(diagnostic$attempt_id %in% ledger_keys), paste(nrow(diagnostic), "diagnostic rows"))
  for (name in c("fixed", "smooth", "predictive", "diagnostic")) {
    object <- get(name)
    object_keys <- unique(object$attempt_id)
    add(paste0(name, "_method_replicate_cardinality"), all(retained_keys %in% object_keys), paste(length(object_keys), "output attempts including failures vs", length(retained_keys), "retained attempts"))
  }
  add("one_diagnostic_row_per_method_replicate", !anyDuplicated(diagnostic$attempt_id), paste(nrow(diagnostic), "diagnostic rows; retry-aware attempt identity"))
  retained_rows <- unique(ledger[ledger$retained, c("attempt_id", "study_id", "scenario_id", "margin_family", "method", "replicate")])
  contract_results <- lapply(seq_len(nrow(retained_rows)), function(i) {
    row <- retained_rows[i, ]
    contract <- jss_recovery_metric_contract(row$margin_family, row$method)
    select <- function(x) x$attempt_id == row$attempt_id
    if (is.null(contract)) return(c(fixed = FALSE, smooth = FALSE, predictive = FALSE, diagnostic = FALSE))
    fixed_i <- fixed[select(fixed), , drop = FALSE]
    smooth_i <- smooth[select(smooth), , drop = FALSE]
    predictive_i <- predictive[select(predictive), , drop = FALSE]
    diagnostic_i <- diagnostic[select(diagnostic), , drop = FALSE]
    jss_recovery_attempt_output_contract(fixed_i, smooth_i, predictive_i, diagnostic_i, contract)
  })
  contract_matrix <- if (length(contract_results)) do.call(rbind, contract_results) else matrix(FALSE, nrow = 1L, ncol = 4L, dimnames = list(NULL, c("fixed", "smooth", "predictive", "diagnostic")))
  for (metric_name in colnames(contract_matrix)) {
    add(paste0(metric_name, "_exact_schema_cardinality"), all(contract_matrix[, metric_name]), paste(nrow(retained_rows), "retained method/replicate cells"))
  }
  if (all(c("has_fixed_output", "has_smooth_output", "has_predictive_output", "has_diagnostic_output") %in% names(ledger))) {
    retained_complete <- with(ledger, !retained | (has_fixed_output & has_smooth_output & has_predictive_output & has_diagnostic_output))
    add("retained_attempts_have_complete_outputs", all(retained_complete), paste(sum(ledger$retained), "retained attempts"))
  }
  out <- do.call(rbind, checks)
  rownames(out) <- NULL
  out
}

jss_recovery_input_provenance <- function(paths, conf.level = 0.95, repo_root = normalizePath(".", winslash = "/", mustWork = TRUE)) {
  normalized <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  relative <- ifelse(startsWith(tolower(normalized), paste0(tolower(root), "/")), substring(normalized, nchar(root) + 2L), paste(basename(dirname(normalized)), basename(normalized), sep = "/"))
  data.frame(
    schema_version = jss_recovery_schema_version(),
    producer_sha256 = jss_recovery_producer_sha256(repo_root),
    analysis_config = paste0("conf.level=", format(conf.level, digits = 17)),
    input_id = names(paths),
    path = relative,
    sha256 = unname(jss_recovery_sha256(normalized)),
    bytes = unname(file.info(normalized)$size),
    stringsAsFactors = FALSE
  )
}

jss_recovery_production_provenance <- function(repo_root, ledger) {
  description <- read.dcf(file.path(repo_root, "DESCRIPTION"))
  git <- jss_recovery_git_identity(repo_root)
  data.frame(
    key = c("schema_version", "producer_sha256", "bundle_status", "phase1_contract_version", "package_path", "package_version", "package_source_sha256", "git_sha", "git_state", "tracked_worktree_sha256", "r_version", "platform", "os"),
    value = c(jss_recovery_schema_version(), jss_recovery_producer_sha256(repo_root), jss_recovery_bundle_status(ledger), jss_recovery_phase1_contract_version(),
      ".", description[1, "Version"], jss_recovery_package_source_sha256(repo_root),
      git[["git_sha"]], git[["git_state"]], git[["tracked_worktree_sha256"]], R.version.string, R.version$platform, paste(Sys.info()[c("sysname", "release", "machine")], collapse = "|")),
    stringsAsFactors = FALSE
  )
}

jss_build_main_recovery_evidence <- function(input_root, output_dir, resume = TRUE, conf.level = 0.95,
                                             bcpe_dir = file.path(input_root, "bcpe-t"),
                                             nbi_dir = file.path(input_root, "nbi-clayton"),
                                             repo_root = normalizePath(".", winslash = "/", mustWork = TRUE)) {
  paths <- c(
    bcpe_log = file.path(bcpe_dir, "fit_run_log.csv"),
    bcpe_fixed = file.path(bcpe_dir, "fixed_effects_by_rep.csv"),
    bcpe_smooth = file.path(bcpe_dir, "smooth_integrated_metrics.csv"),
    bcpe_predictive = file.path(bcpe_dir, "predictive_scores_by_rep.csv"),
    bcpe_diagnostic = file.path(bcpe_dir, "joint_distribution_metrics_by_rep.csv"),
    nbi_log = file.path(nbi_dir, "nbi_sigma_compare_logs.csv"),
    nbi_fixed = file.path(nbi_dir, "fixed_effects_by_rep.csv"),
    nbi_smooth = file.path(nbi_dir, "smooth_integrated_metrics.csv"),
    nbi_predictive = file.path(nbi_dir, "predictive_scores_by_rep.csv"),
    nbi_diagnostic = file.path(nbi_dir, "joint_distribution_metrics_by_rep.csv")
  )
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing main-recovery input(s): ", paste(names(missing), collapse = ", "), call. = FALSE)
  control_paths <- c(
    producer = file.path(repo_root, "paper", "R", "main-recovery-evidence.R"),
    evidence_runner = file.path(repo_root, "paper", "scripts", "final-simulations", "main-recovery", "run_main_recovery_evidence.R"),
    bcpe_runner = file.path(repo_root, "paper", "scripts", "final-simulations", "bcpe-t", "simulation_bcpe_t_gamlss_comparison.R"),
    nbi_runner = file.path(repo_root, "paper", "scripts", "final-simulations", "nbi-clayton", "compare_gamlss_ours_nbi_sigma_smooth.R"),
    package_description = file.path(repo_root, "DESCRIPTION"),
    bcpe_settings = file.path(bcpe_dir, "simulation_settings.csv"),
    nbi_settings = file.path(nbi_dir, "nbi_sigma_compare_settings.csv")
  )
  package_sources <- jss_recovery_package_source_files(repo_root)
  names(package_sources) <- paste0("package_source_", seq_along(package_sources), "_", basename(package_sources))
  control_paths <- c(control_paths, package_sources)
  provenance_paths <- c(paths, control_paths[file.exists(control_paths)])
  provenance <- jss_recovery_input_provenance(provenance_paths, conf.level, repo_root)
  expected_outputs <- file.path(output_dir, c(
    "attempt_metadata.csv", "design_table.csv", "attempt_status_summary.csv", "failure_reason_summary.csv",
    "fixed_parameter_recovery.csv", "smooth_recovery.csv", "runtime_summary.csv", "predictive_metrics.csv",
    "diagnostic_metrics.csv", "weak_t_copula_shape_recovery.csv", "failure_inclusive_sensitivity.csv",
    "paired_method_differences.csv", "fixed_by_attempt.csv", "smooth_by_attempt.csv", "predictive_by_attempt.csv", "diagnostic_by_attempt.csv", "runner_settings_identity.csv", "runner_settings_bcpe.csv", "runner_settings_nbi.csv", "control_source_manifest.csv",
    "bundle_status.csv", "production_provenance.csv", "evidence_validation.csv", "input_provenance.csv", "output_manifest.csv", "bundle_checkpoint.csv"
  ))
  if (isTRUE(resume) && jss_recovery_bundle_is_current(output_dir, expected_outputs, provenance)) {
    return(structure(expected_outputs, resumed = TRUE))
  }

  read <- function(path) utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  ledger <- jss_recovery_attempt_ledger(read(paths[["bcpe_log"]]), read(paths[["nbi_log"]]))
  identity_fields <- c("runner_settings_signature", "runner_settings_sha256", "runner_sha256", "package_source_path", "package_version", "package_source_sha256")
  settings_by_margin <- lapply(c(BCPE = control_paths[["bcpe_settings"]], NBI = control_paths[["nbi_settings"]]), function(path) if (file.exists(path)) read(path) else data.frame())
  settings_identity <- do.call(rbind, lapply(names(settings_by_margin), function(margin) {
    source <- settings_by_margin[[margin]]
    row <- data.frame(margin_family = margin, study_id = unique(ledger$study_id[ledger$margin_family == margin])[[1]], stringsAsFactors = FALSE)
    if (nrow(source) > 1L) stop("Expected at most one runner-settings row for ", margin, ".", call. = FALSE)
    for (name in identity_fields) row[[name]] <- if (nrow(source)) as.character(jss_recovery_column(source, name, "legacy_unrecorded")[[1]]) else unique(as.character(ledger[[name]][ledger$margin_family == margin]))[[1]]
    row
  }))
  production <- ledger$evidence_status == "post_phase1_production"
  if (any(production)) {
    jss_recovery_require_clean_checkout(jss_recovery_git_identity(repo_root))
    runner_path <- c(BCPE = control_paths[["bcpe_runner"]], NBI = control_paths[["nbi_runner"]])
    for (margin in names(settings_by_margin)) {
      rows <- ledger$margin_family == margin & production
      settings_row <- settings_by_margin[[margin]]
      if (nrow(settings_row) != 1L || !all(identity_fields %in% names(settings_row))) stop("Runner settings identity is incomplete for ", margin, ".", call. = FALSE)
      if (any(vapply(identity_fields, function(name) !identical(unique(as.character(ledger[[name]][rows])), as.character(settings_row[[name]][[1]])), logical(1)))) stop("Attempt/settings identity mismatch for ", margin, ".", call. = FALSE)
      canonical_signature <- tryCatch(jss_recovery_settings_signature(settings_row, margin), error = function(e) NA_character_)
      if (is.na(canonical_signature) || !identical(canonical_signature, as.character(settings_row$runner_settings_signature[[1]])) ||
          !identical(digest::digest(canonical_signature, algo = "sha256", serialize = FALSE), as.character(settings_row$runner_settings_sha256[[1]]))) stop("Canonical runner-settings signature mismatch for ", margin, ".", call. = FALSE)
      git_identity <- jss_recovery_git_identity(repo_root)
      if (!all(c("git_sha", "git_state") %in% names(settings_row)) || !identical(as.character(settings_row$git_sha[[1]]), unname(git_identity[["git_sha"]])) ||
          !identical(as.character(settings_row$git_state[[1]]), unname(git_identity[["git_state"]]))) stop("Runner-settings Git provenance mismatch for ", margin, ".", call. = FALSE)
      if (!identical(unique(as.character(ledger$runner_sha256[rows])), unname(jss_recovery_sha256(runner_path[[margin]])))) stop("Runner source SHA mismatch for ", margin, ".", call. = FALSE)
    }
    actual_package_version <- as.character(read.dcf(file.path(repo_root, "DESCRIPTION"))[1, "Version"])
    if (!identical(unique(as.character(ledger$package_source_sha256[production])), jss_recovery_package_source_sha256(repo_root)) ||
        anyNA(jss_recovery_portable_source_identity(ledger$package_source_path[production])) ||
        !identical(unique(as.character(ledger$package_version[production])), actual_package_version)) stop("Checked-out package source identity does not match production attempts.", call. = FALSE)
  }
  fixed <- jss_recovery_adapt_fixed(read(paths[["bcpe_fixed"]]), read(paths[["nbi_fixed"]]), ledger)
  smooth <- jss_recovery_adapt_smooth(read(paths[["bcpe_smooth"]]), read(paths[["nbi_smooth"]]), ledger)
  study_ids <- c(jss_recovery_study_identity(ledger, "BCPE")$study_id, jss_recovery_study_identity(ledger, "NBI")$study_id)
  predictive <- jss_recovery_adapt_by_rep(
    read(paths[["bcpe_predictive"]]), read(paths[["nbi_predictive"]]),
    study_ids, ledger
  )
  diagnostic <- jss_recovery_adapt_by_rep(
    read(paths[["bcpe_diagnostic"]]), read(paths[["nbi_diagnostic"]]),
    study_ids, ledger
  )
  ledger <- jss_recovery_reconcile_retention(ledger, fixed, smooth, predictive, diagnostic)
  validation <- jss_recovery_validate(ledger, fixed, smooth, predictive, diagnostic)
  if (any(!validation$pass)) {
    stop("Main-recovery evidence validation failed: ", paste(validation$check_id[!validation$pass], collapse = ", "), call. = FALSE)
  }

  outputs <- list(
    attempt_metadata = ledger,
    design_table = jss_recovery_design_table(ledger),
    attempt_status_summary = jss_recovery_counts(ledger),
    failure_reason_summary = jss_recovery_failure_summary(ledger, conf.level),
    fixed_parameter_recovery = jss_recovery_fixed_summary(fixed, ledger, conf.level),
    smooth_recovery = jss_recovery_smooth_summary(smooth, ledger, conf.level),
    runtime_summary = jss_recovery_runtime_summary(ledger, conf.level),
    predictive_metrics = jss_recovery_predictive_summary(predictive, ledger, conf.level),
    diagnostic_metrics = jss_recovery_long_metric_summary(
      diagnostic, ledger,
      c("rosenblatt_ks", "rosenblatt_cvm", "abs_rosenblatt_lag1_cor", "abs_rosenblatt_normal_lag1_cor", "rosenblatt_mean_abs_time_mean", "rosenblatt_normal_mean_abs_time_mean"),
      conf.level
    ),
    weak_t_copula_shape_recovery = jss_recovery_t_shape_summary(fixed, ledger, conf.level),
    failure_inclusive_sensitivity = jss_recovery_failure_sensitivity(ledger, fixed, smooth, predictive, conf.level),
    paired_method_differences = jss_recovery_paired_method_differences(ledger, fixed, smooth, predictive, diagnostic, conf.level),
    fixed_by_attempt = fixed,
    smooth_by_attempt = smooth,
    predictive_by_attempt = predictive,
    diagnostic_by_attempt = diagnostic,
    runner_settings_identity = settings_identity,
    runner_settings_bcpe = settings_by_margin$BCPE,
    runner_settings_nbi = settings_by_margin$NBI,
    control_source_manifest = jss_recovery_control_source_identity(repo_root),
    bundle_status = data.frame(
      schema_version = jss_recovery_schema_version(),
      status = if (identical(jss_recovery_bundle_status(ledger), "authoritative_post_phase1")) "candidate_post_phase1_unapproved" else "legacy_reconciliation_not_authoritative",
      publication_eligible = FALSE,
      reason = if (identical(jss_recovery_bundle_status(ledger), "authoritative_post_phase1")) "candidate bundle requires checkout-external detached Ed25519 attestation from the pinned production key before publication use" else "contains legacy pre-Phase1 or route-mismatched evidence; reconciliation only",
      stringsAsFactors = FALSE
    ),
    production_provenance = jss_recovery_production_provenance(repo_root, ledger),
    evidence_validation = validation
  )
  outputs$input_provenance <- jss_recovery_frame_identity(outputs[c(
    "attempt_metadata", "fixed_by_attempt", "smooth_by_attempt", "predictive_by_attempt",
    "diagnostic_by_attempt", "runner_settings_identity", "runner_settings_bcpe", "runner_settings_nbi"
  )])
  paths_written <- jss_recovery_install_bundle(outputs, output_dir, provenance, repo_root)
  structure(unname(paths_written), resumed = FALSE)
}

jss_main_recovery_artifact_names <- function() {
  c(
    "attempt_metadata.csv", "design_table.csv", "attempt_status_summary.csv",
    "failure_reason_summary.csv", "fixed_parameter_recovery.csv",
    "smooth_recovery.csv", "runtime_summary.csv", "predictive_metrics.csv",
    "diagnostic_metrics.csv", "weak_t_copula_shape_recovery.csv",
    "failure_inclusive_sensitivity.csv", "paired_method_differences.csv",
    "fixed_by_attempt.csv", "smooth_by_attempt.csv", "predictive_by_attempt.csv", "diagnostic_by_attempt.csv", "runner_settings_identity.csv",
    "runner_settings_bcpe.csv", "runner_settings_nbi.csv", "control_source_manifest.csv",
    "bundle_status.csv",
    "production_provenance.csv", "evidence_validation.csv",
    "input_provenance.csv", "output_manifest.csv", "bundle_checkpoint.csv"
  )
}

jss_recovery_frames_equal <- function(actual, expected, tolerance = 1e-10) {
  if (!identical(names(actual), names(expected)) || nrow(actual) != nrow(expected)) return(FALSE)
  if (!nrow(actual)) return(TRUE)
  order_frame <- function(x) {
    key <- do.call(paste, c(lapply(x, function(column) ifelse(is.na(column), "<NA>", format(column, digits = 17, scientific = FALSE, trim = TRUE))), sep = "\r"))
    x[order(key), , drop = FALSE]
  }
  actual <- order_frame(actual); expected <- order_frame(expected)
  all(vapply(names(actual), function(name) {
    if (is.numeric(actual[[name]]) && is.numeric(expected[[name]])) {
      same_na <- identical(is.na(actual[[name]]), is.na(expected[[name]]))
      same_inf <- identical(is.infinite(actual[[name]]), is.infinite(expected[[name]])) && identical(sign(actual[[name]][is.infinite(actual[[name]])]), sign(expected[[name]][is.infinite(expected[[name]])]))
      finite <- is.finite(actual[[name]]) & is.finite(expected[[name]])
      same_na && same_inf && all(abs(actual[[name]][finite] - expected[[name]][finite]) <= tolerance * pmax(1, abs(expected[[name]][finite])))
    } else identical(as.character(actual[[name]]), as.character(expected[[name]]))
  }, logical(1)))
}

jss_recovery_pinned_public_key <- function() {
  key <- as.raw(c(66, 91, 71, 233, 21, 247, 172, 45, 215, 202, 170, 0, 64,
    43, 83, 206, 23, 50, 48, 154, 25, 217, 178, 37, 252, 59, 158, 195, 237, 0, 31, 216))
  fingerprint <- digest::digest(key, "sha256", serialize = FALSE)
  if (!identical(fingerprint, "cb73c05cede55bfd56357b1780b90c1bf254413d82765c7a682fc9db3a0d8587")) stop("Pinned main-recovery Ed25519 public-key fingerprint is inconsistent.", call. = FALSE)
  key
}

jss_recovery_parse_rfc3339_utc <- function(x, label = "timestamp") {
  if (!is.character(x) || anyNA(x) || any(!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$", x))) stop(label, " is not an RFC3339 UTC instant.", call. = FALSE)
  parsed <- as.POSIXct(strptime(x, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  if (anyNA(parsed) || !identical(format(parsed, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), x)) stop(label, " is not a real RFC3339 UTC instant.", call. = FALSE)
  parsed
}

jss_recovery_validate_approval_order <- function(approved_at_utc, execution_timestamps) {
  approved <- jss_recovery_parse_rfc3339_utc(approved_at_utc, "Main-recovery approval timestamp")
  executed <- jss_recovery_parse_rfc3339_utc(execution_timestamps, "Main-recovery execution timestamp")
  if (approved <= max(executed)) stop("Main-recovery approval must occur after every execution.", call. = FALSE)
  invisible(TRUE)
}

jss_recovery_bundle_identity <- function(path) {
  artifact_names <- jss_main_recovery_artifact_names()
  files <- file.path(path, artifact_names)
  if (any(!file.exists(files))) stop("Cannot identify an incomplete main-recovery bundle.", call. = FALSE)
  hashes <- unname(jss_recovery_sha256(files))
  bundle_sha <- digest::digest(paste(artifact_names, hashes, sep = "\t", collapse = "\n"), "sha256", serialize = FALSE)
  ledger <- utils::read.csv(file.path(path, "attempt_metadata.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  design <- utils::read.csv(file.path(path, "design_table.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  if (!"execution_completed_at_utc" %in% names(ledger)) stop("Main-recovery attempts lack execution completion instants.", call. = FALSE)
  execution <- ledger[order(ledger$attempt_id), c("attempt_id", "execution_completed_at_utc"), drop = FALSE]
  jss_recovery_parse_rfc3339_utc(execution$execution_completed_at_utc, "Main-recovery execution timestamp")
  settings <- data.frame(
    margin_family = c("BCPE", "NBI"),
    file = c("runner_settings_bcpe.csv", "runner_settings_nbi.csv"),
    sha256 = unname(jss_recovery_sha256(file.path(path, c("runner_settings_bcpe.csv", "runner_settings_nbi.csv")))),
    stringsAsFactors = FALSE
  )
  package_sha <- unique(as.character(ledger$package_source_sha256)); producer <- utils::read.csv(file.path(path, "production_provenance.csv"), stringsAsFactors = FALSE)
  producer_sha <- producer$value[producer$key == "producer_sha256"]
  planned_groups <- unique(ledger[c("study_id", "scenario_id", "method", "target_replicates")])
  registered <- list(n_subjects = 500L, n_time = 4L, target_replicates = 100L,
    planned_method_replicate_cells = 400L, scenario_count = 2L, method_count = 4L)
  if (length(package_sha) != 1L || length(producer_sha) != 1L || nrow(design) != 2L ||
      !all(design$n_subjects == registered$n_subjects & design$n_time == registered$n_time & design$target_replicates == registered$target_replicates) ||
      nrow(planned_groups) != registered$method_count ||
      !all(planned_groups$target_replicates == registered$target_replicates) ||
      sum(planned_groups$target_replicates) != registered$planned_method_replicate_cells) stop("Main-recovery bundle identity does not satisfy the exact registered n/T/R design.", call. = FALSE)
  list(
    schema_version = 1L, study = "main-recovery", bundle_sha256 = bundle_sha,
    output_manifest_sha256 = unname(jss_recovery_sha256(file.path(path, "output_manifest.csv"))),
    bundle_checkpoint_sha256 = unname(jss_recovery_sha256(file.path(path, "bundle_checkpoint.csv"))),
    runner_settings_manifest = settings,
    control_source_manifest_sha256 = unname(jss_recovery_sha256(file.path(path, "control_source_manifest.csv"))),
    package_source_sha256 = package_sha, producer_sha256 = producer_sha,
    registered_configuration = registered, execution_manifest = execution
  )
}

jss_recovery_bundle_trust_hash <- function(path) {
  jss_recovery_bundle_identity(path)$bundle_sha256
}

jss_recovery_require_external_attestation <- function(identity, root, attestation_path, signature_path) {
  if (!requireNamespace("sodium", quietly = TRUE)) stop("sodium is required for main-recovery promotion verification.", call. = FALSE)
  if (!is.character(attestation_path) || length(attestation_path) != 1L || !nzchar(attestation_path) ||
      !is.character(signature_path) || length(signature_path) != 1L || !nzchar(signature_path)) stop("Main-recovery bundle lacks a detached production attestation/signature.", call. = FALSE)
  root_prefix <- paste0(tolower(normalizePath(root, winslash = "/", mustWork = TRUE)), "/")
  paths <- vapply(c(attestation_path, signature_path), normalizePath, character(1), winslash = "/", mustWork = TRUE)
  if (any(startsWith(tolower(paste0(paths, "/")), root_prefix))) stop("Main-recovery attestation/signature must be external to the editable checkout.", call. = FALSE)
  message_raw <- readBin(paths[[1L]], "raw", n = file.info(paths[[1L]])$size)
  signature_raw <- readBin(paths[[2L]], "raw", n = file.info(paths[[2L]])$size)
  if (!isTRUE(tryCatch(sodium::sig_verify(message_raw, signature_raw, jss_recovery_pinned_public_key()), error = function(e) FALSE))) stop("Main-recovery attestation lacks a valid detached production signature from the pinned key.", call. = FALSE)
  attestation <- tryCatch(unserialize(message_raw), error = function(e) NULL)
  expected_names <- c(names(identity), "approved_at_utc", "approver")
  if (!is.list(attestation) || !identical(names(attestation), expected_names) ||
      !identical(attestation[names(identity)], identity) || !is.character(attestation$approver) || length(attestation$approver) != 1L || !nzchar(trimws(attestation$approver))) stop("Main-recovery detached attestation does not bind the exact canonical bundle identity.", call. = FALSE)
  jss_recovery_validate_approval_order(attestation$approved_at_utc, identity$execution_manifest$execution_completed_at_utc)
  attestation
}

jss_main_recovery_validate_public_bundle <- function(path, conf.level = 0.95,
    repo_root = normalizePath(".", winslash = "/", mustWork = TRUE),
    attestation_path = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MAIN_RECOVERY_ATTESTATION", unset = ""),
    signature_path = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MAIN_RECOVERY_ATTESTATION_SIGNATURE", unset = ""),
    require_attestation = TRUE) {
  required <- file.path(path, jss_main_recovery_artifact_names())
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    stop("Tracked main-recovery evidence is incomplete: ", paste(basename(missing), collapse = ", "), call. = FALSE)
  }
  read <- function(name) utils::read.csv(file.path(path, name), stringsAsFactors = FALSE, check.names = FALSE)
  manifest <- read("output_manifest.csv")
  checkpoint <- read("bundle_checkpoint.csv")
  input_provenance <- read("input_provenance.csv")
  expected_inputs <- jss_recovery_public_input_map()
  if (!all(c("input_id", "path", "sha256", "bytes") %in% names(input_provenance)) || anyDuplicated(input_provenance$input_id) ||
      any(!grepl("^[0-9a-f]{64}$", input_provenance$sha256)) || any(dirname(input_provenance$path) != ".")) {
    stop("Main-recovery bundled input provenance is malformed.", call. = FALSE)
  }
  if (!identical(as.character(input_provenance$input_id), names(expected_inputs)) ||
      !identical(as.character(input_provenance$path), unname(expected_inputs))) {
    stop("Main-recovery bundled input provenance does not exactly match the canonical normalized-input allowlist.", call. = FALSE)
  }
  input_paths <- file.path(path, input_provenance$path)
  input_ok <- file.exists(input_paths) & input_provenance$sha256 == unname(jss_recovery_sha256(input_paths)) & input_provenance$bytes == unname(file.info(input_paths)$size)
  if (any(!input_ok)) stop("Main-recovery bundled input provenance hash mismatch.", call. = FALSE)
  protected <- setdiff(jss_main_recovery_artifact_names(), c("output_manifest.csv", "bundle_checkpoint.csv"))
  if (!identical(sort(as.character(manifest$artifact)), sort(protected)) || anyDuplicated(manifest$artifact)) {
    stop("Main-recovery output manifest does not exactly cover the publication allowlist.", call. = FALSE)
  }
  if (!all(manifest$schema_version == jss_recovery_schema_version()) || any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) stop("Main-recovery output manifest has invalid schema or SHA-256 fields.", call. = FALSE)
  manifest_paths <- file.path(path, manifest$artifact)
  hash_ok <- manifest$sha256 == unname(jss_recovery_sha256(manifest_paths)) & manifest$bytes == unname(file.info(manifest_paths)$size)
  if (any(!hash_ok)) stop("Main-recovery publication artifact hash mismatch: ", paste(manifest$artifact[!hash_ok], collapse = ", "), call. = FALSE)
  input_digest <- digest::digest(paste(input_provenance$input_id, input_provenance$sha256, sep = "=", collapse = "\n"), algo = "sha256", serialize = FALSE)
  if (nrow(checkpoint) != 1L || checkpoint$schema_version != jss_recovery_schema_version() ||
      checkpoint$producer_sha256 != jss_recovery_producer_sha256(repo_root) || checkpoint$input_provenance_sha256 != input_digest ||
      checkpoint$output_manifest_sha256 != unname(jss_recovery_sha256(file.path(path, "output_manifest.csv"))) ||
      checkpoint$artifact_count != nrow(manifest)) {
    stop("Main-recovery portable checkpoint is invalid.", call. = FALSE)
  }
  recorded_control_sources <- read("control_source_manifest.csv")
  current_control_sources <- jss_recovery_control_source_identity(repo_root)
  if (!jss_recovery_frames_equal(recorded_control_sources, current_control_sources)) {
    stop("Main-recovery control-source manifest does not match the current validator, runners, target/gate, allowlist, and detached-signing controls.", call. = FALSE)
  }

  ledger_input <- read("attempt_metadata.csv")
  recomputed_convergence <- jss_recovery_recompute_convergence(ledger_input)
  if (!identical(as.logical(ledger_input$converged), as.logical(recomputed_convergence$converged)) ||
      !identical(as.logical(ledger_input$convergence_eligible), as.logical(recomputed_convergence$convergence_eligible)) ||
      !identical(as.character(ledger_input$convergence_status), as.character(recomputed_convergence$convergence_status))) {
    stop("Main-recovery bundled convergence booleans/status contradict strict raw convergence evidence.", call. = FALSE)
  }
  ledger_input$converged <- recomputed_convergence$converged
  ledger_input$convergence_eligible <- recomputed_convergence$convergence_eligible
  ledger_input$convergence_status <- recomputed_convergence$convergence_status
  fixed <- read("fixed_by_attempt.csv"); smooth <- read("smooth_by_attempt.csv")
  predictive <- read("predictive_by_attempt.csv"); diagnostic <- read("diagnostic_by_attempt.csv")
  settings_identity <- read("runner_settings_identity.csv")
  settings_payload <- list(BCPE = read("runner_settings_bcpe.csv"), NBI = read("runner_settings_nbi.csv"))
  git_identity <- jss_recovery_git_identity(repo_root)
  jss_recovery_require_clean_checkout(git_identity)
  identity_fields <- c("runner_settings_signature", "runner_settings_sha256", "runner_sha256", "package_source_path", "package_version", "package_source_sha256")
  if (nrow(settings_identity) != 2L || !all(c("margin_family", "study_id", identity_fields) %in% names(settings_identity))) stop("Main-recovery runner-settings identity artifact is malformed.", call. = FALSE)
  for (margin in c("BCPE", "NBI")) {
    setting <- settings_identity[settings_identity$margin_family == margin, , drop = FALSE]
    attempts <- ledger_input[ledger_input$margin_family == margin, , drop = FALSE]
    if (nrow(setting) != 1L || any(vapply(identity_fields, function(name) !identical(unique(as.character(attempts[[name]])), as.character(setting[[name]][[1]])), logical(1)))) stop("Main-recovery attempt/settings identity mismatch for ", margin, ".", call. = FALSE)
    payload <- settings_payload[[margin]]
    signature <- tryCatch(jss_recovery_settings_signature(payload, margin), error = function(e) NA_character_)
    if (nrow(payload) != 1L || is.na(signature) || !identical(signature, as.character(setting$runner_settings_signature[[1]])) ||
        !identical(digest::digest(signature, algo = "sha256", serialize = FALSE), as.character(setting$runner_settings_sha256[[1]]))) {
      stop("Main-recovery canonical settings payload mismatch for ", margin, ".", call. = FALSE)
    }
    if (!all(c("git_sha", "git_state") %in% names(payload)) || !identical(as.character(payload$git_sha[[1]]), unname(git_identity[["git_sha"]])) ||
        !identical(as.character(payload$git_state[[1]]), unname(git_identity[["git_state"]]))) stop("Main-recovery Git provenance mismatch for ", margin, ".", call. = FALSE)
  }
  sha_fields <- c("runner_settings_sha256", "runner_sha256", "package_source_sha256")
  if (any(vapply(sha_fields, function(name) any(!grepl("^[0-9a-f]{64}$", as.character(ledger_input[[name]]))), logical(1)))) stop("Main-recovery source identity contains a noncanonical SHA-256 value.", call. = FALSE)
  recomputed_settings_sha <- vapply(as.character(ledger_input$runner_settings_signature), digest::digest, character(1), algo = "sha256", serialize = FALSE)
  if (any(as.character(ledger_input$runner_settings_sha256) != recomputed_settings_sha)) stop("Main-recovery canonical settings signature/SHA mismatch.", call. = FALSE)
  runner_paths <- c(BCPE = file.path(repo_root, "paper", "scripts", "final-simulations", "bcpe-t", "simulation_bcpe_t_gamlss_comparison.R"),
    NBI = file.path(repo_root, "paper", "scripts", "final-simulations", "nbi-clayton", "compare_gamlss_ours_nbi_sigma_smooth.R"))
  for (margin in names(runner_paths)) {
    actual <- unname(jss_recovery_sha256(runner_paths[[margin]]))
    if (!identical(unique(as.character(ledger_input$runner_sha256[ledger_input$margin_family == margin])), actual)) stop("Main-recovery runner source SHA mismatch for ", margin, ".", call. = FALSE)
  }
  actual_package_version <- as.character(read.dcf(file.path(repo_root, "DESCRIPTION"))[1, "Version"])
  if (!identical(unique(as.character(ledger_input$package_source_sha256)), jss_recovery_package_source_sha256(repo_root)) ||
      anyNA(jss_recovery_portable_source_identity(ledger_input$package_source_path)) ||
      !identical(unique(as.character(ledger_input$package_version)), actual_package_version)) stop("Main-recovery checkout package identity mismatch.", call. = FALSE)
  ledger <- jss_recovery_reconcile_retention(ledger_input, fixed, smooth, predictive, diagnostic)
  validation <- jss_recovery_validate(ledger, fixed, smooth, predictive, diagnostic)
  if (any(!validation$pass)) stop("Main-recovery recomputed validation failed: ", paste(validation$check_id[!validation$pass], collapse = ", "), call. = FALSE)
  if (!identical(jss_recovery_bundle_status(ledger), "authoritative_post_phase1")) stop("Main-recovery authority could not be recomputed from attempts.", call. = FALSE)
  cell_count <- nrow(unique(ledger[c("study_id", "scenario_id", "method", "replicate")]))
  if (cell_count != 400L || nrow(ledger) < cell_count) stop("Main-recovery planned-cell/attempt-row cardinality is invalid.", call. = FALSE)

  expected <- list(
    attempt_metadata = ledger,
    design_table = jss_recovery_design_table(ledger), attempt_status_summary = jss_recovery_counts(ledger),
    failure_reason_summary = jss_recovery_failure_summary(ledger, conf.level),
    fixed_parameter_recovery = jss_recovery_fixed_summary(fixed, ledger, conf.level),
    smooth_recovery = jss_recovery_smooth_summary(smooth, ledger, conf.level),
    runtime_summary = jss_recovery_runtime_summary(ledger, conf.level),
    predictive_metrics = jss_recovery_predictive_summary(predictive, ledger, conf.level),
    diagnostic_metrics = jss_recovery_long_metric_summary(diagnostic, ledger,
      c("rosenblatt_ks", "rosenblatt_cvm", "abs_rosenblatt_lag1_cor", "abs_rosenblatt_normal_lag1_cor", "rosenblatt_mean_abs_time_mean", "rosenblatt_normal_mean_abs_time_mean"), conf.level),
    weak_t_copula_shape_recovery = jss_recovery_t_shape_summary(fixed, ledger, conf.level),
    failure_inclusive_sensitivity = jss_recovery_failure_sensitivity(ledger, fixed, smooth, predictive, conf.level),
    paired_method_differences = jss_recovery_paired_method_differences(ledger, fixed, smooth, predictive, diagnostic, conf.level),
    evidence_validation = validation
  )
  for (name in names(expected)) {
    actual <- read(paste0(name, ".csv"))
    if (!jss_recovery_frames_equal(actual, expected[[name]])) stop("Main-recovery recomputation mismatch: ", name, call. = FALSE)
  }
  status <- read("bundle_status.csv")
  expected_candidate_status <- if (identical(jss_recovery_bundle_status(ledger), "authoritative_post_phase1")) "candidate_post_phase1_unapproved" else "legacy_reconciliation_not_authoritative"
  if (nrow(status) != 1L || isTRUE(as.logical(status$publication_eligible)) || status$status != expected_candidate_status) stop("Main-recovery candidate status is inconsistent with recomputed authority and external-approval policy.", call. = FALSE)
  provenance <- read("production_provenance.csv")
  recorded_producer <- provenance$value[provenance$key == "producer_sha256"]
  if (length(recorded_producer) != 1L || !identical(recorded_producer, jss_recovery_producer_sha256(repo_root))) stop("Main-recovery public producer/control-source identity is inconsistent.", call. = FALSE)
  recorded_source <- provenance$value[provenance$key == "package_source_sha256"]
  if (length(recorded_source) != 1L || !identical(recorded_source, unique(ledger$package_source_sha256))) stop("Main-recovery package-source identity is inconsistent.", call. = FALSE)
  recorded_git_sha <- provenance$value[provenance$key == "git_sha"]
  recorded_git_state <- provenance$value[provenance$key == "git_state"]
  recorded_worktree_sha <- provenance$value[provenance$key == "tracked_worktree_sha256"]
  recorded_package_path <- provenance$value[provenance$key == "package_path"]
  if (!identical(recorded_git_sha, unname(git_identity[["git_sha"]])) || !identical(recorded_git_state, unname(git_identity[["git_state"]])) ||
      !identical(recorded_worktree_sha, unname(git_identity[["tracked_worktree_sha256"]])) ||
      length(recorded_package_path) != 1L || is.na(jss_recovery_portable_source_identity(recorded_package_path))) {
    stop("Main-recovery public Git/worktree or portable checkout provenance does not match the current checkout.", call. = FALSE)
  }
  identity <- jss_recovery_bundle_identity(path)
  if (isTRUE(require_attestation)) jss_recovery_require_external_attestation(identity, repo_root, attestation_path, signature_path)
  attr(required, "candidate_trust_sha256") <- identity$bundle_sha256
  attr(required, "bundle_identity") <- identity
  invisible(required)
}

jss_recovery_write_manuscript_assets <- function(bundle_dir, settings) {
  status <- utils::read.csv(file.path(bundle_dir, "attempt_status_summary.csv"), stringsAsFactors = FALSE)
  fixed <- utils::read.csv(file.path(bundle_dir, "fixed_parameter_recovery.csv"), stringsAsFactors = FALSE)
  smooth <- utils::read.csv(file.path(bundle_dir, "smooth_recovery.csv"), stringsAsFactors = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for main-recovery manuscript figures.", call. = FALSE)
  write_table <- function(margin, path, label) {
    study <- unique(status$study_id[grepl(tolower(margin), tolower(status$study_id), fixed = TRUE)])
    rows <- status[status$study_id %in% study, ]
    body <- paste(sprintf("%s & %d & %d & %d & %.3f \\\\", rows$method, rows$planned_method_replicate_cells, rows$actual_attempt_rows, rows$retained, rows$retention_rate), collapse = "\n")
    lines <- c("\\begin{table}[t]", "\\centering", "\\begin{tabular}{lrrrr}", "Method & Planned cells & Attempt rows & Retained & Retention rate \\\\ \\hline", body, "\\end{tabular}",
      paste0("\\caption{Main-recovery lifecycle counts derived from the authoritative attempt ledger for ", margin, ". Nonconverged and incomplete attempts are excluded from ordinary summaries and retained in failure sensitivity.}"), paste0("\\label{", label, "}"), "\\end{table}")
    writeLines(lines, path, useBytes = TRUE)
  }
  write_recovery_table <- function(margin, path, label) {
    study <- unique(status$study_id[grepl(tolower(margin), tolower(status$study_id), fixed = TRUE)])
    rows <- fixed[fixed$study_id %in% study, ]
    rows <- stats::aggregate(cbind(abs_bias = abs(bias), rmse = rmse, coverage = coverage) ~ method, rows, mean, na.rm = TRUE)
    body <- paste(sprintf("%s & %.4f & %.4f & %.3f \\\\", rows$method, rows$abs_bias, rows$rmse, rows$coverage), collapse = "\n")
    lines <- c("\\begin{table}[t]", "\\centering", "\\begin{tabular}{lrrr}", "Method & Mean absolute bias & Mean RMSE & Mean coverage \\\\ \\hline", body, "\\end{tabular}",
      paste0("\\caption{Fixed-parameter recovery for ", margin, " from converged, complete, publication-retained attempts in the authoritative main-recovery bundle. Full estimand-specific MCSEs and intervals are in the public evidence table.}"), paste0("\\label{", label, "}"), "\\end{table}")
    writeLines(lines, path, useBytes = TRUE)
  }
  table_paths <- file.path(settings$out_dir, c("paper_simulation_bcpe_t_fit_characteristics.tex", "paper_simulation_bcpe_t_fixed_parameter_bias_rmse.tex", "paper_simulation_nbi_clayton_highsignal_fit_characteristics.tex", "paper_simulation_nbi_clayton_highsignal_fixed_parameter_bias_rmse.tex"))
  write_table("bcpe", table_paths[[1]], "tab:bcpe-fit")
  write_recovery_table("bcpe", table_paths[[2]], "tab:bcpe-recovery")
  write_table("nbi", table_paths[[3]], "tab:nbi-fit")
  write_recovery_table("nbi", table_paths[[4]], "tab:nbi-recovery")
  figure_paths <- file.path(settings$out_dir, c("paper_simulation_bcpe_t_fixed_effect_recovery.png", "paper_simulation_bcpe_t_smooth_effect_recovery.png", "paper_simulation_nbi_clayton_highsignal_fixed_effect_recovery.png", "paper_simulation_nbi_clayton_highsignal_smooth_effect_recovery.png"))
  for (margin in c("BCPE", "NBI")) {
    study <- unique(status$study_id[grepl(tolower(margin), tolower(status$study_id), fixed = TRUE)])
    fx <- fixed[fixed$study_id %in% study, ]
    sm <- smooth[smooth$study_id %in% study, ]
    p_fixed <- ggplot2::ggplot(fx, ggplot2::aes(x = term, y = bias, colour = method, shape = method)) + ggplot2::geom_hline(yintercept = 0, linetype = 2) + ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.4)) + ggplot2::facet_wrap(~parameter, scales = "free_y") + ggplot2::theme_bw() + ggplot2::labs(x = NULL, y = "Bias", subtitle = "Converged, complete, publication-retained attempts only")
    p_smooth <- ggplot2::ggplot(sm, ggplot2::aes(x = parameter, y = irmse, fill = method)) + ggplot2::geom_col(position = "dodge") + ggplot2::theme_bw() + ggplot2::labs(x = NULL, y = "Integrated RMSE", subtitle = "Authoritative main-recovery attempt set")
    index <- if (margin == "BCPE") 1:2 else 3:4
    ggplot2::ggsave(figure_paths[index[[1]]], p_fixed, width = 9, height = 5.5, dpi = 180)
    ggplot2::ggsave(figure_paths[index[[2]]], p_smooth, width = 8, height = 5, dpi = 180)
  }
  list(tables = table_paths, figures = figure_paths)
}

jss_recovery_copy_exact <- function(source, destination) {
  if (length(source) != length(destination) || any(!file.exists(source))) stop("Main-recovery immutable staging source is incomplete.", call. = FALSE)
  invisible(lapply(unique(dirname(destination)), dir.create, recursive = TRUE, showWarnings = FALSE))
  if (!all(file.copy(source, destination, overwrite = FALSE))) stop("Could not create immutable main-recovery staged copy.", call. = FALSE)
  destination
}

jss_recovery_stage_approved_bundle <- function(source_dir, destination_dir, repo_root, attestation_path, signature_path) {
  names <- jss_main_recovery_artifact_names(); original <- file.path(source_dir, names)
  before <- unname(jss_recovery_sha256(original))
  first <- jss_main_recovery_validate_public_bundle(source_dir, repo_root = repo_root,
    attestation_path = attestation_path, signature_path = signature_path, require_attestation = TRUE)
  staging <- tempfile("main-recovery-approved-stage-", tmpdir = dirname(destination_dir)); dir.create(staging)
  on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)
  staged <- jss_recovery_copy_exact(original, file.path(staging, names))
  second <- jss_main_recovery_validate_public_bundle(staging, repo_root = repo_root,
    attestation_path = attestation_path, signature_path = signature_path, require_attestation = TRUE)
  after <- unname(jss_recovery_sha256(original))
  if (!identical(before, after) || !identical(attr(first, "bundle_identity"), attr(second, "bundle_identity"))) stop("Main-recovery source changed during immutable staged validation.", call. = FALSE)
  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
  destination <- file.path(destination_dir, names)
  if (!all(file.copy(staged, destination, overwrite = TRUE))) stop("Could not install validated main-recovery staged copy.", call. = FALSE)
  third <- jss_main_recovery_validate_public_bundle(destination_dir, repo_root = repo_root,
    attestation_path = attestation_path, signature_path = signature_path, require_attestation = TRUE)
  if (!identical(attr(second, "bundle_identity"), attr(third, "bundle_identity")) || !identical(after, unname(jss_recovery_sha256(original)))) stop("Main-recovery source or destination changed after staged validation/re-attestation.", call. = FALSE)
  destination
}

jss_run_phase2_main_recovery <- function(settings) {
  module_id <- "phase2-main-recovery"
  output_dir <- file.path(settings$data_dir, "main-recovery")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (identical(settings$profile, "paper")) {
    tracked <- file.path(settings$public_data_dir, "main-recovery")
    attestation <- if (!is.null(settings$main_recovery_attestation)) settings$main_recovery_attestation else Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MAIN_RECOVERY_ATTESTATION", unset = "")
    signature <- if (!is.null(settings$main_recovery_attestation_signature)) settings$main_recovery_attestation_signature else Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MAIN_RECOVERY_ATTESTATION_SIGNATURE", unset = "")
    paths <- jss_recovery_stage_approved_bundle(tracked, output_dir, settings$root, attestation, signature)
  } else if (identical(settings$profile, "full")) {
    input_root <- file.path(settings$out_dir, "main-recovery")
    bcpe_raw <- file.path(input_root, "bcpe-t"); nbi_raw <- file.path(input_root, "nbi-clayton")
    jss_run_script(file.path(settings$root, "paper", "scripts", "final-simulations", "bcpe-t", "simulation_bcpe_t_gamlss_comparison.R"),
      c(OUT_DIR = bcpe_raw, N_FITS = "100", N_CORES = as.character(settings$workers), SAVE_FITS = "0", COMPUTE_PREDICTIVE_SCORES = "1", VARIOGRAM_P_VALUES = "0.5,2", MAX_ELAPSED_SEC = "180", GAMLSS_LONGITUDINAL_SOURCE_ROOT = settings$root), settings$root)
    jss_run_script(file.path(settings$root, "paper", "scripts", "final-simulations", "nbi-clayton", "compare_gamlss_ours_nbi_sigma_smooth.R"),
      c(NBI_COMPARE_OUT_DIR = nbi_raw, NBI_COMPARE_REPS = "100", NBI_COMPARE_RESUME = "TRUE", NBI_COMPARE_SAVE_FITS = "FALSE", NBI_COMPARE_ENGINES = "gamlss|ours_rs_joint", NBI_COMPARE_VARIOGRAM_P_VALUES = "0.5,2", GAMLSS_LONGITUDINAL_SOURCE_ROOT = settings$root), settings$root)
    paths <- jss_build_main_recovery_evidence(
      input_root = input_root,
      output_dir = output_dir,
      resume = TRUE,
      repo_root = settings$root
    )
    validated <- jss_main_recovery_validate_public_bundle(output_dir, repo_root = settings$root, require_attestation = FALSE)
    jss_recovery_atomic_csv(data.frame(schema_version = jss_recovery_schema_version(), study = "main-recovery",
      status = "candidate_pending_detached_attestation", bundle_sha256 = attr(validated, "candidate_trust_sha256"), stringsAsFactors = FALSE),
      file.path(settings$out_dir, "main-recovery-candidate-identity.csv"))
    return(list(
      module_id = module_id, status = "candidate_pending_external_approval", data = unname(paths),
      tables = character(), figures = character(),
      notes = "Fresh main-recovery candidate generated and validated structurally; manuscript assets require a checkout-external detached Ed25519 attestation from the pinned production key."
    ))
  } else {
    return(list(
      module_id = module_id, status = "not_run_in_smoke", data = character(),
      tables = character(), figures = character(),
      notes = "The main recovery Monte Carlo is excluded from the smoke profile."
    ))
  }
  assets <- jss_recovery_write_manuscript_assets(output_dir, settings)
  list(
    module_id = module_id,
    status = "current",
    data = unname(paths),
    tables = assets$tables,
    figures = assets$figures,
    notes = "Validated authoritative n=500, T=4, R=100 BCPE/t and NBI/Clayton attempt-level evidence; exactly 400 planned method/replicate cells and a retry-aware dynamic attempt ledger."
  )
}
