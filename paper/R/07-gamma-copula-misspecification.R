jss_misspec_copulas <- function() {
  # Keep the evidence design inside the package's locked Phase 1 capability
  # contract.  The Gamma margin is currently validated end to end with the
  # Gaussian and Clayton copulas only; using the package-wide copula list here
  # would turn unsupported preflight errors into purported selection evidence.
  c("N", "C")
}

jss_misspec_capability_contract <- function(repo_root) {
  registry_path <- file.path(repo_root, "R", "capability-registry.R")
  if (!file.exists(registry_path)) {
    stop("Module 07 cannot locate the checked-out Phase 1 capability registry.", call. = FALSE)
  }
  registry <- new.env(parent = baseenv())
  sys.source(registry_path, envir = registry)
  version <- registry$.gl_capability_registry_version()
  gamma_copulas <- as.character(
    registry$.gl_capability_margin_specs()$GA$compatible_copulas
  )
  if (!identical(version, "2026.2") ||
      !identical(gamma_copulas, jss_misspec_copulas())) {
    stop(
      "Module 07's Gamma copula-selection design no longer matches the locked Phase 1 capability registry; review and version the evidence design before running.",
      call. = FALSE
    )
  }
  list(registry_version = version, copulas = gamma_copulas)
}

jss_misspec_tau_levels <- function() {
  data.frame(
    tau_label = c("moderate", "high"),
    target_tau = c(0.25, 0.55),
    stringsAsFactors = FALSE
  )
}

jss_misspec_stage <- function(profile) {
  requested <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MISSPEC_STAGE", unset = "")
  if (nzchar(requested)) {
    return(match.arg(requested, c("smoke", "pilot", "full")))
  }
  if (identical(profile, "expanded")) "pilot" else "smoke"
}

jss_misspec_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "DESCRIPTION")) &&
        file.exists(file.path(current, "paper", "R", "07-gamma-copula-misspecification.R"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) stop("Could not locate the Module 07 repository root.", call. = FALSE)
    current <- parent
  }
}

jss_misspec_config <- function(settings, stage = jss_misspec_stage(settings$profile)) {
  stage_explicit <- !missing(stage)
  stage <- match.arg(stage, c("smoke", "pilot", "full"))
  reps_override <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MISSPEC_REPS", unset = "")
  reps <- switch(stage, smoke = 1L, pilot = 10L, full = 100L)
  if (nzchar(reps_override)) {
    requested_reps <- suppressWarnings(as.integer(reps_override))
    if (is.na(requested_reps) || requested_reps < 1L) {
      stop("GAMLSS_LONGITUDINAL_JSS_MISSPEC_REPS must be a positive integer.", call. = FALSE)
    }
    if (identical(stage, "full") && requested_reps != 100L) {
      stop("Module 07 full configuration is immutable and requires exactly 100 replicates.", call. = FALSE)
    }
    if (!identical(stage, "full") && !isTRUE(stage_explicit)) {
      stop("Module 07 replicate overrides require an explicit smoke or pilot nonpublication stage.", call. = FALSE)
    }
    reps <- if (identical(stage, "full")) 100L else requested_reps
  }

  timeout_raw <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MISSPEC_TIMEOUT_SEC", unset = "")
  timeout_default <- if (identical(stage, "smoke")) 20 else 600
  timeout <- if (nzchar(timeout_raw)) suppressWarnings(as.numeric(timeout_raw)) else timeout_default
  repo_root <- if (!is.null(settings$root) && nzchar(settings$root)) settings$root else jss_misspec_repo_root()
  capability <- jss_misspec_capability_contract(repo_root)
  source_path <- normalizePath(
    file.path(repo_root, "paper", "R", "07-gamma-copula-misspecification.R"),
    winslash = "/", mustWork = TRUE
  )
  config <- list(
    stage = stage,
    margin_family = "GA",
    copulas = capability$copulas,
    capability_registry_version = capability$registry_version,
    tau_levels = if (identical(stage, "smoke")) {
      jss_misspec_tau_levels()[1L, , drop = FALSE]
    } else {
      jss_misspec_tau_levels()
    },
    sample_sizes = if (identical(stage, "smoke")) 20L else c(50L, 150L, 500L),
    times = 1:4,
    reps = reps,
    margin_mu = 2,
    margin_sigma = 0.5,
    max_outer_iter = if (identical(stage, "smoke")) 3L else 100L,
    max_inner_iter = if (identical(stage, "smoke")) 3L else 100L,
    max_elapsed_sec = timeout,
    compute_vcov = !identical(stage, "smoke"),
    seed = settings$seed + 700000L,
    repo_root = normalizePath(repo_root, winslash = "/", mustWork = TRUE),
    source_path = source_path,
    producer_sha256 = jss_misspec_sha256_file(source_path),
    package_source_sha256 = jss_misspec_package_source_sha256(repo_root)
  )
  config$configuration_sha256 <- jss_misspec_sha256_object(
    jss_misspec_registered_config_identity(config)
  )
  jss_misspec_validate_config(config, require_callr = !identical(stage, "smoke"))
  config
}

jss_misspec_validate_config <- function(config, require_callr = TRUE, verify_source = TRUE) {
  capability <- jss_misspec_capability_contract(config$repo_root)
  if (!identical(config$margin_family, "GA") ||
      !identical(config$capability_registry_version, capability$registry_version) ||
      !identical(as.character(config$copulas), capability$copulas)) {
    stop(
      "Module 07 configuration diverges from the locked Gamma capability contract.",
      call. = FALSE
    )
  }
  if (identical(config$stage, "full") &&
      (length(config$reps) != 1L || !identical(as.integer(config$reps), 100L))) {
    stop("Module 07 full configuration is immutable and requires exactly 100 replicates.", call. = FALSE)
  }
  valid_timeout <- length(config$max_elapsed_sec) == 1L &&
    is.finite(config$max_elapsed_sec) && config$max_elapsed_sec > 0
  if (!valid_timeout) stop("Module 07 requires a finite positive registered kill timeout.", call. = FALSE)
  if (isTRUE(require_callr) && !requireNamespace("callr", quietly = TRUE)) {
    stop("Module 07 paper/full execution requires package 'callr' for killable fit boundaries.", call. = FALSE)
  }
  if (!file.exists(config$source_path) || !file.exists(file.path(config$repo_root, "DESCRIPTION"))) {
    stop("Module 07 checked-out source/config paths are invalid.", call. = FALSE)
  }
  hashes <- c(
    producer_sha256 = config$producer_sha256,
    package_source_sha256 = config$package_source_sha256,
    configuration_sha256 = config$configuration_sha256
  )
  if (!is.character(hashes) || length(hashes) != 3L ||
      any(!grepl("^[0-9a-f]{64}$", hashes))) {
    stop("Module 07 configuration lacks valid checked-out source identities.", call. = FALSE)
  }
  if (isTRUE(verify_source)) {
    current_configuration <- jss_misspec_sha256_object(
      jss_misspec_registered_config_identity(config)
    )
    if (!identical(config$producer_sha256, jss_misspec_sha256_file(config$source_path)) ||
        !identical(config$package_source_sha256, jss_misspec_package_source_sha256(config$repo_root)) ||
        !identical(config$configuration_sha256, current_configuration)) {
      stop("Module 07 configuration is stale relative to the checked-out producer/package source.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

jss_misspec_validate_public_full_bundle <- function(results, grid, config) {
  if (!identical(config$stage, "full") || !identical(as.integer(config$reps), 100L)) {
    stop("Module 07 public integration requires the immutable full R=100 configuration.", call. = FALSE)
  }
  if (!is.data.frame(results) || !nrow(results) || !"stage" %in% names(results) ||
      any(is.na(results$stage)) || any(as.character(results$stage) != "full")) {
    stop("Module 07 public input is mixed-stage or non-full and is ineligible; filtering is forbidden.", call. = FALSE)
  }
  expected_ids <- as.character(grid$fit_id)
  actual_ids <- as.character(results$fit_id)
  if (nrow(results) != nrow(grid) || anyDuplicated(actual_ids) ||
      !identical(sort(actual_ids), sort(expected_ids))) {
    stop("Module 07 public input is not the exact complete full R=100 grid.", call. = FALSE)
  }
  issues <- jss_misspec_result_contract_issues(results, grid = grid, config = config)
  if (length(issues)) {
    stop("Module 07 public full bundle failed authoritative binding: ", paste(issues, collapse = "; "), call. = FALSE)
  }
  invisible(TRUE)
}

jss_misspec_paths <- function(settings) {
  module_id <- "07-gamma-copula-misspecification"
  list(
    results = file.path(settings$data_dir, paste0(module_id, "-results.csv")),
    grid = file.path(settings$data_dir, paste0(module_id, "-grid.csv")),
    summary = file.path(settings$tables_dir, paste0(module_id, "-summary.csv")),
    selection = file.path(settings$tables_dir, paste0(module_id, "-selection.csv")),
    selection_attempts = file.path(settings$data_dir, paste0(module_id, "-selection-attempts.csv")),
    selection_confusion = file.path(settings$tables_dir, paste0(module_id, "-selection-confusion.csv")),
    selection_failures = file.path(settings$tables_dir, paste0(module_id, "-selection-failures.csv")),
    paired_effects = file.path(settings$tables_dir, paste0(module_id, "-paired-effects.csv")),
    warning_audit = file.path(settings$tables_dir, paste0(module_id, "-warning-audit.csv")),
    execution_manifest = file.path(settings$data_dir, paste0(module_id, "-execution-manifest.csv")),
    review = file.path(settings$tables_dir, paste0(module_id, "-review-gate.csv")),
    heatmap = file.path(settings$figures_dir, paste0(module_id, "-delta-heatmap.png")),
    margin_rmse_heatmap = file.path(settings$figures_dir, paste0(module_id, "-margin-rmse-heatmap.png")),
    tau_error_heatmap = file.path(settings$figures_dir, paste0(module_id, "-tau-error-heatmap.png")),
    paper_summary_heatmap = file.path(settings$figures_dir, paste0(module_id, "-paper-summary-heatmap.png")),
    convergence = file.path(settings$figures_dir, paste0(module_id, "-convergence.png")),
    checkpoints = file.path(settings$data_dir, paste0(module_id, "-checkpoints-ga-nc-v2")),
    smoke_checkpoints = file.path(settings$data_dir, paste0(module_id, "-smoke-checkpoints"))
  )
}

jss_misspec_grid <- function(config) {
  base <- expand.grid(
    generating_copula = config$copulas,
    fitted_copula = config$copulas,
    tau_label = config$tau_levels$tau_label,
    n_subject = config$sample_sizes,
    rep = seq_len(config$reps),
    stringsAsFactors = FALSE
  )
  tau_map <- stats::setNames(config$tau_levels$target_tau, config$tau_levels$tau_label)
  base$target_tau <- unname(tau_map[base$tau_label])
  base$margin_family <- config$margin_family
  base$capability_registry_version <- config$capability_registry_version
  base$producer_sha256 <- config$producer_sha256
  base$package_source_sha256 <- config$package_source_sha256
  base$configuration_sha256 <- config$configuration_sha256
  base$n_time <- length(config$times)
  base$seed <- mapply(
    jss_misspec_seed,
    base$generating_copula,
    base$fitted_copula,
    base$tau_label,
    base$n_subject,
    base$rep,
    MoreArgs = list(base_seed = config$seed),
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
  base$dataset_seed <- mapply(
    jss_misspec_dataset_seed,
    base$generating_copula,
    base$tau_label,
    base$n_subject,
    base$rep,
    MoreArgs = list(base_seed = config$seed),
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
  base$fit_id <- mapply(
    jss_misspec_fit_id,
    base$generating_copula,
    base$fitted_copula,
    base$tau_label,
    base$n_subject,
    base$rep,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
  base
}

jss_misspec_seed <- function(generating_copula, fitted_copula, tau_label, n_subject, rep, base_seed) {
  copulas <- jss_misspec_copulas()
  gen_i <- match(generating_copula, copulas)
  fit_i <- match(fitted_copula, copulas)
  tau_i <- match(tau_label, jss_misspec_tau_levels()$tau_label)
  as.integer(base_seed + rep * 10000L + n_subject * 10L + gen_i * 100L + fit_i * 7L + tau_i)
}

jss_misspec_dataset_seed <- function(generating_copula, tau_label, n_subject, rep, base_seed) {
  jss_misspec_seed(generating_copula, generating_copula, tau_label, n_subject, rep, base_seed)
}

jss_misspec_fit_id <- function(generating_copula, fitted_copula, tau_label, n_subject, rep) {
  paste(
    paste0("gen-", generating_copula),
    paste0("fit-", fitted_copula),
    paste0("tau-", tau_label),
    paste0("n-", n_subject),
    sprintf("rep-%03d", as.integer(rep)),
    sep = "__"
  )
}

jss_misspec_checkpoint_dir <- function(paths, stage) {
  if (identical(stage, "smoke")) paths$smoke_checkpoints else paths$checkpoints
}

jss_misspec_checkpoint_path <- function(checkpoint_dir, fit_id) {
  file.path(checkpoint_dir, paste0(gsub("[^A-Za-z0-9_-]+", "_", fit_id), ".csv"))
}

jss_misspec_pending_grid <- function(grid, checkpoint_dir, config = NULL) {
  paths <- jss_misspec_checkpoint_path(checkpoint_dir, grid$fit_id)
  grid$checkpoint_path <- paths
  valid <- vapply(seq_len(nrow(grid)), function(i) {
    jss_misspec_checkpoint_valid(paths[[i]], grid[i, , drop = FALSE], config = config)
  }, logical(1L))
  rejected <- which(file.exists(paths) & !valid)
  if (length(rejected)) {
    quarantine_dir <- file.path(checkpoint_dir, "quarantine")
    dir.create(quarantine_dir, recursive = TRUE, showWarnings = FALSE)
    ledger <- lapply(rejected, function(i) {
      reason <- paste(jss_misspec_checkpoint_issues(paths[[i]], grid[i, , drop = FALSE], config), collapse = " | ")
      destination <- file.path(
        quarantine_dir,
        paste0(tools::file_path_sans_ext(basename(paths[[i]])), "-", format(Sys.time(), "%Y%m%d%H%M%OS3"), ".rejected.csv")
      )
      if (!file.rename(paths[[i]], destination)) stop("Could not quarantine invalid Module 07 checkpoint: ", paths[[i]], call. = FALSE)
      data.frame(
        fit_id = grid$fit_id[[i]], rejected_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        reason = reason, quarantined_path = normalizePath(destination, winslash = "/", mustWork = TRUE),
        stringsAsFactors = FALSE
      )
    })
    ledger <- do.call(rbind, ledger)
    ledger_path <- file.path(checkpoint_dir, "checkpoint-rejections-ledger.csv")
    prior <- if (file.exists(ledger_path)) utils::read.csv(ledger_path, stringsAsFactors = FALSE) else data.frame()
    jss_misspec_write_csv_atomic(jss_misspec_row_bind(list(prior, ledger)), ledger_path)
  }
  grid$checkpoint_exists <- valid
  grid[!valid, , drop = FALSE]
}

jss_misspec_row_bind <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(data.frame())
  }
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(all_names, names(row))
    for (name in missing) row[[name]] <- rep(NA, nrow(row))
    row[all_names]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

jss_misspec_is_png <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 8L) return(FALSE)
  identical(readBin(path, what = "raw", n = 8L), as.raw(c(137, 80, 78, 71, 13, 10, 26, 10)))
}

jss_misspec_sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (requireNamespace("openssl", quietly = TRUE)) {
    connection <- file(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    return(paste0(as.character(openssl::sha256(connection)), collapse = ""))
  }
  if (requireNamespace("digest", quietly = TRUE)) return(digest::digest(file = path, algo = "sha256", serialize = FALSE))
  stop("Package 'openssl' or 'digest' is required for SHA-256 output verification.", call. = FALSE)
}

jss_misspec_sha256_object <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("Package 'digest' is required for Module 07 identity checks.", call. = FALSE)
  digest::digest(x, algo = "sha256", serialize = TRUE)
}

jss_misspec_registered_config_identity <- function(config) {
  list(
    stage = config$stage, margin_family = config$margin_family,
    copulas = as.character(config$copulas),
    capability_registry_version = as.character(config$capability_registry_version),
    tau_levels = config$tau_levels, sample_sizes = as.integer(config$sample_sizes),
    times = as.integer(config$times), reps = as.integer(config$reps),
    margin_mu = config$margin_mu, margin_sigma = config$margin_sigma,
    max_outer_iter = as.integer(config$max_outer_iter),
    max_inner_iter = as.integer(config$max_inner_iter),
    max_elapsed_sec = as.numeric(config$max_elapsed_sec),
    compute_vcov = isTRUE(config$compute_vcov), seed = as.integer(config$seed),
    producer_sha256 = config$producer_sha256,
    package_source_sha256 = config$package_source_sha256
  )
}

jss_misspec_attestation_path <- function() {
  Sys.getenv("GAMLSS_LONGITUDINAL_JSS_COPULA_ATTESTATION", unset = "")
}

jss_misspec_signature_path <- function() {
  Sys.getenv("GAMLSS_LONGITUDINAL_JSS_COPULA_ATTESTATION_SIGNATURE", unset = "")
}

jss_misspec_pinned_public_key <- function() {
  as.raw(c(66, 91, 71, 233, 21, 247, 172, 45, 215, 202, 170, 0, 64,
    43, 83, 206, 23, 50, 48, 154, 25, 217, 178, 37, 252, 59, 158, 195,
    237, 0, 31, 216))
}

jss_misspec_parse_rfc3339_utc <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$", x)) {
    stop("Module 07 approval timestamp is not RFC3339 UTC.", call. = FALSE)
  }
  parsed <- as.POSIXct(strptime(x, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  if (is.na(parsed) || !identical(format(parsed, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), x)) {
    stop("Module 07 approval timestamp is not a real RFC3339 UTC instant.", call. = FALSE)
  }
  parsed
}

jss_misspec_package_source_sha256 <- function(root) {
  files <- sort(c(file.path(root, "DESCRIPTION"),
    list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE, recursive = TRUE)))
  files <- files[file.exists(files)]
  rel <- substring(normalizePath(files, winslash = "/", mustWork = TRUE),
    nchar(normalizePath(root, winslash = "/", mustWork = TRUE)) + 2L)
  jss_misspec_sha256_object(data.frame(
    file = rel, sha256 = vapply(files, jss_misspec_sha256_file, character(1L)),
    stringsAsFactors = FALSE
  ))
}

jss_misspec_load_checkout <- function(repo_root, expected_package_source_sha256) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required for checked-out Module 07 execution.", call. = FALSE)
  }
  pkgload::load_all(repo_root, quiet = TRUE, export_all = TRUE, helpers = FALSE)
  expected_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  namespace_path <- normalizePath(
    getNamespaceInfo(asNamespace("gamlss.longitudinal"), "path"),
    winslash = "/", mustWork = TRUE
  )
  actual_hash <- jss_misspec_package_source_sha256(repo_root)
  if (!identical(namespace_path, expected_root) ||
      !identical(actual_hash, expected_package_source_sha256)) {
    stop("Module 07 execution did not load and attest the immutable checked-out package source.", call. = FALSE)
  }
  c(namespace_path = namespace_path, package_source_sha256 = actual_hash)
}

jss_misspec_evidence_paths <- function(results_path) {
  results_path <- normalizePath(results_path, winslash = "/", mustWork = FALSE)
  base <- basename(results_path)
  suffix <- if (identical(base, "results.csv")) {
    function(name) file.path(dirname(results_path), paste0(name, ".csv"))
  } else if (endsWith(base, "-results.csv")) {
    stem <- substring(base, 1L, nchar(base) - nchar("-results.csv"))
    function(name) file.path(dirname(results_path), paste0(stem, "-", name, ".csv"))
  } else {
    stop("Module 07 results path is not a canonical candidate or generated filename.", call. = FALSE)
  }
  c(
    results = results_path,
    paired_effects = suffix("paired-effects"), warning_audit = suffix("warning-audit"),
    execution_manifest = suffix("execution-manifest"),
    selection_attempts = suffix("selection-attempts"), selection = suffix("selection")
  )
}

jss_misspec_assert_csv_reconstruction <- function(reconstructed, path, label) {
  temporary <- tempfile("module07-reconstruction-", fileext = ".csv")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  utils::write.csv(reconstructed, temporary, row.names = FALSE)
  if (!identical(jss_misspec_sha256_file(temporary), jss_misspec_sha256_file(path))) {
    stop("Module 07 ", label, " is not an exact reconstruction from approved attempts.", call. = FALSE)
  }
  invisible(TRUE)
}

jss_misspec_validate_claim_registry <- function(root) {
  path <- file.path(root, "paper", "phase2-claims.csv")
  claims <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  claims <- claims[claims$study_id == "copula-misspecification", , drop = FALSE]
  paired <- claims[claims$claim_type == "paired_effect", , drop = FALSE]
  selection <- claims[claims$claim_type == "selection_rate", , drop = FALSE]
  required <- c("claim_id", "tex_anchor", "metric", "lhs_method", "rhs_method",
    "attempt_artifact_id", "effect_artifact_id", "source_row_key", "wording_strength")
  valid <- all(required %in% names(claims)) && !anyDuplicated(claims$claim_id) &&
    nrow(paired) >= 2L && nrow(selection) >= 2L &&
    all(paired$attempt_artifact_id == "p2_copula_results") &&
    all(paired$effect_artifact_id == "p2_copula_paired_effects") &&
    all(paired$lhs_method == "misspecified_copula") &&
    all(paired$rhs_method == "correctly_specified_copula") &&
    all(paired$metric %in% c("joint_loglik", "margin_param_rmse", "tau_abs_error")) &&
    all(grepl("generating_copula=", paired$source_row_key, fixed = TRUE)) &&
    all(grepl("fitted_copula=", paired$source_row_key, fixed = TRUE)) &&
    all(grepl("metric=", paired$source_row_key, fixed = TRUE)) &&
    all(selection$attempt_artifact_id == "p2_copula_selection_attempts") &&
    all(selection$effect_artifact_id == "p2_copula_selection") &&
    all(selection$metric == "aic_correct_selection_rate") &&
    all(grepl("generating_copula=", selection$source_row_key, fixed = TRUE)) &&
    all(claims$tex_anchor == "sec:sim_sensitivity_copula")
  if (!isTRUE(valid)) {
    stop("Fatal Module 07 claim registry lacks explicit attempt/effect mappings for off-diagonal effects and selection.", call. = FALSE)
  }
  invisible(claims)
}

jss_misspec_validate_evidence_bundle <- function(results_path, config) {
  jss_misspec_validate_claim_registry(config$repo_root)
  paths <- jss_misspec_evidence_paths(results_path)
  if (any(!file.exists(paths))) {
    stop("Module 07 candidate evidence bundle is incomplete: ", paste(names(paths)[!file.exists(paths)], collapse = ", "), call. = FALSE)
  }
  results <- jss_misspec_upgrade_result_contract(utils::read.csv(paths[["results"]], stringsAsFactors = FALSE, check.names = FALSE))
  grid <- jss_misspec_grid(config)
  jss_misspec_validate_public_full_bundle(results, grid, config)
  derived <- jss_misspec_add_deltas(results)
  paired <- jss_misspec_paired_effects(derived)
  warnings <- jss_misspec_warning_audit(derived)
  attempts <- jss_misspec_selection_attempts(derived)
  selection <- jss_misspec_selection(derived)
  jss_misspec_assert_csv_reconstruction(paired, paths[["paired_effects"]], "paired-effect artifact")
  jss_misspec_assert_csv_reconstruction(warnings, paths[["warning_audit"]], "warning audit")
  jss_misspec_assert_csv_reconstruction(attempts, paths[["selection_attempts"]], "selection-attempt artifact")
  jss_misspec_assert_csv_reconstruction(selection, paths[["selection"]], "selection artifact")
  jss_misspec_validate_warning_audit(warnings)
  required_effect <- c(
    "metric", "lhs_method", "rhs_method", "n_attempted", "n_retained",
    "estimate", "mcse", "ci_lower", "ci_upper", "confidence_level", "ci_method", "seed_pairing"
  )
  if (!all(required_effect %in% names(paired)) || !nrow(paired) ||
      any(paired$n_attempted != config$reps) || any(paired$n_retained < 2L) ||
      any(!is.finite(paired$estimate) | !is.finite(paired$mcse) |
        !is.finite(paired$ci_lower) | !is.finite(paired$ci_upper)) ||
      any(paired$ci_lower > paired$estimate | paired$ci_upper < paired$estimate) ||
      any(paired$seed_pairing != "registered_dataset_seed_exact")) {
    stop("Fatal Module 07 paired-effect gate lacks complete same-dataset uncertainty evidence.", call. = FALSE)
  }
  manifest <- utils::read.csv(paths[["execution_manifest"]], stringsAsFactors = FALSE, check.names = FALSE)
  required_manifest <- c(
    "execution_schema_version", "study", "stage", "execution_started_at_utc", "execution_completed_at_utc",
    "git_sha", "git_state", "r_version", "platform", "os", "package_version",
    "dependency_versions", "dependency_versions_sha256", "workers_requested", "base_seed",
    "seed_min", "seed_max", "dataset_seed_count", "fit_seed_count", "planned_fit_rows", "result_rows",
    "producer_sha256", "package_source_sha256", "configuration_sha256", "results_sha256",
    "paired_effects_sha256", "warning_audit_sha256", "selection_attempts_sha256", "selection_sha256",
    "expected_warning_events", "nonconvergence_warning_events", "unexpected_warning_events"
  )
  git <- jss_misspec_git_identity(config$repo_root)
  hash_fields <- c("dependency_versions_sha256", "producer_sha256", "package_source_sha256",
    "configuration_sha256", "results_sha256", "paired_effects_sha256", "warning_audit_sha256",
    "selection_attempts_sha256", "selection_sha256")
  timestamp_ok <- function(x) is.character(x) && length(x) == 1L &&
    grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$", x)
  if (!identical(names(manifest), required_manifest) || nrow(manifest) != 1L ||
      manifest$execution_schema_version[[1L]] != 1L || manifest$study[[1L]] != "copula-misspecification" ||
      manifest$stage[[1L]] != "full" || manifest$git_state[[1L]] != "clean" ||
      !jss_misspec_git_commit_exists(config$repo_root, manifest$git_sha[[1L]]) || git$state != "clean" ||
      !timestamp_ok(manifest$execution_started_at_utc[[1L]]) || !timestamp_ok(manifest$execution_completed_at_utc[[1L]]) ||
      as.POSIXct(manifest$execution_completed_at_utc[[1L]], format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC") <
        as.POSIXct(manifest$execution_started_at_utc[[1L]], format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC") ||
      manifest$r_version[[1L]] != R.version.string || manifest$platform[[1L]] != R.version$platform ||
      manifest$package_version[[1L]] != as.character(utils::packageVersion("gamlss.longitudinal")) ||
      manifest$dependency_versions[[1L]] != jss_misspec_dependency_versions() ||
      manifest$dependency_versions_sha256[[1L]] != jss_misspec_sha256_object(manifest$dependency_versions[[1L]]) ||
      manifest$workers_requested[[1L]] < 1L || manifest$base_seed[[1L]] != config$seed ||
      manifest$seed_min[[1L]] != min(grid$seed) || manifest$seed_max[[1L]] != max(grid$seed) ||
      manifest$dataset_seed_count[[1L]] != length(unique(grid$dataset_seed)) ||
      manifest$fit_seed_count[[1L]] != length(unique(grid$seed)) ||
      manifest$planned_fit_rows[[1L]] != nrow(grid) || manifest$result_rows[[1L]] != nrow(results) ||
      manifest$producer_sha256[[1L]] != config$producer_sha256 ||
      manifest$package_source_sha256[[1L]] != config$package_source_sha256 ||
      manifest$configuration_sha256[[1L]] != config$configuration_sha256 ||
      manifest$results_sha256[[1L]] != jss_misspec_sha256_file(paths[["results"]]) ||
      manifest$paired_effects_sha256[[1L]] != jss_misspec_sha256_file(paths[["paired_effects"]]) ||
      manifest$warning_audit_sha256[[1L]] != jss_misspec_sha256_file(paths[["warning_audit"]]) ||
      manifest$selection_attempts_sha256[[1L]] != jss_misspec_sha256_file(paths[["selection_attempts"]]) ||
      manifest$selection_sha256[[1L]] != jss_misspec_sha256_file(paths[["selection"]]) ||
      manifest$expected_warning_events[[1L]] != warnings$n_warning_events[warnings$warning_classification == "expected_numerical_fallback"] ||
      manifest$nonconvergence_warning_events[[1L]] != warnings$n_warning_events[warnings$warning_classification == "expected_nonconvergence_signal"] ||
      manifest$unexpected_warning_events[[1L]] != 0L ||
      any(!grepl("^[0-9a-f]{64}$", unlist(manifest[hash_fields], use.names = FALSE)))) {
    stop("Module 07 execution manifest is stale, incomplete, or not bound to generation-time source/runtime evidence.", call. = FALSE)
  }
  list(paths = paths, results = results, paired_effects = paired, warning_audit = warnings,
    selection_attempts = attempts, selection = selection, execution_manifest = manifest)
}

jss_misspec_candidate_identity <- function(results_path, config) {
  if (!identical(config$stage, "full") || !identical(as.integer(config$reps), 100L)) {
    stop("Module 07 candidate identity requires the immutable full R=100 configuration.", call. = FALSE)
  }
  bundle <- jss_misspec_validate_evidence_bundle(results_path, config)
  manifest_entries <- paste(names(bundle$paths), vapply(bundle$paths, jss_misspec_sha256_file, character(1L)), sep = "\t")
  data.frame(
    approval_schema_version = 3L, profile = "paper", status = "candidate_pending_independent_approval",
    results_sha256 = jss_misspec_sha256_file(results_path),
    results_rows = as.integer(nrow(bundle$results)),
    bundle_sha256 = jss_misspec_sha256_object(manifest_entries),
    execution_manifest_sha256 = jss_misspec_sha256_file(bundle$paths[["execution_manifest"]]),
    configuration_sha256 = jss_misspec_sha256_object(jss_misspec_registered_config_identity(config)),
    producer_sha256 = jss_misspec_sha256_file(config$source_path),
    package_source_sha256 = jss_misspec_package_source_sha256(config$repo_root),
    approved_at_utc = "", approver = "", stringsAsFactors = FALSE
  )
}

jss_misspec_verify_signed_approval <- function(identity, config,
    attestation_path = jss_misspec_attestation_path(),
    signature_path = jss_misspec_signature_path()) {
  if (!requireNamespace("sodium", quietly = TRUE)) stop("sodium is required for Module 07 approval verification.", call. = FALSE)
  if (!nzchar(attestation_path) || !nzchar(signature_path)) {
    stop("Module 07 full R=100 bundle lacks a detached production approval signature.", call. = FALSE)
  }
  root <- paste0(tolower(normalizePath(config$repo_root, winslash = "/", mustWork = TRUE)), "/")
  paths <- vapply(c(attestation_path, signature_path), normalizePath, character(1L),
    winslash = "/", mustWork = TRUE)
  if (any(startsWith(tolower(paste0(paths, "/")), root))) {
    stop("Module 07 approval attestation/signature must be external to the checkout.", call. = FALSE)
  }
  message_raw <- readBin(paths[[1L]], "raw", n = file.info(paths[[1L]])$size)
  signature_raw <- readBin(paths[[2L]], "raw", n = file.info(paths[[2L]])$size)
  if (!isTRUE(tryCatch(sodium::sig_verify(message_raw, signature_raw,
      jss_misspec_pinned_public_key()), error = function(e) FALSE))) {
    stop("Module 07 approval lacks a valid detached production signature.", call. = FALSE)
  }
  x <- tryCatch(unserialize(message_raw), error = function(e) NULL)
  expected <- c(
    "schema_version", "study", "results_sha256", "results_rows",
    "bundle_sha256", "execution_manifest_sha256", "configuration_sha256",
    "producer_sha256", "package_source_sha256", "approved_at_utc", "approver"
  )
  if (!is.list(x) || !identical(names(x), expected) || !identical(x$schema_version, 3L) ||
      !identical(x$study, "copula-misspecification") ||
      !identical(x$results_sha256, identity$results_sha256[[1L]]) ||
      !identical(as.integer(x$results_rows), as.integer(nrow(jss_misspec_grid(config)))) ||
      !identical(as.integer(x$results_rows), identity$results_rows[[1L]]) ||
      !identical(x$bundle_sha256, identity$bundle_sha256[[1L]]) ||
      !identical(x$execution_manifest_sha256, identity$execution_manifest_sha256[[1L]]) ||
      !identical(x$configuration_sha256, identity$configuration_sha256[[1L]]) ||
      !identical(x$producer_sha256, identity$producer_sha256[[1L]]) ||
      !identical(x$package_source_sha256, identity$package_source_sha256[[1L]]) ||
      !is.character(x$approver) || length(x$approver) != 1L || !nzchar(trimws(x$approver))) {
    stop("Module 07 signed approval does not bind the exact full bundle/source/provenance.", call. = FALSE)
  }
  jss_misspec_parse_rfc3339_utc(x$approved_at_utc)
  invisible(x)
}

jss_misspec_validate_approved_public_bundle <- function(
    results_path, grid, config,
    attestation_path = jss_misspec_attestation_path(),
    signature_path = jss_misspec_signature_path()) {
  source_paths <- jss_misspec_evidence_paths(results_path)
  if (any(!file.exists(source_paths))) stop("Module 07 approved evidence bundle is incomplete.", call. = FALSE)
  before <- vapply(source_paths, jss_misspec_sha256_file, character(1L))
  snapshot_dir <- tempfile("module07-approved-snapshot-"); dir.create(snapshot_dir)
  on.exit(unlink(snapshot_dir, recursive = TRUE, force = TRUE), add = TRUE)
  snapshot_paths <- file.path(snapshot_dir, basename(source_paths)); names(snapshot_paths) <- names(source_paths)
  if (!all(file.copy(source_paths, snapshot_paths, overwrite = FALSE)) ||
      !identical(before, vapply(source_paths, jss_misspec_sha256_file, character(1L))) ||
      !identical(before, vapply(snapshot_paths, jss_misspec_sha256_file, character(1L)))) {
    stop("Module 07 source changed while creating its immutable validation snapshot.", call. = FALSE)
  }
  snapshot <- snapshot_paths[["results"]]
  bundle <- jss_misspec_validate_evidence_bundle(snapshot, config)
  results <- bundle$results
  identity <- jss_misspec_candidate_identity(snapshot, config)
  approval <- jss_misspec_verify_signed_approval(identity, config, attestation_path, signature_path)
  if (!identical(before, vapply(source_paths, jss_misspec_sha256_file, character(1L))) ||
      !identical(before, vapply(snapshot_paths, jss_misspec_sha256_file, character(1L)))) {
    stop("Module 07 source or immutable snapshot changed during approval validation.", call. = FALSE)
  }
  identity$status <- "approved"
  identity$approved_at_utc <- approval$approved_at_utc
  identity$approver <- approval$approver
  invisible(c(identity[1L, ], list(
    results = results, source_sha256 = before[["results"]], source_bundle_sha256 = before,
    attestation_path = normalizePath(attestation_path, winslash = "/", mustWork = TRUE),
    signature_path = normalizePath(signature_path, winslash = "/", mustWork = TRUE)
  )))
}

jss_misspec_revalidate_approved_source <- function(results_path, approval, config) {
  paths <- jss_misspec_evidence_paths(results_path)
  current <- vapply(paths, jss_misspec_sha256_file, character(1L))
  if (!identical(current, approval$source_bundle_sha256)) {
    stop("Module 07 approved source changed after immutable validation.", call. = FALSE)
  }
  identity_names <- names(jss_misspec_candidate_identity(results_path, config))
  identity <- as.data.frame(approval[intersect(names(approval), identity_names)], stringsAsFactors = FALSE)
  jss_misspec_verify_signed_approval(identity, config, approval$attestation_path, approval$signature_path)
  invisible(TRUE)
}

jss_misspec_replace_file_atomic <- function(temporary, path) {
  if (!file.exists(path)) {
    if (!file.rename(temporary, path)) stop("Could not atomically install CSV: ", path, call. = FALSE)
    return(invisible(path))
  }
  backup <- tempfile(paste0(".", basename(path), "-previous-"), tmpdir = dirname(path))
  if (.Platform$OS.type == "windows") {
    literal <- function(x) paste0("'", gsub("'", "''", normalizePath(x, winslash = "/", mustWork = FALSE), fixed = TRUE), "'")
    command <- paste0("[System.IO.File]::Replace(", literal(temporary), ",", literal(path), ",", literal(backup), ",$true)")
    output <- suppressWarnings(system2("powershell", c("-NoProfile", "-NonInteractive", "-Command", shQuote(command)), stdout = TRUE, stderr = TRUE))
    exit_status <- attr(output, "status")
    if (is.null(exit_status)) exit_status <- 0L
    if (!identical(as.integer(exit_status), 0L) || !file.exists(path)) {
      stop("Windows atomic Module 07 CSV replacement failed: ", paste(output, collapse = " | "), call. = FALSE)
    }
  } else if (!file.rename(temporary, path)) {
    stop("Could not atomically replace Module 07 CSV: ", path, call. = FALSE)
  }
  if (file.exists(backup)) unlink(backup, force = TRUE)
  invisible(path)
}

jss_misspec_write_csv_atomic <- function(x, path) {
  if (!identical(tolower(tools::file_ext(path)), "csv")) stop("Refusing non-CSV output path: ", path, call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  utils::write.csv(x, temporary, row.names = FALSE)
  if (jss_misspec_is_png(temporary)) stop("CSV serializer produced PNG bytes: ", path, call. = FALSE)
  expected_hash <- jss_misspec_sha256_file(temporary)
  jss_misspec_replace_file_atomic(temporary, path)
  if (jss_misspec_is_png(path)) stop("Results CSV contains PNG bytes: ", path, call. = FALSE)
  installed <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(jss_misspec_sha256_file(path), expected_hash) ||
      !identical(names(installed), names(x)) || nrow(installed) != nrow(x)) {
    stop("Installed Module 07 CSV failed hash/schema re-verification: ", path, call. = FALSE)
  }
  invisible(path)
}

jss_misspec_registered_grid_fields <- function() {
  c(
    "fit_id", "margin_family", "capability_registry_version",
    "producer_sha256", "package_source_sha256", "configuration_sha256",
    "generating_copula", "fitted_copula", "tau_label", "target_tau",
    "n_subject", "n_time", "rep", "seed", "dataset_seed"
  )
}

jss_misspec_required_result_fields <- function() {
  c(
    jss_misspec_registered_grid_fields(), "stage", "correctly_specified",
    "attempted", "success", "converged", "retained", "stop_reason",
    "failure_type", "error", "warnings", "elapsed_sec"
  )
}

jss_misspec_seed_base_candidates <- function(results) {
  as.integer(results$seed) - mapply(
    jss_misspec_seed,
    results$generating_copula, results$fitted_copula, results$tau_label,
    results$n_subject, results$rep,
    MoreArgs = list(base_seed = 0L), SIMPLIFY = TRUE, USE.NAMES = FALSE
  )
}

jss_misspec_upgrade_result_contract <- function(results) {
  if (!nrow(results)) return(results)
  if (!"attempted" %in% names(results)) results$attempted <- TRUE
  if (!"retained" %in% names(results)) results$retained <- results$success %in% TRUE & results$converged %in% TRUE
  if (!"stop_reason" %in% names(results)) {
    results$stop_reason <- ifelse(results$retained %in% TRUE, "converged", as.character(results$failure_type))
  }
  if (!"dataset_seed" %in% names(results) && all(c(
    "generating_copula", "fitted_copula", "tau_label", "n_subject", "rep", "seed"
  ) %in% names(results))) {
    candidates <- jss_misspec_seed_base_candidates(results)
    if (length(unique(candidates)) == 1L) {
      results$dataset_seed <- mapply(
        jss_misspec_dataset_seed,
        results$generating_copula, results$tau_label, results$n_subject, results$rep,
        MoreArgs = list(base_seed = candidates[[1L]]), SIMPLIFY = TRUE, USE.NAMES = FALSE
      )
    } else {
      results$dataset_seed <- NA_integer_
    }
  }
  results
}

jss_misspec_result_contract_issues <- function(results, grid = NULL, config = NULL) {
  issues <- character()
  required <- jss_misspec_required_result_fields()
  missing <- setdiff(required, names(results))
  if (!nrow(results)) issues <- c(issues, "no result rows")
  if (length(missing)) return(c(issues, paste0("missing fields: ", paste(missing, collapse = ", "))))
  expected_id <- mapply(
    jss_misspec_fit_id,
    results$generating_copula, results$fitted_copula, results$tau_label,
    results$n_subject, results$rep, SIMPLIFY = TRUE, USE.NAMES = FALSE
  )
  if (!identical(as.character(results$fit_id), as.character(expected_id))) issues <- c(issues, "fit_id is not canonical for registered grid fields")
  if (anyDuplicated(results$fit_id)) issues <- c(issues, "fit_id is not unique")
  expected_correct <- as.character(results$generating_copula) == as.character(results$fitted_copula)
  if (any(is.na(results$correctly_specified)) || any((results$correctly_specified %in% TRUE) != expected_correct)) {
    issues <- c(issues, "correctly_specified contradicts generating/fitted copula")
  }
  if (!is.null(config) && (any(is.na(results$stage)) || any(as.character(results$stage) != as.character(config$stage)))) {
    issues <- c(issues, "stage does not match the authoritative registered configuration")
  }
  if (any(!as.character(results$generating_copula) %in% jss_misspec_copulas()) ||
      any(!as.character(results$fitted_copula) %in% jss_misspec_copulas()) ||
      any(!as.character(results$tau_label) %in% jss_misspec_tau_levels()$tau_label)) {
    issues <- c(issues, "copula or tau-level value is outside the registered grid")
  }
  if (any(is.na(results$margin_family)) || any(as.character(results$margin_family) != "GA") ||
      any(is.na(results$capability_registry_version)) ||
      any(as.character(results$capability_registry_version) != "2026.2")) {
    issues <- c(issues, "margin family or capability-registry version is outside the registered design")
  }
  if (!is.null(config) && (
      any(is.na(results$producer_sha256)) ||
      any(as.character(results$producer_sha256) != config$producer_sha256) ||
      any(is.na(results$package_source_sha256)) ||
      any(as.character(results$package_source_sha256) != config$package_source_sha256) ||
      any(is.na(results$configuration_sha256)) ||
      any(as.character(results$configuration_sha256) != config$configuration_sha256))) {
    issues <- c(issues, "result/checkpoint source or configuration fingerprint is stale")
  }
  tau_map <- stats::setNames(jss_misspec_tau_levels()$target_tau, jss_misspec_tau_levels()$tau_label)
  expected_tau <- unname(tau_map[as.character(results$tau_label)])
  actual_tau <- suppressWarnings(as.numeric(results$target_tau))
  if (any(!is.finite(actual_tau)) || any(actual_tau != expected_tau)) {
    issues <- c(issues, "target_tau contradicts the registered tau level")
  }
  base_candidates <- if (!is.null(config)) rep(as.integer(config$seed), nrow(results)) else jss_misspec_seed_base_candidates(results)
  if (any(!is.finite(base_candidates)) || length(unique(base_candidates)) != 1L) issues <- c(issues, "fit seeds do not share one registered base seed")
  base_seed <- if (length(base_candidates)) base_candidates[[1L]] else NA_integer_
  expected_fit_seed <- mapply(
    jss_misspec_seed,
    results$generating_copula, results$fitted_copula, results$tau_label,
    results$n_subject, results$rep,
    MoreArgs = list(base_seed = base_seed), SIMPLIFY = TRUE, USE.NAMES = FALSE
  )
  expected_dataset_seed <- mapply(
    jss_misspec_dataset_seed,
    results$generating_copula, results$tau_label, results$n_subject, results$rep,
    MoreArgs = list(base_seed = base_seed), SIMPLIFY = TRUE, USE.NAMES = FALSE
  )
  if (any(as.integer(results$seed) != expected_fit_seed, na.rm = TRUE) || anyNA(results$seed)) issues <- c(issues, "fit seed mismatch")
  if (any(as.integer(results$dataset_seed) != expected_dataset_seed, na.rm = TRUE) || anyNA(results$dataset_seed)) issues <- c(issues, "dataset seed mismatch")
  attempted <- results$attempted %in% TRUE
  success <- results$success %in% TRUE
  converged <- results$converged %in% TRUE
  retained <- results$retained %in% TRUE
  if (!all(attempted) || any(success != retained) || any(retained & !converged)) issues <- c(issues, "attempt/success/convergence/retention semantics contradict")
  if (any(retained & as.character(results$stop_reason) != "converged") ||
      any(!retained & (!nzchar(trimws(as.character(results$stop_reason))) | is.na(results$stop_reason)))) {
    issues <- c(issues, "stop reason semantics contradict")
  }
  failure_type <- as.character(results$failure_type)
  error_text <- as.character(results$error)
  successful_lifecycle <- retained & failure_type == "none" &
    (is.na(error_text) | !nzchar(trimws(error_text))) & converged &
    as.character(results$stop_reason) == "converged"
  failed_lifecycle <- !retained & !converged & !is.na(failure_type) &
    nzchar(trimws(failure_type)) & failure_type != "none" &
    !is.na(error_text) & nzchar(trimws(error_text)) &
    as.character(results$stop_reason) == failure_type
  if (any(retained & !successful_lifecycle) || any(!retained & !failed_lifecycle)) {
    issues <- c(issues, "success/failure lifecycle truth table contradicts failure_type, error, convergence, or stop_reason")
  }
  if (any(!is.finite(suppressWarnings(as.numeric(results$elapsed_sec))) | suppressWarnings(as.numeric(results$elapsed_sec)) < 0)) {
    issues <- c(issues, "elapsed_sec must be finite and nonnegative")
  }
  accuracy <- intersect(c(
    "marginal_loglik", "copula_loglik", "joint_loglik", "aic_joint", "bic_joint",
    "fitted_mu", "fitted_sigma", "mu_bias", "sigma_bias", "margin_param_rmse",
    "fitted_copula_tau", "tau_abs_error", grep("^benchmark_", names(results), value = TRUE)
  ), names(results))
  if (any(!retained) && length(accuracy) && any(vapply(
    results[!retained, accuracy, drop = FALSE],
    function(x) any(is.finite(suppressWarnings(as.numeric(x)))), logical(1L)
  ))) issues <- c(issues, "nonretained fit contributes accuracy evidence")
  if (!is.null(grid)) {
    fields <- jss_misspec_registered_grid_fields()
    if (!all(fields %in% names(grid))) {
      issues <- c(issues, "grid lacks registered binding fields")
    } else {
      left <- grid[order(grid$fit_id), fields, drop = FALSE]
      right <- results[order(results$fit_id), fields, drop = FALSE]
      rownames(left) <- rownames(right) <- NULL
      for (name in fields) {
        if (name %in% c("target_tau")) {
          if (!identical(as.numeric(left[[name]]), as.numeric(right[[name]]))) issues <- c(issues, paste0("grid binding mismatch: ", name))
        } else if (!identical(as.character(left[[name]]), as.character(right[[name]]))) {
          issues <- c(issues, paste0("grid binding mismatch: ", name))
        }
      }
      if (nrow(left) != nrow(right)) issues <- c(issues, "grid/result cardinality mismatch")
    }
  }
  unique(issues)
}

jss_misspec_checkpoint_issues <- function(path, grid_row, config = NULL) {
  if (!file.exists(path)) return("missing")
  if (jss_misspec_is_png(path)) return("PNG bytes in CSV checkpoint")
  row <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) e)
  if (inherits(row, "error")) return(paste0("read_error: ", conditionMessage(row)))
  if (nrow(row) != 1L) return(paste0("expected one row, found ", nrow(row)))
  jss_misspec_result_contract_issues(row, grid = grid_row, config = config)
}

jss_misspec_checkpoint_valid <- function(path, grid_row, config = NULL) {
  length(jss_misspec_checkpoint_issues(path, grid_row, config)) == 0L
}

jss_misspec_binding_review_gate <- function(results, grid = NULL, paths, context, config = NULL) {
  results <- jss_misspec_upgrade_result_contract(results)
  required <- jss_misspec_required_result_fields()
  installed <- tryCatch(suppressWarnings(utils::read.csv(paths$results, stringsAsFactors = FALSE, check.names = FALSE)), error = function(e) e)
  semantic_issues <- jss_misspec_result_contract_issues(results, grid = grid, config = config)
  rows <- list(
    data.frame(check = "results_schema", status = if (nrow(results) > 0L && all(required %in% names(results))) "pass" else "fail", detail = paste("rows", nrow(results)), stringsAsFactors = FALSE),
    data.frame(check = "results_unique", status = if ("fit_id" %in% names(results) && !anyDuplicated(results$fit_id)) "pass" else "fail", detail = "fit_id must be unique", stringsAsFactors = FALSE),
    data.frame(check = "results_path_is_csv", status = if (file.exists(paths$results) && !jss_misspec_is_png(paths$results)) "pass" else "fail", detail = paths$results, stringsAsFactors = FALSE),
    data.frame(check = "png_paths_are_distinct", status = if (!paths$results %in% unname(unlist(paths[grepl("heatmap|convergence", names(paths))]))) "pass" else "fail", detail = "results CSV cannot alias a figure path", stringsAsFactors = FALSE),
    data.frame(check = "semantic_binding", status = if (!length(semantic_issues)) "pass" else "fail", detail = if (length(semantic_issues)) paste(semantic_issues, collapse = " | ") else "canonical fields, seeds, and attempt semantics verified", stringsAsFactors = FALSE),
    data.frame(
      check = "installed_hash_schema_binding",
      status = if (!inherits(installed, "error") && identical(names(installed), names(results)) && nrow(installed) == nrow(results) &&
        identical(jss_misspec_sha256_file(paths$results), local({
          tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
          utils::write.csv(results, tmp, row.names = FALSE); jss_misspec_sha256_file(tmp)
        }))) "pass" else "fail",
      detail = "installed result was reread and SHA-256/schema compared", stringsAsFactors = FALSE
    )
  )
  if (!is.null(grid)) rows[[length(rows) + 1L]] <- data.frame(
    check = "grid_binding", status = if (!any(grepl("grid binding|cardinality", semantic_issues))) "pass" else "fail",
    detail = paste("exact registered-field join; grid", nrow(grid), "results", nrow(results)), stringsAsFactors = FALSE
  )
  review <- do.call(rbind, rows)
  review$context <- context
  if (any(review$status != "pass")) {
    stop("Fatal Module 07 binding review gate failed: ", paste(review$check[review$status != "pass"], collapse = ", "), call. = FALSE)
  }
  review
}

jss_misspec_simulate <- function(generating_copula, target_tau, n_subject, config, seed) {
  gamlss.longitudinal::simulate_longitudinal_dataset(
    n = n_subject,
    times = config$times,
    margin_dist = gamlss.dist::GA(mu.link = "log", sigma.link = "log"),
    copula_dist = generating_copula,
    margin_params = list(mu = config$margin_mu, sigma = config$margin_sigma),
    copula_params = list(tau = target_tau),
    seed = seed,
    include_truth = TRUE
  )
}

jss_misspec_fit_one_in_process <- function(dat, row, config) {
  warnings <- character(0)
  start <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(
    tryCatch(
      gamlss.longitudinal::gamlss_longitudinal(
        dataset = dat,
        margin_dist = gamlss.dist::GA(mu.link = "log", sigma.link = "log"),
        copula_dist = row$fitted_copula,
        time_var = "time",
        subject_var = "subject",
        mu.formula = response ~ 1,
        sigma.formula = ~1,
        nu.formula = ~1,
        tau.formula = ~1,
        theta.formula = ~1,
        zeta.formula = ~1,
        method = "RS",
        optimizer_control = gamlss.longitudinal::gamlss_longitudinal_control(
          max_outer_iter = config$max_outer_iter,
          max_elapsed_sec = config$max_elapsed_sec,
          rs = list(max_inner_iter = config$max_inner_iter)
        ),
        compute_vcov = config$compute_vcov,
        verbose = 0
      ),
      error = function(e) e
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  elapsed <- proc.time()[["elapsed"]] - start

  jss_misspec_result_row(
    fit = fit,
    dat = dat,
    row = row,
    config = config,
    elapsed = elapsed,
    warnings = warnings
  )
}

jss_misspec_killable_call <- function(fun_name, args, config) {
  jss_misspec_validate_config(config, require_callr = TRUE, verify_source = FALSE)
  tryCatch(
    callr::r(
      function(fun_name, args, config) {
        setwd(config$repo_root)
        source(config$source_path, local = .GlobalEnv)
        if (!identical(jss_misspec_sha256_file(config$source_path), config$producer_sha256)) {
          stop("Module 07 killable subprocess producer source changed after configuration.")
        }
        jss_misspec_load_checkout(config$repo_root, config$package_source_sha256)
        do.call(get(fun_name, envir = .GlobalEnv), args)
      },
      args = list(fun_name = fun_name, args = args, config = config),
      timeout = as.numeric(config$max_elapsed_sec),
      spinner = FALSE
    ),
    error = function(e) e
  )
}

jss_misspec_timeout_probe <- function(seconds) {
  Sys.sleep(seconds)
  TRUE
}

jss_misspec_fit_one <- function(dat, row, config) {
  jss_misspec_validate_config(config, require_callr = TRUE, verify_source = FALSE)
  start <- proc.time()[["elapsed"]]
  child_config <- config
  child_config$max_elapsed_sec <- max(0.1, as.numeric(config$max_elapsed_sec) * 0.9)
  result <- jss_misspec_killable_call(
    "jss_misspec_fit_one_in_process",
    list(dat = dat, row = row, config = child_config),
    config = config
  )
  if (!inherits(result, "error")) return(result)
  elapsed <- proc.time()[["elapsed"]] - start
  message <- conditionMessage(result)
  timeout <- grepl("timed out|timeout", message, ignore.case = TRUE)
  failure <- simpleError(if (timeout) {
    paste0("Killable Module 07 fit timed out after ", config$max_elapsed_sec, " seconds.")
  } else {
    paste0("Killable Module 07 subprocess failed: ", message)
  })
  out <- jss_misspec_result_row(
    fit = failure, dat = dat, row = row, config = config,
    elapsed = elapsed, warnings = character()
  )
  if (timeout) {
    out$failure_type <- "timeout"
    out$stop_reason <- "timeout"
    out$error <- conditionMessage(failure)
  }
  out
}

jss_misspec_result_row <- function(fit, dat, row, config, elapsed, warnings) {
  is_fit <- inherits(fit, "gamlss.longitudinal")
  loglik <- if (is_fit) fit$calc_lik_out_end$log_lik else c(marginal = NA_real_, copula = NA_real_, joint = NA_real_)
  estimates <- if (is_fit) gamlss.longitudinal:::.coverage_natural_estimates(fit) else stats::setNames(numeric(0), character(0))
  fitted_tau <- if (is_fit && "theta" %in% names(estimates)) {
    zeta <- if ("zeta" %in% names(estimates)) estimates[["zeta"]] else 0
    suppressWarnings(gamlss.longitudinal:::.copula_par_to_tau(row$fitted_copula, estimates[["theta"]], zeta))
  } else {
    NA_real_
  }
  k <- if (is_fit) length(fit$par) else NA_integer_
  n_obs <- sum(is.finite(dat$response))
  joint_loglik <- as.numeric(loglik["joint"])
  aic_joint <- if (is.finite(joint_loglik) && is.finite(k)) -2 * joint_loglik + 2 * k else NA_real_
  bic_joint <- if (is.finite(joint_loglik) && is.finite(k)) -2 * joint_loglik + log(n_obs) * k else NA_real_
  truth_metrics <- if (is_fit) {
    gamlss.longitudinal:::.coverage_benchmark_gamlss_metrics(dat, fit, "GA")
  } else {
    gamlss.longitudinal:::.coverage_benchmark_gamlss_metrics(dat, NULL, "GA")
  }
  ci <- jss_misspec_margin_ci(fit)

  converged <- is_fit && !is.null(fit$convergence$converged) && isTRUE(fit$convergence$converged)
  success <- is_fit && converged && all(is.finite(as.numeric(loglik)))
  retained <- success
  failure_type <- if (retained) "none" else if (inherits(fit, "error")) "error" else if (is_fit && !converged) "nonconverged" else "nonfinite_or_no_fit"
  error_text <- if (retained) {
    NA_character_
  } else if (inherits(fit, "error")) {
    conditionMessage(fit)
  } else if (is_fit && !converged) {
    "Optimizer did not converge."
  } else {
    "Fit was unavailable or produced non-finite likelihood values."
  }
  out <- data.frame(
    fit_id = row$fit_id,
    margin_family = row$margin_family,
    capability_registry_version = row$capability_registry_version,
    producer_sha256 = row$producer_sha256,
    package_source_sha256 = row$package_source_sha256,
    configuration_sha256 = row$configuration_sha256,
    stage = config$stage,
    generating_copula = row$generating_copula,
    fitted_copula = row$fitted_copula,
    correctly_specified = identical(row$generating_copula, row$fitted_copula),
    tau_label = row$tau_label,
    target_tau = row$target_tau,
    n_subject = row$n_subject,
    n_time = row$n_time,
    rep = row$rep,
    seed = row$seed,
    dataset_seed = row$dataset_seed,
    attempted = TRUE,
    success = success,
    converged = converged,
    retained = retained,
    stop_reason = if (retained) "converged" else failure_type,
    failure_type = failure_type,
    error = error_text,
    warnings = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
    marginal_loglik = as.numeric(loglik["marginal"]),
    copula_loglik = as.numeric(loglik["copula"]),
    joint_loglik = joint_loglik,
    aic_joint = aic_joint,
    bic_joint = bic_joint,
    elapsed_sec = elapsed,
    n_parameters = k,
    fitted_mu = if ("mu" %in% names(estimates)) unname(estimates[["mu"]]) else NA_real_,
    fitted_sigma = if ("sigma" %in% names(estimates)) unname(estimates[["sigma"]]) else NA_real_,
    mu_bias = if ("mu" %in% names(estimates)) unname(estimates[["mu"]] - config$margin_mu) else NA_real_,
    sigma_bias = if ("sigma" %in% names(estimates)) unname(estimates[["sigma"]] - config$margin_sigma) else NA_real_,
    margin_param_rmse = jss_misspec_margin_rmse(estimates, config),
    mu_covered_95 = ci$mu_covered_95,
    sigma_covered_95 = ci$sigma_covered_95,
    fitted_copula_tau = as.numeric(fitted_tau)[1],
    tau_abs_error = abs(as.numeric(fitted_tau)[1] - row$target_tau),
    benchmark_mae = unname(truth_metrics[["benchmark_mae"]]),
    benchmark_rmse = unname(truth_metrics[["benchmark_rmse"]]),
    benchmark_mean_bias = unname(truth_metrics[["benchmark_mean_bias"]]),
    benchmark_mean_mae = unname(truth_metrics[["benchmark_mean_mae"]]),
    benchmark_mean_rmse = unname(truth_metrics[["benchmark_mean_rmse"]]),
    benchmark_q90_mae = unname(truth_metrics[["benchmark_q90_mae"]]),
    benchmark_neg_log_score = unname(truth_metrics[["benchmark_neg_log_score"]]),
    benchmark_upper_tail_error_90 = unname(truth_metrics[["benchmark_upper_tail_error_90"]]),
    benchmark_interval_coverage_95 = unname(truth_metrics[["benchmark_interval_coverage_95"]]),
    benchmark_interval_width_95 = unname(truth_metrics[["benchmark_interval_width_95"]]),
    benchmark_pit_ks_p_value = unname(truth_metrics[["benchmark_pit_ks_p_value"]]),
    benchmark_pit_mean_abs_error = unname(truth_metrics[["benchmark_pit_mean_abs_error"]]),
    benchmark_tail_error_lower_05 = unname(truth_metrics[["benchmark_tail_error_lower_05"]]),
    benchmark_tail_error_upper_05 = unname(truth_metrics[["benchmark_tail_error_upper_05"]]),
    stringsAsFactors = FALSE
  )
  if (!retained) {
    accuracy <- intersect(c(
      "marginal_loglik", "copula_loglik", "joint_loglik", "aic_joint", "bic_joint",
      "n_parameters", "fitted_mu", "fitted_sigma", "mu_bias", "sigma_bias",
      "margin_param_rmse", "mu_covered_95", "sigma_covered_95", "fitted_copula_tau",
      "tau_abs_error", grep("^benchmark_", names(out), value = TRUE)
    ), names(out))
    for (name in accuracy) out[[name]] <- NA
  }
  out
}

jss_misspec_margin_rmse <- function(estimates, config) {
  err <- c(
    if ("mu" %in% names(estimates)) estimates[["mu"]] - config$margin_mu else NA_real_,
    if ("sigma" %in% names(estimates)) estimates[["sigma"]] - config$margin_sigma else NA_real_
  )
  if (!any(is.finite(err))) NA_real_ else sqrt(mean(err[is.finite(err)]^2))
}

jss_misspec_margin_ci <- function(fit) {
  empty <- list(mu_covered_95 = NA, sigma_covered_95 = NA)
  if (!inherits(fit, "gamlss.longitudinal")) {
    return(empty)
  }
  ci <- tryCatch(
    stats::confint(fit, parm = c("mu.intercept", "sigma.intercept"), level = 0.95),
    error = function(e) NULL
  )
  if (is.null(ci)) {
    return(empty)
  }
  empty$mu_covered_95 <- isTRUE(ci["mu.intercept", 1] <= log(2) && ci["mu.intercept", 2] >= log(2))
  empty$sigma_covered_95 <- isTRUE(ci["sigma.intercept", 1] <= log(0.5) && ci["sigma.intercept", 2] >= log(0.5))
  empty
}

jss_misspec_atomic_csv <- function(x, path) jss_misspec_write_csv_atomic(x, path)

jss_misspec_run_scenario_checkpoints <- function(idx, pending, config) {
  scenario <- pending[idx[1L], , drop = FALSE]
  dat_seed <- as.integer(scenario$dataset_seed)
  dat <- jss_misspec_simulate(
    generating_copula = scenario$generating_copula,
    target_tau = scenario$target_tau,
    n_subject = scenario$n_subject,
    config = config,
    seed = dat_seed
  )
  for (j in idx) {
    fit_row <- pending[j, , drop = FALSE]
    out <- jss_misspec_fit_one(dat, fit_row, config)
    jss_misspec_atomic_csv(out, fit_row$checkpoint_path)
    if (!jss_misspec_checkpoint_valid(fit_row$checkpoint_path, fit_row, config)) {
      stop("Installed Module 07 checkpoint failed semantic revalidation: ", fit_row$fit_id, call. = FALSE)
    }
  }
  length(idx)
}

jss_misspec_prepare_workers <- function(cl, config,
    source_env = environment(jss_misspec_run_checkpoints)) {
  functions <- ls(pattern = "^jss_misspec_", envir = source_env)
  parallel::clusterExport(cl, functions, envir = source_env)
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(gamlss.dist))
    suppressPackageStartupMessages(library(VineCopula))
    NULL
  })
  attestations <- parallel::clusterCall(
    cl,
    function(repo_root, expected_hash) {
      jss_misspec_load_checkout(repo_root, expected_hash)
    },
    config$repo_root,
    config$package_source_sha256
  )
  expected <- c(
    namespace_path = normalizePath(config$repo_root, winslash = "/", mustWork = TRUE),
    package_source_sha256 = config$package_source_sha256
  )
  if (!all(vapply(attestations, identical, logical(1L), y = expected))) {
    stop("One or more Module 07 PSOCK workers failed checked-out source attestation.", call. = FALSE)
  }
  invisible(attestations)
}

jss_misspec_run_checkpoints <- function(grid, config, checkpoint_dir, workers = 1L) {
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  pending <- jss_misspec_pending_grid(grid, checkpoint_dir, config = config)
  if (!nrow(pending)) {
    return(invisible(0L))
  }

  scenario_cols <- c("generating_copula", "tau_label", "target_tau", "n_subject", "n_time", "rep")
  scenario_key <- interaction(pending[scenario_cols], drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(pending)), scenario_key)
  workers <- min(max(1L, as.integer(workers)), length(groups))
  if (workers == 1L) {
    return(invisible(sum(vapply(groups, jss_misspec_run_scenario_checkpoints, integer(1), pending = pending, config = config))))
  }
  worker_log <- file.path(dirname(checkpoint_dir), paste0(basename(checkpoint_dir), "-workers.log"))
  cl <- parallel::makePSOCKcluster(workers, outfile = worker_log, setup_strategy = "parallel")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  jss_misspec_prepare_workers(cl, config)
  completed <- parallel::parLapplyLB(
    cl, groups, jss_misspec_run_scenario_checkpoints,
    pending = pending, config = config
  )
  post_attestations <- parallel::clusterCall(
    cl,
    function(repo_root, expected_hash) {
      jss_misspec_load_checkout(repo_root, expected_hash)
    },
    config$repo_root,
    config$package_source_sha256
  )
  expected_attestation <- c(
    namespace_path = normalizePath(config$repo_root, winslash = "/", mustWork = TRUE),
    package_source_sha256 = config$package_source_sha256
  )
  if (!all(vapply(post_attestations, identical, logical(1L), y = expected_attestation))) {
    stop("One or more Module 07 PSOCK workers failed post-task source attestation.", call. = FALSE)
  }
  invisible(sum(unlist(completed, use.names = FALSE)))
}

jss_misspec_read_checkpoints <- function(checkpoint_dir, grid = NULL, config = NULL) {
  if (!dir.exists(checkpoint_dir)) {
    return(data.frame())
  }
  files <- list.files(checkpoint_dir, pattern = "[.]csv$", full.names = TRUE)
  files <- files[basename(files) != "checkpoint-rejections-ledger.csv"]
  if (!length(files)) {
    return(data.frame())
  }
  character_cols <- c(
    "fit_id", "stage", "generating_copula", "fitted_copula", "tau_label",
    "failure_type", "error", "warnings"
  )
  rows <- lapply(files, function(path) {
    header <- names(utils::read.csv(path, nrows = 0, stringsAsFactors = FALSE))
    classes <- stats::setNames(rep(NA_character_, length(header)), header)
    classes[intersect(character_cols, header)] <- "character"
    utils::read.csv(path, stringsAsFactors = FALSE, colClasses = classes)
  })
  out <- jss_misspec_row_bind(rows)
  if (!is.null(grid)) {
    expected_paths <- jss_misspec_checkpoint_path(checkpoint_dir, grid$fit_id)
    invalid <- vapply(seq_len(nrow(grid)), function(i) {
      !jss_misspec_checkpoint_valid(expected_paths[[i]], grid[i, , drop = FALSE], config)
    }, logical(1L))
    if (any(invalid)) stop("Module 07 checkpoint set contains invalid or missing rows: ", paste(grid$fit_id[invalid], collapse = ", "), call. = FALSE)
    out <- out[match(grid$fit_id, out$fit_id), , drop = FALSE]
  }
  out
}

jss_misspec_add_deltas <- function(results) {
  if (!nrow(results)) {
    return(results)
  }
  derived_names <- grepl("^(correct_|delta_)", names(results))
  if (any(derived_names)) results <- results[, !derived_names, drop = FALSE]
  keys <- c("generating_copula", "tau_label", "n_subject", "rep")
  correct <- results[results$correctly_specified %in% TRUE, , drop = FALSE]
  keep <- c(keys, "joint_loglik", "copula_loglik", "tau_abs_error", "margin_param_rmse", "elapsed_sec")
  correct <- correct[keep]
  names(correct)[-(seq_along(keys))] <- paste0("correct_", names(correct)[-(seq_along(keys))])
  merged <- merge(results, correct, by = keys, all.x = TRUE, sort = FALSE)
  merged$delta_joint_loglik_vs_correct <- merged$joint_loglik - merged$correct_joint_loglik
  merged$delta_copula_loglik_vs_correct <- merged$copula_loglik - merged$correct_copula_loglik
  merged$delta_tau_abs_error_vs_correct <- merged$tau_abs_error - merged$correct_tau_abs_error
  merged$delta_margin_rmse_vs_correct <- merged$margin_param_rmse - merged$correct_margin_param_rmse
  merged$delta_elapsed_sec_vs_correct <- merged$elapsed_sec - merged$correct_elapsed_sec
  merged
}

jss_misspec_summary <- function(results) {
  if (!nrow(results)) {
    return(data.frame())
  }
  group_cols <- c("generating_copula", "fitted_copula", "tau_label", "target_tau", "n_subject")
  metric_cols <- c(
    "success", "converged", "joint_loglik", "copula_loglik", "aic_joint", "bic_joint",
    "elapsed_sec", "mu_bias", "sigma_bias", "margin_param_rmse", "tau_abs_error",
    "benchmark_neg_log_score", "benchmark_pit_mean_abs_error",
    "delta_joint_loglik_vs_correct", "delta_copula_loglik_vs_correct",
    "delta_tau_abs_error_vs_correct", "delta_margin_rmse_vs_correct",
    "delta_elapsed_sec_vs_correct"
  )
  split_key <- interaction(results[group_cols], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(seq_len(nrow(results)), split_key), function(idx) {
    row <- results[idx[1L], group_cols, drop = FALSE]
    row$n_attempted <- length(idx)
    row$success_rate <- mean(results$success[idx] %in% TRUE)
    row$convergence_rate <- mean(results$converged[idx] %in% TRUE)
    for (metric in setdiff(metric_cols, c("success", "converged"))) {
      row[[paste0("mean_", metric)]] <- mean(results[[metric]][idx], na.rm = TRUE)
    }
    row
  })
  jss_misspec_row_bind(rows)
}

jss_misspec_paired_effects <- function(results, level = 0.95) {
  if (!nrow(results)) return(data.frame())
  required <- c(
    "generating_copula", "fitted_copula", "tau_label", "target_tau",
    "n_subject", "rep", "dataset_seed", "retained", "joint_loglik",
    "margin_param_rmse", "tau_abs_error"
  )
  missing <- setdiff(required, names(results))
  if (length(missing)) stop("Module 07 paired effects lack fields: ", paste(missing, collapse = ", "), call. = FALSE)
  metrics <- c("joint_loglik", "margin_param_rmse", "tau_abs_error")
  keys <- c("generating_copula", "tau_label", "target_tau", "n_subject", "rep", "dataset_seed")
  correct <- results[results$generating_copula == results$fitted_copula,
    c(keys, "retained", metrics), drop = FALSE]
  wrong <- results[results$generating_copula != results$fitted_copula,
    c(keys, "fitted_copula", "retained", metrics), drop = FALSE]
  names(correct)[names(correct) == "retained"] <- "rhs_retained"
  names(wrong)[names(wrong) == "retained"] <- "lhs_retained"
  names(correct)[match(metrics, names(correct))] <- paste0("rhs_", metrics)
  names(wrong)[match(metrics, names(wrong))] <- paste0("lhs_", metrics)
  pairs <- merge(wrong, correct, by = keys, all = TRUE, sort = FALSE)
  if (nrow(pairs) != nrow(wrong) || anyDuplicated(pairs[c(keys, "fitted_copula")])) {
    stop("Module 07 paired-effect join is not one-to-one on generated datasets.", call. = FALSE)
  }
  group_cols <- c("generating_copula", "fitted_copula", "tau_label", "target_tau", "n_subject")
  group_key <- interaction(pairs[group_cols], drop = TRUE, lex.order = TRUE)
  rows <- list(); k <- 1L
  for (idx in split(seq_len(nrow(pairs)), group_key)) {
    for (metric in metrics) {
      lhs <- suppressWarnings(as.numeric(pairs[[paste0("lhs_", metric)]][idx]))
      rhs <- suppressWarnings(as.numeric(pairs[[paste0("rhs_", metric)]][idx]))
      retained <- pairs$lhs_retained[idx] %in% TRUE & pairs$rhs_retained[idx] %in% TRUE &
        is.finite(lhs) & is.finite(rhs)
      difference <- lhs[retained] - rhs[retained]
      n_retained <- length(difference)
      estimate <- if (n_retained) mean(difference) else NA_real_
      mcse <- if (n_retained > 1L) stats::sd(difference) / sqrt(n_retained) else NA_real_
      critical <- if (n_retained > 1L) stats::qt(1 - (1 - level) / 2, df = n_retained - 1L) else NA_real_
      rows[[k]] <- cbind(
        pairs[idx[[1L]], group_cols, drop = FALSE],
        data.frame(
          metric = metric,
          lhs_method = "misspecified_copula",
          rhs_method = "correctly_specified_copula",
          contrast = "misspecified_minus_correct_same_dataset",
          n_attempted = length(idx), n_retained = n_retained,
          estimate = estimate, mcse = mcse,
          ci_lower = estimate - critical * mcse,
          ci_upper = estimate + critical * mcse,
          confidence_level = level, ci_method = "paired_t",
          seed_pairing = "registered_dataset_seed_exact",
          stringsAsFactors = FALSE
        )
      )
      k <- k + 1L
    }
  }
  jss_misspec_row_bind(rows)
}

jss_misspec_warning_class <- function(message) {
  message <- trimws(as.character(message))
  fallback <- grepl(
    paste(c(
      "fall(ing)? back to (a )?(numerical|finite[- ]difference)",
      "(analytical|analytic).*(hessian|derivative).*(unavailable|failed|singular|invert)",
      "(hessian|information matrix).*(singular|not positive definite|could not be inverted)",
      "variance[- ]covariance.*(unavailable|could not be computed)",
      "could not precompute variance[- ]covariance"
    ), collapse = "|"),
    message, ignore.case = TRUE, perl = TRUE
  )
  nonconvergence <- grepl(
    "without satisfying the optimizer convergence contract|model (stopped|did not converge)|max_outer_iter",
    message, ignore.case = TRUE, perl = TRUE
  )
  ifelse(fallback, "expected_numerical_fallback",
    ifelse(nonconvergence, "expected_nonconvergence_signal", "unexpected"))
}

jss_misspec_warning_audit <- function(results, level = 0.95) {
  if (!all(c("fit_id", "warnings", "retained") %in% names(results))) {
    stop("Module 07 warning audit lacks fit_id, warnings, or retained fields.", call. = FALSE)
  }
  denominator <- nrow(results)
  warning_text <- as.character(results$warnings)
  has_warning <- !is.na(warning_text) & nzchar(trimws(warning_text))
  records <- lapply(which(has_warning), function(i) {
    messages <- trimws(strsplit(warning_text[[i]], " | ", fixed = TRUE)[[1L]])
    messages <- unique(messages[nzchar(messages)])
    data.frame(
      fit_id = results$fit_id[[i]], message = messages,
      classification = jss_misspec_warning_class(messages),
      retained = results$retained[[i]] %in% TRUE,
      stringsAsFactors = FALSE
    )
  })
  records <- jss_misspec_row_bind(records)
  classes <- c("expected_numerical_fallback", "expected_nonconvergence_signal", "unexpected")
  rows <- lapply(classes, function(classification) {
    sub <- if (nrow(records)) records[records$classification == classification, , drop = FALSE] else records
    affected <- if (nrow(sub)) unique(sub$fit_id) else character()
    n_affected <- length(affected)
    interval <- jss_misspec_wilson_interval(n_affected, denominator, level = level)
    data.frame(
      warning_classification = classification,
      denominator_attempts = denominator,
      n_warning_events = nrow(sub),
      n_affected_attempts = n_affected,
      warning_rate = n_affected / denominator,
      warning_rate_mcse = sqrt((n_affected / denominator) * (1 - n_affected / denominator) / denominator),
      warning_ci_lower = interval[["lower"]], warning_ci_upper = interval[["upper"]],
      confidence_level = level, ci_method = "Wilson",
      n_retained_affected = if (nrow(sub)) length(unique(sub$fit_id[sub$retained %in% TRUE])) else 0L,
      retention_policy = if (classification == "unexpected") "fatal_zero_tolerance" else if (classification == "expected_nonconvergence_signal") "must_be_nonretained" else "retain_only_if_fit_lifecycle_is_retained",
      consequence = if (classification == "unexpected") {
        "Any event fails the publication gate; classify and rerun before approval."
      } else if (classification == "expected_nonconvergence_signal") {
        "The warning is disclosed and the associated nonconverged fit must be excluded from effect estimation."
      } else {
        "Numerical fallback is disclosed; the fit contributes only when its ordinary convergence/retention contract passes."
      },
      messages = if (nrow(sub)) paste(sort(unique(sub$message)), collapse = " || ") else "<none>",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

jss_misspec_validate_warning_audit <- function(audit) {
  required <- c(
    "warning_classification", "denominator_attempts", "n_warning_events",
    "n_affected_attempts", "warning_rate", "warning_rate_mcse",
    "warning_ci_lower", "warning_ci_upper", "confidence_level", "ci_method",
    "n_retained_affected", "retention_policy", "consequence", "messages"
  )
  expected_classes <- c("expected_numerical_fallback", "expected_nonconvergence_signal", "unexpected")
  valid <- identical(names(audit), required) && nrow(audit) == length(expected_classes) &&
    identical(as.character(audit$warning_classification), expected_classes) &&
    length(unique(audit$denominator_attempts)) == 1L && audit$denominator_attempts[[1L]] > 0L &&
    all(audit$n_warning_events >= audit$n_affected_attempts) &&
    all(audit$n_affected_attempts >= audit$n_retained_affected) &&
    isTRUE(all.equal(audit$warning_rate, audit$n_affected_attempts / audit$denominator_attempts, tolerance = 0)) &&
    all(audit$warning_ci_lower <= audit$warning_rate & audit$warning_ci_upper >= audit$warning_rate) &&
    all(audit$confidence_level == 0.95) && all(audit$ci_method == "Wilson") &&
    all(nzchar(audit$retention_policy)) && all(nzchar(audit$consequence)) && all(nzchar(audit$messages))
  unexpected <- audit[audit$warning_classification == "unexpected", , drop = FALSE]
  nonconvergence <- audit[audit$warning_classification == "expected_nonconvergence_signal", , drop = FALSE]
  if (!isTRUE(valid) || nrow(unexpected) != 1L || unexpected$n_warning_events[[1L]] != 0L ||
      unexpected$n_affected_attempts[[1L]] != 0L || unexpected$retention_policy[[1L]] != "fatal_zero_tolerance" ||
      nrow(nonconvergence) != 1L || nonconvergence$n_retained_affected[[1L]] != 0L) {
    stop("Fatal Module 07 warning gate found malformed evidence or one or more unexpected warnings.", call. = FALSE)
  }
  invisible(TRUE)
}

jss_misspec_dependency_versions <- function() {
  packages <- c("gamlss.longitudinal", "gamlss", "gamlss.dist", "VineCopula", "ggplot2", "callr", "digest", "openssl", "sodium")
  versions <- vapply(packages, function(package) {
    if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else "<not-installed>"
  }, character(1L))
  paste(packages, versions, sep = "=", collapse = ";")
}

jss_misspec_git_identity <- function(root) {
  sha <- tryCatch(system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE), error = function(e) character())
  status <- tryCatch(system2("git", c("-C", shQuote(root), "status", "--porcelain", "--untracked-files=no"), stdout = TRUE, stderr = FALSE), error = function(e) "<unknown>")
  if (length(sha) != 1L || !grepl("^[0-9a-f]{40}$", sha) || identical(status, "<unknown>")) {
    stop("Module 07 could not determine its generation-time Git identity.", call. = FALSE)
  }
  list(sha = sha[[1L]], state = if (length(status)) "dirty" else "clean")
}

jss_misspec_git_commit_exists <- function(root, sha) {
  if (!is.character(sha) || length(sha) != 1L || !grepl("^[0-9a-f]{40}$", sha)) return(FALSE)
  status <- suppressWarnings(system2(
    "git", c("-C", shQuote(root), "cat-file", "-e", paste0(sha, "^{commit}")),
    stdout = FALSE, stderr = FALSE
  ))
  identical(as.integer(status), 0L)
}

jss_misspec_execution_start <- function(config, workers) {
  git <- jss_misspec_git_identity(config$repo_root)
  if (identical(config$stage, "full") && !identical(git$state, "clean")) {
    stop("Module 07 full evidence must start from a clean Git worktree.", call. = FALSE)
  }
  list(
    started_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    git_sha = git$sha, git_state = git$state,
    r_version = R.version.string, platform = R.version$platform,
    os = paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "),
    dependency_versions = jss_misspec_dependency_versions(),
    workers_requested = as.integer(workers)
  )
}

jss_misspec_execution_manifest <- function(start, config, grid, results, artifacts) {
  required <- c("results", "paired_effects", "warning_audit", "selection_attempts", "selection")
  if (!identical(sort(names(artifacts)), sort(required)) || any(!file.exists(artifacts))) {
    stop("Module 07 execution manifest requires the exact generated evidence artifacts.", call. = FALSE)
  }
  warning_audit <- utils::read.csv(artifacts[["warning_audit"]], stringsAsFactors = FALSE)
  dependency_sha <- jss_misspec_sha256_object(start$dependency_versions)
  data.frame(
    execution_schema_version = 1L, study = "copula-misspecification", stage = config$stage,
    execution_started_at_utc = start$started_at_utc,
    execution_completed_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    git_sha = start$git_sha, git_state = start$git_state,
    r_version = start$r_version, platform = start$platform, os = start$os,
    package_version = as.character(utils::packageVersion("gamlss.longitudinal")),
    dependency_versions = start$dependency_versions,
    dependency_versions_sha256 = dependency_sha,
    workers_requested = start$workers_requested,
    base_seed = as.integer(config$seed), seed_min = min(grid$seed), seed_max = max(grid$seed),
    dataset_seed_count = length(unique(grid$dataset_seed)), fit_seed_count = length(unique(grid$seed)),
    planned_fit_rows = nrow(grid), result_rows = nrow(results),
    producer_sha256 = config$producer_sha256,
    package_source_sha256 = config$package_source_sha256,
    configuration_sha256 = config$configuration_sha256,
    results_sha256 = jss_misspec_sha256_file(artifacts[["results"]]),
    paired_effects_sha256 = jss_misspec_sha256_file(artifacts[["paired_effects"]]),
    warning_audit_sha256 = jss_misspec_sha256_file(artifacts[["warning_audit"]]),
    selection_attempts_sha256 = jss_misspec_sha256_file(artifacts[["selection_attempts"]]),
    selection_sha256 = jss_misspec_sha256_file(artifacts[["selection"]]),
    expected_warning_events = sum(warning_audit$n_warning_events[warning_audit$warning_classification == "expected_numerical_fallback"]),
    nonconvergence_warning_events = sum(warning_audit$n_warning_events[warning_audit$warning_classification == "expected_nonconvergence_signal"]),
    unexpected_warning_events = sum(warning_audit$n_warning_events[warning_audit$warning_classification == "unexpected"]),
    stringsAsFactors = FALSE
  )
}

jss_misspec_wilson_interval <- function(events, attempts, level = 0.95) {
  if (!is.finite(events) || !is.finite(attempts) || attempts < 1L || events < 0L || events > attempts) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- events / attempts
  denominator <- 1 + z^2 / attempts
  centre <- (p + z^2 / (2 * attempts)) / denominator
  half_width <- z * sqrt(p * (1 - p) / attempts + z^2 / (4 * attempts^2)) / denominator
  c(lower = max(0, centre - half_width), upper = min(1, centre + half_width))
}

jss_misspec_bootstrap_median_interval <- function(x, reps = 2000L, seed = 20260901L, level = 0.95) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(c(median = NA_real_, bootstrap_se = NA_real_, lower = NA_real_, upper = NA_real_))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  boot <- replicate(as.integer(reps), stats::median(sample(x, length(x), replace = TRUE)))
  alpha <- (1 - level) / 2
  c(
    median = stats::median(x),
    bootstrap_se = if (length(boot) > 1L) stats::sd(boot) else NA_real_,
    lower = as.numeric(stats::quantile(boot, alpha, names = FALSE, type = 6)),
    upper = as.numeric(stats::quantile(boot, 1 - alpha, names = FALSE, type = 6))
  )
}

jss_misspec_selection_attempts <- function(results, close_aic_threshold = 2) {
  if (!nrow(results)) {
    return(data.frame())
  }
  required <- c(
    "generating_copula", "fitted_copula", "tau_label", "target_tau",
    "n_subject", "rep", "success", "converged", "aic_joint", "bic_joint",
    "joint_loglik", "failure_type"
  )
  missing <- setdiff(required, names(results))
  if (length(missing)) {
    stop("Copula selection attempt data are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  candidates <- jss_misspec_copulas()
  keys <- c("generating_copula", "tau_label", "target_tau", "n_subject", "rep")
  split_key <- interaction(results[keys], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(seq_len(nrow(results)), split_key), function(idx) {
    dat <- results[idx, , drop = FALSE]
    finite_criteria <- is.finite(dat$aic_joint) & is.finite(dat$bic_joint) & is.finite(dat$joint_loglik)
    usable_flag <- dat$success %in% TRUE & dat$converged %in% TRUE & finite_criteria
    usable <- dat[usable_flag, , drop = FALSE]
    base <- dat[1L, keys, drop = FALSE]
    candidate_counts <- table(factor(dat$fitted_copula, levels = candidates))
    base$n_candidates_expected <- length(candidates)
    base$n_candidate_rows <- nrow(dat)
    base$n_candidates_attempted <- sum(candidate_counts > 0L)
    base$n_candidates_usable <- sum(vapply(candidates, function(candidate) {
      any(dat$fitted_copula == candidate & usable_flag)
    }, logical(1)))
    base$n_candidates_failed <- base$n_candidates_expected - base$n_candidates_usable
    base$all_candidates_present_once <- all(candidate_counts == 1L)
    base$selection_eligible <- isTRUE(base$all_candidates_present_once) &&
      identical(base$n_candidates_usable, base$n_candidates_expected)
    missing_candidates <- candidates[candidate_counts == 0L]
    duplicate_candidates <- candidates[candidate_counts > 1L]
    base$missing_candidates <- paste(missing_candidates, collapse = ",")
    base$duplicate_candidates <- paste(duplicate_candidates, collapse = ",")
    failed <- union(dat$fitted_copula[!usable_flag], c(missing_candidates, duplicate_candidates))
    base$failed_candidates <- if (length(failed)) paste(sort(unique(failed)), collapse = ",") else ""
    unusable_idx <- which(!usable_flag)
    failure_detail <- if (length(unusable_idx)) {
      paste0(
        dat$fitted_copula[unusable_idx], "=",
        ifelse(!is.na(dat$failure_type[unusable_idx]) & nzchar(dat$failure_type[unusable_idx]), dat$failure_type[unusable_idx], "nonconverged_or_nonfinite")
      )
    } else {
      character()
    }
    structural_detail <- c(
      if (length(missing_candidates)) paste0("missing_candidate=", missing_candidates) else character(),
      if (length(duplicate_candidates)) paste0("duplicate_candidate=", duplicate_candidates, "(n=", as.integer(candidate_counts[duplicate_candidates]), ")") else character()
    )
    failure_detail <- c(failure_detail, structural_detail)
    base$failure_reasons <- if (length(failure_detail)) paste(sort(unique(failure_detail)), collapse = ";") else ""
    base$best_aic_copula <- NA_character_
    base$best_bic_copula <- NA_character_
    base$best_loglik_copula <- NA_character_
    base$aic_gap_to_runner_up <- NA_real_
    base$close_aic <- NA
    if (isTRUE(base$selection_eligible)) {
      base$best_aic_copula <- jss_misspec_best_copula(usable, "aic_joint", decreasing = FALSE)
      base$best_bic_copula <- jss_misspec_best_copula(usable, "bic_joint", decreasing = FALSE)
      base$best_loglik_copula <- jss_misspec_best_copula(usable, "joint_loglik", decreasing = TRUE)
      ordered_aic <- sort(as.numeric(usable$aic_joint), na.last = NA)
      base$aic_gap_to_runner_up <- if (length(ordered_aic) >= 2L) ordered_aic[[2L]] - ordered_aic[[1L]] else NA_real_
      base$close_aic <- is.finite(base$aic_gap_to_runner_up) && base$aic_gap_to_runner_up <= close_aic_threshold
    }
    base
  })
  jss_misspec_row_bind(rows)
}

jss_misspec_selection_confusion <- function(attempts, criterion = "aic") {
  criterion <- match.arg(criterion, c("aic", "bic", "loglik"))
  if (!nrow(attempts)) return(data.frame())
  selected_col <- paste0("best_", criterion, "_copula")
  group_cols <- c("generating_copula", "tau_label", "target_tau", "n_subject")
  outcomes <- c(jss_misspec_copulas(), "<no_selection>")
  rows <- list()
  k <- 1L
  groups <- unique(attempts[group_cols])
  for (i in seq_len(nrow(groups))) {
    idx <- rep(TRUE, nrow(attempts))
    for (col in group_cols) idx <- idx & attempts[[col]] == groups[[col]][i]
    sub <- attempts[idx, , drop = FALSE]
    selected <- as.character(sub[[selected_col]])
    selected[is.na(selected) | !nzchar(selected)] <- "<no_selection>"
    for (outcome in outcomes) {
      count <- sum(selected == outcome)
      interval <- jss_misspec_wilson_interval(count, nrow(sub))
      rows[[k]] <- cbind(
        groups[i, , drop = FALSE],
        data.frame(
          criterion = criterion,
          selected_copula = outcome,
          n_attempted = nrow(sub),
          n_selected = count,
          selection_rate = count / nrow(sub),
          selection_mcse = sqrt((count / nrow(sub)) * (1 - count / nrow(sub)) / nrow(sub)),
          selection_ci_lower = interval[["lower"]],
          selection_ci_upper = interval[["upper"]],
          stringsAsFactors = FALSE
        )
      )
      k <- k + 1L
    }
  }
  jss_misspec_row_bind(rows)
}

jss_misspec_selection_failures <- function(attempts) {
  if (!nrow(attempts)) return(data.frame())
  group_cols <- c("generating_copula", "tau_label", "target_tau", "n_subject")
  split_key <- interaction(attempts[group_cols], drop = TRUE, lex.order = TRUE)
  out <- lapply(split(seq_len(nrow(attempts)), split_key), function(idx) {
    row <- attempts[idx[1L], group_cols, drop = FALSE]
    sub <- attempts[idx, , drop = FALSE]
    row$n_attempted <- nrow(sub)
    row$n_selection_eligible <- sum(sub$selection_eligible %in% TRUE)
    row$n_no_selection <- sum(!(sub$selection_eligible %in% TRUE))
    row$n_with_any_candidate_failure <- sum(sub$n_candidates_failed > 0L)
    row$n_with_duplicate_or_missing_candidate <- sum(!(sub$all_candidates_present_once %in% TRUE))
    row$failure_rate <- row$n_no_selection / row$n_attempted
    interval <- jss_misspec_wilson_interval(row$n_no_selection, row$n_attempted)
    row$failure_mcse <- sqrt(row$failure_rate * (1 - row$failure_rate) / row$n_attempted)
    row$failure_ci_lower <- interval[["lower"]]
    row$failure_ci_upper <- interval[["upper"]]
    row$failure_reasons <- paste(sort(unique(sub$failure_reasons[nzchar(sub$failure_reasons)])), collapse = " | ")
    row
  })
  jss_misspec_row_bind(out)
}

jss_misspec_selection <- function(results, close_aic_threshold = 2) {
  attempts <- jss_misspec_selection_attempts(results, close_aic_threshold = close_aic_threshold)
  if (!nrow(attempts)) return(data.frame())
  group_cols <- c("generating_copula", "tau_label", "target_tau", "n_subject")
  split_key <- interaction(attempts[group_cols], drop = TRUE, lex.order = TRUE)
  out <- lapply(split(seq_len(nrow(attempts)), split_key), function(idx) {
    row <- attempts[idx[1L], group_cols, drop = FALSE]
    sub <- attempts[idx, , drop = FALSE]
    row$n_reps <- nrow(sub)
    row$n_attempted <- nrow(sub)
    row$n_selection_eligible <- sum(sub$selection_eligible %in% TRUE)
    for (criterion in c("aic", "bic", "loglik")) {
      selected <- sub[[paste0("best_", criterion, "_copula")]]
      correct <- !is.na(selected) & selected == sub$generating_copula
      count <- sum(correct)
      rate <- count / nrow(sub)
      interval <- jss_misspec_wilson_interval(count, nrow(sub))
      row[[paste0(criterion, "_n_correct")]] <- count
      row[[paste0(criterion, "_correct_selection_rate")]] <- rate
      row[[paste0(criterion, "_correct_selection_mcse")]] <- sqrt(rate * (1 - rate) / nrow(sub))
      row[[paste0(criterion, "_correct_selection_ci_lower")]] <- interval[["lower"]]
      row[[paste0(criterion, "_correct_selection_ci_upper")]] <- interval[["upper"]]
    }
    row$n_close_aic <- sum(sub$close_aic %in% TRUE)
    row$close_aic_rate <- row$n_close_aic / nrow(sub)
    close_interval <- jss_misspec_wilson_interval(row$n_close_aic, nrow(sub))
    row$close_aic_mcse <- sqrt(row$close_aic_rate * (1 - row$close_aic_rate) / nrow(sub))
    row$close_aic_ci_lower <- close_interval[["lower"]]
    row$close_aic_ci_upper <- close_interval[["upper"]]
    median_interval <- jss_misspec_bootstrap_median_interval(
      sub$aic_gap_to_runner_up,
      seed = 20260901L + as.integer(row$n_subject) + 100L * match(row$generating_copula, jss_misspec_copulas()) + 10L * match(row$tau_label, jss_misspec_tau_levels()$tau_label)
    )
    row$median_aic_gap_to_runner_up <- median_interval[["median"]]
    row$median_aic_gap_bootstrap_se <- median_interval[["bootstrap_se"]]
    row$median_aic_gap_ci_lower <- median_interval[["lower"]]
    row$median_aic_gap_ci_upper <- median_interval[["upper"]]
    row$median_aic_gap_bootstrap_reps <- 2000L
    row$close_aic_threshold <- close_aic_threshold
    row
  })
  jss_misspec_row_bind(out)
}

jss_misspec_best_copula <- function(results, metric, decreasing = FALSE) {
  value <- results[[metric]]
  ok <- is.finite(value)
  if (!any(ok)) {
    return(NA_character_)
  }
  usable <- results[ok, , drop = FALSE]
  value <- value[ok]
  if (isTRUE(decreasing)) {
    usable$fitted_copula[which.max(value)]
  } else {
    usable$fitted_copula[which.min(value)]
  }
}

jss_misspec_review_gate <- function(results, grid, paths) {
  expected <- nrow(grid)
  actual <- nrow(results)
  combo_cols <- c("generating_copula", "fitted_copula", "tau_label", "n_subject")
  expected_combos <- unique(grid[combo_cols])
  actual_combos <- if (nrow(results)) unique(results[combo_cols]) else results[combo_cols]

  high_clayton <- results[
    results$tau_label == "high" &
      (results$generating_copula == "C" | results$fitted_copula == "C"),
    ,
    drop = FALSE
  ]
  correct <- results[results$correctly_specified %in% TRUE, , drop = FALSE]
  wrong <- results[!(results$correctly_specified %in% TRUE), , drop = FALSE]
  correct_success <- if (nrow(correct)) mean(correct$success %in% TRUE) else NA_real_
  wrong_success <- if (nrow(wrong)) mean(wrong$success %in% TRUE) else NA_real_

  data.frame(
    check = c(
      "row_count",
      "full_combination_coverage",
      "high_tau_clayton_convergence_recorded",
      "finite_information_criteria",
      "correct_copula_not_systematically_worse",
      "selection_attempt_totals_reconcile",
      "summary_artifacts_rendered"
    ),
    status = c(
      if (actual == expected) "pass" else "review",
      if (nrow(actual_combos) == nrow(expected_combos)) "pass" else "review",
      if (!nrow(high_clayton) || all(!is.na(high_clayton$converged))) "pass" else "review",
      if (!nrow(results) || mean(is.finite(results$aic_joint) & is.finite(results$bic_joint), na.rm = TRUE) > 0.5) "pass" else "review",
      if (is.na(correct_success) || is.na(wrong_success) || correct_success + 0.05 >= wrong_success) "pass" else "review",
      {
        attempts <- jss_misspec_selection_attempts(results)
        expected_attempts <- length(unique(grid$fit_id)) / length(jss_misspec_copulas())
        if (nrow(attempts) == expected_attempts && sum(attempts$n_candidate_rows) == actual) "pass" else "review"
      },
      if (all(file.exists(c(
        paths$summary, paths$selection, paths$selection_attempts,
        paths$selection_confusion, paths$selection_failures, paths$paired_effects,
        paths$warning_audit, paths$execution_manifest,
        paths$heatmap, paths$convergence
      )))) "pass" else "review"
    ),
    detail = c(
      paste(actual, "of", expected, "fit rows available"),
      paste(nrow(actual_combos), "of", nrow(expected_combos), "scenario combinations represented"),
      paste(nrow(high_clayton), "registered high-tau rows involving Clayton have convergence status"),
      paste("finite AIC/BIC rate:", round(mean(is.finite(results$aic_joint) & is.finite(results$bic_joint), na.rm = TRUE), 3)),
      paste("correct success:", round(correct_success, 3), "wrong success:", round(wrong_success, 3)),
      paste("fit rows:", actual, "selection attempts:", nrow(jss_misspec_selection_attempts(results))),
      "summary, paired effects, warning audit, execution manifest, selection, heatmap, and convergence artifacts checked"
    ),
    stringsAsFactors = FALSE
  )
}

jss_misspec_write_heatmap <- function(summary, path) {
  if (!nrow(summary)) {
    return(jss_misspec_write_empty_plot(path, "Gamma copula mis-specification heatmap"))
  }
  p <- ggplot2::ggplot(
    summary,
    ggplot2::aes(x = fitted_copula, y = generating_copula, fill = mean_delta_joint_loglik_vs_correct)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::facet_grid(tau_label ~ n_subject, labeller = ggplot2::label_both) +
    ggplot2::scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0, na.value = "grey85") +
    ggplot2::labs(x = "Fitted copula", y = "Generating copula", fill = "Delta joint logLik") +
    ggplot2::theme_minimal(base_size = 10)
  ggplot2::ggsave(path, p, width = 10, height = 6, dpi = 320, bg = "white")
  path
}

jss_misspec_write_metric_heatmap <- function(summary, path, metric, title, fill_label) {
  if (!nrow(summary)) {
    return(jss_misspec_write_empty_plot(path, title))
  }
  p <- ggplot2::ggplot(
    summary,
    ggplot2::aes(x = fitted_copula, y = generating_copula, fill = .data[[metric]])
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::facet_grid(tau_label ~ n_subject, labeller = ggplot2::label_both) +
    ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08519c", na.value = "grey85") +
    ggplot2::labs(title = title, x = "Fitted copula", y = "Generating copula", fill = fill_label) +
    ggplot2::theme_minimal(base_size = 10)
  ggplot2::ggsave(path, p, width = 10, height = 6, dpi = 320, bg = "white")
  path
}

jss_misspec_summary_panel <- function(plot_data, metric, title, fill_label, scale) {
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = fitted_copula, y = generating_copula, fill = .data[[metric]])
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    scale +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, x = "Fitted copula", y = "Generating copula", fill = fill_label) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "bottom",
      legend.key.width = grid::unit(1.25, "cm"),
      panel.grid = ggplot2::element_blank()
    )
}

jss_misspec_figure_metadata <- function(summary, n_subject = 500L, tau_label = "high", reps = 100L) {
  selected <- summary[summary$n_subject == n_subject & summary$tau_label == tau_label, , drop = FALSE]
  tau <- unique(as.numeric(selected$target_tau))
  if (length(tau) != 1L || !is.finite(tau)) {
    tau <- unname(stats::setNames(jss_misspec_tau_levels()$target_tau, jss_misspec_tau_levels()$tau_label)[tau_label])
  }
  list(
    n_subject = as.integer(n_subject), target_tau = as.numeric(tau), reps = as.integer(reps),
    visible_label = sprintf("Gamma margins: n = %d subjects; Kendall's tau = %.2f; R = %d datasets per generating copula",
      as.integer(n_subject), as.numeric(tau), as.integer(reps)),
    tau_panel_title = "Mean absolute Kendall tau error"
  )
}

jss_misspec_write_paper_summary_heatmap <- function(summary, path, n_subject = 500L, tau_label = "high", reps = 100L) {
  metadata <- jss_misspec_figure_metadata(summary, n_subject, tau_label, reps)
  title <- metadata$visible_label
  if (!nrow(summary)) {
    return(jss_misspec_write_empty_plot(path, title))
  }

  plot_data <- summary[
    summary$n_subject == n_subject & summary$tau_label == tau_label,
    ,
    drop = FALSE
  ]
  if (!nrow(plot_data)) {
    return(jss_misspec_write_empty_plot(path, title))
  }

  copulas <- jss_misspec_copulas()
  plot_data$generating_copula <- factor(plot_data$generating_copula, levels = rev(copulas))
  plot_data$fitted_copula <- factor(plot_data$fitted_copula, levels = copulas)

  likelihood_panel <- jss_misspec_summary_panel(
    plot_data,
    metric = "mean_delta_joint_loglik_vs_correct",
    title = "Delta joint log-likelihood",
    fill_label = "Delta logLik",
    scale = ggplot2::scale_fill_gradient2(
      low = "#b2182b",
      mid = "white",
      high = "#2166ac",
      midpoint = 0,
      na.value = "grey85"
    )
  )
  margin_panel <- jss_misspec_summary_panel(
    plot_data,
    metric = "mean_margin_param_rmse",
    title = "Margin parameter RMSE",
    fill_label = "RMSE",
    scale = ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08519c", na.value = "grey85")
  )
  tau_panel <- jss_misspec_summary_panel(
    plot_data,
    metric = "mean_tau_abs_error",
    title = metadata$tau_panel_title,
    fill_label = "Abs. error",
    scale = ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08519c", na.value = "grey85")
  )

  grDevices::png(path, width = 11, height = 4.6, units = "in", res = 320, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    2, 3, heights = grid::unit(c(0.35, 4.25), "in")
  )))
  grid::grid.text(title, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:3),
    gp = grid::gpar(fontface = "bold", fontsize = 11))
  print(likelihood_panel, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(margin_panel, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  print(tau_panel, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 3))
  grid::popViewport()
  path
}

jss_misspec_write_convergence_plot <- function(summary, path) {
  if (!nrow(summary)) {
    return(jss_misspec_write_empty_plot(path, "Gamma copula mis-specification convergence"))
  }
  p <- ggplot2::ggplot(
    summary,
    ggplot2::aes(x = fitted_copula, y = convergence_rate, fill = fitted_copula)
  ) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::facet_grid(generating_copula + tau_label ~ n_subject, labeller = ggplot2::label_both) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Fitted copula", y = "Convergence rate") +
    ggplot2::theme_minimal(base_size = 9)
  ggplot2::ggsave(path, p, width = 11, height = 8, dpi = 320, bg = "white")
  path
}

jss_misspec_write_empty_plot <- function(path, title) {
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_blank() +
    ggplot2::labs(title = title, subtitle = "No result rows available yet.") +
    ggplot2::theme_void(base_size = 11)
  ggplot2::ggsave(path, p, width = 7, height = 4, dpi = 220, bg = "white")
  path
}

jss_run_07_gamma_copula_misspecification <- function(settings, stage = jss_misspec_stage(settings$profile)) {
  stage_explicit <- !missing(stage)
  if (nzchar(Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MISSPEC_REPS", unset = "")) &&
      !isTRUE(stage_explicit) && !identical(stage, "full")) {
    stop("Module 07 replicate overrides require an explicit smoke or pilot nonpublication stage.", call. = FALSE)
  }
  paths <- jss_misspec_paths(settings)
  config <- jss_misspec_config(settings, stage = stage)
  jss_misspec_load_checkout(config$repo_root, config$package_source_sha256)
  execution_start <- jss_misspec_execution_start(config, settings$workers)
  checkpoint_dir <- jss_misspec_checkpoint_dir(paths, config$stage)
  grid <- jss_misspec_grid(config)

  jss_misspec_write_csv_atomic(grid, paths$grid)
  jss_misspec_run_checkpoints(grid, config, checkpoint_dir, workers = settings$workers)

  results <- jss_misspec_read_checkpoints(checkpoint_dir, grid = grid, config = config)
  results <- jss_misspec_add_deltas(results)
  summary <- jss_misspec_summary(results)
  selection_attempts <- jss_misspec_selection_attempts(results)
  selection <- jss_misspec_selection(results)
  selection_confusion <- jss_misspec_selection_confusion(selection_attempts, criterion = "aic")
  selection_failures <- jss_misspec_selection_failures(selection_attempts)
  paired_effects <- jss_misspec_paired_effects(results)
  warning_audit <- jss_misspec_warning_audit(results)

  jss_misspec_write_csv_atomic(results, paths$results)
  installed_results <- utils::read.csv(paths$results, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(names(installed_results), names(results)) || nrow(installed_results) != nrow(results) ||
      !identical(jss_misspec_sha256_file(paths$results), local({
        tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
        utils::write.csv(results, tmp, row.names = FALSE); jss_misspec_sha256_file(tmp)
      }))) stop("Fatal Module 07 installed-result verification failed.", call. = FALSE)
  jss_misspec_write_csv_atomic(summary, paths$summary)
  jss_misspec_write_csv_atomic(selection, paths$selection)
  jss_misspec_write_csv_atomic(selection_attempts, paths$selection_attempts)
  jss_misspec_write_csv_atomic(selection_confusion, paths$selection_confusion)
  jss_misspec_write_csv_atomic(selection_failures, paths$selection_failures)
  jss_misspec_write_csv_atomic(paired_effects, paths$paired_effects)
  jss_misspec_write_csv_atomic(warning_audit, paths$warning_audit)
  jss_misspec_write_heatmap(summary, paths$heatmap)
  jss_misspec_write_metric_heatmap(
    summary,
    paths$margin_rmse_heatmap,
    metric = "mean_margin_param_rmse",
    title = "Gamma margin parameter RMSE by fitted copula",
    fill_label = "Margin RMSE"
  )
  jss_misspec_write_metric_heatmap(
    summary,
    paths$tau_error_heatmap,
    metric = "mean_tau_abs_error",
    title = "Kendall tau absolute error by fitted copula",
    fill_label = "Tau abs. error"
  )
  jss_misspec_write_paper_summary_heatmap(summary, paths$paper_summary_heatmap)
  jss_misspec_write_convergence_plot(summary, paths$convergence)
  execution_manifest <- jss_misspec_execution_manifest(
    execution_start, config, grid, installed_results,
    c(results = paths$results, paired_effects = paths$paired_effects,
      warning_audit = paths$warning_audit, selection_attempts = paths$selection_attempts,
      selection = paths$selection)
  )
  jss_misspec_write_csv_atomic(execution_manifest, paths$execution_manifest)
  review <- jss_misspec_review_gate(results, grid, paths)
  binding_review <- jss_misspec_binding_review_gate(installed_results, grid, paths, context = paste0("full:", config$stage), config = config)
  review <- rbind(review, binding_review[names(review)])
  jss_misspec_write_csv_atomic(review, paths$review)
  if (identical(config$stage, "full") && any(review$status != "pass")) {
    stop("Fatal Module 07 full review gate failed: ", paste(review$check[review$status != "pass"], collapse = ", "), call. = FALSE)
  }
  if (identical(config$stage, "full")) jss_misspec_validate_evidence_bundle(paths$results, config)

  list(
    module_id = "07-gamma-copula-misspecification",
    title = "Gamma copula mis-specification simulation",
    status = "current",
    stage = config$stage,
    data = c(paths$grid, paths$results, paths$selection_attempts, paths$execution_manifest),
    tables = c(paths$summary, paths$selection, paths$selection_confusion, paths$selection_failures,
      paths$paired_effects, paths$warning_audit, paths$review),
    figures = c(
      paths$heatmap,
      paths$margin_rmse_heatmap,
      paths$tau_error_heatmap,
      paths$paper_summary_heatmap,
      paths$convergence
    ),
    notes = paste("Checkpointed", config$stage, "run with", nrow(grid), "planned fits.")
  )
}
