jss_missing_mechanism_registry <- function() {
  data.frame(
    missing_mechanism = c(
      "monotone_dropout",
      "time_dependent_intermittent_mar",
      "mar"
    ),
    display_name = c(
      "subject-level monotone dropout MAR",
      "time-dependent intermittent MAR",
      "covariate-dependent intermittent MAR"
    ),
    pattern = c("monotone_dropout", "intermittent", "intermittent"),
    analysis_role = c("headline", "sensitivity", "sensitivity"),
    stringsAsFactors = FALSE
  )
}

jss_missing_validate_mechanisms <- function(x) {
  x <- trimws(as.character(x))
  if (any(x == "time_mar")) {
    stop(
      "'time_mar' was renamed to 'time_dependent_intermittent_mar'; update MISSING_MECHANISMS.",
      call. = FALSE
    )
  }
  allowed <- jss_missing_mechanism_registry()$missing_mechanism
  if (!length(x) || any(!nzchar(x)) || any(!x %in% allowed)) {
    stop("Missingness mechanisms must be drawn from: ", paste(allowed, collapse = ", "), call. = FALSE)
  }
  x
}

jss_missing_finite_scalar_or_na <- function(x, integer = FALSE) {
  if (!is.logical(integer) || length(integer) != 1L || is.na(integer)) {
    stop("integer must be one non-missing logical value.", call. = FALSE)
  }
  missing_value <- if (isTRUE(integer)) NA_integer_ else NA_real_
  if (length(x) != 1L || !(is.integer(x) || is.double(x)) || !is.finite(x[[1L]])) {
    return(missing_value)
  }
  value <- as.numeric(x[[1L]])
  if (!isTRUE(integer)) return(value)
  if (value != floor(value) || abs(value) > .Machine$integer.max) {
    return(NA_integer_)
  }
  as.integer(value)
}

jss_missing_calibrate_intercept <- function(mean_rate, target_rate) {
  if (target_rate <= 0) return(-Inf)
  stats::uniroot(
    function(a) mean_rate(a) - target_rate,
    interval = c(-35, 35)
  )$root
}

jss_missing_subject_scale <- function(x) {
  out <- as.numeric(scale(x))
  out[!is.finite(out)] <- 0
  out
}

jss_missing_monotone_probabilities <- function(subjects, d, target_rate) {
  if (d < 2L) stop("Monotone dropout requires at least two scheduled visits.", call. = FALSE)
  maximum_rate <- (d - 1) / d
  if (target_rate >= maximum_rate) {
    stop(
      "Target monotone-dropout rate must be below ", signif(maximum_rate, 4),
      " when there are ", d, " visits.",
      call. = FALSE
    )
  }
  transition <- seq_len(d - 1L)
  time_effect <- as.numeric(scale(transition))
  if (any(!is.finite(time_effect))) time_effect[] <- 0
  subject_effect <- 0.70 * jss_missing_subject_scale(subjects$x1) +
    0.50 * subjects$x2 - 0.40 * jss_missing_subject_scale(subjects$s1)
  linear <- outer(subject_effect, rep(1, d - 1L)) +
    outer(rep(1, nrow(subjects)), 1.10 * time_effect)
  expected_rate <- function(intercept) {
    hazard <- stats::plogis(intercept + linear)
    survival <- matrix(1, nrow = nrow(subjects), ncol = d)
    for (j in 2:d) survival[, j] <- survival[, j - 1L] * (1 - hazard[, j - 1L])
    mean(1 - survival)
  }
  intercept <- jss_missing_calibrate_intercept(expected_rate, target_rate)
  hazard <- stats::plogis(intercept + linear)
  cumulative_missing <- matrix(0, nrow = nrow(subjects), ncol = d)
  survival <- rep(1, nrow(subjects))
  for (j in 2:d) {
    survival <- survival * (1 - hazard[, j - 1L])
    cumulative_missing[, j] <- 1 - survival
  }
  list(hazard = hazard, cumulative_missing = cumulative_missing)
}

jss_missing_summarize <- function(dat, mechanism, target_rate) {
  dat <- dat[order(dat$id, dat$time), , drop = FALSE]
  patterns <- split(is.finite(dat$response), dat$id)
  has_interior_gap <- vapply(patterns, function(observed) {
    later_observed <- rev(cummax(as.integer(rev(observed))) > 0)
    any(!observed & later_observed)
  }, logical(1))
  is_complete <- vapply(patterns, all, logical(1))
  is_monotone <- vapply(patterns, function(observed) {
    missing <- which(!observed)
    !length(missing) || all(!observed[min(missing):length(observed)])
  }, logical(1))
  pair_groups <- split(seq_len(nrow(dat)), dat$id)
  total_pairs <- sum(pmax(lengths(pair_groups) - 1L, 0L))
  complete_pairs <- sum(vapply(pair_groups, function(idx) {
    if (length(idx) < 2L) return(0L)
    sum(is.finite(dat$response[idx[-length(idx)]]) & is.finite(dat$response[idx[-1L]]))
  }, integer(1)))
  registry <- jss_missing_mechanism_registry()
  meta <- registry[registry$missing_mechanism == mechanism, , drop = FALSE]
  data.frame(
    missing_mechanism = mechanism,
    missingness_label = meta$display_name,
    missingness_pattern = meta$pattern,
    analysis_role = meta$analysis_role,
    target_missing_rate = target_rate,
    observed_missing_rate = mean(!is.finite(dat$response)),
    n_rows = nrow(dat),
    n_observed_rows = sum(is.finite(dat$response)),
    n_subjects = length(patterns),
    n_complete_subjects = sum(is_complete),
    n_dropout_subjects = if (identical(meta$pattern[[1L]], "monotone_dropout")) {
      sum(!is_complete & is_monotone)
    } else {
      0L
    },
    n_monotone_incomplete_subjects = sum(!is_complete & is_monotone),
    n_subjects_with_interior_gaps = sum(has_interior_gap),
    no_observations_after_dropout = all(is_monotone),
    total_adjacent_pairs = total_pairs,
    complete_adjacent_pairs = complete_pairs,
    complete_adjacent_pair_rate = if (total_pairs) complete_pairs / total_pairs else NA_real_,
    stringsAsFactors = FALSE
  )
}

jss_missing_reconstruct_pattern_summary <- function(pattern, target_rate) {
  required <- c("scenario", "n", "d", "rep", "missing_mechanism", "id", "time",
    "response_observed")
  if (!is.data.frame(pattern) || !nrow(pattern) || !all(required %in% names(pattern)) ||
      length(unique(pattern$scenario)) != 1L || length(unique(pattern$rep)) != 1L ||
      length(unique(pattern$missing_mechanism)) != 1L) {
    stop("Missingness pattern task is not uniquely keyed.", call. = FALSE)
  }
  dat <- pattern[c("id", "time")]
  dat$response <- ifelse(as.logical(pattern$response_observed), 0, NA_real_)
  out <- jss_missing_summarize(dat, unique(pattern$missing_mechanism), target_rate)
  out$scenario <- unique(pattern$scenario)
  out$n <- as.integer(unique(pattern$n))
  out$d <- as.integer(unique(pattern$d))
  out$rep <- as.integer(unique(pattern$rep))
  out[c("scenario", "n", "d", "rep", setdiff(names(out), c("scenario", "n", "d", "rep")))]
}

jss_missing_apply <- function(dat, rate, mechanism, seed) {
  mechanism <- jss_missing_validate_mechanisms(mechanism)
  if (length(mechanism) != 1L) stop("Apply exactly one missingness mechanism at a time.", call. = FALSE)
  dat <- dat[order(dat$id, dat$time), , drop = FALSE]
  if (rate <= 0) {
    dat$missing_probability <- 0
    dat$response_observed <- TRUE
    attr(dat, "missingness_summary") <- jss_missing_summarize(dat, mechanism, rate)
    return(dat)
  }
  set.seed(seed)
  if (identical(mechanism, "monotone_dropout")) {
    first <- !duplicated(dat$id)
    subjects <- dat[first, c("id", "x1", "x2", "s1"), drop = FALSE]
    times <- sort(unique(dat$time))
    probs <- jss_missing_monotone_probabilities(subjects, length(times), rate)
    observed <- matrix(TRUE, nrow = nrow(subjects), ncol = length(times))
    for (i in seq_len(nrow(subjects))) {
      dropped <- FALSE
      for (j in seq_len(length(times) - 1L)) {
        if (!dropped && stats::runif(1) < probs$hazard[i, j]) dropped <- TRUE
        if (dropped) observed[i, j + 1L] <- FALSE
      }
    }
    row_i <- match(dat$id, subjects$id)
    col_i <- match(dat$time, times)
    dat$missing_probability <- probs$cumulative_missing[cbind(row_i, col_i)]
    dat$response_observed <- observed[cbind(row_i, col_i)]
  } else {
    linear <- 0.70 * jss_missing_subject_scale(dat$x1) + 0.50 * dat$x2 -
      0.40 * jss_missing_subject_scale(dat$s1)
    if (identical(mechanism, "time_dependent_intermittent_mar")) {
      linear <- linear + 1.35 * jss_missing_subject_scale(dat$t)
    }
    intercept <- jss_missing_calibrate_intercept(
      function(a) mean(stats::plogis(a + linear)),
      rate
    )
    dat$missing_probability <- stats::plogis(intercept + linear)
    dat$response_observed <- stats::runif(nrow(dat)) >= dat$missing_probability
  }
  dat$response[!dat$response_observed] <- NA_real_
  attr(dat, "missingness_summary") <- jss_missing_summarize(dat, mechanism, rate)
  dat
}

jss_missing_short_hash <- function(x) {
  text <- paste(capture.output(dput(x)), collapse = "")
  ints <- utf8ToInt(enc2utf8(text))
  checksum <- if (length(ints)) sum((seq_along(ints) %% 1009L + 1L) * ints) %% 2147483629 else 0
  paste0(nchar(text), "-", format(checksum, scientific = FALSE, trim = TRUE))
}

jss_missing_sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required for missingness checkpoint integrity.", call. = FALSE)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

jss_missing_content_sha256 <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required for missingness checkpoint integrity.", call. = FALSE)
  digest::digest(x, algo = "sha256", serialize = TRUE)
}

jss_missing_producer_sha256 <- function(producer_sources) {
  sources <- sort(normalizePath(producer_sources, winslash = "/", mustWork = TRUE))
  digest::digest(
    paste(basename(sources), vapply(sources, jss_missing_sha256_file, character(1)), collapse = "\n"),
    algo = "sha256", serialize = FALSE
  )
}

jss_missing_portable_task_sha256 <- function(runs, fixed, smooth, missingness, missingness_pattern = NULL) {
  encode_frame <- function(x) {
    if (is.null(x)) return("<NULL>")
    x <- x[, setdiff(names(x), "public_payload_sha256"), drop = FALSE]
    encode <- function(value) {
      if (is.logical(value)) return(ifelse(is.na(value), "NA", ifelse(value, "TRUE", "FALSE")))
      if (is.numeric(value)) return(ifelse(is.na(value), "NA", sprintf("%.15g", value)))
      value <- as.character(value); value[is.na(value) | !nzchar(value)] <- "NA"; enc2utf8(value)
    }
    x <- x[, sort(names(x)), drop = FALSE]
    rows <- if (nrow(x)) apply(as.data.frame(lapply(x, encode), stringsAsFactors = FALSE), 1L,
      paste, collapse = "\r") else character()
    paste(sort(rows), collapse = "\n")
  }
  text <- paste(vapply(list(runs = runs, fixed = fixed, smooth = smooth,
    missingness = missingness, missingness_pattern = missingness_pattern),
    encode_frame, character(1)), collapse = "\f")
  digest::digest(text, algo = "sha256", serialize = FALSE)
}

jss_missing_checkout_identity <- function(project_root) {
  root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  files <- sort(normalizePath(c(
    list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE),
    file.path(root, c("DESCRIPTION", "NAMESPACE"))
  ), winslash = "/", mustWork = TRUE))
  rel <- substring(files, nchar(root) + 2L)
  manifest <- paste(rel, vapply(files, jss_missing_sha256_file, character(1)), sep = "\t", collapse = "\n")
  desc <- read.dcf(file.path(root, "DESCRIPTION"))
  list(package = "gamlss.longitudinal", version = unname(desc[1L, "Version"]),
    checkout_path = root, source_sha256 = digest::digest(manifest, "sha256", serialize = FALSE),
    fingerprint_scope = "sorted R/*.R + DESCRIPTION + NAMESPACE", source_file_count = length(files),
    load_strategy = "pkgload_checkout")
}

jss_missing_verify_checkout <- function(project_root, expected, load_checkout = TRUE) {
  if (isTRUE(load_checkout)) {
    if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required for verified missingness execution.", call. = FALSE)
    pkgload::load_all(project_root, quiet = TRUE, export_all = TRUE, helpers = FALSE)
  }
  actual <- jss_missing_checkout_identity(project_root)
  namespace_path <- normalizePath(getNamespaceInfo(asNamespace(expected$package), "path"), winslash = "/", mustWork = TRUE)
  valid <- identical(actual$source_sha256, expected$source_sha256) &&
    identical(actual$version, expected$version) && identical(namespace_path, expected$checkout_path)
  if (!valid) stop("Missingness worker package/source identity mismatch.", call. = FALSE)
  c(expected, list(verified = TRUE, verified_package_path = namespace_path,
    verified_source_sha256 = actual$source_sha256, verified_version = actual$version))
}

jss_missing_reverify_sources <- function(project_root, expected_package_identity,
                                         producer_sources, expected_producer_sha256) {
  actual_package <- jss_missing_checkout_identity(project_root)
  actual_producer <- jss_missing_producer_sha256(producer_sources)
  hex <- "^[0-9a-f]{64}$"
  if (!grepl(hex, expected_package_identity$source_sha256) ||
      !grepl(hex, expected_producer_sha256) ||
      !identical(actual_package$source_sha256, expected_package_identity$source_sha256) ||
      !identical(actual_package$version, expected_package_identity$version) ||
      !identical(actual_producer, expected_producer_sha256)) {
    stop("Missingness checkout or producer source changed after attestation.", call. = FALSE)
  }
  invisible(list(package = actual_package, producer_sha256 = actual_producer))
}

jss_missing_runtime_identity <- function() {
  ext <- extSoftVersion()
  portable_library <- function(value) {
    if (is.na(value) || !nzchar(value)) return("R-default-or-unreported")
    if (grepl("^([A-Za-z]:[/\\\\]|/)", value)) basename(value) else value
  }
  blas <- if ("BLAS" %in% names(ext)) unname(ext[["BLAS"]]) else NA_character_
  lapack <- if ("LAPACK" %in% names(ext)) unname(ext[["LAPACK"]]) else NA_character_
  list(timestamp_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), pid = Sys.getpid(),
    host = unname(Sys.info()[["nodename"]]),
    os = paste(unname(Sys.info()[c("sysname", "release", "machine")]), collapse = " "),
    platform = R.version$platform,
    r_version = R.version.string, rng_kind = paste(RNGkind(), collapse = "/"),
    blas = portable_library(blas), lapack = portable_library(lapack),
    rlibs_user_sha256 = digest::digest(gsub("\\\\", "/", Sys.getenv("R_LIBS_USER", unset = "")),
      "sha256", serialize = FALSE),
    libpaths_sha256 = digest::digest(paste(gsub("\\\\", "/", .libPaths()), collapse = ";"),
      "sha256", serialize = FALSE))
}

jss_missing_estimand_registry <- function() {
  data.frame(
    component = c("fixed_intercept", "fixed_nonintercept", "smooth_curve"),
    estimand = c(
      "population coefficient plus smooth mean under Uniform(0,1) s1",
      "full-data population generating coefficient",
      "population smooth centered under the fixed Uniform(0,1) s1 reference grid"
    ),
    observed_data_projection = FALSE,
    reference_distribution = c("s1 ~ Uniform(0,1)", "generating model", "fixed 10001-point s1 grid over [0,1]"),
    stringsAsFactors = FALSE
  )
}

jss_missing_population_smooth_mean <- function(fun) {
  mean(fun(seq(0, 1, length.out = 10001L)))
}

jss_missing_seed_offset <- function(mechanism) {
  offsets <- c(
    mar = 150000L,
    monotone_dropout = 250000L,
    time_dependent_intermittent_mar = 350000L
  )
  value <- unname(offsets[[as.character(mechanism)[[1L]]]])
  if (is.null(value)) stop("No registered seed offset for missingness mechanism.", call. = FALSE)
  as.integer(value)
}

jss_missing_checkpoint_spec <- function(task, configuration) {
  list(
    schema_version = 4L,
    scenario_id = as.integer(task$scenario_id[[1L]]),
    scenario = sprintf(
      "n%d_d%d_%s_miss%02d",
      as.integer(task$n[[1L]]), as.integer(task$d[[1L]]),
      as.character(task$missing_mechanism[[1L]]),
      round(100 * as.numeric(task$missing_rate[[1L]]))
    ),
    n = as.integer(task$n[[1L]]),
    d = as.integer(task$d[[1L]]),
    replicate = as.integer(task$rep[[1L]]),
    missing_mechanism = as.character(task$missing_mechanism[[1L]]),
    missing_rate = as.numeric(task$missing_rate[[1L]]),
    simulation_seed = 100000L + as.integer(task$rep[[1L]]),
    missingness_seed = 100000L + as.integer(task$rep[[1L]]) +
      jss_missing_seed_offset(task$missing_mechanism[[1L]]) +
      as.integer(round(1000 * as.numeric(task$missing_rate[[1L]]))),
    configuration = configuration,
    rng_kind = RNGkind()
  )
}

jss_missing_checkpoint_content <- function(result) {
  result[c("runs", "fixed", "smooth", "joint", "predictive", "missingness", "missingness_pattern")]
}

jss_missing_checkpoint_archive_record <- function(result) {
  if (!is.list(result) || !is.list(result$checkpoint_spec) ||
      !is.character(result$checkpoint_configuration_key) ||
      length(result$checkpoint_configuration_key) != 1L) {
    stop("Missingness checkpoint archive payload is malformed.", call. = FALSE)
  }
  spec <- result$checkpoint_spec
  checkpoint_name <- sprintf("%s_scenario%d_%s_miss%02d_n%d_d%d_rep%d.rds",
    result$checkpoint_configuration_key, spec$scenario_id, spec$missing_mechanism,
    round(100 * spec$missing_rate), spec$n, spec$d, spec$replicate)
  data.frame(scenario = spec$scenario, rep = spec$replicate,
    checkpoint = gsub("\\\\", "/", file.path("rep_results", checkpoint_name)),
    checkpoint_content_sha256 = jss_missing_content_sha256(jss_missing_checkpoint_content(result)),
    public_payload_sha256 = jss_missing_portable_task_sha256(result$runs, result$fixed,
      result$smooth, result$missingness, result$missingness_pattern), stringsAsFactors = FALSE)
}

jss_missing_eligible <- function(runs) {
  if (!is.data.frame(runs) || !all(c("success", "converged") %in% names(runs))) {
    return(rep(FALSE, if (is.data.frame(runs)) nrow(runs) else 0L))
  }
  runs$success %in% TRUE & runs$converged %in% TRUE
}

jss_missing_filter_payload <- function(payload, runs) {
  if (is.null(payload)) return(NULL)
  if (!is.data.frame(payload)) stop("Missingness checkpoint payload must be a data frame.", call. = FALSE)
  required <- c("scenario", "model", "rep")
  if (!all(required %in% names(payload)) || !all(required %in% names(runs))) {
    stop("Missingness checkpoint payload lacks scenario/model/rep keys.", call. = FALSE)
  }
  eligible_runs <- runs[jss_missing_eligible(runs), , drop = FALSE]
  make_key <- function(x) do.call(paste, c(lapply(x[required], as.character), sep = "\r"))
  payload_key <- make_key(payload)
  run_key <- make_key(eligible_runs)
  matched <- match(payload_key, run_key)
  keep <- !is.na(matched)
  metadata <- intersect(
    c("n", "d", "missing_mechanism", "target_missing_rate"),
    intersect(names(payload), names(eligible_runs))
  )
  for (field in metadata) {
    comparable <- keep
    if (identical(field, "target_missing_rate")) {
      comparable[keep] <- abs(
        as.numeric(payload[[field]][keep]) - as.numeric(eligible_runs[[field]][matched[keep]])
      ) < 1e-12
    } else {
      comparable[keep] <- as.character(payload[[field]][keep]) ==
        as.character(eligible_runs[[field]][matched[keep]])
    }
    comparable[is.na(comparable)] <- FALSE
    keep <- keep & comparable
  }
  payload[keep, , drop = FALSE]
}

jss_missing_checkpoint_valid <- function(result, task, configuration = list()) {
  expected <- jss_missing_checkpoint_spec(task, configuration)
  result_schema <- c("checkpoint_schema_version", "checkpoint_configuration_key", "checkpoint_spec",
    "fixed", "smooth", "joint", "predictive", "runs", "missingness", "missingness_pattern",
    "checkpoint_provenance", "checkpoint_content_sha256")
  if (!is.list(result) || !setequal(names(result), result_schema) ||
      !identical(result$checkpoint_schema_version, 4L) ||
      !is.character(result$checkpoint_configuration_key) || length(result$checkpoint_configuration_key) != 1L ||
      !identical(result$checkpoint_spec, expected) || !is.list(result$checkpoint_provenance) ||
      !is.character(result$checkpoint_content_sha256) || length(result$checkpoint_content_sha256) != 1L ||
      !identical(jss_missing_content_sha256(jss_missing_checkpoint_content(result)), result$checkpoint_content_sha256)) return(FALSE)
  runs <- result$runs
  payload_names <- c("fixed", "smooth", "joint", "predictive")
  required_runs <- c(
    "scenario", "n", "d", "rep", "missing_mechanism", "target_missing_rate",
    "simulation_seed", "missingness_seed", "model", "success", "converged", "retained",
    "stop_reason", "failure_type", "logLik", "df", "elapsed_sec", "error"
  )
  base_valid <- is.data.frame(runs) && nrow(runs) == 2L &&
    all(required_runs %in% names(runs)) &&
    all(vapply(runs[c("n", "d", "rep", "simulation_seed", "missingness_seed")],
      is.integer, logical(1))) &&
    all(vapply(runs[c("scenario", "missing_mechanism", "model", "stop_reason",
      "failure_type", "error")], is.character, logical(1))) &&
    all(vapply(runs[c("success", "converged", "retained")], is.logical, logical(1))) &&
    is.numeric(runs$logLik) && is.numeric(runs$df) && is.numeric(runs$elapsed_sec) &&
    !anyNA(runs[c("success", "converged", "retained")]) &&
    identical(runs$retained, jss_missing_eligible(runs)) &&
    all(runs$scenario == expected$scenario) &&
    all(runs$n == expected$n) && all(runs$d == expected$d) &&
    all(runs$rep == expected$replicate) &&
    all(runs$missing_mechanism == task$missing_mechanism) &&
    all(abs(runs$target_missing_rate - task$missing_rate) < 1e-12) &&
    all(runs$simulation_seed == expected$simulation_seed) &&
    all(runs$missingness_seed == expected$missingness_seed) &&
    identical(sort(as.character(runs$model)), c("gamlss.longitudinal", "gamlss2")) &&
    !any(!(runs$success %in% TRUE) & runs$converged %in% TRUE) &&
    all(runs$success %in% TRUE | (is.na(runs$logLik) & is.na(runs$df))) &&
    all(is.finite(runs$elapsed_sec) & runs$elapsed_sec >= 0 & runs$elapsed_sec <= 1e12)
  if (!isTRUE(base_valid) || !all(c(payload_names, "missingness_pattern") %in% names(result))) return(FALSE)
  provenance <- result$checkpoint_provenance
  configured_identity <- expected$configuration$package_identity
  runtime_fields <- c("timestamp_utc", "pid", "host", "os", "platform", "r_version", "rng_kind",
    "blas", "lapack", "rlibs_user_sha256", "libpaths_sha256")
  identity_fields <- if (is.null(configured_identity)) character() else
    c("package_identity_verified", "package_source_sha256", "producer_sha256_verified")
  if (!setequal(names(provenance), c(runtime_fields, identity_fields)) ||
      !is.integer(provenance$pid) || provenance$pid < 1L ||
      !all(vapply(provenance[setdiff(runtime_fields, "pid")], is.character, logical(1))) ||
      any(!nzchar(unlist(provenance[setdiff(runtime_fields, "pid")], use.names = FALSE))) ||
      !grepl("^[0-9a-f]{64}$", provenance$rlibs_user_sha256) ||
      !grepl("^[0-9a-f]{64}$", provenance$libpaths_sha256) ||
      !identical(provenance$r_version, R.version.string) ||
      !identical(provenance$rng_kind, "Mersenne-Twister/Inversion/Rejection") ||
      !identical(provenance$platform, R.version$platform)) return(FALSE)
  parsed_time <- as.POSIXct(strptime(provenance$timestamp_utc, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  if (is.na(parsed_time) || !identical(format(parsed_time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      provenance$timestamp_utc)) return(FALSE)
  observed <- jss_missing_runtime_identity()
  observed_fields <- c("host", "os", "platform", "r_version", "rng_kind", "blas", "lapack",
    "rlibs_user_sha256", "libpaths_sha256")
  if (any(vapply(observed_fields, function(field)
      !identical(provenance[[field]], observed[[field]]), logical(1)))) return(FALSE)
  if (!is.null(configured_identity) && (!isTRUE(provenance$package_identity_verified) ||
      !identical(provenance$package_source_sha256, configured_identity$source_sha256) ||
      !identical(provenance$producer_sha256_verified, expected$configuration$producer_sha256))) return(FALSE)
  miss <- result$missingness
  required_missing <- c("scenario", "n", "d", "rep", "missing_mechanism",
    "missingness_label", "missingness_pattern", "analysis_role", "target_missing_rate",
    "observed_missing_rate", "n_rows", "n_observed_rows", "n_subjects", "n_complete_subjects",
    "n_dropout_subjects", "n_monotone_incomplete_subjects", "n_subjects_with_interior_gaps",
    "no_observations_after_dropout", "total_adjacent_pairs", "complete_adjacent_pairs",
    "complete_adjacent_pair_rate")
  if (!is.data.frame(miss) || nrow(miss) != 1L || !all(required_missing %in% names(miss)) ||
      !is.character(miss$scenario) || !is.character(miss$missing_mechanism) ||
      !is.character(miss$missingness_label) || !is.character(miss$missingness_pattern) ||
      !is.character(miss$analysis_role) ||
      !is.integer(miss$n) || !is.integer(miss$d) || !is.integer(miss$rep) ||
      !all(vapply(miss[c("n_rows", "n_observed_rows", "n_subjects", "n_complete_subjects",
        "n_dropout_subjects", "n_monotone_incomplete_subjects", "n_subjects_with_interior_gaps",
        "total_adjacent_pairs", "complete_adjacent_pairs")], is.integer, logical(1))) ||
      !is.logical(miss$no_observations_after_dropout) || anyNA(miss$no_observations_after_dropout) ||
      !identical(as.character(miss$scenario), expected$scenario) ||
      miss$n[[1L]] != expected$n || miss$d[[1L]] != expected$d || miss$rep[[1L]] != expected$replicate ||
      !identical(as.character(miss$missing_mechanism), expected$missing_mechanism) ||
      abs(miss$target_missing_rate[[1L]] - expected$missing_rate) >= 1e-12 ||
      miss$observed_missing_rate[[1L]] < 0 || miss$observed_missing_rate[[1L]] > 1 ||
      miss$complete_adjacent_pair_rate[[1L]] < 0 || miss$complete_adjacent_pair_rate[[1L]] > 1 ||
      (identical(expected$missing_mechanism, "monotone_dropout") &&
        !isTRUE(miss$no_observations_after_dropout[[1L]]))) return(FALSE)
  pattern <- result$missingness_pattern
  required_pattern <- c("scenario", "n", "d", "rep", "missing_mechanism", "id", "time", "response_observed")
  if (!is.data.frame(pattern) || !nrow(pattern) || !all(required_pattern %in% names(pattern)) ||
      !is.character(pattern$scenario) || !is.character(pattern$missing_mechanism) ||
      !all(vapply(pattern[c("n", "d", "rep", "id")], is.integer, logical(1))) ||
      !is.numeric(pattern$time) || !is.logical(pattern$response_observed) || anyNA(pattern$response_observed) ||
      any(pattern$scenario != expected$scenario) || any(pattern$n != expected$n) || any(pattern$d != expected$d) ||
      any(pattern$rep != expected$replicate) || any(pattern$missing_mechanism != expected$missing_mechanism) ||
      any(!is.finite(pattern$time)) || nrow(pattern) != expected$n * expected$d ||
      anyDuplicated(pattern[c("id", "time")]) || length(unique(pattern$id)) != expected$n ||
      any(table(pattern$id) != expected$d) ||
      ("observed_missing_rate" %in% names(miss) &&
        abs(mean(!pattern$response_observed) - miss$observed_missing_rate[[1L]]) > 1e-12)) return(FALSE)
  if (identical(expected$missing_mechanism, "monotone_dropout")) {
    monotone_ok <- vapply(split(pattern[order(pattern$id, pattern$time), ], pattern$id), function(x) {
      missing <- which(!x$response_observed)
      !length(missing) || all(!x$response_observed[min(missing):nrow(x)])
    }, logical(1))
    if (!all(monotone_ok)) return(FALSE)
  }
  reconstructed_miss <- tryCatch(
    jss_missing_reconstruct_pattern_summary(pattern, expected$missing_rate),
    error = function(e) NULL)
  if (is.null(reconstructed_miss)) return(FALSE)
  audit_fields <- c("scenario", "n", "d", "rep", "missing_mechanism", "missingness_label",
    "missingness_pattern", "analysis_role", "target_missing_rate", "observed_missing_rate",
    "n_rows", "n_observed_rows", "n_subjects", "n_complete_subjects", "n_dropout_subjects",
    "n_monotone_incomplete_subjects", "n_subjects_with_interior_gaps",
    "no_observations_after_dropout", "total_adjacent_pairs", "complete_adjacent_pairs",
    "complete_adjacent_pair_rate")
  if (!all(vapply(audit_fields, function(field) isTRUE(all.equal(miss[[field]],
      reconstructed_miss[[field]], tolerance = 1e-12, check.attributes = FALSE)), logical(1)))) return(FALSE)
  meta <- jss_missing_mechanism_registry()
  meta <- meta[meta$missing_mechanism == expected$missing_mechanism, , drop = FALSE]
  if (nrow(meta) != 1L || miss$missingness_label[[1L]] != meta$display_name[[1L]] ||
      miss$missingness_pattern[[1L]] != meta$pattern[[1L]] ||
      miss$analysis_role[[1L]] != meta$analysis_role[[1L]] ||
      miss$n_rows[[1L]] != nrow(pattern) || miss$n_observed_rows[[1L]] != sum(pattern$response_observed) ||
      miss$n_subjects[[1L]] != length(unique(pattern$id))) return(FALSE)
  eligible <- jss_missing_eligible(runs)
  allowed_stops <- c("converged", "relative_deviance_tolerance", "max_iterations", "max_stall",
    "objective_deterioration", "invalid_likelihood", "numerical_failure", "time_limit",
    "outer_iteration_limit_or_invalid_loglik", "fit_error")
  expected_failure <- ifelse(eligible, "none",
    ifelse(runs$success %in% TRUE, paste0("optimizer_nonconvergence:", runs$stop_reason), "fit_error"))
  if (any(is.na(runs$stop_reason)) || any(!runs$stop_reason %in% allowed_stops) ||
      !identical(as.character(runs$failure_type), as.character(expected_failure)) ||
      any(eligible & runs$stop_reason %in% c("fit_error", "invalid_likelihood", "numerical_failure")) ||
      any(!runs$success & runs$stop_reason != "fit_error") ||
      any(!runs$success & (is.na(runs$error) | !nzchar(runs$error))) ||
      any(runs$success & !is.na(runs$error) & nzchar(runs$error))) return(FALSE)
  numeric_fields <- names(runs)[vapply(runs, is.numeric, logical(1))]
  integer_diagnostics <- intersect(c("outer_iterations", "best_raw_loglik_iteration",
    "complete_adjacent_pairs", "n_dropout_subjects", "n_subjects_with_interior_gaps"), names(runs))
  if (length(integer_diagnostics) && !all(vapply(runs[integer_diagnostics], is.integer, logical(1)))) return(FALSE)
  rate_fields <- intersect(c("target_missing_rate", "observed_missing_rate",
    "complete_adjacent_pair_rate"), names(runs))
  if (any(vapply(runs[rate_fields], function(x)
      any(!is.na(x) & (x < 0 | x > 1)), logical(1)))) return(FALSE)
  if (any(eligible & (!is.finite(runs$logLik) | !is.finite(runs$df)))) return(FALSE)
  if (any(vapply(numeric_fields, function(nm) any(!is.na(runs[[nm]]) &
      (!is.finite(runs[[nm]]) | abs(runs[[nm]]) > 1e12)), logical(1)))) return(FALSE)
  all(vapply(payload_names, function(name) {
    payload <- result[[name]]
    if (is.null(payload)) return(TRUE)
    required_keys <- c("scenario", "model", "rep")
    if (!all(required_keys %in% names(payload)) || !is.character(payload$scenario) ||
        !is.character(payload$model) || !is.integer(payload$rep)) return(FALSE)
    required_numeric <- switch(name,
      fixed = c("estimate", "true_value"),
      smooth = c("s1", "smooth_hat", "smooth_true"),
      joint = character(), predictive = character())
    if (length(required_numeric) && (!all(required_numeric %in% names(payload)) ||
        !all(vapply(payload[required_numeric], is.numeric, logical(1))))) return(FALSE)
    integer_keys <- intersect(c("n", "d", "simulation_seed", "missingness_seed"), names(payload))
    if (length(integer_keys) && !all(vapply(payload[integer_keys], is.integer, logical(1)))) return(FALSE)
    numeric_fields <- names(payload)[vapply(payload, is.numeric, logical(1))]
    if (any(vapply(numeric_fields, function(field) any(!is.na(payload[[field]]) &
        (!is.finite(payload[[field]]) | abs(payload[[field]]) > 1e12)), logical(1)))) return(FALSE)
    bounded_rates <- intersect(c("target_missing_rate", "observed_missing_rate",
      "complete_adjacent_pair_rate", "s1"), names(payload))
    if (any(vapply(bounded_rates, function(field) any(!is.na(payload[[field]]) &
        (payload[[field]] < 0 | payload[[field]] > 1)), logical(1)))) return(FALSE)
    nonnegative <- intersect(c("std_error", "rmse", "irmse", "variogram_score", "elapsed_sec"), names(payload))
    if (any(vapply(nonnegative, function(field) any(!is.na(payload[[field]]) & payload[[field]] < 0), logical(1)))) return(FALSE)
    filtered <- tryCatch(jss_missing_filter_payload(payload, runs), error = function(e) NULL)
    !is.null(filtered) && nrow(filtered) == nrow(payload)
  }, logical(1)))
}
