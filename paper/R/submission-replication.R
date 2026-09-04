jss_submission_schema_version <- function() "jss-submission-inputs-2"

jss_submission_input_registry <- function(root = getwd()) {
  path <- file.path(root, "paper", "submission-input-schemas.csv")
  if (!file.exists(path)) stop("Submission input-schema registry is missing: ", path, call. = FALSE)
  registry <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(
    "schema_version", "study", "relative_path", "min_rows", "exact_rows",
    "key_columns", "required_columns", "numeric_columns", "rate_columns"
  )
  if (!identical(names(registry), required) || !nrow(registry) ||
      anyDuplicated(registry$study) ||
      any(registry$schema_version != jss_submission_schema_version())) {
    stop("Submission input-schema registry is malformed or has an unsupported version.", call. = FALSE)
  }
  registry
}

jss_submission_split_fields <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) return(character())
  strsplit(x, "|", fixed = TRUE)[[1L]]
}

jss_submission_read_registered_input <- function(study, path, root = getwd()) {
  registry <- jss_submission_input_registry(root)
  row <- registry[registry$study == study, , drop = FALSE]
  if (nrow(row) != 1L) stop("No unique submission schema is registered for study: ", study, call. = FALSE)
  if (!file.exists(path)) stop("Registered submission input is missing: ", path, call. = FALSE)
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- jss_submission_split_fields(row$required_columns)
  keys <- jss_submission_split_fields(row$key_columns)
  numeric <- jss_submission_split_fields(row$numeric_columns)
  rates <- jss_submission_split_fields(row$rate_columns)
  missing <- setdiff(unique(c(required, keys, numeric, rates)), names(x))
  if (length(missing)) {
    stop("Registered submission input lacks columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  exact_rows <- suppressWarnings(as.integer(row$exact_rows))
  min_rows <- suppressWarnings(as.integer(row$min_rows))
  if ((!is.na(exact_rows) && nrow(x) != exact_rows) ||
      (is.na(exact_rows) && (is.na(min_rows) || nrow(x) < min_rows))) {
    stop("Registered submission input has an unexpected row count for ", study, ".", call. = FALSE)
  }
  if (length(keys) && anyDuplicated(x[keys])) {
    stop("Registered submission input has duplicate keys for ", study, ".", call. = FALSE)
  }
  for (field in numeric) {
    value <- suppressWarnings(as.numeric(x[[field]]))
    if (any(is.na(value) & !is.na(x[[field]])) || any(!is.na(value) & !is.finite(value))) {
      stop("Registered numeric column is invalid for ", study, ": ", field, call. = FALSE)
    }
  }
  for (field in rates) {
    value <- suppressWarnings(as.numeric(x[[field]]))
    if (any(!is.na(value) & (value < 0 | value > 1))) {
      stop("Registered rate column is outside [0, 1] for ", study, ": ", field, call. = FALSE)
    }
  }
  success_field <- if ("success" %in% names(x)) "success" else if ("fit_returned" %in% names(x)) "fit_returned" else NA_character_
  lifecycle <- unique(c(intersect("attempted", names(x)), success_field,
    intersect(c("converged", "retained"), names(x))))
  lifecycle <- lifecycle[!is.na(lifecycle)]
  if (!is.na(success_field) && all(c("converged", "retained") %in% lifecycle)) {
    logical_value <- lapply(x[lifecycle], function(value) {
      if (is.logical(value)) return(value)
      normalized <- tolower(trimws(as.character(value)))
      if (any(!normalized %in% c("true", "false", "1", "0"))) {
        stop("Attempt lifecycle columns must be logical for ", study, ".", call. = FALSE)
      }
      normalized %in% c("true", "1")
    })
    if (("attempted" %in% lifecycle && any(!logical_value$attempted)) ||
        any(logical_value$retained & !logical_value$converged) ||
        any(logical_value$converged & !logical_value[[success_field]]) ||
        ("attempted" %in% lifecycle && any(logical_value[[success_field]] & !logical_value$attempted))) {
      stop("Attempt lifecycle does not reconcile for ", study, ".", call. = FALSE)
    }
    failure_field <- intersect(c("failure_reason", "failure_type", "exclusion_reason"), names(x))
    if (length(failure_field)) {
      reason_ok <- Reduce(`|`, lapply(x[failure_field], function(value) {
        reason <- trimws(as.character(value)); reason[is.na(reason)] <- ""
        nzchar(reason) & tolower(reason) != "none"
      }))
      failed <- !logical_value$retained
      if (any(failed & !reason_ok)) {
        stop("Unretained attempts require an explicit failure reason for ", study, ".", call. = FALSE)
      }
    }
  }
  attr(x, "submission_schema") <- list(
    version = row$schema_version[[1L]], study = study, rows = nrow(x), keys = keys
  )
  x
}

jss_submission_require_files <- function(paths, label) {
  paths <- unname(paths[nzchar(paths)])
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop(label, " is incomplete: ", paste(basename(missing), collapse = ", "), call. = FALSE)
  small <- paths[file.info(paths)$isdir %in% FALSE & file.info(paths)$size < 20L]
  if (length(small)) stop(label, " contains an empty/truncated artifact: ", paste(basename(small), collapse = ", "), call. = FALSE)
  invisible(paths)
}

jss_submission_validate_figures <- function(paths) {
  paths <- unname(paths[nzchar(paths)])
  if (length(paths) && any(!vapply(paths, jss_png_nonblank, logical(1L)))) {
    stop("Submission figure is missing, blank, or has invalid dimensions.", call. = FALSE)
  }
  invisible(TRUE)
}

jss_submission_validate_main_recovery <- function(path, root = getwd()) {
  files <- file.path(path, jss_main_recovery_artifact_names())
  jss_submission_require_files(files, "Main-recovery bundle")
  jss_main_recovery_validate_public_bundle(
    path, repo_root = root, require_attestation = FALSE
  )
  attempts <- jss_submission_read_registered_input(
    "main-recovery", file.path(path, "attempt_metadata.csv"), root
  )
  if (!setequal(unique(attempts$study_id),
      c("bcpe_t_main_recovery", "nbi_clayton_main_recovery")) ||
      !setequal(unique(attempts$margin_family), c("BCPE", "NBI")) ||
      !setequal(unique(attempts$copula_code), c("t", "C")) ||
      !all(as.integer(attempts$n_subjects) == 500L) ||
      !all(as.integer(attempts$n_time) == 4L)) {
    stop("Main-recovery design metadata do not match the registered n=500, T=4 studies.", call. = FALSE)
  }
  attr(attempts, "submission_schema")
}

jss_submission_validate_optimizer <- function(data_dir, tables_dir, figures_dir, root = getwd()) {
  files <- jss_phase2_optimizer_files()
  paths <- c(file.path(data_dir, files$data), file.path(tables_dir, files$tables),
    file.path(figures_dir, files$figures))
  jss_submission_require_files(paths, "Optimizer benchmark bundle")
  attempts <- jss_submission_read_registered_input(
    "optimizer-benchmark",
    file.path(data_dir, "03-joint-vs-separate-optimization-results.csv"), root
  )
  if (!setequal(unique(attempts$method), c("rs_separate", "rs_joint"))) {
    stop("Optimizer benchmark must contain exactly the registered RS methods.", call. = FALSE)
  }
  pair_key <- interaction(attempts[c("case_id", "joint_review_rep")],
    drop = TRUE, lex.order = TRUE, sep = "\r")
  paired <- split(seq_len(nrow(attempts)), pair_key, drop = TRUE)
  if (any(vapply(paired, function(i) {
      length(i) != 2L || !setequal(attempts$method[i], c("rs_separate", "rs_joint")) ||
        length(unique(attempts$paired_seed[i])) != 1L
    }, logical(1L)))) {
    stop("Optimizer benchmark rows do not reconcile into paired methods and seeds.", call. = FALSE)
  }
  jss_submission_validate_figures(file.path(figures_dir, files$figures))
  attr(attempts, "submission_schema")
}

jss_submission_validate_missingness <- function(path, root = getwd()) {
  artifacts <- file.path(path, jss_missingness_public_artifacts())
  jss_submission_require_files(artifacts, "Missingness bundle")
  attempts <- jss_submission_read_registered_input(
    "missingness", file.path(path, "fit_run_log.csv"), root
  )
  if (!setequal(unique(attempts$missing_mechanism),
      c("monotone_dropout", "time_dependent_intermittent_mar"))) {
    stop("Missingness inputs do not contain the registered mechanisms.", call. = FALSE)
  }
  if (!setequal(unique(attempts$model), c("gamlss2", "gamlss.longitudinal"))) {
    stop("Missingness inputs do not contain exactly the registered model pair.", call. = FALSE)
  }
  cell <- interaction(attempts[c("scenario", "rep")], drop = TRUE,
    lex.order = TRUE, sep = "\r")
  if (any(vapply(split(seq_len(nrow(attempts)), cell, drop = TRUE), function(i)
      length(i) != 2L || !setequal(attempts$model[i], c("gamlss2", "gamlss.longitudinal")), logical(1L)))) {
    stop("Missingness rows do not reconcile into paired model attempts.", call. = FALSE)
  }
  attr(attempts, "submission_schema")
}

jss_submission_validate_copula <- function(path, root = getwd()) {
  results <- jss_submission_read_registered_input("copula-misspecification", path, root)
  if (!setequal(unique(results$fitted_copula), c("C", "N")) ||
      !setequal(unique(results$generating_copula), c("C", "N")) ||
      any(as.integer(results$n_time) != 4L)) {
    stop("Copula-misspecification inputs do not match the registered fit grid.", call. = FALSE)
  }
  dataset <- interaction(results[c(
    "generating_copula", "tau_label", "n_subject", "rep"
  )], drop = TRUE, lex.order = TRUE, sep = "\r")
  if (any(vapply(split(seq_len(nrow(results)), dataset, drop = TRUE), function(i)
      length(i) != 2L || !setequal(results$fitted_copula[i], c("C", "N")), logical(1L)))) {
    stop("Copula-misspecification rows do not reconcile into paired candidate fits.", call. = FALSE)
  }
  results
}

jss_submission_validate_scaling <- function(design, attempts, summary, contrasts, hardware,
    root = getwd()) {
  jss_submission_require_files(c(design, attempts, summary, contrasts, hardware), "Fit-scaling bundle")
  value <- jss_submission_read_registered_input("fit-scaling", attempts, root)
  attr(value, "submission_schema")
}

jss_submission_stage_main_recovery <- function(source_dir, destination_dir, root = getwd()) {
  jss_submission_validate_main_recovery(source_dir, root)
  names <- jss_main_recovery_artifact_names()
  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
  destination <- file.path(destination_dir, names)
  if (!all(file.copy(file.path(source_dir, names), destination, overwrite = TRUE))) {
    stop("Could not stage the versioned main-recovery bundle.", call. = FALSE)
  }
  jss_submission_validate_main_recovery(destination_dir, root)
  destination
}
