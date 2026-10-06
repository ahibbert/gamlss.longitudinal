#!/usr/bin/env Rscript

# Read-only audit for the checkpointed four-route ACTG backward-selection run.
# It is safe to run while the build is active: selection tables are written
# atomically, and incomplete routes are reported rather than treated as final.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path(
    getwd(), "paper", "reviews", "2026-09-05-medical-candidates",
    "audit-zero-inflated-four-route-build.R"
  ), mustWork = TRUE)
}
review_dir <- dirname(script_path)
package_root <- normalizePath(file.path(review_dir, "..", "..", ".."), mustWork = TRUE)

argument_value <- function(name, default = "") {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

results_prefix <- argument_value("results-prefix")
reference_prefix <- argument_value("reference-prefix")
limit_sec <- suppressWarnings(as.numeric(argument_value("limit-sec", "1200")))
output_path <- argument_value("output")
results_root_argument <- argument_value("results-root")
require_complete <- tolower(argument_value("require-complete", "false")) %in%
  c("true", "1", "yes")

if (!nzchar(results_prefix)) stop("--results-prefix is required.", call. = FALSE)
if (!is.finite(limit_sec) || limit_sec <= 0) {
  stop("--limit-sec must be a positive finite number.", call. = FALSE)
}

routes <- c("zisichel-g", "zisichel-t", "zibnb-g", "zibnb-t")
results_root <- if (nzchar(results_root_argument)) {
  normalizePath(results_root_argument, mustWork = TRUE)
} else {
  file.path(review_dir, "results")
}
package_source_files <- c(
  file.path(package_root, "DESCRIPTION"),
  list.files(file.path(package_root, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
)
current_source_fingerprint <- tools::md5sum(package_source_files)

canonical_fingerprint <- function(value) {
  if (is.null(value) || !length(value)) return(character())
  keys <- sub("^.*[/\\\\](DESCRIPTION|R[/\\\\].*)$", "\\1", names(value))
  keys <- gsub("\\\\", "/", keys)
  stats::setNames(unname(value), keys)[order(keys)]
}
current_source_fingerprint_canonical <- canonical_fingerprint(
  current_source_fingerprint
)

empty_rows <- function() {
  data.frame(
    route = character(), stage = character(), table_path = character(),
    step = integer(), parameter = character(), dropped_term = character(),
    change_type = character(), candidate_formula = character(),
    logLik = numeric(), AIC = numeric(), BIC = numeric(), EDF = numeric(),
    converged = logical(), accepted = logical(), elapsed_sec = numeric(),
    outer_iterations = integer(), stop_reason = character(), warnings = character(),
    error = character(), checkpoint = character(), stringsAsFactors = FALSE
  )
}

logical_column <- function(x) {
  if (is.logical(x)) return(x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1")
}

selection_tables <- function(prefix, route) {
  route_dir <- file.path(results_root, paste0(prefix, "-", route))
  if (!dir.exists(route_dir)) return(empty_rows())
  paths <- sort(list.files(
    route_dir,
    pattern = "-model-selection-table\\.csv$",
    full.names = TRUE
  ))
  if (!length(paths)) return(empty_rows())
  rows <- lapply(paths, function(path) {
    value <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    value$route <- route
    value$stage <- sub("-model-selection-table\\.csv$", "", basename(path))
    value$table_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
    value$converged <- logical_column(value$converged)
    value$accepted <- logical_column(value$accepted)
    value
  })
  common <- Reduce(intersect, lapply(rows, names))
  do.call(rbind, lapply(rows, function(value) value[common]))
}

route_completion <- function(prefix, route) {
  route_dir <- file.path(results_root, paste0(prefix, "-", route))
  path <- file.path(route_dir, "SELECTION-COMPLETE.rds")
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(error) NULL)
}

manifest_fingerprint_ok <- function(prefix, route) {
  route_dir <- file.path(results_root, paste0(prefix, "-", route))
  paths <- list.files(route_dir, pattern = "-run-manifest\\.rds$", full.names = TRUE)
  if (!length(paths)) return(NA)
  all(vapply(paths, function(path) {
    manifest <- tryCatch(readRDS(path), error = function(error) NULL)
    is.list(manifest) && identical(
      canonical_fingerprint(manifest$package_source_fingerprint),
      current_source_fingerprint_canonical
    )
  }, logical(1L)))
}

key_columns <- c(
  "stage", "step", "parameter", "dropped_term", "change_type",
  "candidate_formula", "mu_formula", "sigma_formula", "nu_formula",
  "tau_formula", "theta_formula", "zeta_formula"
)

row_key <- function(value) {
  columns <- intersect(key_columns, names(value))
  if (!length(columns)) return(rep("", nrow(value)))
  do.call(paste, c(lapply(value[columns], function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "<NA>"
    x
  }), sep = "\r"))
}

compare_reference <- function(current, reference) {
  if (!nrow(current) || !nrow(reference)) {
    return(list(
      overlap = 0L, missing_reference = nrow(current), decision_mismatches = NA_integer_,
      max_abs_logLik = NA_real_, max_abs_AIC = NA_real_,
      warning_mismatches = NA_integer_
    ))
  }
  current_key <- row_key(current)
  reference_key <- row_key(reference)
  matched <- match(current_key, reference_key)
  overlap <- !is.na(matched)
  if (!any(overlap)) {
    return(list(
      overlap = 0L, missing_reference = nrow(current), decision_mismatches = NA_integer_,
      max_abs_logLik = NA_real_, max_abs_AIC = NA_real_,
      warning_mismatches = NA_integer_
    ))
  }
  old <- reference[matched[overlap], , drop = FALSE]
  new <- current[overlap, , drop = FALSE]
  decision_mismatches <- sum(
    old$converged != new$converged | old$accepted != new$accepted,
    na.rm = TRUE
  )
  list(
    overlap = sum(overlap),
    missing_reference = sum(!overlap),
    decision_mismatches = decision_mismatches,
    max_abs_logLik = max(abs(old$logLik - new$logLik), na.rm = TRUE),
    max_abs_AIC = max(abs(old$AIC - new$AIC), na.rm = TRUE),
    warning_mismatches = sum(
      trimws(old$warnings) != trimws(new$warnings),
      na.rm = TRUE
    )
  )
}

unexpected_warning_count <- function(warnings, converged) {
  if (!length(warnings)) return(0L)
  sum(vapply(seq_along(warnings), function(index) {
    value <- warnings[[index]]
    if (is.na(value) || !nzchar(trimws(value))) return(FALSE)
    events <- trimws(strsplit(value, " | ", fixed = TRUE)[[1L]])
    any(vapply(events, function(event) {
      expected <- grepl("segmented likelihood", event, fixed = TRUE) ||
        grepl("Analytical Hessian vcov failed", event, fixed = TRUE) ||
        grepl("fitted_gradient_too_large", event, fixed = TRUE) ||
        (!isTRUE(converged[[index]]) &&
          grepl("without satisfying the optimizer convergence contract", event, fixed = TRUE))
      !expected
    }, logical(1L)))
  }, logical(1L)))
}

formula_signature <- function(formulas) {
  if (!is.list(formulas)) return(character())
  vapply(formulas, function(formula) {
    if (is.null(formula)) return("<NULL>")
    paste(deparse(formula, width.cutoff = 500L), collapse = " ")
  }, character(1L))
}

summary_rows <- lapply(routes, function(route) {
  current <- selection_tables(results_prefix, route)
  reference <- if (nzchar(reference_prefix)) {
    selection_tables(reference_prefix, route)
  } else {
    empty_rows()
  }
  comparison <- compare_reference(current, reference)
  completion <- route_completion(results_prefix, route)
  reference_completion <- if (nzchar(reference_prefix)) {
    route_completion(reference_prefix, route)
  } else {
    NULL
  }
  complete <- is.list(completion) && isTRUE(completion$stable)
  selected_formulas_match <- if (complete && is.list(reference_completion) &&
      isTRUE(reference_completion$stable)) {
    identical(
      formula_signature(completion$selected_formulas),
      formula_signature(reference_completion$selected_formulas)
    )
  } else {
    NA
  }
  finite_elapsed <- current$elapsed_sec[is.finite(current$elapsed_sec)]
  warning_values <- if ("warnings" %in% names(current)) {
    current$warnings
  } else {
    character(nrow(current))
  }
  warning_values[is.na(warning_values)] <- ""
  warning_rows <- nzchar(trimws(warning_values))
  data.frame(
    route = route,
    complete = complete,
    stages = length(unique(current$stage)),
    candidates = nrow(current),
    converged = sum(current$converged %in% TRUE),
    nonconverged = sum(!(current$converged %in% TRUE)),
    max_elapsed_sec = if (length(finite_elapsed)) max(finite_elapsed) else NA_real_,
    runtime_pass = length(finite_elapsed) == nrow(current) &&
      all(finite_elapsed < limit_sec),
    warning_rows = sum(warning_rows),
    unexpected_warning_rows = unexpected_warning_count(
      warning_values, current$converged
    ),
    fingerprint_pass = manifest_fingerprint_ok(results_prefix, route),
    reference_overlap = comparison$overlap,
    reference_missing = comparison$missing_reference,
    decision_mismatches = comparison$decision_mismatches,
    max_abs_logLik = comparison$max_abs_logLik,
    max_abs_AIC = comparison$max_abs_AIC,
    warning_mismatches = comparison$warning_mismatches,
    selected_formulas_match = selected_formulas_match,
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)

print(summary, row.names = FALSE)

if (nzchar(output_path)) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("actg-four-route-audit-", tmpdir = dirname(output_path), fileext = ".csv")
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(summary, temporary, row.names = FALSE, na = "")
  if (!file.copy(temporary, output_path, overwrite = TRUE)) {
    stop("Could not write audit output: ", output_path, call. = FALSE)
  }
}

failures <- character()
present <- summary$candidates > 0L
if (any(present & !summary$runtime_pass)) failures <- c(failures, "runtime")
if (any(present & summary$fingerprint_pass %in% FALSE)) failures <- c(failures, "fingerprint")
if (any(summary$unexpected_warning_rows > 0L)) failures <- c(failures, "warnings")
if (any(summary$decision_mismatches > 0L, na.rm = TRUE)) failures <- c(failures, "decisions")
if (any(summary$max_abs_logLik > 1e-8, na.rm = TRUE) ||
    any(summary$max_abs_AIC > 1e-8, na.rm = TRUE)) {
  failures <- c(failures, "metrics")
}
if (any(summary$selected_formulas_match %in% FALSE)) failures <- c(failures, "formulas")
if (require_complete && any(!summary$complete)) failures <- c(failures, "incomplete")

if (length(failures)) {
  stop(
    "ACTG four-route audit failed: ", paste(unique(failures), collapse = ", "),
    call. = FALSE
  )
}
