#!/usr/bin/env Rscript

# Full, checkpointed ZIBNB/Gumbel build for the ACTG data. This mirrors the
# original staged AIC workflow, adds the ZIBNB tau parameter, and then cycles
# over every parameter until the selected formulas stabilize.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path(
    getwd(), "paper", "reviews", "2026-09-05-medical-candidates",
    "run-zibnb-full-build.R"
  ), mustWork = TRUE)
}
package_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
screen_script <- file.path(dirname(script_path), "run-screening-v2.R")

argument_value <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

workers <- as.integer(argument_value("workers", "4"))
max_cycles <- as.integer(argument_value("max-cycles", "3"))
max_elapsed_sec <- as.numeric(argument_value("max-elapsed-sec", "7200"))
if (!is.finite(workers) || workers < 1L) stop("--workers must be positive.")
if (!is.finite(max_cycles) || max_cycles < 1L) stop("--max-cycles must be positive.")

results_subdir <- "zibnb-full-build"
results_dir <- file.path(
  package_root, "paper", "reviews", "2026-09-05-medical-candidates",
  "results", results_subdir
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(results_dir, "zibnb-full-build-stage-log.csv")
progress_path <- file.path(results_dir, "zibnb-full-build-progress.rds")
stage_rows <- if (file.exists(log_path)) {
  utils::read.csv(log_path, stringsAsFactors = FALSE)
} else {
  data.frame(
    stage = character(), started_at = character(), completed_at = character(),
    status = integer(), formula_path = character(), fit_path = character(),
    stringsAsFactors = FALSE
  )
}

atomic_save_rds <- function(value, path) {
  temporary <- tempfile("zibnb-progress-", tmpdir = dirname(path), fileext = ".rds")
  saveRDS(value, temporary, compress = "gzip")
  copied <- file.copy(temporary, path, overwrite = TRUE)
  unlink(temporary)
  if (!copied) stop("Could not save ", path)
}

write_log <- function() {
  temporary <- tempfile("zibnb-stage-log-", tmpdir = dirname(log_path), fileext = ".csv")
  utils::write.csv(stage_rows, temporary, row.names = FALSE, na = "")
  copied <- file.copy(temporary, log_path, overwrite = TRUE)
  unlink(temporary)
  if (!copied) stop("Could not save ", log_path)
}

formula_signature <- function(path) {
  formulas <- readRDS(path)
  vapply(formulas, function(x) paste(deparse(x, width.cutoff = 500L), collapse = " "), character(1L))
}

run_screen <- function(label, extra_args, stem) {
  formula_path <- file.path(results_dir, paste0(stem, "-final-formulas.rds"))
  fit_path <- file.path(results_dir, paste0(stem, "-final-fit.rds"))
  manifest_path <- file.path(results_dir, paste0(stem, "-run-manifest.rds"))
  if (file.exists(formula_path) && file.exists(fit_path) && file.exists(manifest_path)) {
    manifest <- tryCatch(readRDS(manifest_path), error = function(error) NULL)
    if (is.list(manifest) && identical(manifest$margin_family, "ZIBNB")) {
      cat("Reusing completed stage ", label, ".\n", sep = "")
      return(list(formula_path = formula_path, fit_path = fit_path, reused = TRUE))
    }
  }

  started_at <- Sys.time()
  cat("\nStarting ", label, " at ", format(started_at), "\n", sep = "")
  args <- c(
    screen_script,
    "--margin=ZIBNB",
    "--criterion=AIC",
    "--hierarchy=strong",
    paste0("--workers=", workers),
    paste0("--max-elapsed-sec=", max_elapsed_sec),
    paste0("--results-subdir=", results_subdir),
    paste0("--output-stem=", stem),
    "--resume=true",
    extra_args
  )
  stdout_path <- file.path(results_dir, paste0(stem, ".stdout.log"))
  stderr_path <- file.path(results_dir, paste0(stem, ".stderr.log"))
  status <- system2(
    "Rscript", vapply(args, shQuote, character(1L)),
    stdout = stdout_path, stderr = stderr_path
  )
  completed_at <- Sys.time()
  stage_rows <<- rbind(stage_rows, data.frame(
    stage = label, started_at = format(started_at, tz = "Australia/Sydney"),
    completed_at = format(completed_at, tz = "Australia/Sydney"),
    status = status, formula_path = formula_path, fit_path = fit_path,
    stringsAsFactors = FALSE
  ))
  write_log()
  atomic_save_rds(list(
    updated_at = completed_at, active_stage = label, status = status,
    formula_path = formula_path, fit_path = fit_path
  ), progress_path)
  if (!identical(status, 0L) || !file.exists(formula_path) || !file.exists(fit_path)) {
    stop("Stage ", label, " failed; inspect ", stderr_path)
  }
  list(formula_path = formula_path, fit_path = fit_path, reused = FALSE)
}

initial_stages <- c("mu-smooth", "sigma", "nu", "tau", "theta")
previous <- NULL
for (stage in initial_stages) {
  stem <- paste0("actg-zibnb-backward-selection-", stage)
  value <- run_screen(stage, paste0("--stage=", stage), stem)
  previous <- value
}

stable <- FALSE
cycle_summary <- list()
parameters <- c("mu", "sigma", "nu", "tau", "theta")
for (cycle in seq_len(max_cycles)) {
  cycle_start <- formula_signature(previous$formula_path)
  for (parameter in parameters) {
    stem <- sprintf("actg-zibnb-refinement-cycle-%02d-%s", cycle, parameter)
    label <- sprintf("refinement cycle %02d: %s", cycle, parameter)
    previous <- run_screen(label, c(
      "--stage=refine",
      paste0("--parameters=", parameter),
      paste0("--input-formulas=", normalizePath(previous$formula_path, winslash = "/"))
    ), stem)
  }
  cycle_end <- formula_signature(previous$formula_path)
  stable <- identical(cycle_start, cycle_end)
  cycle_summary[[cycle]] <- data.frame(
    cycle = cycle, stable = stable,
    started_signature = paste(cycle_start, collapse = " | "),
    completed_signature = paste(cycle_end, collapse = " | "),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    do.call(rbind, cycle_summary),
    file.path(results_dir, "zibnb-refinement-cycle-summary.csv"),
    row.names = FALSE
  )
  if (stable) break
}

final_formulas <- readRDS(previous$formula_path)
final_fit <- readRDS(previous$fit_path)
atomic_save_rds(final_formulas, file.path(results_dir, "zibnb-selected-formulas.rds"))
atomic_save_rds(final_fit, file.path(results_dir, "zibnb-selected-fit.rds"))
completion <- list(
  completed_at = Sys.time(), stable = stable,
  cycles_completed = length(cycle_summary), selected_formulas = final_formulas,
  selected_fit_path = file.path(results_dir, "zibnb-selected-fit.rds"),
  selection_log = log_path
)
atomic_save_rds(completion, file.path(results_dir, "ZIBNB-SELECTION-COMPLETE.rds"))
atomic_save_rds(list(
  updated_at = Sys.time(), active_stage = "selection complete", status = 0L,
  stable = stable, cycles_completed = length(cycle_summary)
), progress_path)

cat("\nZIBNB staged selection completed.\n")
print(final_formulas)
print(final_fit$model_selection)
cat("Stable after refinement: ", stable, "\n", sep = "")

run_followup <- function(label, script, args, stem) {
  cat("\nStarting ", label, " at ", format(Sys.time()), "\n", sep = "")
  stdout_path <- file.path(results_dir, paste0(stem, ".stdout.log"))
  stderr_path <- file.path(results_dir, paste0(stem, ".stderr.log"))
  status <- system2(
    "Rscript", vapply(c(script, args), shQuote, character(1L)),
    stdout = stdout_path, stderr = stderr_path
  )
  if (!identical(status, 0L)) {
    stop(label, " failed; inspect ", stderr_path)
  }
  invisible(TRUE)
}

targeted_script <- file.path(dirname(script_path), "run-actg-time1-targeted-screen.R")
targeted_dir <- file.path(results_dir, "time1-targeted-screen")
targeted_completion <- file.path(targeted_dir, "RUN-COMPLETE.rds")
if (!file.exists(targeted_completion)) {
  run_followup(
    "ZIBNB visit-1 targeted screen", targeted_script,
    c(
      "--margin=ZIBNB",
      paste0("--base-fit=", normalizePath(file.path(results_dir, "zibnb-selected-fit.rds"), winslash = "/")),
      paste0("--base-formulas=", normalizePath(file.path(results_dir, "zibnb-selected-formulas.rds"), winslash = "/")),
      paste0("--source-results-dir=", normalizePath(results_dir, winslash = "/")),
      paste0("--output-dir=", normalizePath(targeted_dir, winslash = "/", mustWork = FALSE)),
      paste0("--workers=", workers),
      paste0("--max-elapsed-sec=", max_elapsed_sec)
    ),
    "zibnb-time1-targeted-screen"
  )
}

review_script <- file.path(dirname(script_path), "build-final-review.R")
review_dir <- file.path(results_dir, "zibnb-final-model-review")
review_completion <- file.path(review_dir, "REVIEW-BUNDLE-COMPLETE.rds")
targeted_fit_path <- file.path(targeted_dir, "time1-revised-candidate-fit.rds")
if (!file.exists(review_completion)) {
  run_followup(
    "ZIBNB final diagnostics and publication bundle", review_script,
    c(
      paste0("--fit=", normalizePath(targeted_fit_path, winslash = "/")),
      paste0("--source-results-dir=", normalizePath(results_dir, winslash = "/")),
      paste0("--output-dir=", normalizePath(review_dir, winslash = "/", mustWork = FALSE)),
      "--pit-seed=20260917"
    ),
    "zibnb-final-review"
  )
}

atomic_save_rds(list(
  completed_at = Sys.time(), selection = file.path(results_dir, "ZIBNB-SELECTION-COMPLETE.rds"),
  targeted_screen = targeted_completion, final_fit = targeted_fit_path,
  review_bundle = review_dir
), file.path(results_dir, "ZIBNB-FULL-BUILD-COMPLETE.rds"))
cat("\nComplete ZIBNB review bundle: ", review_dir, "\n", sep = "")
