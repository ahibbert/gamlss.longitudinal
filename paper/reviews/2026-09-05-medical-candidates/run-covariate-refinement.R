#!/usr/bin/env Rscript

# Revisit all four parameter formulas after the initial sequential ACTG search.
# Each coordinate is reset to the full smooth-aware formula, backward-selected
# conditional on the latest other formulas, and complete sweeps repeat until
# the selected formula set is stable.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path(
    getwd(), "paper", "reviews", "2026-09-05-medical-candidates",
    "run-covariate-refinement.R"
  ), mustWork = TRUE)
}
package_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
screening_script <- file.path(dirname(script_path), "run-screening-v2.R")
review_script <- file.path(dirname(script_path), "build-final-review.R")
results_dir <- file.path(dirname(script_path), "results")

argument_value <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

argument_flag <- function(name, default) {
  value <- tolower(argument_value(name, if (default) "true" else "false"))
  if (!value %in% c("true", "false", "1", "0", "yes", "no")) {
    stop("--", name, " must be true or false.")
  }
  value %in% c("true", "1", "yes")
}

workers <- as.integer(argument_value("workers", "4"))
smooth_k <- as.integer(argument_value("smooth-k", "8"))
max_cycles <- as.integer(argument_value("max-cycles", "5"))
dry_run <- argument_flag("dry-run", FALSE)
if (!is.finite(workers) || workers < 1L) stop("--workers must be positive.")
if (!is.finite(smooth_k) || smooth_k < 3L) stop("--smooth-k must be at least 3.")
if (!is.finite(max_cycles) || max_cycles < 1L) stop("--max-cycles must be positive.")

formula_text <- function(formula) paste(deparse(formula, width.cutoff = 500L), collapse = " ")
formula_signature <- function(formulas) vapply(
  formulas[c("mu", "sigma", "nu", "theta")], formula_text, character(1L)
)

source_files <- c(
  file.path(package_root, "DESCRIPTION"),
  list.files(file.path(package_root, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
)
source_fingerprint <- tools::md5sum(source_files)

manifest_path <- function(stem) file.path(results_dir, paste0(stem, "-run-manifest.rds"))
formula_path <- function(stem) file.path(results_dir, paste0(stem, "-final-formulas.rds"))
fit_path <- function(stem) file.path(results_dir, paste0(stem, "-final-fit.rds"))

valid_manifest <- function(stem, parameter) {
  path <- manifest_path(stem)
  if (!file.exists(path) || !file.exists(formula_path(stem)) || !file.exists(fit_path(stem))) return(FALSE)
  manifest <- tryCatch(readRDS(path), error = function(error) NULL)
  is.list(manifest) && identical(manifest$stage, "refine") &&
    identical(manifest$parameters, parameter) &&
    identical(manifest$package_source_fingerprint, source_fingerprint) &&
    is.data.frame(manifest$model_table) && nrow(manifest$model_table) > 0L &&
    all(manifest$model_table$converged)
}

run_coordinate <- function(cycle, parameter, input_path, attempts = 3L) {
  stem <- sprintf("actg-covariate-refinement-cycle-%02d-%s", cycle, parameter)
  if (valid_manifest(stem, parameter)) {
    cat("Reusing completed ", stem, ".\n", sep = "")
    return(stem)
  }

  rscript <- file.path(R.home("bin"), "Rscript.exe")
  for (attempt in seq_len(attempts)) {
    cat(
      "\nStarting ", stem, " (attempt ", attempt, ") at ",
      format(Sys.time()), "...\n", sep = ""
    )
    status <- system2(
      rscript,
      c(
        shQuote(screening_script),
        "--stage=refine",
        paste0("--parameters=", parameter),
        paste0("--input-formulas=", shQuote(input_path)),
        paste0("--output-stem=", stem),
        paste0("--smooth-k=", smooth_k),
        "--subjects=all",
        paste0("--workers=", workers),
        "--criterion=AIC",
        "--hierarchy=strong",
        "--resume=true",
        "--max-elapsed-sec=7200"
      ),
      stdout = "",
      stderr = ""
    )
    if (identical(status, 0L) && valid_manifest(stem, parameter)) return(stem)
    cat(stem, " was incomplete; retrying failed models.\n", sep = "")
  }
  stop(stem, " did not complete after ", attempts, " attempts.")
}

initial_formulas_path <- file.path(
  results_dir, "actg-backward-selection-theta-final-formulas.rds"
)
if (!file.exists(initial_formulas_path)) {
  stop("Initial theta-stage formulas were not found: ", initial_formulas_path)
}

if (dry_run) {
  cat(
    "ACTG coordinate-refinement dry run\n",
    "  initial formulas: ", initial_formulas_path, "\n",
    "  order: mu -> sigma -> nu -> theta\n",
    "  maximum cycles: ", max_cycles, "\n",
    "  workers: ", workers, "\n",
    "  smooth k: ", smooth_k, "\n",
    sep = ""
  )
  print(readRDS(initial_formulas_path))
  quit(save = "no", status = 0L)
}

parameters <- c("mu", "sigma", "nu", "theta")
current_path <- initial_formulas_path
cycle_rows <- list()
stable <- FALSE
last_stem <- NULL

for (cycle in seq_len(max_cycles)) {
  cycle_start <- readRDS(current_path)
  cycle_start_signature <- formula_signature(cycle_start)
  cat("\n=== Covariate refinement cycle ", cycle, " ===\n", sep = "")

  for (parameter in parameters) {
    input_signature <- formula_signature(readRDS(current_path))
    stem <- run_coordinate(cycle, parameter, current_path)
    manifest <- readRDS(manifest_path(stem))
    selected <- readRDS(formula_path(stem))
    output_signature <- formula_signature(selected)
    cycle_rows[[length(cycle_rows) + 1L]] <- data.frame(
      cycle = cycle,
      parameter = parameter,
      stem = stem,
      input_formula = input_signature[[parameter]],
      selected_formula = output_signature[[parameter]],
      formula_changed = !identical(input_signature[[parameter]], output_signature[[parameter]]),
      AIC = tail(manifest$accepted_path$AIC, 1L),
      BIC = tail(manifest$accepted_path$BIC, 1L),
      EDF = tail(manifest$accepted_path$EDF, 1L),
      candidates = nrow(manifest$model_table),
      stringsAsFactors = FALSE
    )
    current_path <- formula_path(stem)
    last_stem <- stem
  }

  cycle_end_signature <- formula_signature(readRDS(current_path))
  stable <- identical(cycle_start_signature, cycle_end_signature)
  cycle_table <- do.call(rbind, cycle_rows)
  utils::write.csv(
    cycle_table,
    file.path(results_dir, "actg-covariate-refinement-cycle-summary.csv"),
    row.names = FALSE,
    na = ""
  )
  saveRDS(
    list(
      updated_at = Sys.time(), stable = stable, completed_cycles = cycle,
      current_formulas = readRDS(current_path), cycle_table = cycle_table,
      source_fingerprint = source_fingerprint
    ),
    file.path(results_dir, "actg-covariate-refinement-progress.rds")
  )
  cat("Cycle ", cycle, if (stable) " is stable.\n" else " changed at least one formula.\n", sep = "")
  if (stable) break
}

if (!stable) stop("Covariate refinement did not stabilize within ", max_cycles, " cycles.")

final_formulas <- readRDS(current_path)
final_fit_source <- fit_path(last_stem)
final_fit_path <- file.path(results_dir, "actg-covariate-refinement-final-fit.rds")
final_formulas_path <- file.path(results_dir, "actg-covariate-refinement-final-formulas.rds")
if (!file.copy(final_fit_source, final_fit_path, overwrite = TRUE)) {
  stop("Could not copy the stable final fit.")
}
saveRDS(final_formulas, final_formulas_path)

final_fit <- readRDS(final_fit_path)
final_manifest_path <- file.path(results_dir, "actg-covariate-refinement-run-manifest.rds")
saveRDS(list(
  completed_at = Sys.time(), stable = stable,
  completed_cycles = max(cycle_rows[[length(cycle_rows)]]$cycle),
  selected_formulas = final_formulas,
  final_fit_path = final_fit_path,
  final_AIC = AIC(final_fit),
  final_BIC = BIC(final_fit),
  cycle_table = do.call(rbind, cycle_rows),
  source_fingerprint = source_fingerprint
), final_manifest_path)

review_output <- file.path(results_dir, "actg-final-model-review-refined")
review_completion <- file.path(review_output, "REVIEW-BUNDLE-COMPLETE.rds")
cat("\nBuilding refined final diagnostics and publication tables...\n")
review_status <- system2(
  file.path(R.home("bin"), "Rscript.exe"),
  c(
    shQuote(review_script),
    paste0("--fit=", shQuote(final_fit_path)),
    paste0("--output-dir=", shQuote(review_output)),
    "--pit-seed=20260909"
  ),
  stdout = "",
  stderr = ""
)
if (!identical(review_status, 0L) || !file.exists(review_completion)) {
  stop("Refined final review bundle did not complete successfully.")
}

cat("\nStable refinement complete.\n")
print(final_formulas)
cat(
  "Final AIC: ", format(AIC(final_fit), digits = 12),
  "\nFinal BIC: ", format(BIC(final_fit), digits = 12),
  "\nReview folder: ", review_output, "\n",
  sep = ""
)
