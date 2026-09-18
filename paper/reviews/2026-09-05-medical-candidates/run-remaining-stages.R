#!/usr/bin/env Rscript

# Continue the checkpointed ACTG selection stages after the smooth-aware mu
# run. This process waits for mu rather than starting overlapping model fits.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path(
    getwd(), "paper", "reviews", "2026-09-05-medical-candidates",
    "run-remaining-stages.R"
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

workers <- as.integer(argument_value("workers", "4"))
smooth_k <- as.integer(argument_value("smooth-k", "8"))
wait_hours <- as.numeric(argument_value("wait-hours", "24"))
if (!is.finite(workers) || workers < 1L) stop("--workers must be positive.")
if (!is.finite(smooth_k) || smooth_k < 3L) stop("--smooth-k must be at least 3.")
if (!is.finite(wait_hours) || wait_hours <= 0) stop("--wait-hours must be positive.")

manifest_path <- function(stage) {
  file.path(
    results_dir,
    paste0("actg-backward-selection-", stage, "-run-manifest.rds")
  )
}

valid_manifest <- function(stage) {
  path <- manifest_path(stage)
  if (!file.exists(path)) return(FALSE)
  manifest <- tryCatch(readRDS(path), error = function(error) NULL)
  is.list(manifest) && identical(manifest$stage, stage) &&
    is.data.frame(manifest$model_table) && nrow(manifest$model_table) > 0L &&
    all(manifest$model_table$converged)
}

run_stage <- function(stage, attempts = 3L) {
  rscript <- file.path(R.home("bin"), "Rscript.exe")
  for (attempt in seq_len(attempts)) {
    cat(
      "\nStarting stage ", stage, " (attempt ", attempt, ") at ",
      format(Sys.time()), "...\n", sep = ""
    )
    status <- system2(
      rscript,
      c(
        shQuote(screening_script),
        paste0("--stage=", stage),
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
    if (identical(status, 0L) && valid_manifest(stage)) return(invisible(TRUE))
    cat("Stage ", stage, " was incomplete; retrying failed models.\n", sep = "")
  }
  stop("Stage '", stage, "' did not complete successfully after ", attempts, " attempts.")
}

cat("Waiting for the smooth-aware mu stage to complete...\n")
started <- Sys.time()
while (!valid_manifest("mu-smooth")) {
  if (as.numeric(difftime(Sys.time(), started, units = "hours")) >= wait_hours) {
    stop("Timed out waiting for the smooth-aware mu stage.")
  }
  if (file.exists(manifest_path("mu-smooth"))) {
    run_stage("mu-smooth")
  } else {
    Sys.sleep(60)
  }
}
cat("Smooth-aware mu stage complete.\n")

for (stage in c("sigma", "nu", "theta")) {
  if (valid_manifest(stage)) {
    cat("Stage ", stage, " is already complete; reusing it.\n", sep = "")
    next
  }
  run_stage(stage)
  cat("Completed stage ", stage, " at ", format(Sys.time()), ".\n", sep = "")
}

completion_path <- file.path(
  results_dir, "actg-backward-selection-all-stages-complete.rds"
)
saveRDS(list(
  completed_at = Sys.time(),
  stages = c("mu-smooth", "sigma", "nu", "theta"),
  manifests = setNames(
    lapply(c("mu-smooth", "sigma", "nu", "theta"), manifest_path),
    c("mu-smooth", "sigma", "nu", "theta")
  )
), completion_path)
cat("\nAll staged searches completed: ", completion_path, "\n", sep = "")

review_output_dir <- file.path(results_dir, "actg-final-model-review")
review_completion <- file.path(review_output_dir, "REVIEW-BUNDLE-COMPLETE.rds")
if (!file.exists(review_completion)) {
  cat("\nBuilding final diagnostics and publication tables...\n")
  review_status <- system2(
    file.path(R.home("bin"), "Rscript.exe"),
    c(
      shQuote(review_script),
      paste0("--fit=", shQuote(file.path(results_dir, "actg-backward-selection-theta-final-fit.rds"))),
      paste0("--output-dir=", shQuote(review_output_dir)),
      "--pit-seed=20260909"
    ),
    stdout = "",
    stderr = ""
  )
  if (!identical(review_status, 0L) || !file.exists(review_completion)) {
    stop("Final review bundle did not complete successfully.")
  }
} else {
  cat("Final review bundle is already complete; reusing it.\n")
}
cat("Final review bundle: ", review_output_dir, "\n", sep = "")
