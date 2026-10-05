#!/usr/bin/env Rscript

# Route-aware replacement for the ad-hoc ACTG joint-screen launcher retained
# in the original paper checkout's temp directory.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: run-live-joint-screen-candidate.R MARGIN COPULA OUTPUT")
}

margin_family <- toupper(args[[1L]])
copula_family <- if (tolower(args[[2L]]) == "t") "t" else toupper(args[[2L]])
output_path <- args[[3L]]

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]),
                             winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."),
                      winslash = "/", mustWork = TRUE)
setwd(root)
suppressPackageStartupMessages(pkgload::load_all(root, quiet = TRUE))

workflow <- new.env(parent = globalenv())
sys.source(file.path(root, "paper", "R", "05-actg-workflow.R"), envir = workflow)
data <- workflow$jss_actg_load_data()
control <- workflow$jss_actg_optimizer_control(
  max_elapsed_sec = 1200,
  margin_family = margin_family,
  copula_family = copula_family,
  screening = TRUE
)

result <- select_joint_distribution(
  data = data,
  time_var = "week_factor",
  subject_var = "subject",
  response_var = "cd4",
  margin_families = margin_family,
  copula_families = copula_family,
  time_intercepts = TRUE,
  candidate_timeout = 1250,
  isolate_candidates = "always",
  progress = TRUE,
  checkpoint_dir = paste0(tools::file_path_sans_ext(output_path), "-checkpoints"),
  resume = TRUE,
  fit_args = list(
    missingness = "segment",
    method = "RS",
    optimization = "joint",
    optimizer_control = control,
    compute_vcov = FALSE,
    verbose = 0
  )
)

saveRDS(result, output_path, version = 3L)
print(result)
