#!/usr/bin/env Rscript

# Export completed author-workflow checkpoints into the small reviewer-facing
# supplementary bundle. Paths are inputs, never published evidence fields.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: run-03-optimizer-publication-lean.R CHECKPOINT_DIR DATA_DIR TABLES_DIR FIGURES_DIR",
    call. = FALSE
  )
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot locate the exporter script.", call. = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "03-optimizer-publication-lean.R"), local = TRUE)

result <- jss_optimizer_lean_publish(
  checkpoint_dir = args[[1L]],
  data_dir = args[[2L]],
  tables_dir = args[[3L]],
  figures_dir = args[[4L]]
)
message(
  "Lean optimizer supplementary bundle validated: ",
  length(result$data), " data, ", length(result$tables), " tables, ",
  length(result$figures), " figure."
)
