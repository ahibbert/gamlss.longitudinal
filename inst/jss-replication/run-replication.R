#!/usr/bin/env Rscript

set.seed(20260528)

out_dir <- file.path("results", "jss-replication")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

methods <- c("gamlss", "rs_separate", "rs_joint", "cg")

coverage_results <- gamlss.longitudinal::run_coverage_simulations(
  families = c("NO", "GA", "PO", "NBI"),
  copulas = c("N", "C"),
  methods = methods,
  designs = c("intercept", "covariate"),
  n = 60L,
  times = 1:3,
  seed = 20260528,
  max_outer_iter = 8L,
  max_inner_iter = 8L,
  max_elapsed_sec = 30,
  output_dir = out_dir,
  write_results = TRUE
)

utils::write.csv(
  coverage_results,
  file.path(out_dir, "coverage_results_latest.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(utils::sessionInfo()),
  file.path(out_dir, "session_info.txt")
)

message("Wrote JSS replication outputs to: ", normalizePath(out_dir))
