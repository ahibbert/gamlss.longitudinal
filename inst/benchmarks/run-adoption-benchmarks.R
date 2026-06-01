# Opt-in adoption benchmark runner for gamlss.longitudinal.
#
# Run from an installed package:
#   Rscript -e "source(system.file('benchmarks', 'run-adoption-benchmarks.R', package = 'gamlss.longitudinal'))"
#
# Environment controls:
#   GAMLSS_LONGITUDINAL_ADOPTION_REPS=100
#   GAMLSS_LONGITUDINAL_ADOPTION_SCENARIOS=gaussian_heteroskedastic,gamma_positive
#   GAMLSS_LONGITUDINAL_ADOPTION_METHODS=rs_separate,gee,glmm,gam
#   GAMLSS_LONGITUDINAL_ADOPTION_METRICS=benchmark_mean_rmse,elapsed_sec
#   GAMLSS_LONGITUDINAL_ADOPTION_MAX_ELAPSED_SEC=60
#   GAMLSS_LONGITUDINAL_ADOPTION_OUTPUT_DIR=results/adoption_benchmarks
#   GAMLSS_LONGITUDINAL_ADOPTION_SEED=20260528

suppressPackageStartupMessages(library(gamlss.longitudinal))

.benchmark_env <- function(name, default = "") {
  value <- Sys.getenv(name, unset = default)
  value <- trimws(value)
  if (!nzchar(value)) default else value
}

.benchmark_env_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(.benchmark_env(name, as.character(default))))
  if (length(value) == 0L || !is.finite(value)) default else value[1L]
}

.benchmark_env_numeric <- function(name, default) {
  value <- suppressWarnings(as.numeric(.benchmark_env(name, as.character(default))))
  if (length(value) == 0L || !is.finite(value)) default else value[1L]
}

.benchmark_env_vector <- function(name, default = NULL) {
  value <- .benchmark_env(name, "")
  if (!nzchar(value)) return(default)
  out <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  out[nzchar(out)]
}

.benchmark_flatten <- function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  for (nm in names(x)) {
    if (is.list(x[[nm]]) && !is.data.frame(x[[nm]])) {
      x[[nm]] <- vapply(x[[nm]], function(value) paste(value, collapse = ";"), character(1))
    }
  }
  x
}

.benchmark_write_csv <- function(x, path) {
  utils::write.csv(.benchmark_flatten(x), path, row.names = FALSE, na = "")
  invisible(path)
}

reps <- .benchmark_env_integer("GAMLSS_LONGITUDINAL_ADOPTION_REPS", 100L)
seed <- .benchmark_env_integer("GAMLSS_LONGITUDINAL_ADOPTION_SEED", 20260528L)
max_elapsed_sec <- .benchmark_env_numeric("GAMLSS_LONGITUDINAL_ADOPTION_MAX_ELAPSED_SEC", 60)
output_dir <- .benchmark_env(
  "GAMLSS_LONGITUDINAL_ADOPTION_OUTPUT_DIR",
  file.path("results", "adoption_benchmarks")
)
scenario_names <- .benchmark_env_vector("GAMLSS_LONGITUDINAL_ADOPTION_SCENARIOS", NULL)
methods <- .benchmark_env_vector(
  "GAMLSS_LONGITUDINAL_ADOPTION_METHODS",
  c("rs_separate", "gee", "glmm", "gam")
)
metrics <- .benchmark_env_vector("GAMLSS_LONGITUDINAL_ADOPTION_METRICS", NULL)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

scenarios <- adoption_benchmark_scenarios(scenario_names)
status <- benchmark_comparator_status()

cat("gamlss.longitudinal adoption benchmark\n")
cat("Output directory:", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n")
cat("Scenarios:", paste(scenarios$scenario, collapse = ", "), "\n")
cat("Methods:", paste(methods, collapse = ", "), "\n")
cat("Replicates per scenario:", reps, "\n")
cat("Max elapsed seconds per fit:", max_elapsed_sec, "\n\n")

cat("Comparator availability:\n")
print(status, row.names = FALSE)
cat("\n")

bench <- run_adoption_benchmarks(
  scenarios = scenarios,
  reps = reps,
  methods = methods,
  metrics = metrics,
  seed = seed,
  max_elapsed_sec = max_elapsed_sec,
  write_results = FALSE,
  output_dir = output_dir
)

results_path <- file.path(output_dir, "adoption_benchmark_results.csv")
summary_path <- file.path(output_dir, "adoption_benchmark_summary.csv")
case_results_path <- file.path(output_dir, "adoption_benchmark_case_results.csv")
scenarios_path <- file.path(output_dir, "adoption_benchmark_scenarios.csv")
status_path <- file.path(output_dir, "adoption_benchmark_comparator_status.csv")
report_path <- file.path(output_dir, "adoption_benchmark_report.md")
object_path <- file.path(output_dir, "adoption_benchmark_object.rds")

.benchmark_write_csv(bench$results, results_path)
.benchmark_write_csv(bench$summary$summary, summary_path)
.benchmark_write_csv(bench$summary$case_results, case_results_path)
.benchmark_write_csv(bench$scenarios, scenarios_path)
.benchmark_write_csv(status, status_path)
write_benchmark_report(bench, path = report_path, comparator_status = status)
saveRDS(bench, object_path)

cat("Benchmark complete.\n")
cat("Results:", normalizePath(results_path, winslash = "/", mustWork = FALSE), "\n")
cat("Summary:", normalizePath(summary_path, winslash = "/", mustWork = FALSE), "\n")
cat("Case results:", normalizePath(case_results_path, winslash = "/", mustWork = FALSE), "\n")
cat("Report:", normalizePath(report_path, winslash = "/", mustWork = FALSE), "\n")
cat("Benchmark object:", normalizePath(object_path, winslash = "/", mustWork = FALSE), "\n\n")

print(bench)
