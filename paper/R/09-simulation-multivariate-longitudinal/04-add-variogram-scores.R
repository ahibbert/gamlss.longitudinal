source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))

mvt_load_package()
run_dir <- mvt_read_run_dir()

message(
  "Variogram scores are computed during fitting for methods with an implemented predictive simulator. ",
  "This script validates that the output exists and refreshes variogram summaries."
)

path <- file.path(run_dir, "variogram_scores_by_rep.csv")
if (!file.exists(path)) {
  stop("Missing variogram output: ", path, call. = FALSE)
}

mvt_summarise_results(run_dir)
message("Variogram summaries refreshed in: ", run_dir)
