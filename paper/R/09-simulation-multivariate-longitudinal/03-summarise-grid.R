source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))

mvt_load_package()
run_dir <- mvt_read_run_dir()
mvt_summarise_results(run_dir)
message("Summaries written to: ", run_dir)
