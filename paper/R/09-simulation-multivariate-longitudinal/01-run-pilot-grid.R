source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))

source_used <- mvt_load_package()
mvt_require_runtime_packages(include_gamcopula = TRUE)

dir.create(mvt_output_root, recursive = TRUE, showWarnings = FALSE)

reps <- mvt_rep_ids(mvt_env_int("GAMLSS_LONGITUDINAL_MVT_REPS", 5L))
time_names <- mvt_env_vector("GAMLSS_LONGITUDINAL_MVT_TIMEPOINTS", c("t5", "t20"))
family_names <- mvt_env_vector("GAMLSS_LONGITUDINAL_MVT_FAMILIES", c("gaussian", "gamma", "binomial"))
dependence_names <- mvt_env_vector(
  "GAMLSS_LONGITUDINAL_MVT_DEPENDENCE",
  c("external_exchangeable", "external_ar1", "native_covariate_dependent_adjacent")
)
checkpoint_every <- mvt_env_int("GAMLSS_LONGITUDINAL_MVT_CHECKPOINT_EVERY", 5L)
seed_base <- mvt_env_int("GAMLSS_LONGITUDINAL_MVT_SEED", 20260818L)

grid <- mvt_expand_grid(
  time_names = time_names,
  family_names = family_names,
  dependence_names = dependence_names,
  reps = reps,
  include_special = FALSE,
  include_appendix = FALSE
)

run_dir <- mvt_env("GAMLSS_LONGITUDINAL_MVT_OUTPUT_DIR", mvt_default_run_dir("pilot"))
message("Package source: ", source_used)
message("Writing pilot run to: ", run_dir)

mvt_run_grid(
  grid = grid,
  run_dir = run_dir,
  seed_base = seed_base,
  checkpoint_every = checkpoint_every,
  require_gamcopula = TRUE
)

mvt_summarise_results(run_dir)
message("Pilot run complete: ", run_dir)
