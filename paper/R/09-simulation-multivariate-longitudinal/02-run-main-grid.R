source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))

source_used <- mvt_load_package()
mvt_require_runtime_packages(include_gamcopula = TRUE)

dir.create(mvt_output_root, recursive = TRUE, showWarnings = FALSE)

reps <- mvt_rep_ids(mvt_env_int("GAMLSS_LONGITUDINAL_MVT_REPS", 100L))
main_scope <- mvt_main_scope()
time_names <- mvt_env_vector("GAMLSS_LONGITUDINAL_MVT_TIMEPOINTS", mvt_default_main_time_names(main_scope))
family_names <- mvt_env_vector("GAMLSS_LONGITUDINAL_MVT_FAMILIES", c("gaussian", "poisson", "gamma", "binomial"))
dependence_names <- mvt_env_vector(
  "GAMLSS_LONGITUDINAL_MVT_DEPENDENCE",
  c(
    "external_exchangeable",
    "external_ar1",
    "native_time_varying_adjacent",
    "native_covariate_dependent_adjacent"
  )
)
include_special <- mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_INCLUDE_SPECIAL", FALSE)
if (isTRUE(include_special) && !"gg_continuous" %in% family_names) {
  family_names <- c(family_names, "gg_continuous")
}
checkpoint_every <- mvt_env_int("GAMLSS_LONGITUDINAL_MVT_CHECKPOINT_EVERY", 10L)
seed_base <- mvt_env_int("GAMLSS_LONGITUDINAL_MVT_SEED", 20260818L)

grid <- mvt_expand_grid(
  time_names = time_names,
  family_names = family_names,
  dependence_names = dependence_names,
  reps = reps,
  include_special = include_special,
  include_appendix = FALSE
)

run_dir <- mvt_env("GAMLSS_LONGITUDINAL_MVT_OUTPUT_DIR", mvt_default_run_dir("main"))
message("Package source: ", source_used)
message("Main scope: ", main_scope)
message("Writing main run to: ", run_dir)

mvt_run_grid(
  grid = grid,
  run_dir = run_dir,
  seed_base = seed_base,
  checkpoint_every = checkpoint_every,
  require_gamcopula = TRUE
)

mvt_summarise_results(run_dir)
message("Main run complete: ", run_dir)
