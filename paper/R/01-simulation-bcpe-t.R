jss_run_script <- function(script, env, root) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)
  old <- Sys.getenv(names(env), unset = NA_character_)
  on.exit({
    for (i in seq_along(env)) {
      if (is.na(old[[i]])) Sys.unsetenv(names(env)[[i]]) else do.call(Sys.setenv, stats::setNames(list(old[[i]]), names(env)[[i]]))
    }
  }, add = TRUE)
  for (i in seq_along(env)) do.call(Sys.setenv, stats::setNames(list(env[[i]]), names(env)[[i]]))
  withCallingHandlers(
    sys.source(script, envir = new.env(parent = globalenv()), chdir = FALSE),
    warning = function(w) {
      expected <- grepl("max_outer_iter|not converged|BiCopHfunc|aes_string|deprecated", conditionMessage(w), ignore.case = TRUE)
      if (expected) invokeRestart("muffleWarning")
    }
  )
}

jss_run_01_bcpe_t <- function(settings) {
  input <- file.path(settings$public_data_dir, "bcpe-t")
  if (identical(settings$profile, "full")) {
    input <- file.path(settings$data_dir, "bcpe-t-full")
    jss_run_script(
      file.path(settings$root, "paper", "scripts", "final-simulations", "bcpe-t", "simulation_bcpe_t_gamlss_comparison.R"),
      c(OUT_DIR = input, N_FITS = "100", N_CORES = as.character(settings$workers), SAVE_FITS = "0", COMPUTE_PREDICTIVE_SCORES = "1", VARIOGRAM_P_VALUES = "0.5,2", MAX_ELAPSED_SEC = "180"),
      settings$root
    )
  }
  jss_run_script(
    file.path(settings$root, "paper", "scripts", "final-simulations", "bcpe-t", "make_bcpe_t_paper_outputs.R"),
    c(BCPE_T_PAPER_RS_JOINT_DIR = input, BCPE_T_PAPER_COMPARISON_DIR = settings$out_dir),
    settings$root
  )
  list(
    module_id = "01-simulation-bcpe-t", status = "regenerated", data = character(),
    tables = file.path(settings$out_dir, c("paper_simulation_bcpe_t_fit_characteristics.tex", "paper_simulation_bcpe_t_fixed_parameter_bias_rmse.tex")),
    figures = file.path(settings$out_dir, c("paper_simulation_bcpe_t_fixed_effect_recovery.png", "paper_simulation_bcpe_t_smooth_effect_recovery.png"))
  )
}
