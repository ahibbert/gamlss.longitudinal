jss_run_02_delaporte_clayton <- function(settings) {
  input <- file.path(settings$public_data_dir, "nbi-clayton")
  if (identical(settings$profile, "full")) {
    input <- file.path(settings$data_dir, "nbi-clayton-full")
    jss_run_script(
      file.path(settings$root, "paper", "scripts", "final-simulations", "nbi-clayton", "compare_gamlss_ours_nbi_sigma_smooth.R"),
      c(NBI_COMPARE_OUT_DIR = input, NBI_COMPARE_REPS = "100", NBI_COMPARE_RESUME = "TRUE", NBI_COMPARE_SAVE_FITS = "FALSE", NBI_COMPARE_ENGINES = "gamlss|ours_rs_joint", NBI_COMPARE_VARIOGRAM_P_VALUES = "0.5,2"),
      settings$root
    )
  }
  jss_run_script(
    file.path(settings$root, "paper", "scripts", "final-simulations", "nbi-clayton", "make_nbi_paper_outputs.R"),
    c(NBI_PAPER_INPUT_DIR = input, NBI_PAPER_OUTPUT_DIR = settings$out_dir),
    settings$root
  )
  list(
    module_id = "02-simulation-delaporte-clayton", status = "regenerated", data = character(),
    tables = file.path(settings$out_dir, c("paper_simulation_nbi_clayton_highsignal_fit_characteristics.tex", "paper_simulation_nbi_clayton_highsignal_fixed_parameter_bias_rmse.tex")),
    figures = file.path(settings$out_dir, c("paper_simulation_nbi_clayton_highsignal_fixed_effect_recovery.png", "paper_simulation_nbi_clayton_highsignal_smooth_effect_recovery.png"))
  )
}
