jss_run_02_delaporte_clayton <- function(settings) {
  jss_write_stub_module(
    settings = settings,
    module_id = "02-simulation-delaporte-clayton",
    title = "Delaporte margin with Clayton copula simulation",
    family = "DEL",
    copula = "C",
    focus = "discrete simulation, count response recovery, Clayton dependence recovery",
    planned_outputs = paste(
      "Simulation data summary, parameter recovery table, count calibration diagnostics,",
      "and at least one recovery figure for the discrete JSS simulation."
    )
  )
}
