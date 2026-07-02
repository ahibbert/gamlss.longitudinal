jss_run_01_bcpe_t <- function(settings) {
  jss_write_stub_module(
    settings = settings,
    module_id = "01-simulation-bcpe-t",
    title = "BCPE margin with t-copula simulation",
    family = "BCPE",
    copula = "t",
    focus = "continuous simulation, parameter recovery, fitted distribution diagnostics",
    planned_outputs = paste(
      "Simulation data summary, parameter recovery table, convergence diagnostics,",
      "and at least one fitted-versus-truth figure for the continuous JSS simulation."
    )
  )
}
