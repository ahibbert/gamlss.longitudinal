jss_run_05_lipid_application <- function(settings) {
  jss_write_stub_module(
    settings = settings,
    module_id = "05-application-lipid",
    title = "LIPID clinical-trial application",
    family = "continuous",
    copula = "to be selected",
    focus = "secure continuous-response clinical-trial application",
    planned_outputs = paste(
      "Application data audit, fitted model summary table, diagnostic table,",
      "and at least one figure for fitted clinical-trial response trajectories."
    ),
    external_envvar = "GAMLSS_LONGITUDINAL_LIPID_DATA",
    external_data_note = "Future secure-environment input; data are not expected in this repository."
  )
}
