jss_run_04_missingness_dropout <- function(settings) {
  jss_write_stub_module(
    settings = settings,
    module_id = "04-missingness-dropout-sensitivity",
    title = "Missingness and dropout sensitivity",
    family = "mixed",
    copula = "mixed",
    focus = "sensitivity to intermittent missing visits and dropout-like panels",
    planned_outputs = paste(
      "Sensitivity data summary, bias/convergence table by missingness pattern,",
      "and at least one figure showing robustness across missingness/dropout settings."
    )
  )
}
