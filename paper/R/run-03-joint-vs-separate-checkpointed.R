# Standalone entry point for the registered optimizer benchmark. Checkpoint and
# resume behavior lives in 03-joint-vs-separate-optimization.R so the targets
# module and this direct production runner share exactly the same contract.

Sys.setenv(
  GAMLSS_LONGITUDINAL_JSS_PROFILE = Sys.getenv(
    "GAMLSS_LONGITUDINAL_JSS_PROFILE",
    unset = "full"
  ),
  GAMLSS_LONGITUDINAL_JSS_JVS_RESUME = Sys.getenv(
    "GAMLSS_LONGITUDINAL_JSS_JVS_RESUME",
    unset = "1"
  )
)

suppressPackageStartupMessages(library(gamlss.longitudinal))

source(file.path("paper", "R", "replication-helpers.R"))
source(file.path("paper", "R", "03-joint-vs-separate-optimization.R"))

settings <- jss_settings()
result <- jss_run_03_joint_vs_separate(settings)
cat(
  "Completed registered paired optimizer benchmark. Resume checkpoints: ",
  jss_joint_checkpoint_dir(settings),
  "\n",
  sep = ""
)
invisible(result)
