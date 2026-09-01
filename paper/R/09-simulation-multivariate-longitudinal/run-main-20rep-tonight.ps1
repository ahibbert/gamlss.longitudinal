$ErrorActionPreference = "Stop"

$RepoRoot = "C:\Users\Aydin\OneDrive - The University of Sydney (Students)\gamlss.longitudinal"
$RunDir = "results/jss-exploratory/09-simulation-multivariate-longitudinal/main_t20_reps100_simplified_core"
$AbsRunDir = Join-Path $RepoRoot $RunDir
$LogDir = Join-Path $AbsRunDir "scheduled_logs"

New-Item -ItemType Directory -Path $LogDir -Force | Out-Null
$LogPath = Join-Path $LogDir ("main_t20_reps020_" + (Get-Date -Format "yyyyMMdd_HHmmss") + ".log")

Set-Location $RepoRoot

$env:GAMLSS_LONGITUDINAL_MVT_SOURCE = "local"
$env:GAMLSS_LONGITUDINAL_MVT_OUTPUT_DIR = $RunDir
$env:GAMLSS_LONGITUDINAL_MVT_REPS = "20"
$env:GAMLSS_LONGITUDINAL_MVT_TIMEPOINTS = "t20"
$env:GAMLSS_LONGITUDINAL_MVT_FAMILIES = "gaussian,poisson,gamma,binomial"
$env:GAMLSS_LONGITUDINAL_MVT_DEPENDENCE = "external_exchangeable,external_ar1,native_time_varying_adjacent,native_covariate_dependent_adjacent"
$env:GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm,glmm,gee_independence,gee_exchangeable,gee_ar1,gamlss.longitudinal,gamCopula_markov,gamCopula_vine_simplified"
$env:GAMLSS_LONGITUDINAL_MVT_RESUME = "true"
$env:GAMLSS_LONGITUDINAL_MVT_CHECKPOINT_EVERY = "1"
$env:GAMLSS_LONGITUDINAL_MVT_PRIMARY_MAX_ELAPSED_SEC = "7200"
$env:GAMLSS_LONGITUDINAL_MVT_GAMCOPULA_VINE_TIMEOUT_SEC = "7200"
$env:GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM = "20"

try {
  "Started at $(Get-Date -Format o)" | Tee-Object -FilePath $LogPath
  Rscript "paper/R/09-simulation-multivariate-longitudinal/02-run-main-grid.R" 2>&1 | Tee-Object -FilePath $LogPath -Append

  $env:GAMLSS_LONGITUDINAL_MVT_RUN_DIR = $RunDir
  Rscript "paper/R/09-simulation-multivariate-longitudinal/03-summarise-grid.R" 2>&1 | Tee-Object -FilePath $LogPath -Append
  Rscript "paper/R/09-simulation-multivariate-longitudinal/05-write-paper-tables.R" 2>&1 | Tee-Object -FilePath $LogPath -Append
  Rscript "paper/R/09-simulation-multivariate-longitudinal/06-make-diagnostic-plots.R" 2>&1 | Tee-Object -FilePath $LogPath -Append
  Rscript "paper/R/09-simulation-multivariate-longitudinal/07-review-audit.R" 2>&1 | Tee-Object -FilePath $LogPath -Append

  "Finished at $(Get-Date -Format o)" | Tee-Object -FilePath $LogPath -Append
} catch {
  "Failed at $(Get-Date -Format o)" | Tee-Object -FilePath $LogPath -Append
  $_ | Out-String | Tee-Object -FilePath $LogPath -Append
  exit 1
}
