$ErrorActionPreference = "Stop"

$RepoRoot = "C:\Users\Aydin\OneDrive - The University of Sydney (Students)\gamlss.longitudinal"
$RunDir = "results/jss-exploratory/09-simulation-multivariate-longitudinal/main_t20_reps100_simplified_core_safe_vario"
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
$env:GAMLSS_LONGITUDINAL_MVT_WORKERS = "1"
$env:GAMLSS_LONGITUDINAL_MVT_CHECKPOINT_EVERY = "1"
$env:GAMLSS_LONGITUDINAL_MVT_PRIMARY_MAX_ELAPSED_SEC = "7200"
$env:GAMLSS_LONGITUDINAL_MVT_GAMLSS_TIMEOUT_SEC = "420"
$env:GAMLSS_LONGITUDINAL_MVT_GAMCOPULA_MARKOV_TIMEOUT_SEC = "7200"
$env:GAMLSS_LONGITUDINAL_MVT_GAMCOPULA_VINE_TIMEOUT_SEC = "7200"
$env:GAMLSS_LONGITUDINAL_MVT_GEE_TIMEOUT_SEC = "300"
$env:GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_TIMEOUT_SEC = "180"
$env:GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM = "20"

function Invoke-RScriptLogged {
  param([Parameter(Mandatory = $true)][string]$ScriptPath)
  $PreviousPreference = $ErrorActionPreference
  try {
    # Windows PowerShell promotes native stderr (including harmless R package
    # startup messages) to ErrorRecord objects under Stop. Allow the process to
    # finish and gate on its actual exit code instead.
    $ErrorActionPreference = "Continue"
    & Rscript $ScriptPath 2>&1 | Tee-Object -FilePath $LogPath -Append
    $ExitCode = $LASTEXITCODE
  } finally {
    $ErrorActionPreference = $PreviousPreference
  }
  if ($ExitCode -ne 0) {
    throw "Rscript failed with exit code $ExitCode`: $ScriptPath"
  }
}

try {
  "Started at $(Get-Date -Format o)" | Tee-Object -FilePath $LogPath
  Invoke-RScriptLogged "paper/R/09-simulation-multivariate-longitudinal/02-run-main-grid.R"

  $env:GAMLSS_LONGITUDINAL_MVT_RUN_DIR = $RunDir
  Invoke-RScriptLogged "paper/R/09-simulation-multivariate-longitudinal/03-summarise-grid.R"
  Invoke-RScriptLogged "paper/R/09-simulation-multivariate-longitudinal/05-write-paper-tables.R"
  Invoke-RScriptLogged "paper/R/09-simulation-multivariate-longitudinal/06-make-diagnostic-plots.R"
  Invoke-RScriptLogged "paper/R/09-simulation-multivariate-longitudinal/07-review-audit.R"

  "Finished at $(Get-Date -Format o)" | Tee-Object -FilePath $LogPath -Append
} catch {
  "Failed at $(Get-Date -Format o)" | Tee-Object -FilePath $LogPath -Append
  $_ | Out-String | Tee-Object -FilePath $LogPath -Append
  exit 1
}
