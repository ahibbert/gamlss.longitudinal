param(
  [int]$Reps = 100,
  [string]$OutDir = "results\jss-exploratory\01-continuous-bcpe-t\bcpe_t_current_defaults_rep100_rs_joint_p05_p2_fits",
  [string]$ComparisonDir = "results\jss-exploratory\01-continuous-bcpe-t\bcpe_t_current_defaults_rep100_comparison_p05_p2_fits",
  [int]$PredictiveNSim = 200,
  [int]$Cores = 3,
  [string]$VariogramPValues = "0.5|2"
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path (Join-Path $PSScriptRoot "..\..\..\..")).Path
Set-Location $repo

$resolvedOutDir = Join-Path $repo $OutDir
$resolvedComparisonDir = Join-Path $repo $ComparisonDir
New-Item -ItemType Directory -Force -Path $resolvedOutDir | Out-Null
New-Item -ItemType Directory -Force -Path $resolvedComparisonDir | Out-Null

$runner = Join-Path $PSScriptRoot "simulation_bcpe_t_gamlss_comparison.R"
$paper = Join-Path $PSScriptRoot "make_bcpe_t_paper_outputs.R"
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$logFile = Join-Path $resolvedOutDir ("run_reps{0}_{1}.log" -f $Reps, $stamp)

$env:OUT_DIR = $resolvedOutDir
$env:N_FITS = [string]$Reps
$env:REP_IDS = (1..$Reps) -join ","
$env:N_CORES = [string]$Cores
$env:COMPUTE_SE = "1"
$env:SAVE_FITS = "1"
$env:COMPUTE_PREDICTIVE_SCORES = "1"
$env:PREDICTIVE_NSIM = [string]$PredictiveNSim
$env:VARIOGRAM_P_VALUES = $VariogramPValues
$env:OPT_METHOD = "RS"
$env:INCLUDE_DLCOPDPAR = "1"
$env:WARM_START_JOINT = "1"
$env:VERBOSE_FITS = "0"

$env:BCPE_T_PAPER_RS_JOINT_DIR = $resolvedOutDir
$env:BCPE_T_PAPER_COMPARISON_DIR = $resolvedComparisonDir

"[$(Get-Date -Format o)] Starting BCPE/t RS-joint SE/diagnostics run: reps=$Reps cores=$Cores p=$VariogramPValues out=$resolvedOutDir" |
  Tee-Object -FilePath $logFile

$ErrorActionPreference = "Continue"
& Rscript $runner 2>&1 | Tee-Object -FilePath $logFile -Append
$runnerExitCode = $LASTEXITCODE
$ErrorActionPreference = "Stop"
if ($runnerExitCode -ne 0) {
  throw "BCPE/t runner failed with exit code $runnerExitCode"
}

"[$(Get-Date -Format o)] Generating BCPE/t paper outputs from $resolvedOutDir" |
  Tee-Object -FilePath $logFile -Append

$ErrorActionPreference = "Continue"
& Rscript $paper 2>&1 | Tee-Object -FilePath $logFile -Append
$paperExitCode = $LASTEXITCODE
$ErrorActionPreference = "Stop"
if ($paperExitCode -ne 0) {
  throw "BCPE/t paper output generation failed with exit code $paperExitCode"
}

"[$(Get-Date -Format o)] Completed BCPE/t SE/diagnostics run" |
  Tee-Object -FilePath $logFile -Append
