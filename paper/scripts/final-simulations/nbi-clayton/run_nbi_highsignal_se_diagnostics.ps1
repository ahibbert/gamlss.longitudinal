param(
  [int]$Reps = 100,
  [switch]$Resume,
  [string]$OutDir = "results\jss-exploratory\02-discrete-delaporte-clayton\nbi_highsignal_n500_rep100_se_diagnostics",
  [int]$PredictiveNSim = 100,
  [string]$VariogramPValues = "0.5|2",
  [switch]$SaveFits
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path (Join-Path $PSScriptRoot "..\..\..\..")).Path
Set-Location $repo

$resolvedOutDir = Join-Path $repo $OutDir
New-Item -ItemType Directory -Force -Path $resolvedOutDir | Out-Null

$runner = Join-Path $PSScriptRoot "compare_gamlss_ours_nbi_sigma_smooth.R"
$paper = Join-Path $PSScriptRoot "make_nbi_paper_outputs.R"
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$logFile = Join-Path $resolvedOutDir ("run_reps{0}_{1}.log" -f $Reps, $stamp)

$env:NBI_COMPARE_OUT_DIR = $resolvedOutDir
$env:NBI_COMPARE_REPS = [string]$Reps
$env:NBI_COMPARE_RESUME = if ($Resume) { "TRUE" } else { "FALSE" }
$env:NBI_COMPARE_ENGINES = "gamlss|ours_rs_joint"
$env:NBI_COMPARE_COMPUTE_SE = "TRUE"
$env:NBI_COMPARE_VCOV_METHOD = "analytical"
$env:NBI_COMPARE_COMPUTE_PREDICTIVE = "TRUE"
$env:NBI_COMPARE_PREDICTIVE_NSIM = [string]$PredictiveNSim
$env:NBI_COMPARE_VARIOGRAM_P_VALUES = $VariogramPValues
$env:NBI_COMPARE_SAVE_FITS = if ($SaveFits) { "TRUE" } else { "FALSE" }
$env:NBI_PAPER_COMPARISON_DIR = $resolvedOutDir

"[$(Get-Date -Format o)] Starting NBI SE/diagnostics run: reps=$Reps resume=$($Resume.IsPresent) p=$VariogramPValues save_fits=$($SaveFits.IsPresent) out=$resolvedOutDir" |
  Tee-Object -FilePath $logFile

$ErrorActionPreference = "Continue"
& Rscript $runner 2>&1 | Tee-Object -FilePath $logFile -Append
$runnerExitCode = $LASTEXITCODE
$ErrorActionPreference = "Stop"
if ($runnerExitCode -ne 0) {
  throw "NBI comparison runner failed with exit code $runnerExitCode"
}

"[$(Get-Date -Format o)] Generating paper outputs from $resolvedOutDir" |
  Tee-Object -FilePath $logFile -Append

$ErrorActionPreference = "Continue"
& Rscript $paper 2>&1 | Tee-Object -FilePath $logFile -Append
$paperExitCode = $LASTEXITCODE
$ErrorActionPreference = "Stop"
if ($paperExitCode -ne 0) {
  throw "NBI paper output generation failed with exit code $paperExitCode"
}

"[$(Get-Date -Format o)] Completed NBI SE/diagnostics run" |
  Tee-Object -FilePath $logFile -Append
