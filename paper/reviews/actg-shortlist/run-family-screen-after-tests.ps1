$ErrorActionPreference = 'Stop'

$ProjectRoot = 'C:\Users\Aydin\OneDrive - The University of Sydney (Students)\gamlss.longitudinal'
$ResultsDir = Join-Path $ProjectRoot 'paper\reviews\actg-shortlist\final-family-copula-screen'
$ScreenScript = Join-Path $ProjectRoot 'paper\reviews\actg-shortlist\run-final-family-copula-screen.R'
$SupervisorLog = Join-Path $ResultsDir 'supervisor.log'
$ScreenLog = Join-Path $ResultsDir 'screen.log'
$ExitStatus = Join-Path $ResultsDir 'exit-status.txt'

New-Item -ItemType Directory -Path $ResultsDir -Force | Out-Null
Set-Location -LiteralPath $ProjectRoot

function Write-SupervisorLog([string]$Message) {
    $line = '{0:o} {1}' -f (Get-Date), $Message
    Add-Content -LiteralPath $SupervisorLog -Value $line
}

Write-SupervisorLog 'Supervisor started; waiting for package test processes to finish.'
while ($true) {
    $tests = @(Get-CimInstance Win32_Process | Where-Object {
        $_.Name -eq 'Rscript.exe' -and $_.CommandLine -match 'testthat::test'
    })
    if ($tests.Count -eq 0) { break }
    Write-SupervisorLog ("Waiting; {0} package test process(es) remain." -f $tests.Count)
    Start-Sleep -Seconds 30
}

$Rscript = (Get-Command Rscript.exe -ErrorAction Stop).Source
Write-SupervisorLog 'Package tests cleared; starting checkpointed family/copula screen.'
& $Rscript $ScreenScript '--workers=3' *>> $ScreenLog
$code = $LASTEXITCODE
Set-Content -LiteralPath $ExitStatus -Value $code
Write-SupervisorLog ("Family/copula screen exited with code {0}." -f $code)
exit $code
