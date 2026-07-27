# Standard Model Benchmarking

This directory contains the module 08 runner for the JSS correlation
misspecification sensitivity analysis. It regenerates the standard-model
benchmark source run used by
`08-simulation-sensitivity-correlation-misspecification.R`, including the
subject-cluster sandwich estimator tables.

Outputs from new runs are written to:

```text
results/jss-exploratory/08-simulation-sensitivity-correlation-misspecification/standard-model-benchmarking/
```

The archived run used for the current paper-facing tables remains under
`results/jss-exploratory/05-standard-model-benchmarking/...` for provenance.
Set `GAMLSS_LONGITUDINAL_JSS_CORR_MISSPEC_RUN_DIR` to point module 08 at a
different completed run.

## Scripts

- `00-benchmark-setup.R`: shared simulation, fitting, package-source, GEE, and
  summary helpers.
- `02-run-rs-joint-standard-model-grid.R`: repeated simulation runner for the
  base `gamlss.longitudinal`, `glm`, and `geepack` comparisons.
- `03-summarise-rs-joint-standard-model-grid.R`: benchmark, coefficient,
  dependence, complexity, and status summaries.
- `05-add-gee-unstructured-rows.R`: resume-safe helper for adding `geepack`
  unstructured working-correlation rows to an existing combined run.
- `08-recompute-all-pair-dependence.R`: recomputes dependence recovery over all
  within-subject pairs rather than adjacent pairs only.
- `11-run-t20-sandwich-grid.R`: refits the T=20 grid and replaces coefficient
  uncertainty with subject-cluster sandwich standard errors.
- `07-write-story-tables.R`: writes the paper-facing scenario, main T=20, and
  appendix T=50 table fragments.

## Stable Package Source

By default, the scripts use the newest root package tarball when present, then
fall back to the installed package. Override with:

```powershell
$env:GAMLSS_LONGITUDINAL_BENCHMARK_SOURCE = "installed" # tarball, installed, github, or local
```

For a pinned GitHub source:

```powershell
$env:GAMLSS_LONGITUDINAL_BENCHMARK_SOURCE = "github"
$env:GAMLSS_LONGITUDINAL_BENCHMARK_GITHUB_REPO = "ahibbert/gamlss.longitudinal"
$env:GAMLSS_LONGITUDINAL_BENCHMARK_GITHUB_REF = "<commit-sha>"
```

## Regenerate A Source Run

Small pilot:

```powershell
$env:GAMLSS_LONGITUDINAL_BENCHMARK_SOURCE = "installed"
$env:GAMLSS_LONGITUDINAL_BENCHMARK_FAMILIES = "gaussian,gamma,binary"
$env:GAMLSS_LONGITUDINAL_BENCHMARK_SCENARIOS = "external_ar1_moderate"
$env:GAMLSS_LONGITUDINAL_BENCHMARK_TIMEPOINTS = "20"
$env:GAMLSS_LONGITUDINAL_BENCHMARK_REPS = "2"
$env:GAMLSS_LONGITUDINAL_BENCHMARK_N = "20"
Rscript "paper/R/08-simulation-sensitivity-correlation-misspecification/standard-model-benchmarking/02-run-rs-joint-standard-model-grid.R"
Rscript "paper/R/08-simulation-sensitivity-correlation-misspecification/standard-model-benchmarking/03-summarise-rs-joint-standard-model-grid.R"
```

Paper-facing T=20/T=50 grid:

```powershell
$env:GAMLSS_LONGITUDINAL_BENCHMARK_SOURCE = "installed"
$env:GAMLSS_LONGITUDINAL_BENCHMARK_FAMILIES = "gaussian,gamma,binary"
$env:GAMLSS_LONGITUDINAL_BENCHMARK_TIMEPOINTS = "20,50"
$env:GAMLSS_LONGITUDINAL_BENCHMARK_REPS = "20"
$env:GAMLSS_LONGITUDINAL_BENCHMARK_N = "120"
$env:GAMLSS_LONGITUDINAL_BENCHMARK_PRIMARY_TIMEOUT_SEC = "300"
$env:GAMLSS_LONGITUDINAL_BENCHMARK_GEE_TIMEOUT_SEC = "30"
$env:GAMLSS_LONGITUDINAL_BENCHMARK_GEE_UNSTRUCTURED_TIMEOUT_SEC = "30"
Rscript "paper/R/08-simulation-sensitivity-correlation-misspecification/standard-model-benchmarking/02-run-rs-joint-standard-model-grid.R"
Rscript "paper/R/08-simulation-sensitivity-correlation-misspecification/standard-model-benchmarking/03-summarise-rs-joint-standard-model-grid.R"
```

The scripts write `latest_run_dir.txt` in the output root. Use that directory
as `GAMLSS_LONGITUDINAL_BENCHMARK_REPORT_DIR` for the follow-up steps.

## Add Sandwich And Paper Tables

```powershell
$env:GAMLSS_LONGITUDINAL_BENCHMARK_REPORT_DIR = "results/jss-exploratory/08-simulation-sensitivity-correlation-misspecification/standard-model-benchmarking/run_YYYYMMDD_HHMMSS"
Rscript "paper/R/08-simulation-sensitivity-correlation-misspecification/standard-model-benchmarking/08-recompute-all-pair-dependence.R"
Rscript "paper/R/08-simulation-sensitivity-correlation-misspecification/standard-model-benchmarking/11-run-t20-sandwich-grid.R"
Rscript "paper/R/08-simulation-sensitivity-correlation-misspecification/standard-model-benchmarking/07-write-story-tables.R"
```

If the source run does not already include unstructured GEE rows, add them
before the all-pair and sandwich steps:

```powershell
Rscript "paper/R/08-simulation-sensitivity-correlation-misspecification/standard-model-benchmarking/05-add-gee-unstructured-rows.R"
```

By default, the story tables include Gaussian, Gamma, and binary margins. Poisson
is intentionally excluded from the current paper-facing table set.

## Main Outputs

The source run directory contains per-replicate CSVs, summaries, status files,
and diagnostic figures. The final paper-facing artifacts are written under:

```text
<run_dir>/sandwich_t20_grid/
<run_dir>/story_tables/
```

Module 08 then copies the selected table fragments and validation summaries into
`results/jss-replication/<profile>/`.
