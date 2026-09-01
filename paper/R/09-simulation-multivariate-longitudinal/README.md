# Multivariate Longitudinal Simulation Study

This module extends the JSS workflow to multivariate longitudinal settings with
`T = 5, 20, 50`, total-row controlled sample sizes, standard longitudinal
comparators, `gamCopula`, and `gamlss.longitudinal`.

Outputs are written to:

```text
results/jss-exploratory/09-simulation-multivariate-longitudinal/
```

The module is not run by CRAN checks. It requires `gamCopula` at runtime because
the two-stage conditional copula comparator is part of the planned study.

## Dependencies

Install optional paper-workflow packages before running:

```r
install.packages(c(
  "gamCopula", "VineCopula", "mvtnorm", "geepack", "lme4",
  "callr", "scoringRules"
))
```

Use local package source when developing:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SOURCE = "local"
```

Use installed package source for a stable replication run:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SOURCE = "installed"
```

Runs resume by default when the target output directory already contains
`*_by_rep.csv` or `*_checkpoint.csv` files. Disable this when intentionally
rerunning a directory:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_RESUME = "false"
```

Limit methods for smoke or debugging runs with `GAMLSS_LONGITUDINAL_MVT_COMPARATORS`.
Allowed values are `glm`, `glmm`, `gee_independence`, `gee_exchangeable`,
`gee_ar1`, `gee_unstructured`, `gamlss.longitudinal`, `gamCopula_markov`,
`gamCopula_vine_simplified`, `gamCopula_vine`, and the appendix-only
`glmm_slope` sensitivity. The legacy shortcut `gamCopula` maps to
`gamCopula_markov`, while `gee` expands to all GEE working correlations.

The main workflow uses `gamCopula_vine_simplified` to keep routine runs
tractable. Use `gamCopula_vine` for targeted full-vine sensitivity runs,
especially the covariate-dependent adjacent dependence scenario.

## Tiny Smoke Run

This uses one small case and writes all expected output files if optional
dependencies are installed.

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SOURCE = "local"
$env:GAMLSS_LONGITUDINAL_MVT_TIMEPOINTS = "t5"
$env:GAMLSS_LONGITUDINAL_MVT_FAMILIES = "gaussian"
$env:GAMLSS_LONGITUDINAL_MVT_DEPENDENCE = "external_exchangeable"
$env:GAMLSS_LONGITUDINAL_MVT_REPS = "1"
$env:GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = "10"
$env:GAMLSS_LONGITUDINAL_MVT_PRIMARY_MAX_ELAPSED_SEC = "30"
$env:GAMLSS_LONGITUDINAL_MVT_GEE_UNSTRUCTURED_TIMEOUT_SEC = "10"
$env:GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM = "2"
$env:GAMLSS_LONGITUDINAL_MVT_RESUME = "false"
Rscript "paper/R/09-simulation-multivariate-longitudinal/01-run-pilot-grid.R"
```

## Review Smoke Run

This exercises all four main margins and both generator tracks on a very small
dataset, then writes summaries, tables, figures, and the review audit. It is the
recommended pre-pilot check after editing the workflow.

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SOURCE = "local"
$env:GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = "8"
$env:GAMLSS_LONGITUDINAL_MVT_REPS = "1"
$env:GAMLSS_LONGITUDINAL_MVT_TIMEPOINTS = "t5"
$env:GAMLSS_LONGITUDINAL_MVT_FAMILIES = "gaussian,poisson,gamma,binomial"
$env:GAMLSS_LONGITUDINAL_MVT_DEPENDENCE = "external_exchangeable,native_covariate_dependent_adjacent"
$env:GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm,gamCopula_markov,gamCopula_vine_simplified,gamlss.longitudinal"
$env:GAMLSS_LONGITUDINAL_MVT_GEE_UNSTRUCTURED_TIMEOUT_SEC = "5"
$env:GAMLSS_LONGITUDINAL_MVT_PRIMARY_MAX_ELAPSED_SEC = "45"
$env:GAMLSS_LONGITUDINAL_MVT_MAX_OUTER_ITER = "8"
$env:GAMLSS_LONGITUDINAL_MVT_MAX_INNER_ITER = "8"
$env:GAMLSS_LONGITUDINAL_MVT_COMPUTE_VCOV = "true"
$env:GAMLSS_LONGITUDINAL_MVT_SUMMARY_VCOV = "true"
$env:GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM = "1"
$env:GAMLSS_LONGITUDINAL_MVT_RESUME = "false"
$env:GAMLSS_LONGITUDINAL_MVT_OUTPUT_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/review_smoke_allfamilies"
Rscript "paper/R/09-simulation-multivariate-longitudinal/01-run-pilot-grid.R"
$env:GAMLSS_LONGITUDINAL_MVT_RUN_DIR = $env:GAMLSS_LONGITUDINAL_MVT_OUTPUT_DIR
Rscript "paper/R/09-simulation-multivariate-longitudinal/05-write-paper-tables.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/06-make-diagnostic-plots.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/07-review-audit.R"
```

## Full Method Stack Smoke Run

This fits every main comparator on one tiny Gaussian scenario. It is not
publication evidence, but it checks the full comparator stack before longer
runs.

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SOURCE = "local"
$env:GAMLSS_LONGITUDINAL_MVT_OUTPUT_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/full_method_stack_smoke"
$env:GAMLSS_LONGITUDINAL_MVT_TIMEPOINTS = "t5"
$env:GAMLSS_LONGITUDINAL_MVT_FAMILIES = "gaussian"
$env:GAMLSS_LONGITUDINAL_MVT_DEPENDENCE = "external_exchangeable"
$env:GAMLSS_LONGITUDINAL_MVT_REPS = "1"
$env:GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = "8"
$env:GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm,glmm,gee_independence,gee_exchangeable,gee_ar1,gee_unstructured,gamlss.longitudinal,gamCopula_markov,gamCopula_vine_simplified,gamCopula_vine"
$env:GAMLSS_LONGITUDINAL_MVT_GEE_UNSTRUCTURED_TIMEOUT_SEC = "10"
$env:GAMLSS_LONGITUDINAL_MVT_GEE_TIMEOUT_SEC = "10"
$env:GAMLSS_LONGITUDINAL_MVT_PRIMARY_MAX_ELAPSED_SEC = "45"
$env:GAMLSS_LONGITUDINAL_MVT_MAX_OUTER_ITER = "8"
$env:GAMLSS_LONGITUDINAL_MVT_MAX_INNER_ITER = "8"
$env:GAMLSS_LONGITUDINAL_MVT_COMPUTE_VCOV = "true"
$env:GAMLSS_LONGITUDINAL_MVT_SUMMARY_VCOV = "true"
$env:GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM = "2"
$env:GAMLSS_LONGITUDINAL_MVT_RESUME = "false"
Rscript "paper/R/09-simulation-multivariate-longitudinal/01-run-pilot-grid.R"
$env:GAMLSS_LONGITUDINAL_MVT_RUN_DIR = $env:GAMLSS_LONGITUDINAL_MVT_OUTPUT_DIR
Rscript "paper/R/09-simulation-multivariate-longitudinal/05-write-paper-tables.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/06-make-diagnostic-plots.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/07-review-audit.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/09-make-review-bundle.R"
```

## Pilot Grid

## Real-Size Pilot Shard

This uses the planned pilot row counts and scenario grid but only one replicate.
It is useful for checking real-size runtime, convergence warnings, and timeout
classification before launching the default 5-replicate pilot.

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SOURCE = "local"
$env:GAMLSS_LONGITUDINAL_MVT_OUTPUT_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/pilot_real_size_rep1_shard"
$env:GAMLSS_LONGITUDINAL_MVT_REPS = "1"
$env:GAMLSS_LONGITUDINAL_MVT_TIMEPOINTS = "t5,t20"
$env:GAMLSS_LONGITUDINAL_MVT_FAMILIES = "gaussian,gamma,binomial"
$env:GAMLSS_LONGITUDINAL_MVT_DEPENDENCE = "external_exchangeable,external_ar1,native_covariate_dependent_adjacent"
$env:GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm,glmm,gee_independence,gee_exchangeable,gee_ar1,gee_unstructured,gamlss.longitudinal,gamCopula_markov,gamCopula_vine_simplified"
$env:GAMLSS_LONGITUDINAL_MVT_GEE_UNSTRUCTURED_TIMEOUT_SEC = "30"
$env:GAMLSS_LONGITUDINAL_MVT_GEE_TIMEOUT_SEC = "30"
$env:GAMLSS_LONGITUDINAL_MVT_PRIMARY_MAX_ELAPSED_SEC = "180"
$env:GAMLSS_LONGITUDINAL_MVT_MAX_OUTER_ITER = "30"
$env:GAMLSS_LONGITUDINAL_MVT_MAX_INNER_ITER = "30"
$env:GAMLSS_LONGITUDINAL_MVT_COMPUTE_VCOV = "true"
$env:GAMLSS_LONGITUDINAL_MVT_SUMMARY_VCOV = "true"
$env:GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM = "1"
$env:GAMLSS_LONGITUDINAL_MVT_RESUME = "false"
Rscript "paper/R/09-simulation-multivariate-longitudinal/01-run-pilot-grid.R"
$env:GAMLSS_LONGITUDINAL_MVT_RUN_DIR = $env:GAMLSS_LONGITUDINAL_MVT_OUTPUT_DIR
Rscript "paper/R/09-simulation-multivariate-longitudinal/03-summarise-grid.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/07-review-audit.R"
```

Default pilot:

- `T = 5, 20`
- Gaussian, Gamma, Binomial
- exchangeable, AR(1), covariate-dependent adjacent dependence
- 5 replicates

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SOURCE = "local"
$env:GAMLSS_LONGITUDINAL_MVT_REPS = "5"
Rscript "paper/R/09-simulation-multivariate-longitudinal/01-run-pilot-grid.R"
```

## Main Grid

## Main-Core One-Rep Shard

This exercises the core main grid once at the intended `T = 20` row size. It is
useful for diagnosing main-run feasibility before launching the 100-replicate
main grid.

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SOURCE = "local"
$env:GAMLSS_LONGITUDINAL_MVT_OUTPUT_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/main_core_rep1_shard"
$env:GAMLSS_LONGITUDINAL_MVT_REPS = "1"
$env:GAMLSS_LONGITUDINAL_MVT_TIMEPOINTS = "t20"
$env:GAMLSS_LONGITUDINAL_MVT_FAMILIES = "gaussian,poisson,gamma,binomial"
$env:GAMLSS_LONGITUDINAL_MVT_DEPENDENCE = "external_exchangeable,external_ar1,native_time_varying_adjacent,native_covariate_dependent_adjacent"
$env:GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm,glmm,gee_independence,gee_exchangeable,gee_ar1,gee_unstructured,gamlss.longitudinal,gamCopula_markov,gamCopula_vine_simplified"
$env:GAMLSS_LONGITUDINAL_MVT_GEE_UNSTRUCTURED_TIMEOUT_SEC = "30"
$env:GAMLSS_LONGITUDINAL_MVT_GEE_TIMEOUT_SEC = "30"
$env:GAMLSS_LONGITUDINAL_MVT_PRIMARY_MAX_ELAPSED_SEC = "180"
$env:GAMLSS_LONGITUDINAL_MVT_MAX_OUTER_ITER = "30"
$env:GAMLSS_LONGITUDINAL_MVT_MAX_INNER_ITER = "30"
$env:GAMLSS_LONGITUDINAL_MVT_COMPUTE_VCOV = "true"
$env:GAMLSS_LONGITUDINAL_MVT_SUMMARY_VCOV = "true"
$env:GAMLSS_LONGITUDINAL_MVT_VARIOGRAM_NSIM = "1"
$env:GAMLSS_LONGITUDINAL_MVT_RESUME = "false"
Rscript "paper/R/09-simulation-multivariate-longitudinal/02-run-main-grid.R"
$env:GAMLSS_LONGITUDINAL_MVT_RUN_DIR = $env:GAMLSS_LONGITUDINAL_MVT_OUTPUT_DIR
Rscript "paper/R/09-simulation-multivariate-longitudinal/03-summarise-grid.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/07-review-audit.R"
```

Default main grid:

- `T = 20`
- Gaussian, Poisson, Gamma, Binomial
- exchangeable, AR(1), time-varying adjacent, covariate-dependent adjacent
- 100 replicates

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SOURCE = "local"
$env:GAMLSS_LONGITUDINAL_MVT_REPS = "100"
$env:GAMLSS_LONGITUDINAL_MVT_GEE_UNSTRUCTURED_TIMEOUT_SEC = "30"
Rscript "paper/R/09-simulation-multivariate-longitudinal/02-run-main-grid.R"
```

Run the broader appendix grid over `T = 5, 20, 50`:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_MAIN_SCOPE = "appendix"
Rscript "paper/R/09-simulation-multivariate-longitudinal/02-run-main-grid.R"
```

Include the GAMLSS-only generalized gamma scenario:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_INCLUDE_SPECIAL = "true"
Rscript "paper/R/09-simulation-multivariate-longitudinal/02-run-main-grid.R"
```

Run the richer GLMM appendix sensitivity:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm,glmm,glmm_slope,gamCopula_markov,gamCopula_vine_simplified,gamlss.longitudinal"
$env:GAMLSS_LONGITUDINAL_MVT_MAIN_SCOPE = "appendix"
Rscript "paper/R/09-simulation-multivariate-longitudinal/02-run-main-grid.R"
```

## Resumable Main Run

For a 20-replicate overnight run, use the simplified vine for the full core grid
and checkpoint every case:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SOURCE = "local"
$env:GAMLSS_LONGITUDINAL_MVT_OUTPUT_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/main_t20_reps100_simplified_core"
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
Rscript "paper/R/09-simulation-multivariate-longitudinal/02-run-main-grid.R"
```

To continue the same evidence set from 20 to 100 replicates, keep the same
output directory and change only the replicate target:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_REPS = "100"
$env:GAMLSS_LONGITUDINAL_MVT_RESUME = "true"
Rscript "paper/R/09-simulation-multivariate-longitudinal/02-run-main-grid.R"
```

Add targeted full-vine evidence for the covariate-dependent scenarios into that
same run directory:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_RUN_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/main_t20_reps100_simplified_core"
$env:GAMLSS_LONGITUDINAL_MVT_RERUN_METHODS = "gamCopula_vine"
$env:GAMLSS_LONGITUDINAL_MVT_RERUN_DEPENDENCE = "native_covariate_dependent_adjacent"
$env:GAMLSS_LONGITUDINAL_MVT_GAMCOPULA_VINE_TIMEOUT_SEC = "7200"
Rscript "paper/R/09-simulation-multivariate-longitudinal/15-rerun-selected-methods.R"
```

## Follow-Up Outputs

Use the latest completed run:

```powershell
Rscript "paper/R/09-simulation-multivariate-longitudinal/03-summarise-grid.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/04-add-variogram-scores.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/05-write-paper-tables.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/06-make-diagnostic-plots.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/07-review-audit.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/08-publication-readiness-audit.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/09-make-review-bundle.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/11-write-study-protocol.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/12-write-implementation-status.R"
```

Or target a specific run:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_RUN_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/pilot_YYYYMMDD_HHMMSS"
Rscript "paper/R/09-simulation-multivariate-longitudinal/03-summarise-grid.R"
```

## Merge Review Shards

Use `16-merge-run-directories.R` when a review or publication evidence set is
assembled from separate shard directories. The merge rewrites a single
`scenario_grid.csv`, combines all `*_by_rep.csv` files, records source
directories in `run_metadata.csv`, refreshes summaries, and writes an artifact
manifest.

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_MERGE_SOURCE_DIRS = "results/jss-exploratory/09-simulation-multivariate-longitudinal/round2_t20_completed16,results/jss-exploratory/09-simulation-multivariate-longitudinal/round2_t20_native_reps001"
$env:GAMLSS_LONGITUDINAL_MVT_MERGE_TARGET_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/round2_t20_review_merged"
$env:GAMLSS_LONGITUDINAL_MVT_MERGE_NOTE = "Corrected round-2 T=20 review evidence; not a full 100-replicate production run."
$env:GAMLSS_LONGITUDINAL_MVT_MERGE_OVERWRITE = "true"
Rscript "paper/R/09-simulation-multivariate-longitudinal/16-merge-run-directories.R"

$env:GAMLSS_LONGITUDINAL_MVT_RUN_DIR = $env:GAMLSS_LONGITUDINAL_MVT_MERGE_TARGET_DIR
Rscript "paper/R/09-simulation-multivariate-longitudinal/05-write-paper-tables.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/06-make-diagnostic-plots.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/07-review-audit.R"
Rscript "paper/R/09-simulation-multivariate-longitudinal/09-make-review-bundle.R"
```

## Publication Readiness Gate

The run-level audit accepts smoke and filtered runs. Before treating outputs as
paper-facing evidence, run the stricter publication-readiness audit against the
completed pilot and main directories:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_PILOT_RUN_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/pilot_YYYYMMDD_HHMMSS"
$env:GAMLSS_LONGITUDINAL_MVT_MAIN_RUN_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/main_YYYYMMDD_HHMMSS"
$env:GAMLSS_LONGITUDINAL_MVT_APPENDIX_RUN_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/main_appendix_YYYYMMDD_HHMMSS"
Rscript "paper/R/09-simulation-multivariate-longitudinal/08-publication-readiness-audit.R"
```

Optional sensitivity runs can be included with:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SPECIAL_RUN_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/special_YYYYMMDD_HHMMSS"
$env:GAMLSS_LONGITUDINAL_MVT_GLMM_SENSITIVITY_RUN_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/glmm_sensitivity_YYYYMMDD_HHMMSS"
Rscript "paper/R/09-simulation-multivariate-longitudinal/08-publication-readiness-audit.R"
```

## Publication Suite Driver

The suite driver writes a full run plan by default and does not start the long
simulations until explicitly enabled:

```powershell
Rscript "paper/R/09-simulation-multivariate-longitudinal/10-run-publication-suite.R"
```

The default suite plan includes the 5-rep pilot, 100-rep `T = 20` main grid, and
100-rep appendix grid over `T = 5, 20, 50`, with named run folders under a single
suite directory. The dry-run plan reports estimated case counts and maximum row
counts before any models are fit. Execute the planned suite only when ready for
the full compute:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SOURCE = "installed"
$env:GAMLSS_LONGITUDINAL_MVT_RUN_PUBLICATION_SUITE = "true"
Rscript "paper/R/09-simulation-multivariate-longitudinal/10-run-publication-suite.R"
```

Write the full suite preflight without fitting models:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SUITE_ID = "publication_suite_preflight_default"
$env:GAMLSS_LONGITUDINAL_MVT_SUITE_PREFLIGHT_ONLY = "true"
Rscript "paper/R/09-simulation-multivariate-longitudinal/10-run-publication-suite.R"
```

This creates each role directory, writes `scenario_grid.csv`,
`preflight_checks.*`, `run_metadata.csv`, `package_versions.csv`,
`session_info.txt`, and refreshes the suite `README.md`,
`publication_suite_artifacts.csv`, and `publication_suite_preflight.*`.

Optional evidence roles can be added before launching:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SUITE_INCLUDE_SPECIAL = "true"
$env:GAMLSS_LONGITUDINAL_MVT_SUITE_INCLUDE_GLMM_SENSITIVITY = "true"
```

Useful suite controls:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SUITE_ID = "publication_suite_final_01"
$env:GAMLSS_LONGITUDINAL_MVT_SUITE_ROLES = "pilot,main_core,appendix"
$env:GAMLSS_LONGITUDINAL_MVT_SUITE_MAIN_CORE_REPS = "100"
$env:GAMLSS_LONGITUDINAL_MVT_SUITE_APPENDIX_REPS = "100"
```

The suite is strict by default: full execution requires at least the required
`pilot` and `main_core` roles and runs the publication-readiness audit at the
end. For intentionally partial runs, disable the final strict gate:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SUITE_ROLES = "pilot"
$env:GAMLSS_LONGITUDINAL_MVT_SUITE_STRICT_READINESS = "false"
```

The suite folder writes `README.md` and `publication_suite_artifacts.csv` during
dry runs and refreshes them after each completed role, so interrupted runs still
have a review/status index.

## Sharded Main Runs

For long main runs, use explicit replicate IDs to split the 100-replicate grid
across several output directories. Each shard keeps the true replicate labels in
`case_id`.

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SOURCE = "local"
$env:GAMLSS_LONGITUDINAL_MVT_REPS = "100"
$env:GAMLSS_LONGITUDINAL_MVT_REP_IDS = "1,2,3,4,5,6,7,8,9,10"
$env:GAMLSS_LONGITUDINAL_MVT_OUTPUT_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/main_core_shards/shard_001_010"
$env:GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm,glmm,gee_independence,gee_exchangeable,gee_ar1,gamlss.longitudinal,gamCopula_markov,gamCopula_vine_simplified"
Rscript "paper/R/09-simulation-multivariate-longitudinal/02-run-main-grid.R"
```

After each shard completes, add the full-vine rows only for the
covariate-dependent scenarios in that shard:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_RUN_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/main_core_shards/shard_001_010"
$env:GAMLSS_LONGITUDINAL_MVT_RERUN_METHODS = "gamCopula_vine"
$env:GAMLSS_LONGITUDINAL_MVT_RERUN_DEPENDENCE = "native_covariate_dependent_adjacent"
$env:GAMLSS_LONGITUDINAL_MVT_RERUN_CHECKPOINT_EVERY = "1"
$env:GAMLSS_LONGITUDINAL_MVT_GAMCOPULA_VINE_TIMEOUT_SEC = "7200"
Rscript "paper/R/09-simulation-multivariate-longitudinal/15-rerun-selected-methods.R"
```

Merge completed shards into one readiness-compatible run directory:

```powershell
$env:GAMLSS_LONGITUDINAL_MVT_SHARD_DIRS = "results/jss-exploratory/09-simulation-multivariate-longitudinal/main_core_shards/shard_001_010,results/jss-exploratory/09-simulation-multivariate-longitudinal/main_core_shards/shard_011_020"
$env:GAMLSS_LONGITUDINAL_MVT_MERGED_DIR = "results/jss-exploratory/09-simulation-multivariate-longitudinal/main_core_100rep_merged"
Rscript "paper/R/09-simulation-multivariate-longitudinal/13-merge-run-shards.R"
$env:GAMLSS_LONGITUDINAL_MVT_RUN_DIR = $env:GAMLSS_LONGITUDINAL_MVT_MERGED_DIR
Rscript "paper/R/09-simulation-multivariate-longitudinal/07-review-audit.R"
```

## Main Artifacts

Each run writes:

```text
scenario_grid.csv
preflight_checks.csv
preflight_checks.md
artifact_manifest.csv
artifact_manifest.md
fit_status_by_rep.csv
benchmark_results_by_rep.csv
coefficient_results_by_rep.csv
dependence_recovery_by_rep.csv
variogram_scores_by_rep.csv
runtime_by_rep.csv
benchmark_summary.csv
coefficient_summary.csv
dependence_recovery_summary.csv
variogram_summary.csv
case_method_completion_summary.csv
pilot_feasibility_by_method.csv
pilot_feasibility_by_scenario.csv
pilot_feasibility_overall_method.csv
pilot_feasibility.md
run_metadata.csv
package_versions.csv
session_info.txt
review_audit.csv
review_audit.md
paper_tables/
figures/
review_bundle/
```

The module output root also stores:

```text
publication_readiness_audit.csv
publication_readiness_audit.md
implementation_status.csv
implementation_status.md
```

Publication suite folders additionally write:

```text
README.md
publication_suite_plan.csv
publication_suite_plan.md
publication_suite_artifacts.csv
publication_suite_preflight.csv
publication_suite_preflight.md
publication_readiness_audit.csv
publication_readiness_audit.md
```

The module output root additionally writes the generated study protocol:

```text
study_protocol.md
study_protocol/
```

## Notes

- `N` is interpreted as total rows. Subject counts are reduced as `T` grows.
- `study_protocol.md` is generated from the workflow definitions and maps the
  design, comparators, metrics, expected artifacts, readiness evidence, and
  requirement-to-evidence checklist.
- `implementation_status.md` reports whether source/protocol/preflight evidence
  is ready and explicitly marks full model-fitting evidence as incomplete until
  the publication suite has been run.
- `preflight_checks.*` is written before model fitting and records package,
  comparator, row-cap, resume, and timeout checks for long-run reproducibility.
- `artifact_manifest.*` records file paths, sizes, timestamps, and MD5 hashes;
  `review_bundle/README.md` gives a human-facing index for review.
- The primary GLMM comparator is a random-intercept model.
- Unstructured GEE is attempted with strict timeout at `T = 50` and reported as
  feasibility evidence when it fails or times out.
- `fit_status_by_rep.csv` contains every attempted method, including classified
  `ok`, `warning`, `timeout`, and `error` statuses, plus a short classified
  reason for review.
- `case_method_completion_summary.csv` gives one row per scenario replicate with
  the methods completed, warned, timed out, or failed.
- Variogram scores are computed during fitting for simulation-capable
  `gamlss.longitudinal`, `gamCopula_markov`, `gamCopula_vine_simplified`, and targeted
  `gamCopula_vine` rows. Standard marginal/working
  correlation comparators keep explicit empty rows because they do not define a
  fitted joint response simulator in this workflow.
- Missingness and dropout are intentionally excluded from this module.
- GJRM is retained as historical bivariate/trivariate context rather than a
  high-dimensional comparator.
