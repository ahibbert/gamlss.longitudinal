# Reviewer Guide

This guide is the front door for CRAN reviewers, statistical software
reviewers, JSS reviewers, and contributors. It points to the main review paths
without trying to list every helper function.

The exhaustive source-to-test map is
`inst/standards/method-traceability.csv`. It is checked by
`tests/testthat/test-p2-method-traceability.R`, including coverage for every
`R/*.R` source file.

## How To Use This Guide

Most reviewers need to answer four questions first:

- How is the model fit?
- How are likelihood, Hessian, and covariance calculations computed?
- How do prediction, simulation, diagnostics, and plots work?
- How are benchmark, coverage, and paper-replication claims generated?

Use the route sections below for a first pass. Use
`inst/standards/method-traceability.csv` for a complete audit of source files,
tests, paper modules, and standards evidence.

The source files in `R/` use prefix-based names. For example, `model-*` files
contain fitting, prediction, simulation, and fitted-object helpers;
`optimizer-*` files contain RS and CG optimization logic; `likelihood-*`,
`hessian-*`, and `model-vcov*` files contain uncertainty calculations; and
`benchmark-*`/`coverage-*` files contain opt-in reviewer evidence.

## Reviewer Routes

### 1. Model Fitting

Start with the public fit wrapper:

- `R/model-fit.R`: exported `gamlss_longitudinal()`.
- `R/model-fit-entrypoint.R`: sequences fit preparation, optimization, and
  finalization.

The fit follows three phases:

| Phase | Main entry point | What to review |
|---|---|---|
| Prepare | `R/model-fit-workflow.R::.gl_prepare_fit_workflow()` | input validation, formulas, model matrices, starting values |
| Optimize | `R/model-fit-optimizer.R::.gl_run_prepared_fit_optimizer()` | RS/CG dispatch and optimizer state |
| Finalize | `R/model-fit-finalize.R::.gl_finalize_prepared_fit()` | convergence metadata, fitted object, optional vcov cache |

Supporting source areas:

- Input/data policy: `R/model-preprocess*.R`, `R/model-column-policy.R`,
  `R/model-data-shape.R`, `R/model-time.R`, `R/model-formulas.R`,
  `R/model-missing-panels.R`.
- Model matrices and eta: `R/model-matrix*.R`, `R/model-eta.R`.
- Starting values and warm starts: `R/model-starting*.R`,
  `R/model-parameter-transforms.R`, `R/model-warm-start*.R`.
- Result assembly: `R/model-fit-result.R`, `R/model-fit-convergence.R`,
  `R/model-fit-object*.R`, `R/model-fit-reporting.R`.

Optimizer entry points:

- RS: `R/optimizer-rs-runner.R::.gl_run_rs_optimizer()`.
- RS backfitting calculation:
  `R/optimizer-rs-backfitting-iteration.R::.gl_rs_backfitting_iteration()`.
- RS score dispatch: `R/optimizer-rs-score-dispatch.R`,
  `R/optimizer-rs-score-margin.R`, `R/optimizer-rs-score-copula.R`.
- CG: `R/optimizer-cg-runner.R::.gl_run_cg_optimizer()`.
- CG loop/line search: `R/optimizer-cg-outer-loop.R`,
  `R/optimizer-cg-outer-iteration.R`, `R/optimizer-cg-line-search.R`,
  `R/optimizer-cg-acceptance.R`.

Evidence: `tests/testthat/test-p0-model-fit-entrypoint.R`,
`tests/testthat/test-p0-model-fit-workflow.R`,
`tests/testthat/test-p0-model-optimizer-dispatch.R`,
`tests/testthat/test-p0-optimizer-rs-runner.R`,
`tests/testthat/test-p0-optimizer-cg-runner.R`, and
`tests/testthat/test-p1-core.R`.

### 2. Likelihood, Hessian, And Vcov

Start with likelihood evaluation:

- `R/likelihood-evaluation.R`: joint likelihood entry point,
  including `calc_likelihood_minimal()`.
- `R/likelihood-evaluation-margin.R`: margin density, CDF, and derivative
  evaluation.
- `R/likelihood-evaluation-copula.R` and
  `R/likelihood-evaluation-copula-update.R`: pair-level copula likelihood
  evaluation.
- `R/likelihood-copula-rectangle-*.R`: discrete rectangle likelihood pieces.
- `R/likelihood-gradient.R` and `R/likelihood-scores.R`: analytical score
  paths.

Then review uncertainty:

- `R/hessian-*.R`: analytical and semi-analytical Hessian setup, derivatives,
  and block assembly.
- `R/model-vcov.R`: exported `vcov.gamlss.longitudinal()`.
- `R/model-vcov-setup.R`, `R/model-vcov-hessian.R`,
  `R/model-vcov-solve.R`, `R/model-vcov-result.R`, `R/model-vcov-smooth.R`:
  covariance setup, solving, fallback handling, and result shaping.

Evidence: `tests/testthat/test-p0-likelihood-*.R`,
`tests/testthat/test-p0-discrete-rectangle-likelihood.R`,
`tests/testthat/test-p0-hessian-*.R`,
`tests/testthat/test-p0-model-vcov.R`, and
`tests/testthat/test-p0-cg-gradient.R`.

### 3. Prediction, Simulation, Diagnostics, And Plots

Prediction:

- `R/model-predict.R`: `predict.gamlss.longitudinal()`.
- `R/model-newdata*.R`: newdata translation, default columns, factor-level
  alignment, response validation, and model-matrix alignment.
- `R/model-predict-values.R`, `R/model-predict-distributional.R`,
  `R/model-predict-intervals.R`, `R/model-predict-se.R`: evaluation values,
  distributional summaries, intervals, and standard errors.

Simulation:

- `R/model-simulate.R`: `simulate.gamlss.longitudinal()`.
- `R/model-simulate-newdata.R`, `R/model-simulate-newdata-helpers.R`,
  `R/model-simulate-copula.R`: simulation grids, dependence alignment, and
  copula draws.
- `R/simulation-*.R`: helper functions for tests and JSS simulation modules.

Diagnostics and plots:

- `R/diagnostics-*.R`: diagnostic families, scoring, PIT-style diagnostics,
  and split-family checks.
- `R/missingness_diagnostics.R`, `R/missingness-diagnostics-*.R`: missingness
  diagnostic workflow.
- `R/plot-method.R`, `R/plot-terms*.R`, `R/plot-dist*.R`,
  `R/plot-margin-fit*.R`, `R/plot-copula-fit*.R`,
  `R/plot-copula-diagnostics*.R`: dashboard, term, distribution, margin,
  copula, and diagnostic plots.

Evidence: `tests/testthat/test-p0-model-newdata.R`,
`tests/testthat/test-p0-prediction-intervals.R`,
`tests/testthat/test-p0-simulation-helpers.R`,
`tests/testthat/test-p0-diagnostics-namespace.R`,
`tests/testthat/test-p0-missingness-diagnostics.R`, and
`tests/testthat/test-p1-plots-summary.R`.

### 4. Distribution And Copula Backend

Start here:

- `R/copula-family-links.R`: family and link metadata.
- `R/copula-backend-*.R`: CDF, density, h-functions, tau conversion, and
  derivatives.
- `R/copula_selection.R`, `R/copula-selection-*.R`: copula selection.
- `R/joint_distribution_selection.R`, `R/joint-distribution-selection-*.R`:
  joint margin/copula selection.
- `R/link-functions-*.R`: link and inverse-link helpers.

Evidence: `tests/testthat/test-p0-copula-backend.R`,
`tests/testthat/test-p0-copula-selection.R`,
`tests/testthat/test-p0-joint-distribution-selection.R`, and
`tests/testthat/test-p2-srr-distribution-methods.R`.

### 5. Benchmark And Coverage Evidence

These are opt-in reviewer workflows, not the main user-facing fit path.

Benchmark evidence:

- `R/benchmark-adoption-scenarios.R`: predefined scenarios.
- `R/benchmark-comparators.R`: comparator fits and metrics.
- `R/benchmark-report.R`, `R/benchmark-summary.R`: report and summary output.
- `R/select-margin-*.R`: margin screening and `fitDist()` orchestration.
- `R/reporting-model-spec.R`, `R/reporting-table.R`: report helpers.

Coverage evidence:

- `R/coverage-runner.R`: coverage simulation runner.
- `R/coverage-fitters.R`: longitudinal and comparator fitters.
- `R/coverage-estimands.R`: estimands and truth extraction.
- `R/coverage-benchmarks.R`: scenario summaries.
- `R/coverage-report.R`: report generation.

Evidence: `tests/testthat/test-p0-benchmark-*.R`,
`tests/testthat/test-p1-benchmark-comparators.R`,
`tests/testthat/test-p0-coverage-*.R`,
`tests/testthat/test-p1-coverage-simulations.R`, and
`tests/testthat/test-p1-adoption-workflow.R`.

### 6. JSS Paper Replication

The paper workflow is repository-only and excluded from CRAN/source builds.

Start here:

- `paper/README.md`: replication instructions.
- `paper/manifest.csv`: claim-to-artifact map.
- `paper/_targets.R`: pipeline wiring.
- `paper/R/*.R`: numbered replication modules.
- `paper/replicate.R`: smoke and expanded replication entry point.

Planned JSS modules:

1. Continuous simulation: BCPE margin with t copula.
2. Discrete simulation: Delaporte margin with Clayton copula.
3. Joint optimization versus separate optimization.
4. Sensitivity to missingness and dropout.
5. LIPID clinical-trial application.
6. RAND doctor-visits application.

Commands:

```r
source("paper/replicate.R")
Sys.setenv(GAMLSS_LONGITUDINAL_JSS_PROFILE = "expanded")
source("paper/replicate.R")
```

Generated tables, figures, logs, session information, and hashes are written
under `results/jss-replication/<profile>/`. Private application data are not
committed; local paths are provided with `GAMLSS_LONGITUDINAL_LIPID_DATA` and
`GAMLSS_LONGITUDINAL_RAND_DATA`.

Evidence: `tests/testthat/test-p2-jss-replication-skeleton.R`,
`paper/manifest.csv`, and `paper/README.md`.

### 7. CRAN And rOpenSci Review

CRAN-facing checks:

```r
devtools::test(reporter = "summary")
rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"), error_on = "warning")
urlchecker::url_check()
spelling::spell_check_package()
```

Source-package hygiene:

```sh
R CMD build .
tar -tzf gamlss.longitudinal_*.tar.gz
```

The tarball should not include local state (`.RData`, `.Rhistory`,
`.Rprofile`), generated sites/results (`docs/`, `results/`, `examples/`),
paper sources (`paper/`), check directories, `Rplots.pdf`, or root benchmark
metric outputs.

rOpenSci evidence:

- `CONTRIBUTING.md`
- `R/srr-stats-standards.R`
- `inst/standards/ropensci-srr-compliance.md`
- `inst/standards/method-traceability.csv`
- `tests/testthat/test-p2-ropensci-standards.R`
- `tests/testthat/test-p2-srr-input-policy.R`
- `tests/testthat/test-p2-srr-error-map.R`
- `tests/testthat/test-p2-srr-extended.R`

Extended tests:

```sh
GAMLSS_LONGITUDINAL_EXTENDED_TESTS=true Rscript -e "pkgload::load_all(); testthat::test_dir('tests/testthat')"
```

Optional comparator packages are not hard dependencies. For example, `gamlss2`
is used only by opt-in comparator workflows and should skip or report
unavailability when absent.

## Quick Traceability Table

| Review area | Primary source | Main evidence |
|---|---|---|
| Fit workflow | `R/model-fit*.R` | `test-p0-model-fit-*.R`, `test-p1-core.R` |
| Input policy | `R/model-preprocess*.R`, `R/model-column-policy.R` | `test-p0-model-preprocess.R`, `test-p2-srr-input-policy.R` |
| Matrices | `R/model-matrix*.R`, `R/model-eta.R` | `test-p0-model-matrix-bundle.R` |
| Optimizers | `R/optimizer-rs-*.R`, `R/optimizer-cg-*.R` | `test-p0-optimizer-*.R` |
| Likelihood | `R/likelihood-*.R` | likelihood and rectangle tests |
| Hessian/vcov | `R/hessian-*.R`, `R/model-vcov*.R` | Hessian and vcov tests |
| Prediction | `R/model-predict*.R`, `R/model-newdata*.R` | prediction and edge-case tests |
| Simulation | `R/model-simulate*.R`, `R/simulation-*.R` | simulation tests |
| Diagnostics/plots | `R/diagnostics-*.R`, `R/plot-*.R` | diagnostics and plot tests |
| Copula/backend | `R/copula-*.R`, `R/joint-distribution-*.R` | copula and selection tests |
| Benchmarks | `R/benchmark-*.R`, `R/select-margin-*.R` | benchmark/adoption tests |
| Coverage | `R/coverage-*.R` | coverage tests |
| Paper | `paper/_targets.R`, `paper/R/*.R`, `paper/manifest.csv` | JSS skeleton test |
| Standards | `R/srr-stats-standards.R`, `inst/standards/*` | p2 standards tests |

## Reviewer Checklist

- Read the relevant route section above.
- Open the source entry points listed there.
- Check the corresponding tests.
- Use `inst/standards/method-traceability.csv` for complete file coverage.
- Use `inst/standards/ropensci-srr-compliance.md` for standards claims and
  remaining TODOs.
- Use `paper/manifest.csv` for JSS claim-to-artifact mapping.
- Use `cran-comments.md` for final submission check results.

## Known Limitations

Remaining standards work should stay explicit rather than be implied as met.
Current known limitations are tracked in
`inst/standards/ropensci-srr-compliance.md`.

Expected remaining review tasks:

- freeze final JSS manuscript tables and figures;
- replace paper-module stubs with final analyses;
- record final platform-specific CRAN/GitHub check results in
  `cran-comments.md`;
- decide whether legacy helpers can be removed after historical scripts are
  reviewed;
- update benchmark and coverage evidence if the comparison design changes.

## Most Recent Local Verification

Refresh this section immediately before CRAN submission or external review.
The last full submission-readiness checkpoint recorded here was run on
2026-06-14 on Windows 11 x64 with R 4.4.1. It recorded passing documentation,
routine tests, extended tests, source build, CRAN-style check, URL check,
spelling check, new-user smoke test, and default JSS smoke replication. The
CRAN-style check had 0 ERRORs, 0 WARNINGs, and 3 NOTEs.

Since then, the codebase has continued to be reorganized for reviewability.
Use `cran-comments.md` as the authoritative record for final submission checks.
