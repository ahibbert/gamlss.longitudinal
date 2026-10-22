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
tests and standards evidence. Paper artifact traceability is maintained in the
companion repository described in route 6.

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

The manuscript and its replication workflow live in the separate
[paper repository](https://github.com/ahibbert/working-paper-gamlss-long), outside this
software checkout and the installed R package.

Start with its [reviewer guide](https://github.com/ahibbert/working-paper-gamlss-long/blob/main/paper/REVIEWER.md).
This is a release-preparation draft at its intended publication location, not a
claim that the prepared files are already published or that all three profiles
have passed. Use the exact package/paper source pair identified in the final
release record. Follow its clone/archive assembly and environment instructions;
run paper commands from that assembled workspace, not this package checkout.

The three paper profiles have distinct purposes:

- `smoke` checks installation and a small synthetic workflow.
- `paper` rebuilds assets from versioned analysis-ready inputs and refits the
  declared ACTG application analyses.
- `full` recomputes the declared simulation and timing studies and application
  analyses before rebuilding their assets.

Both `paper` and `full` retain the declared historical ACTG selection ledgers;
the final ACTG fits and post-fit analyses are separate computations. Consult the
guide for the precise scope, retained warnings and uncertainty eligibility.

In the paired paper revision, inspect [paper/manifest.csv](https://github.com/ahibbert/working-paper-gamlss-long/blob/main/paper/manifest.csv),
[paper/_targets.R](https://github.com/ahibbert/working-paper-gamlss-long/blob/main/paper/_targets.R),
[paper/replicate.R](https://github.com/ahibbert/working-paper-gamlss-long/blob/main/paper/replicate.R)
and [paper/tests](https://github.com/ahibbert/working-paper-gamlss-long/tree/main/paper/tests).
These branch links help navigation; the release's pinned revisions govern
reproduction. Its [submission checklist](https://github.com/ahibbert/working-paper-gamlss-long/blob/main/paper/jss-submission-checklist.md)
tracks additional manuscript and submission gates.

The package checks in route 7 and the installed new-user smoke script in the
README review the software. Neither substitutes for the three paper-profile
runs, output checks or the independent reviewer rehearsal.

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
| Paper | Companion repository: `paper/_targets.R`, `paper/R/*.R`, `paper/manifest.csv` | Pinned paper tests, profile receipts and output checks; see route 6 |
| Standards | `R/srr-stats-standards.R`, `inst/standards/*` | p2 standards tests |

## Reviewer Checklist

- Read the relevant route section above.
- Open the source entry points listed there.
- Check the corresponding tests.
- Use `inst/standards/method-traceability.csv` for complete file coverage.
- Use `inst/standards/ropensci-srr-compliance.md` for standards claims and
  remaining TODOs.
- Use the paired paper repository's `paper/manifest.csv` for JSS
  claim-to-artifact mapping; follow route 6 for release and rehearsal gates.
- Use `cran-comments.md` for final submission check results.

## Known Limitations

Remaining standards work should stay explicit rather than be implied as met.
Current known limitations are tracked in
`inst/standards/ropensci-srr-compliance.md`.

Expected remaining review tasks:

- freeze final JSS manuscript tables and figures;
- verify the assembled paper workflow end to end for all three profiles;
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

That checkpoint is historical; it does not verify the current split repositories
or the current paper profiles. Since then, the codebase has continued to change.

Use `cran-comments.md` as the authoritative record for final submission checks.
