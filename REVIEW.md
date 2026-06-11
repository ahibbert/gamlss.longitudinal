# Reviewer Guide

This file is the reviewer-facing entry point for `gamlss.longitudinal`. It
collects the package checks, statistical-software standards evidence, and paper
replication workflow in one place.

## Source Map

The `R/` directory is organised by review task:

- `model-*.R`: main fitting, preprocessing, model matrices, prediction,
  simulation, inference, vcov, checks, effects, accessors, and copula summaries.
- `likelihood.R` and `hessian-*.R`: likelihood evaluation, score components,
  and semi-analytical Hessian assembly.
- `copula-backend-*.R`: native and optional VineCopula backend utilities,
  family-specific copula formulas, and backend dispatch.
- `adoption-margin-selection.R`, `copula_selection.R`, and
  `joint_distribution_selection.R`: margin, copula, and joint model-screening
  workflows.
- `diagnostics-*.R`, `plot-*.R`, and `missingness_diagnostics.R`: residual,
  distributional, copula, term-effect, fit-overlay, and missingness diagnostics.
- `coverage-*.R`, `benchmark-*.R`, and `coverage_report.R`: opt-in simulation,
  comparator benchmark, summary, and reporting workflows.
- `simulation_helpers.R` and `legacy-simulation-utils.R`: current simulation
  helpers and retained historical/research-era helpers.

`R/model-fit.R` remains intentionally large because it contains the main
`gamlss_longitudinal()` fitting loop. It has not been decomposed internally in
this cleanup pass to avoid changing optimiser behaviour.

## CRAN And Package Review

The CRAN-facing objective is a clean source package with fast routine examples,
tests, and vignettes. Generated pkgdown, result, check, manuscript, and local
development artefacts are excluded from the source tarball.

Recommended commands:

```r
devtools::test(reporter = "summary")
rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"), error_on = "warning")
urlchecker::url_check()
spelling::spell_check_package()
```

Source package hygiene can be checked with:

```sh
R CMD build .
tar -tzf gamlss.longitudinal_*.tar.gz
```

The tarball should not include local state (`.RData`, `.Rhistory`,
`.Rprofile`), generated sites or results (`docs/`, `results/`, `examples/`),
paper sources (`paper/`), check directories, `Rplots.pdf`, or root benchmark
metric outputs. `cran-comments.md` records the latest local and remote check
results immediately before submission.

Optional comparator packages are handled as optional dependencies. `gamlss2` is
not a hard dependency; methods that require it are opt-in and should skip or
report unavailability when the package is absent.

## Statistical Software Review

The package is prepared against rOpenSci statistical software standards for
general statistical software, regression and supervised learning, and
probability distributions.

Primary evidence:

- `CONTRIBUTING.md`: terminology, input policy, missingness policy, numerical
  assumptions, and extended-test instructions.
- `R/srr-stats-standards.R`: machine-readable `srr` tags.
- `inst/standards/ropensci-srr-compliance.md`: human-readable standards
  crosswalk, including validation/test evidence and remaining TODO items.
- `tests/testthat/`: deterministic unit and integration tests for validation,
  fitting, prediction, simulation, diagnostics, and copula parity.

Routine tests are intended to be CRAN- and CI-friendly. Extended stochastic and
stress tests are opt-in:

```sh
GAMLSS_LONGITUDINAL_EXTENDED_TESTS=true Rscript -e "pkgload::load_all(); testthat::test_dir('tests/testthat')"
```

Current known standards limitations are intentionally left as `Partial/TODO` in
the crosswalk. They mostly concern broader external benchmark datasets, final
manuscript-output freezing, and platform-specific numerical tolerance evidence
after the full CI matrix has run.

## Paper Replication Review

The paper workflow is split into a CRAN-safe installed smoke entry point and a
repository-only full replication workflow.

Recommended smoke workflow:

```r
source(system.file("smoke-tests", "new-user-smoke.R", package = "gamlss.longitudinal"))
source("paper/replicate.R")
```

The default `paper/replicate.R` profile is smoke-sized. The expanded profile is
not CRAN-safe and may take a long time:

```r
Sys.setenv(GAMLSS_LONGITUDINAL_JSS_PROFILE = "expanded")
source("paper/replicate.R")
```

`paper/manifest.csv` is the source claim-to-artifact map. Generated tables,
figures, logs, session information, and output hashes are written under
`results/jss-replication/<profile>/`, which is excluded from source builds.

## Runtime Tiers

- Fast smoke: installation smoke test and targeted examples; expected to run in
  minutes.
- Routine package review: `devtools::test()` and `R CMD check --as-cran`;
  expected to be CI-friendly.
- Extended review: recovery/stress tests and expanded paper replication; opt-in
  because runtime depends on platform and optional comparator availability.

## Current Local Baseline

The latest local cleanup baseline was run on 2026-06-11 on Windows 11 x64 with
R 4.4.1:

- `devtools::test(reporter = "summary")`: pass.
- Extended `testthat::test_dir()` with `GAMLSS_LONGITUDINAL_EXTENDED_TESTS=true`:
  pass, 0 failures, 1 expected warning, 4 CRAN skips, 1099 passing expectations.
- `R CMD check --as-cran --no-manual`: 0 errors, 0 warnings, 2 notes.
- `urlchecker::url_check()`: pass.
- `spelling::spell_check_package()`: pass.
- New-user smoke and default JSS smoke replication: pass.
