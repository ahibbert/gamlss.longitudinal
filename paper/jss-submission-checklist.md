# JSS Submission Checklist

This checklist converts the JSS-first plan into concrete gates. A submission is
not ready until every `Blocking` item is resolved.

## Current Status

The repository already contains a strong package-review foundation:

- `REVIEW.md` routes reviewers through model fitting, likelihood, diagnostics,
  prediction, simulation, benchmarks, and replication.
- `inst/standards/ropensci-srr-compliance.md` maps rOpenSci statistical
  software standards to evidence.
- `paper/replicate.R`, `paper/_targets.R`, and `paper/manifest.csv` define a
  smoke/expanded replication workflow.
- Several validation modules are current, but the draft manuscript still reads
  like a methods preprint and includes private-data and stub-result risks.

## Blocking Gates

| Gate | Status | Required action |
|---|---|---|
| JSS manuscript style | Blocking | Convert the final manuscript to JSS LaTeX style and remove preprint formatting. |
| Main-paper reproducibility | Blocking | Ensure every kept table and figure has a manifest row, generated output, session info, and hash. |
| Primary example data | Blocking | Use public, package-shipped, or simulated data for the main worked example. |
| Private data policy | Blocking | Keep private-data analyses outside the public workflow unless accepted reviewer access instructions are supplied. |
| Stub modules | Blocking | Do not cite module 04 or 05 outputs as main-paper results while manifest rows are `stub`. |
| Unresolved citations | Blocking | Fix all literal question-mark citation artifacts and missing bibliography entries. |
| Main text balance | Blocking | Expand software workflow sections and reduce derivation/simulation dominance. |

## Manuscript Gates

| Gate | Status | Required action |
|---|---|---|
| Title | Ready to revise | Use `gamlss.longitudinal: Longitudinal GAMLSS Models With Copula Dependence in R`. |
| Abstract | Ready to revise | Use the abstract in `paper/jss-manuscript-blueprint.md` as the first replacement draft. |
| Introduction | Ready to revise | Center the software contribution and gap relative to existing R tools. |
| Statistical model | Needs tightening | State exact first-order likelihood assumptions and discrete/missing-panel handling. |
| Software section | Needs expansion | Add interface, fitted object, S3 methods, diagnostics, prediction, simulation, and optimizer controls. |
| Worked example | Needs replacement | Use a reproducible primary example before any private-data example. |
| Validation | Needs condensation | Keep targeted validation; move large grids to supplement. |
| Discussion | Ready to revise | State limits on first-order dependence, runtime, model selection, and private-data examples. |

## Replication Gates

| Gate | Status | Required action |
|---|---|---|
| Smoke profile | Existing | Run from a clean session and record runtime. |
| Expanded profile | Existing | Run once before submission and record runtime, skipped modules, session info, and hashes. |
| Manifest coverage | Partial | Confirm all main-paper outputs are `current`, not `stub`, and are not private-data dependent. |
| Output hashes | Existing | Keep `logs/output_hashes.csv` in generated replication artifacts. |
| Runtime notes | Needs documentation | Add expected smoke and expanded runtimes after final profile runs. |
| External data notes | Needs documentation | State exact data availability for any non-shipped data. |

## Package-Facing Gates

| Gate | Status | Required action |
|---|---|---|
| Main fit API | Stable default | Keep `gamlss_longitudinal()` as the central interface. |
| S3 methods | Existing | Surface stable methods in the paper and reviewer vignette. |
| Workflow vignettes | Existing | Link them from the new JSS start-here vignette and site guide. |
| Smoke test | Existing | Keep `inst/smoke-tests/new-user-smoke.R` as the install check. |
| Non-goals | Needs visibility | State no automatic imputation, no arbitrary higher-order vine fitting, and no unsupplied time-series forecasting. |

## Recommended Main-Paper Result Set

Keep main text compact:

- One implementation capability table.
- One comparison-to-existing-software table.
- One primary worked example screening table.
- One primary worked example final-model/effects table.
- One primary worked example diagnostics figure.
- One continuous validation table or figure.
- One discrete validation table or figure.
- One optimizer guidance table.
- One first-order misspecification or copula-misspecification summary.

Everything else should be supplementary, vignette material, or omitted.

## Pre-Submission Commands

Run these from the repository root before creating the final submission bundle:

```r
devtools::document()
devtools::test(reporter = "summary")
devtools::check(args = "--as-cran")
source("paper/replicate.R")
Sys.setenv(GAMLSS_LONGITUDINAL_JSS_PROFILE = "expanded")
source("paper/replicate.R")
```

Record all warnings, skips, runtime, platform, and optional dependency state.
