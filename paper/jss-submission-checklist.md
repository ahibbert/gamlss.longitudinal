# JSS Submission Checklist

This checklist converts the JSS-first plan into concrete gates. A submission is
not ready until every `Blocking` item is resolved.

## Current Status

- Phase 1 package correctness is complete; see
  `reviews/2026-08-27-jss-round-1/phase1-acceptance-ledger.md`.
- The simplified Phase 2 scope contains five public evidence modules and no
  multivariate Module 09 dependency. Final main-recovery promotion and the
  clean `paper` build remain in progress.
- Phase 3 manuscript revision remains blocking. The canonical manuscript still
  contains obsolete designs/claims and has not yet received the documented
  public/private reproducibility wording.

## Blocking Gates

| Gate | Status | Required action |
|---|---|---|
| JSS manuscript style | Blocking | Convert the final manuscript to JSS LaTeX style and remove preprint formatting. |
| Main-paper reproducibility | Phase 2 in progress | Finish the main-recovery promotion and clean `paper` build; then apply only allowlisted assets to the manuscript clone. |
| Primary example data | Evidence ready; manuscript edit blocking | Use the fully simulated BCPE/t workflow as the reproducible primary example. |
| Evidence modules | Resolved | Five public modules are active; Module 09 is explicitly deferred and no active item is a stub. |
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
| Smoke profile | Passed locally | Repeat once from the final clean commit and record cold/warm runtime. |
| Paper profile | In progress | Complete after authoritative main-recovery inputs are promoted; require zero target warnings/errors. |
| Full profile | Optional author QA | Retain the checkpointed command, but do not make a publication-scale rerun a routine reviewer or Phase 2 gate. |
| Manifest coverage | Asset classification passed | The read-only canonical audit classifies all 40 current manuscript assets; Phase 3 must update obsolete empirical references. |
| Output hashes | Implemented | Keep output, input, and provenance hashes in each generated replication run. |
| Runtime notes | Pending final run | Record final clean smoke and paper runtimes; document full separately as optional. |

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

```text
Rscript paper/replicate.R --profile smoke
Rscript paper/replicate.R --profile paper
Rscript paper/audit-manuscript.R --paper-repo <clean-paper-clone>
Rscript paper/publish-assets.R --paper-repo <clean-paper-clone> --profile paper --dry-run

Rscript -e "devtools::document(); devtools::test(reporter = 'summary')"
R CMD build .
R CMD check --as-cran <source-tarball>
```

Record all warnings, skips, runtime, platform, and optional dependency state.
