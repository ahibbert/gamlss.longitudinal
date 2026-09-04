# JSS Submission Checklist

This checklist converts the JSS-first plan into concrete gates. A submission is
not ready until every `Blocking` item is resolved.

## Current Status

- Phase 1 package correctness is complete; see
  `reviews/2026-08-27-jss-round-1/phase1-acceptance-ledger.md`.
- The simplified Phase 2 scope contains five public evidence modules and no
  multivariate Module 09 dependency. Main-recovery promotion, the strict
  `paper` build, publisher isolation, and compiled-paper visual QA pass.
- Phase 3 manuscript revision remains blocking. The canonical manuscript still
  contains obsolete designs/claims and has not yet received the documented
  public/private reproducibility wording.

## Blocking Gates

| Gate | Status | Required action |
|---|---|---|
| JSS manuscript style | Blocking | Convert the final manuscript to JSS LaTeX style and remove preprint formatting. |
| Main-paper reproducibility | Phase 2 evidence ready | The strict `paper` profile passes and the publisher stages exactly 17 allowlisted renderable assets. Phase 3 must insert/update the corresponding manuscript references. |
| Primary example data | Evidence ready; manuscript edit blocking | Use the fully simulated BCPE/t workflow as the reproducible primary example. |
| Private data policy | Repository resolved; manuscript edit blocking | Keep the withdrawn clinical application outside the public repository and insert the approved reproducibility wording for any separately retained private-data example. |
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
| Smoke profile | Passed | A detached clean-checkout run completed in 1 minute 4 seconds with no target warnings/errors; nine deliberately short nonconverged fits were recorded structurally. |
| Paper profile | Passed | A detached clean checkout passes all 10 Phase 2 gates with 83 exact public paths and no target warnings/errors; canonical text hashes survive Git LF normalization. |
| Full profile | Optional author QA | Retain the checkpointed command, but do not make a publication-scale rerun a routine reviewer or Phase 2 gate. |
| Manifest coverage | Asset classification passed | The read-only canonical audit classifies all 40 current manuscript assets; Phase 3 must update obsolete empirical references. |
| Output hashes | Implemented | Keep output, input, and provenance hashes in each generated replication run. |
| Runtime notes | Evidence recorded | Fresh strict smoke completed in 1 minute 4 seconds; a cached strict paper pipeline took 1 minute 1.6 seconds and the complete detached clean-checkout paper pipeline took 3 minutes 50.6 seconds. Full remains optional. |
| External data notes | Repository resolved; manuscript edit blocking | Public commands have no private-data dependency; state any private-data access exception in the paper. |

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
