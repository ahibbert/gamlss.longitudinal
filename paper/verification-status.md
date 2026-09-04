# Replication verification status

Last updated: 2026-09-05 (Australia/Sydney).

## Gate definition

The simplified Phase 2 gate is the portable public `paper` profile built from
five versioned evidence modules. A fresh publication-scale `full` rerun,
Module 09, detached signatures, and author-held approval keys are not reviewer
or submission gates. The optional `full` command remains available for future
author-side recomputation.

## Completed checks

- Phase 1 passed independently; see
  `reviews/2026-08-27-jss-round-1/phase1-acceptance-ledger.md`.
- The complete close-out source-tree suite passes all ordinary P0/P1/P2 tests.
  An initial run found one repository-boundary failure caused by two
  documentation strings; those strings were corrected and the full suite then
  passed. The only intentional skips are two opt-in extended SRR stress tests.
- A full local Windows R 4.4.1 source-package check completed with status OK:
  zero errors, zero warnings, and zero notes. The later Phase 2 changes affect
  CRAN-excluded paper files and tests, not packaged R code.
- Strict-warning smoke recomputation from a detached clean checkout completed
  in 1 minute 4 seconds with no target warnings/errors. Nine deliberately short
  nonconverged fits are retained as structured `fit-events.csv` rows rather
  than leaked target warnings.
- The optimizer bundle contains three scenarios and 97 paired replicates per
  scenario: 582/582 attempts converged and were retained. Both producer and
  submission validators pass, including the cross-checkout checkpoint test.
- The missingness bundle contains 240 scenario/replicate tasks and 480 fits:
  399 retained, 81 named failures, 298 expected warnings, and zero unexpected
  warnings. Monotone dropout is the headline and intermittent MAR is a
  sensitivity analysis.
- The 2,400-fit copula-misspecification bundle and 62-attempt fit-scaling bundle
  pass their submission schemas and reconstruction checks.
- BCPE/t completed all 200 method attempts successfully and converged. Its
  exact-settings resume accepted 200/200 checkpoints, found zero pending, and
  changed no checkpoint timestamp or hash.
- NBI/Clayton completed all 200 method attempts successfully and converged. Its
  exact-settings resume also accepted 200/200 checkpoints, found zero pending,
  and changed no checkpoint timestamp or hash.
- The promoted main-recovery bundle contains 26 CSV files and 400/400
  successful, converged, retained attempts. All 31 evidence checks pass; the
  public bundle contains no RDS, fitted object, residual, or private asset.
- The strict `paper` profile passes all 10 Phase 2 gates with 83 exact public
  manifest paths and 17 verified claims. `target-events.csv` reports no target
  warnings/errors; the final cached pipeline took 1 minute 1.6 seconds.
- Publisher dry-run/apply in a clean disposable manuscript clone staged exactly
  17 assets and changed only `charts/` and `tables/`. It did not edit, commit,
  or push manuscript source.
- After applying three known manuscript-only compile fixes in that disposable
  clone, the staged paper compiled to a visually inspected 47-page PDF. The
  generated recovery tables and figures are nonblank and unclipped. Two
  unresolved references and eight citation occurrences remain assigned to
  Phase 3.
- The read-only canonical manuscript audit classified 40/40 assets with zero
  unclassified items. A scan of public-derived text inputs found no absolute
  Windows, OneDrive, `/Users`, or `/home` paths.
- Module 09 is absent from active claims, artifact contracts, manifests, and
  target graphs. Its source and tiny unit-test fixtures remain for future work.

## Checks in progress

- Give a context-free independent reviewer only `paper/REVIEWER.md`; fix every
  in-scope ambiguity or failure before closing Phase 2.

## Separate Phase 3 manuscript work

The canonical Phase 2 manuscript audit is expected to fail until authoring is
updated. At the latest run it found nine pending claim markers, three
pending/missing scope anchors, three pending artifacts, and eight known stale
empirical blocks. These are Phase 3 wording, structure, citation, and asset-use
tasks; they do not invalidate the underlying Phase 2 public evidence gate.
