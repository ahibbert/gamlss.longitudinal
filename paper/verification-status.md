# Replication verification status

Last updated: 2026-09-04 (Australia/Sydney).

## Gate definition

The simplified Phase 2 gate is the portable public `paper` profile built from
five versioned evidence modules. A fresh publication-scale `full` rerun,
Module 09, detached signatures, and author-held approval keys are not reviewer
or submission gates. The optional `full` command remains available for future
author-side recomputation.

## Completed checks

- Phase 1 passed independently; see
  `reviews/2026-08-27-jss-round-1/phase1-acceptance-ledger.md`.
- The complete ordinary source-tree suite passes with no test failures or test
  warnings. The only skips are two explicitly opt-in extended SRR stress tests.
- A full local Windows R 4.4.1 source-package check completed with status OK:
  zero errors, zero warnings, and zero notes. The later Phase 2 changes affect
  CRAN-excluded paper files and tests, not packaged R code.
- Strict-warning smoke recomputation completed with no target warnings/errors.
  Three deliberately short fits that reached their optimizer limit are retained
  as structured `fit-events.csv` rows rather than leaked target warnings.
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
- The read-only canonical manuscript audit classified 40/40 assets with zero
  unclassified items. A scan of public-derived text inputs found no absolute
  Windows, OneDrive, `/Users`, or `/home` paths.
- Module 09 is absent from active claims, artifact contracts, manifests, and
  target graphs. Its source and tiny unit-test fixtures remain for future work.

## Checks in progress

- The NBI/Clayton n=500, T=4, R=100 production run is completing from fixed
  seeds and append-only checkpoints. Its completed checkpoints have one source
  hash and settings hash, no duplicate method/replicate keys, and all completed
  attempts are successful, converged, and descriptively complete.
- After NBI/Clayton finishes, the automatic finalizer will build the 26-file
  CSV-only main-recovery bundle and rerun the lean optimizer exporter.
- The final bundle must then be promoted, validated, and exercised through a
  clean `paper` build, a final clean smoke build, publisher dry-run/apply in a
  disposable paper clone, and an asset-staged manuscript compile/visual check.
- A context-free independent audit will run after those checks and every
  in-scope defect it finds will be fixed before this ledger is closed.

## Separate Phase 3 manuscript work

The canonical Phase 2 manuscript audit is expected to fail until authoring is
updated. At the latest run it found nine pending claim markers, three
pending/missing scope anchors, three pending artifacts, and eight known stale
empirical blocks. These are Phase 3 wording, structure, citation, and asset-use
tasks; they do not invalidate the underlying Phase 2 public evidence gate.
