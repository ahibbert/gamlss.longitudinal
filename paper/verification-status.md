# Replication verification status

Last updated: 2026-09-22 (Australia/Sydney).

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
  sensitivity analysis. A generated retention figure now reports every
  retained/attempted denominator and its Wilson interval before the
  conditional recovery figures are interpreted.
- The focused six-by-six copula-misspecification study completed all 3,600
  candidate fits: 100 data sets from each of six generating copulas, with all
  six candidates fitted to every data set, at `n = 150` and Kendall's tau
  0.25. All fits succeeded, converged, and were retained. The run recorded 116
  expected numerical/inference diagnostic warnings, zero nonconvergence
  warnings, zero unexpected warnings, and no failed selection sets. The clean
  Windows R 4.4.1 production run used four workers and completed in 2 hours 4
  minutes 36 seconds. The 62-attempt fit-scaling bundle also passes its
  submission schema and reconstruction checks.
- BCPE/t completed all 200 method attempts successfully and converged. Its
  exact-settings resume accepted 200/200 checkpoints, found zero pending, and
  changed no checkpoint timestamp or hash.
- NBI/Clayton completed all 200 method attempts successfully and converged. Its
  exact-settings resume also accepted 200/200 checkpoints, found zero pending,
  and changed no checkpoint timestamp or hash.
- The promoted main-recovery bundle contains 26 CSV files and 400/400
  successful, converged, retained attempts. All 31 evidence checks pass; the
  public bundle contains no RDS, fitted object, residual, or private asset.
- The strict `paper` profile passes all 10 Phase 2 gates with 84 exact public
  manifest paths and 17 verified claims. `target-events.csv` reports no target
  warnings/errors; the final cached pipeline took 1 minute 1.6 seconds and a
  detached clean-checkout pipeline took 3 minutes 50.6 seconds.
- The first context-free audit exposed main-recovery hashes recorded against
  CRLF working-tree bytes while Git enforces LF. The repaired bundle uses
  newline/BOM-canonicalized text hashes and canonical byte counts. A regression
  converts the complete bundle to LF before validating it, and a newly detached
  checkout passes the ordinary `paper` command.
- Terminal manifest validation now parses, constrains, covers, and recomputes
  every `output_hashes.csv` row. Failed target runs write available run/target
  metadata and a reviewer summary before returning the original error.
- The final context-free reviewer repeated the paper profile from a new
  detached checkout: 19/19 targets and 10/10 gates passed, with zero target or
  fit events. Independent checks matched all 511 input, 88 output, and 11
  provenance hashes. Original-resolution inspection confirmed the shortened
  NBI term labels and missingness facet labels are legible and unclipped; the
  zero-event paper run contains a schema-stable empty `fit-events.csv`.
- Publisher dry-run/apply in a clean disposable manuscript clone staged exactly
  17 assets and changed only `charts/` and `tables/`. It did not edit, commit,
  or push manuscript source.
- After applying three known manuscript-only compile fixes in that disposable
  clone, the staged paper compiled to a visually inspected 47-page PDF. The
  generated recovery tables and figures are nonblank and unclipped. Two
  unresolved references and eight citation occurrences remain assigned to
  Phase 3.
- The read-only canonical manuscript audit classified 40/40 assets with zero
  unclassified items. A scan of public-derived text inputs found no checkout,
  user-home, OneDrive, `/Users`, or private-data paths. Generic R interpreter
  paths are retained as execution provenance and are not input dependencies.
- Module 09 is absent from active claims, artifact contracts, manifests, and
  target graphs. Its source and tiny unit-test fixtures remain for future work.

## Close-out

The simplified Phase 2 evidence and independent reviewer gates are complete.
Module 09 and publication-scale full recomputation remain explicitly optional.

## Separate Phase 3 manuscript work

The canonical Phase 2 manuscript audit is expected to fail until authoring is
updated. At the latest run it found nine pending claim markers, three
pending/missing scope anchors, three pending artifacts, and eight known stale
empirical blocks. These are Phase 3 wording, structure, citation, and asset-use
tasks; they do not invalidate the underlying Phase 2 public evidence gate.

## Current ACTG manuscript audit

The ACTG workflow added after the Phase 2 close-out has now passed a genuinely
fresh Windows R 4.4.1 `paper` run at package commit
`a5510d650ab0faa5e0d1461618d35eca77f87bca`. All 27 targets completed from an
empty store in 39 minutes 6 seconds with zero target warnings/errors and 23/23
passing figure checks. The run installed the checked-out package source and
recorded a clean Git state, 530 input hashes, and 129 output hashes.

The computational result does not yet close the manuscript-publication gate.
The current 49-page ACTG PDF is newer than the TeX held by the canonical paper
repository, which still contains the earlier application. Audit of the closest
available ACTG TeX classified all 45 assets but found eight retired simulation
artifacts still referenced. Seven generated ACTG tables also remain inline,
and 26 of 29 active public manuscript assets do not yet have author-approved
reference hashes.

The asset publisher stages 36 allowlisted files without changing paper source.
A disposable compile with the closest ACTG TeX succeeds after publication, but
the stale canonical bibliography leaves six undefined citations and the draft
still references the retired joint-versus-separate section. The formal
main-recovery producer also remains out of alignment with four recovery tables
and four recovery figures in the current manuscript; this must be resolved
before final asset hashes are approved. See
`reviews/2026-09-22-end-to-end-replication-audit.md` for the itemized close-out
plan.
