# Replication verification status

Last updated: 2026-08-27 (Australia/Sydney).

## Completed locally

- Windows, R 4.4.1: clean-snapshot `smoke` run with lockfile restore and
  checked-out source installation; no manual intervention.
- Windows: strict-warning `smoke` recomputation from invalidated simulation
  targets and strict-warning `paper` build.
- Manifest audit against the canonical paper clone: 41 of 41 manuscript assets
  classified, zero unclassified.
- Asset publisher dry-run and apply in disposable clean paper clones; only
  allowlisted `charts/` and `tables/` assets changed and no commit/push occurred.
- Asset-staged 47-page manuscript compile and visual inspection. Three existing
  TeX/BibTeX blockers were fixed only in the disposable QA clone; the manual
  fixes and remaining undefined citations are listed in
  `overleaf-manuscript-fixes.md`.
- Missingness checkpoint resume test: a completed checkpoint set resumed with
  zero repeated fits.
- Four context-free reviewer smoke passes, with each reported defect or
  ambiguity fixed before the next clean-snapshot run. The final clean run
  completed without intervention in 2 minutes 27.6 seconds: 12 of 12 targets,
  zero target warnings/errors, 24 explicitly reported optimizer-limit events,
  zero execution failures, and an independently matching raw manifest hash.
- The committed snapshot `fd2ffc1` was checked in a detached clean worktree.
  Lockfile restore and checked-out source installation succeeded; `smoke`
  completed 13 of 13 targets and `paper` completed 16 of 16 targets, both with
  zero target warnings/errors and a clean Git provenance record.
- Module 08 checkpoint replacement and restart behavior was exercised with a
  temporary one-case benchmark. The second invocation reused the same run
  directory, skipped the completed case, retained one status row, and created
  no duplicate run directory.

## Submission gates not yet executed

- The `full` public Monte Carlo profile was started from an empty store and
  checkpoint set on 2026-08-27 and is currently running. The first restart was
  an intentional worker-scaling check; valid atomic checkpoints were retained.
  The published 20-replicate missingness study completed in 51 minutes 22
  seconds. The active run is continuing as a hidden 12-worker process with
  stdout/stderr under `results/jss-replication/full/logs/`. End-to-end
  completion and the subsequent no-repeat resume run remain required. During
  this gate, module 08 exposed a restart defect in its legacy driver; the driver
  now reads its completed/checkpoint tables, validates the scenario and family
  grids, skips synchronized completed cases, and uses replace-safe CSV writes.
  The active run directory is recorded so its existing progress remains usable.
- The metric-level comparison is now wired into the full graph and its complete
  mirrored-layout test passed 2,275 comparisons. Acceptance against freshly
  recomputed results remains pending completion of the active full run.
- The Windows/Ubuntu CI matrix is defined in
  `.github/workflows/paper-replication.yaml` for both R 4.4.1 and current R. The
  generic `replication-gates` branch is pushed, but GitHub cannot manually
  dispatch a workflow that is not yet present on the default branch. A pull
  request or merge is therefore required to register and trigger this matrix.
- The exact live Overleaf source has not been modified. Authors must apply the
  documented source/checklist edits and run the publisher against a clean live
  clone before submission.
