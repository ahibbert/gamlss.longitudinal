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

## Submission gates not yet executed

- The `full` public Monte Carlo profile was started from an empty store and
  checkpoint set on 2026-08-27 and is currently running. The first restart was
  an intentional worker-scaling check; valid atomic checkpoints were retained.
  The published 20-replicate missingness study completed in 51 minutes 22
  seconds. The active run is continuing as a hidden 12-worker process with
  stdout/stderr under `results/jss-replication/full/logs/`. End-to-end
  completion and the subsequent no-repeat resume run remain required.
- The metric-level comparison is now wired into the full graph and its complete
  mirrored-layout test passed 2,275 comparisons. Acceptance against freshly
  recomputed results remains pending completion of the active full run.
- The Windows/Ubuntu CI matrix is defined in
  `.github/workflows/paper-replication.yaml` for both R 4.4.1 and current R, but
  remote CI results do not exist until the changes are committed and pushed.
- The exact live Overleaf source has not been modified. Authors must apply the
  documented source/checklist edits and run the publisher against a clean live
  clone before submission.
