# Joint-vs-Separate Optimization Evidence

## Reviewed Sources

- The active producer is `paper/R/03-optimizer-publication-lean.R`.
- The reviewer-ready inputs are versioned under
  `paper/data/public-derived/optimizer-benchmark/`.
- Earlier exploratory scripts remain recoverable from the annotated Git tag
  `archive/pre-submission-exploratory-2026-09-08`.

## Retained Evidence

- Candidate files from balanced tau, stress search, controlled sweep, broad
  copula screen, tau-calibrated pilot, selected report, and Clayton showcase.
- By-rep and summary CSVs for fit quality, coefficient recovery, and
  joint-minus-separate deltas.
- Report figures and selected PDF/TEX reports.

## Captured Insights

- Candidate tables already contain the main review signals:
  `delta_test_log_score_per_obs_mean`, `delta_train_joint_loglik_mean`,
  RMSE deltas for `mu`, `sigma`, `theta`, and `tau`, and replicate counts.
- These summaries are sufficient for the first-stage 7-case review set without
  reopening thousands of cell `.rds` files.
- The formal paper module now mines these retained summaries and writes a
  candidate-selection CSV before running a fresh RS joint/separate comparison.
- The full all-design expanded grid was too slow for an interactive cleanup
  pass. The default expanded module 03 run now uses a practical review grid
  (`NO`, `NBI`; candidate copulas; intercept, covariate, and time-dependence
  designs), with environment variables available for heavier follow-up runs.

## Archive Decision

Keep only the compact reviewer-ready input bundle active. Exploratory scripts,
raw cell files, and smoke runs are outside the submission tree.
