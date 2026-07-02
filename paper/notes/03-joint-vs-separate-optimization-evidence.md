# Joint-vs-Separate Optimization Evidence

## Reviewed Sources

- Scripts moved to `R (testing)/03-joint-vs-separate-optimization/`.
- Latest retained results live under
  `results/jss-exploratory/03-joint-vs-separate-optimization/`.
- Smoke directories and cell-level raw files archived under
  `results/_archive/testing-cleanup-20260612/results/rs_joint_*`.

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

Keep current summary/candidate/report artifacts active. Raw cell files and smoke
runs are archived, not deleted, until the review set is accepted.
