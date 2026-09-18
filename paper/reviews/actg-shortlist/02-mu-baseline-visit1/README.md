# Final ACTG ZISICHEL–Gumbel model review bundle

Generated: 2026-09-10 13:51:23 AEST

## Final formulas

- `mu`: `cd4 ~ treatment + age + s(baseline_cd4, bs = "ps", k = 8) + week_factor + baseline_cd4:visit1`
- `sigma`: `~treatment + baseline_cd4`
- `nu`: `~1`
- `tau`: `~treatment + sex + baseline_cd4`
- `theta`: `~treatment + s(baseline_cd4, bs = "ps", k = 8) + week_factor`

## Start here

- `review-index.html`: compact browsable summary.
- `tables/coefficient-table-publication.csv`: formatted publication coefficient table.
- `tables/coefficient-table-numeric.csv`: unrounded numeric coefficient results.
- `tables/model-fit-statistics.csv`: likelihood, information criteria, sample sizes, and convergence.
- `tables/smooth-terms-edf.csv`: retained smooth terms and effective degrees of freedom.
- `model-selection/all-candidate-models.csv`: every candidate evaluated in all four staged searches.
- `model-selection/accepted-selection-paths.csv`: the accepted greedy AIC path.
- `figures/`: marginal, residual, copula, calibration, and term-effect diagnostics in PNG and PDF.
- `logs/artifact-status.csv`: success, warnings, runtime, and any failure for every output step.

## Interpretation notes

- The primary coefficient uncertainty uses subject-clustered sandwich standard errors. Model-based standard errors are included as a sensitivity comparison.
- Confidence intervals and p-values are conditional on the final selected model and do not account for the preceding model-selection search.
- AIC/BIC use the stated segmented likelihood: observations in separate contiguous segments for a subject are treated as independent.
- Discrete-margin PIT diagnostics use seeded randomized PIT values (seed 20260910).
- Penalized smooths are summarized by EDF and effect plots; basis coefficients are not suitable as a publication coefficient table.
- Diagnostic thresholds are descriptive aids, not automatic goodness-of-fit verdicts.
