# ACTG shortlist full review

Open `index.html` for the comparison dashboard. Each numbered subdirectory is a
self-contained full HTML report with subject-clustered sandwich SEs, confidence
intervals, coefficient tables, diagnostics, and term plots.

- `01-preserved`: current accepted model.
- `02-mu-baseline-visit1`: adds a visit-1 baseline-CD4 effect to `mu`.
- `03-mu-baseline-visit1-nu-visit1`: also adds a visit-1 Sichel-shape effect.

The original accepted review and time-1 investigation outputs were not overwritten.

The selected model can be refitted directly with `fit-final-model.R`. By default,
that script writes `final-model-refit.rds` without overwriting the reviewed model.
