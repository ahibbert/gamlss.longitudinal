# ACTG analysis-ready provenance

The response and covariates are supplied by the public, package-shipped
actg_cd4 dataset. Its source, licence, and deterministic derivation are
documented in inst/DATA-LICENSES.md and data-raw/build-public-candidate-data.R.

The covariate-selection CSV files are a compact record of the completed search
used to choose the prespecified final model. The four `screen-*.csv` files record
the broad margin/copula screens and focused six-model joint shortlist printed in
the manuscript. They contain no fitted objects, residuals, design matrices, or
private data. The standard `paper` workflow verifies and renders these ledgers
but does not repeat the 51 historical selection fits or the slow joint screen.
The `full` profile reruns the broad screens with checkpoint/resume. The final
model, family comparison, diagnostics, prediction, simulation, and
intermittent-gap sensitivity analysis are refitted from the public data in both
profiles.

The ledger was reduced from the repository's historical ACTG review bundles on
2026-09-22. The original full search is recoverable from Git history.

Canonical SHA-256 digests at promotion:

- covariate-selection-accepted.csv: ce7bb682fb2cdf48a6bdf86737b07a2e1d9f601b0025bc4d2473ee5de0c88278
- covariate-selection-targeted.csv: a5756c1cf16d6f6f7dd467da5b0b6e9ca0b033a2179dda1b82c231c14635c5a
- final-formulas.csv: c0fcec6e65ae345d388f633e0115c47957552a3ec2a1f120a981c760d8133758
- screen-margin-pooled.csv: 32ed5f1d0fbc4cc8d49af921243d45f8f65bd9e449ecaeaa2398aa94ed276013
- screen-margin-time.csv: a6135f12cdcf516e6c119480915f5a30920f99fbe09d555d61495ed31d91f6b6
- screen-copula.csv: e4d58a8458e04436e0b142e7108f65872cb870d8de4f5c10d33c1d67f69f4c4b
- screen-joint.csv: 507a24af537c44e3e6f770be1709247a1f335d4628c18c15c86715c47a4de465
