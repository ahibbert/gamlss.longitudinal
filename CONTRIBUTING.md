# Contributing And rOpenSci Standards Notes

This package is being prepared against the rOpenSci statistical software
standards for general statistical software, regression and supervised learning,
and probability distributions. The machine-readable `srr` tags are in
`R/srr-stats-standards.R`; the reviewer crosswalk is in
`inst/standards/ropensci-srr-compliance.md`.

## Package Positioning

`gamlss.longitudinal` is an R implementation of a longitudinal distributional
regression workflow that combines GAMLSS marginal models with adjacent-time
copula dependence. It is not a replacement for `gamlss`, `gamlss2`,
`VineCopula`, `geepack`, `lme4`, or `mgcv`. It builds on those ideas and, where
available, those packages, adding a longitudinal first-order copula layer,
diagnostics, simulation, inference helpers, and opt-in benchmark scaffolds.

## Glossary

- Margin: the per-time-point marginal response distribution, represented by a
  `gamlss.dist` family object.
- Copula: the dependence model connecting adjacent observed time points within
  a subject.
- `theta`: the primary copula dependence parameter. For one-parameter copulas,
  this is the only copula parameter.
- `zeta`: the second copula parameter, currently used for the t-copula degrees
  of freedom.
- Structural missingness: a subject-time combination absent from submitted data.
  The fitter expands these to explicit rows with missing responses.
- Observed data: rows with finite observed responses used in marginal
  diagnostics and regression accessors.
- Expanded data: the full subject-by-time grid after structural missing rows are
  inserted.
- `newdata`: a user-supplied prediction panel. It may contain future rows or
  new subject identifiers, but categorical levels must have appeared in the
  fitting data.
- Forecast panel: a `newdata` panel used for extrapolation or future-time
  prediction. The package treats this as prediction from supplied covariates,
  not as a separate time-series forecasting model.

## Input And Missingness Policy

The main fitting functions expect a long-format data frame with one row per
observed subject-time combination. Scalar controls must be length-one values.
Formula inputs may be formulas or single strings. String-valued option arguments
are case-sensitive unless documented otherwise; constrained options use
`match.arg()` or an equivalent validation path.

Supported predictor columns are numeric, integer, logical, factor, ordered
factor, or character. Character predictors are treated as unordered categorical
variables; convert them to factors before fitting to control level order.
Ordered factors are encoded with treatment contrasts rather than polynomial
contrasts. List-columns and matrix/data-frame columns are rejected.

Time columns may be numeric, integer, numeric-like character, or factor.
Numeric-like character time is converted to numeric with a warning. Categorical
visit labels should be supplied as factors.

Response values may be finite numeric values or `NA`. Response `NA` values
represent missing outcomes. Response `NaN`, `Inf`, and `-Inf` are rejected
before fitting. Predictor values in submitted rows must be observed and finite;
predictor `NA`, `NaN`, `Inf`, and `-Inf` values are rejected before fitting.
Structurally inserted rows may contain missing responses and internally filled
predictor proxies for model-matrix construction.

The package does not perform statistical imputation. If imputation is required,
do it explicitly before fitting and report the imputation model separately.

## Numerical And Distribution Methods

Marginal distribution operations are delegated to `gamlss.dist` family
functions. The package uses the `d*`, `p*`, and `q*` family interfaces for
density, cumulative probability, quantile prediction, simulation, and
diagnostics.

Copula calculations use the native backend by default, with `VineCopula`
available as a parity/reference backend when installed. Gaussian, Clayton,
Frank, Gumbel, Joe, and t-copulas are supported. The t-copula CDF uses
one-dimensional numerical integration of the conditional representation with a
relative tolerance of `1e-7`; inputs are clamped away from 0 and 1 before
calculation to avoid undefined tail evaluations.

Variance-covariance calculations can use analytical/semi-analytical Hessian
blocks or finite-difference Hessians. Analytical paths fall back to numerical
reference calculations when inversion fails or unsupported regions are detected.
Discrete/count margins use rectangle probability calculations where exact
integer response mass is required; count-tail tests should compare CDF/quantile
consistency rather than approximating integrals with unbounded summation.

Floating-point comparisons in package code should use tolerances unless values
are guaranteed to be counts, dimensions, indices, logical flags, or categorical
codes.

## Extended Tests

Routine tests must remain CRAN- and CI-friendly. Long-running recovery,
benchmark, and stress tests are opt-in:

```sh
Rscript -e "source('inst/smoke-tests/new-user-smoke.R')"
Rscript -e "pkgload::load_all(); testthat::test_file('tests/testthat/test-p0-user-facing-estimands.R')"
GAMLSS_LONGITUDINAL_EXTENDED_TESTS=true Rscript -e "pkgload::load_all(); testthat::test_dir('tests/testthat')"
```

Extended tests may cover multi-seed parameter recovery, larger panels,
near-independent and strong dependence, missing visits, benchmark comparators,
and platform-specific numerical tolerance checks. If an extended test requires
optional packages or external artifacts that are unavailable, it should skip
with an informative diagnostic rather than fail.
