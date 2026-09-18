# Correctness test plan

Status: pre-implementation acceptance contract.

No estimator-performance simulation can substitute for these likelihood and
indexing checks. Tests should use small fixed examples, independent reference
calculations, and explicit tolerances.

## A. Vine topology and pooled predictors

### VINE-001: full edge enumeration

For T = 2, 3, 4, and 5, verify exactly T(T-1)/2 edges, T-ell edges in tree ell,
and the expected conditioning sets. Verify no edge joins subjects.

### VINE-002: input-order invariance

Permuting long-format rows must leave topology, likelihood, and fitted
coefficients unchanged after restoring canonical row order.

### VINE-003: irregular time

Use visits at times 0, 1, 4, and 10. Verify scheduled tree membership separately
from elapsed `time_distance` values 1, 3, 6, 4, 9, and 10 on the appropriate
edges.

### VINE-004: generated endpoint covariates

Verify `left_*`, `right_*`, time variables, tree, order distance, and
conditioning size for every edge. Reject a supposedly subject-constant
covariate that varies within subject.

### VINE-005: pooled formula mapping

Construct `~ tree + s(time_distance) + left_x + right_x`; compare its matrix to
an independently constructed edge data frame. Verify one coefficient vector is
used across every higher tree.

## B. Full outcome density

### MATH-001: h-function orientation

For every initially supported copula family, compare both conditional CDF
orientations against numerical derivatives of the copula CDF.

### MATH-002: conditional inverse

Apply h followed by inverse-h at interior and tail probability values. Include
positive and negative dependence wherever the family permits it.

### MATH-003: density normalization

Integrate two-, three-, and four-visit continuous vine densities over the unit
hypercube. Verify total mass one and uniform one-dimensional margins.

### MATH-004: Gaussian oracle

Map Gaussian vine partial correlations to an ordinary correlation matrix.
Compare joint log densities, conditional CDFs, and simulations against
`mvtnorm` for fixed three- and four-visit cases.

### MATH-005: first-order collapse

Set every higher-tree edge to independence. Compare the new complete-data
likelihood, score, predictions, and simulated serial summaries with the current
first-order implementation.

### MATH-006: boundary stability

Evaluate conditional PITs and log densities near 0 and 1 without returning
spurious `NaN`, silent probability clipping, or finite values for invalid
parameters.

## C. Observation process

### OBS-001: formula availability

Accept fully observed exogenous variables and generated recorded-history
variables. Reject raw current response terms and covariates unavailable on
missing scheduled rows.

### OBS-002: logistic matrix and smooth

Compare the observation fixed and smooth matrices with an independent `mgcv`
construction. Verify link inversion and penalty dimensions.

### OBS-003: independence limit

At independence, verify q(s) = pi for a dense grid of s and pi. Do not equate a
native zero parameter with independence unless the family definition supports
it.

### OBS-004: Bernoulli normalization

For every s and pi grid point, verify q(s) + {1-q(s)} = 1, finite log
probabilities, and valid limiting behavior.

### OBS-005: innovation calibration

Numerically integrate q(w) over w in [0,1] and recover pi. Repeat conditionally
over fixed recorded histories generated from the outcome model.

### OBS-006: connection distinction

Under nonzero serial dependence, exhibit a history for which outcome-level and
innovation-level selection probabilities differ. Verify equality in a declared
special case such as longitudinal independence.

### OBS-007: risk-set semantics

Verify Bernoulli factors at every scheduled visit for intermittent missingness,
only through the first dropout for dropout, and according to the supplied
at-risk variable for mixed processes. Reject non-monotone data declared as pure
dropout.

## D. Observed-data likelihood

### LIK-001: no-missing reduction

With every R = 1, the observed-data likelihood must equal the complete-data
outcome density times the observation probabilities.

### LIK-002: exhaustive pattern normalization

For T = 3 and 4, integrate over all outcomes and sum over every allowed
observation pattern. Recover total probability one for intermittent and dropout
risk rules separately.

### LIK-003: margin preservation

After summing over observation patterns, recover the specified full-data
GAMLSS marginal distribution at every visit.

### LIK-004: one-gap reference

For observed-missing-observed data at T = 3, compare the production likelihood
with independent high-accuracy one-dimensional integration on both response and
uniform scales.

### LIK-005: multi-gap reference

For T = 4, cover every one- and two-missing pattern and compare deterministic
quadrature with an independent nested-integration or high-precision reference.

### LIK-006: innovation-dropout reduction

Compare the analytic terminal-dropout contribution with brute-force
integration. Verify that future visits add neither observation factors nor
integration error.

### LIK-007: outcome-dropout integral

Show by reference calculation that outcome-level terminal dropout generally
requires integration and that the innovation shortcut is not applied.

### LIK-008: MAR interior gap

At selection independence, compare the correct interior-gap marginal
likelihood with the current segmented approximation and with brute-force
integration. The correct implementation must not silently reconnect observed
visits or declare separated segments independent.

### LIK-009: integration audit

Repeat every numerical integral at higher resolution. Verify the stored method,
dimension, node count, error diagnostic, and non-finite evaluation count.

## E. Joint objective and fitting

### FIT-001: parameter packing

Round-trip all margin, tree-1, higher-tree, observation, and smooth coefficients
through pack/unpack operations without reordering.

### FIT-002: score reference

Compare every parameter block with high-accuracy central finite differences on
complete, intermittent, and dropout examples.

### FIT-003: one objective

Instrument the fitter to verify that reported estimates optimize the combined
observed-data objective. Sequential fits may initialize but must not replace
the joint fit.

### FIT-004: fixed-sensitivity refitting

Across a tau grid, verify that tau is fixed while every estimated parameter
block is re-optimized. Warm starts must not cause path-dependent final answers.

### FIT-005: small recovery

On deterministic large-sample or nearly noise-free fixtures, recover known
parameters for complete, MAR, innovation-MNAR, and outcome-MNAR cases.

### FIT-006: reproducibility

Fixed data, starting values, and deterministic quadrature must give identical
objective values and estimates across repeated fits.

## F. User-visible demonstrator

### DEMO-001: paired data provenance

Retain a replicate identifier connecting the complete generated data, imposed
missing data, oracle fit, unadjusted fit, and every adjusted fit.

### DEMO-002: visible adjustment

For a predeclared MNAR scenario, report location, scale, and upper-quantile
truth and errors. The report must display rather than hide cases where the
adjustment fails or the wrong connection performs poorly.

### DEMO-003: operational evidence

Report attempted and retained fits, convergence, boundary events, integration
diagnostics, warning classes, and runtime for every method.

### DEMO-004: non-circular stress case

Include at least one generator outside the fitted selection connection so the
demonstrator is not exclusively a fitted-model self-recovery exercise.

## G. Future discrete capability boundary

### DISC-001: explicit rejection

Until implemented, requesting the new full-vine observed-data model with a
discrete margin must fail with a targeted capability message rather than using
a continuous density approximation.

### DISC-002: stable interface

The specification objects, risk-set definitions, and sensitivity semantics
must not assume continuity except inside the likelihood backend. This is the
contract for adding exact discrete support later.

