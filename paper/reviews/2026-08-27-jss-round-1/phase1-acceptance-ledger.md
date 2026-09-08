# Phase 1 acceptance ledger

Status: Phase 1 passed on 2026-09-01 after independent review.

This ledger maps the original Phase 1 gates to executable evidence.
The approved scope decision for segmented likelihoods is retained: AIC and BIC
are available under the stated between-segment independence assumption, while
ordinary model-Hessian inference is unavailable.

| Gate | Implemented contract | Primary acceptance evidence |
|---|---|---|
| JSS-001 likelihood validity | Observation masks are fixed by the data; invalid included contributions invalidate the objective; exact discrete rectangles use stable log-scale tail evaluation and an independent Gaussian fast path. | `test-p0-likelihood-validity.R`, `test-p0-discrete-rectangle-likelihood.R`, `test-p0-likelihood-evaluation-margin.R` |
| JSS-002 Hessian inference | Inference requires named, finite, symmetric curvature; positive-definite full-rank information; bounded scaled condition number; a checked fitted gradient; method agreement when requested; and positive covariance diagonals. Failures use `gamlss_longitudinal_inference_unavailable`. | `test-p0-inference-validation.R`, including equality and just-over-boundary tests for symmetry, rank, condition, gradient, and method agreement; `test-p0-inference-contract-adversarial.R` |
| JSS-003/JSS-018 missingness | Ordinary likelihood is restricted to complete panels and observed prefixes ending in dropout. Intermittent gaps and delayed entry error by default. Explicit `missingness = "segment"` records segment counts and independence, keeps AIC/BIC under that likelihood, and blocks ordinary Hessian inference. Scheduled order, duplicate rows, irregular spacing, adjacency, and left-row pair covariates are documented. | `test-p0-missingness-contract.R`, including the six-row panel; `test-p0-model-preprocess.R`; `test-p0-phase1-documentation.R` |
| JSS-030 capability scope | A versioned registry defines supported homogeneous margin, copula, likelihood, Hessian, and diagnostic routes; unsupported combinations fail before optimization with named conditions. | `test-p0-capability-registry.R`, `test-p2-srr-input-policy.R`, `test-p2-srr-error-map.R` |
| JSS-009/JSS-031 methods and convergence | `logLik()` is scalar and classed with observed `nobs` and joint EDF; AIC/BIC and stored/summary criteria share the observed-response denominator; `vcov()` is coefficient-aligned; nonconvergence gates inference; every optimizer termination route is named. | `test-p0-standard-generics.R`, `test-p0-model-fit-result.R`, `test-p0-optimizer-control-convergence.R`, `test-p0-model-vcov.R` |
| JSS-010 diagnostics | Supported discrete margins use seeded randomized PIT by default. PIT, tail, and lagged residual summaries are descriptive, higher unavailable lags are recorded, and only a user-supplied cutoff creates a review flag. | `test-p0-model-check-thresholds.R`, including a correctly specified Poisson calibration test; `test-p0-phase1-documentation.R` |
| JSS-032 inference scope | A versioned registry covers fixed analytical/numerical Hessian, sandwich, smooth, parametric and cluster bootstrap, Wald/LR, interval, prediction, term-plot, marginal-effect, and publication outputs, including coefficient blocks, conditioning, omissions, approximations, and failure states. Manuscript claims are limited to this conditional scope. Empirical conditional-coverage assessment is explicitly deferred to the attempt-level Phase 2 recovery study and is not claimed as Phase 1 evidence. | `test-p0-inference-contracts.R`, `test-p0-inference-contract-adversarial.R`, `test-p0-user-facing-estimands.R`, `test-p0-phase1-documentation.R` |

## Full-suite result

Command:

```r
devtools::test(stop_on_failure = FALSE)
```

Result on Windows, 2026-09-01: **5,088 passed, 0 failed, 0 warnings, 2
intentional skips**, duration 252.8 seconds. The skipped tests are the two
publication-scale extended SRR checks guarded by
`GAMLSS_LONGITUDINAL_EXTENDED_TESTS=true`; they are not Phase 1 acceptance
tests.

## Expected numerical and console events

| Event | Scope | Disposition |
|---|---|---|
| `boundary (singular) fit` messages from `lme4` | Deliberately small comparator fixtures in `test-p1-benchmark-comparators.R` | Expected third-party comparator message; it is not used as evidence that package inference is available and does not produce a test warning or failure. |
| Smoothing-parameter progress text | Fit/plot integration fixtures | Expected progress output from exercised fitting paths; no warning or failure. |
| Publication-suite dry-run notices | Multivariate replication planning tests | Expected confirmation that publication-scale work is opt-in; no simulation is silently launched. |

Any new warning, failed assertion, unclassified optimizer stop reason, or
unlisted numerical event fails this ledger and requires re-review.

The final independent re-audit rated JSS-001, JSS-002, JSS-003/JSS-018,
JSS-030, JSS-009/JSS-031, JSS-010, and JSS-032 as PASS. Its affected-test rerun
reported 156 passes with no failures, warnings, or skips. The complete suite
above was then rerun after the final audit corrections.
