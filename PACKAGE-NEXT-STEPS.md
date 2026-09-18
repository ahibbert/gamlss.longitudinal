# Package next steps

This is a collection log for requested package changes. The 2026-09-09 items
were implemented together on 2026-09-09; add later requests as new entries.

## Collected 2026-09-09

### NS-001 — More visible progress by default (implemented)

Affected user-facing functions:

- `select_copula()`
- `select_margin()`
- `gamlss_longitudinal()`

Requested behaviour:

- The default call should emit a modest amount of useful progress output so a
  user can tell that a long-running fit or selection is still advancing.
- `verbose = 0` must suppress progress output completely.
- Output should be stage-level or otherwise throttled rather than overwhelming.
- The three functions should use a consistent verbosity convention, while
  preserving any distinct purpose served by `select_margin(trace = ...)`.

Implementation notes:

- All three functions now default to `verbose = 1` and use the shared
  numeric/logical verbosity normalization.
- `select_margin(trace = ...)` continues to control the underlying GAMLSS
  trace independently of package-level progress.
- Level 1 reports stages and, for model fits and copula screening, iterative
  progress; levels 2 and 3 retain more detailed optimizer output.

Acceptance checks:

- Default calls show clear progress on a deliberately slow example.
- Explicit `verbose = 0` is silent for all three entry points, including nested
  marginal fits and optimizer output.
- Existing higher verbosity levels continue to provide increasingly detailed
  diagnostics.
- Documentation and tests describe the shared verbosity levels and the
  relationship between `verbose` and `trace`.

### NS-002 — Reduce overplotting for high-count discrete margins (implemented)

Affected function: `plot_margin_fit()`.

Observed problem:

- Discrete margins currently draw a unit-width bar and fitted probability mass
  at every integer in the observed response range.
- For a high-count ZIPIG response, the many narrow bins obscure the overall
  observed-versus-fitted trend. The supplied example has four time facets
  (8, 16, 24, and 32) and a long sparse tail extending to roughly 500.

Requested behaviour:

- Put a sensible cap on the number of displayed bins for discrete responses.
- When the discrete support exceeds that cap, aggregate adjacent response
  values into interpretable bins for both the observations and fitted
  probabilities.
- Keep unit-width integer bins for smaller discrete ranges.
- Do not change continuous-margin histogram behaviour.

Implementation notes:

- The existing `bins` argument is now the discrete maximum as well as the
  continuous histogram count.
- Discrete facets share equal-width integer boundaries based on their combined
  displayed support.
- Observed proportions and fitted probability masses are summed over the same
  bins without truncating the sparse tail.

Acceptance checks:

- A ZIPIG example like the supplied plot has a readable number of bars and a
  visible fitted trend.
- Observed and fitted probabilities are aggregated over exactly the same bin
  boundaries.
- Displayed fitted mass is conserved apart from any explicitly documented tail
  treatment.
- Low-range discrete examples retain one bar per integer.
- Faceted plots remain comparable across time.

### NS-003 — Do not recompute inference when calling `summary()` (implemented)

Affected flow: `gamlss_longitudinal()` fit finalization and
`summary.gamlss.longitudinal()`.

Observed problem:

- Calling `summary(fit)` after a fit appears to recalculate an expensive result
  that should already have been computed and stored during fitting.
- Suspected regression area: recent segmented-likelihood warning/inference
  routing changes.
- `summary.gamlss.longitudinal()` continues to call `.resolve_vcov()` when
  `include_vcov = TRUE`; compatible requests now resolve to fit-time cached
  output for both ordinary and segmented fits.

Reproduction supplied by the requester:

```r
fit_actg <- gamlss_longitudinal(
  data = actg_patient_sample,
  time_var = "week_factor",
  subject_var = "subject",
  margin_dist = ZIPIG(),
  copula_dist = "G",
  mu.formula = cd4 ~ treatment + s(age) + sex + baseline_cd4 +
    week_factor + treatment * week_factor,
  verbose = 3
)

summary(fit_actg)
```

Requested behaviour:

- Compute and retain the standard summary/inference ingredients as part of the
  fit when they were requested during fitting.
- `summary()` should reuse those stored results and return promptly rather than
  repeating the expensive calculation.
- Segmented-likelihood warnings and covariance-method routing must not cause a
  previously completed calculation to be discarded or repeated.

Implementation notes:

- Fit finalization no longer marks all segmented model-Hessian requests as
  unavailable before attempting inference. It lets the existing inference
  router compute subject-clustered sandwich covariance and caches that result
  under the originally requested method.
- The general segmented-fit warning already explains this routing, so the
  redundant fit-finalization routing warning is muffled. Direct later requests
  for an incompatible method retain their warning.

Acceptance checks:

- Instrumented tests show no second covariance/Hessian/sandwich calculation for
  `summary()` after a fit created with inference enabled.
- Cached results are used only when compatible with the requested summary
  options (`numderiv`, `vcov_method`, confidence level, and extra arguments).
- A fit made with inference disabled has an explicit, documented on-demand
  summary path rather than silently pretending a cache exists.
- Ordinary and segmented-likelihood fits both preserve their intended warnings,
  covariance method, and inference metadata.
- Repeated compatible `summary()` calls are deterministic and do not refit.

### NS-004 — Report model-based and sandwich standard errors together (implemented)

Affected flow: sandwich covariance output and `summary.gamlss.longitudinal()`.

Requested behaviour:

- Follow the useful part of the `gee` presentation by showing model-based and
  subject-clustered robust standard errors alongside one another whenever
  sandwich covariance is calculated.
- Keep the robust covariance as the primary inferential result: tests,
  p-values, confidence intervals, and the existing `std_error` field must use
  the robust standard error.
- Label the model-based values as assumption-dependent comparison values, not
  as a second equally preferred inferential result.
- Reuse the inverse-Hessian sandwich bread rather than calculating another
  Hessian.

Implementation notes:

- Detailed sandwich covariance output retains `model_based` and
  `cluster_robust` covariance and standard-error entries; `overall` remains the
  cluster-robust result for compatibility.
- Summary coefficient data include `model_std_error` and `robust_std_error`,
  while `std_error` remains the robust value.
- Console summaries show `model_se` and `robust_se` columns and explicitly say
  that tests and p-values use the robust values.
