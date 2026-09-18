# ACTG final margin sensitivity screen

Generated 2026-09-16 from the full ACTG dataset. The accepted model and its
review bundle were not overwritten.

## Decision

Retain the selected ZISICHEL--Gumbel model.

ZIPIG is ruled out under the final covariate structure: it is worse by 64.94 AIC
points and 53.28 BIC points, and its visit-specific PIT and worm diagnostics are
also worse.

ZIBNB has a modest likelihood/information-criterion advantage, improving AIC by
5.96 and BIC by 6.61 points. However, its calibration is materially inferior at
all four visits. At visit 1, its PIT SD is 0.2591 rather than 0.2745 for ZISICHEL
(uniform target 0.2887), CvM is 2.168 rather than 0.495, and worm RMS is 0.1649
rather than 0.0851. Its PIT mean is below 0.48 at every visit, showing a systematic
location bias not present to the same degree under ZISICHEL. The diagnostic loss
is too large to accept for a six-point AIC/BIC improvement.

## Gumbel comparison

| Margin | AIC | BIC | Delta AIC | Delta BIC | Visit-1 PIT SD | Visit-1 CvM | Visit-1 worm RMS |
|---|---:|---:|---:|---:|---:|---:|---:|
| ZIBNB | 26561.21 | 26767.90 | -5.96 | -6.61 | 0.2591 | 2.1684 | 0.1649 |
| **ZISICHEL** | 26567.17 | 26774.51 | 0.00 | 0.00 | **0.2745** | **0.4945** | **0.0851** |
| ZIPIG | 26632.11 | 26827.79 | +64.94 | +53.28 | 0.2576 | 1.7110 | 0.1557 |

All three Gumbel models converged. The same selected covariate structure was used,
with parameters mapped by role. In ZIPIG, `nu` is the zero-inflation parameter;
the selected ZISICHEL `tau` predictor was therefore transferred to ZIPIG `nu`.
ZIBNB and ZISICHEL both use `nu` for count shape and `tau` for zero inflation, so
their formulas transfer directly.

## Student-t sensitivity status

Student-t was not accepted as a fitted alternative:

- ZIPIG--Student-t repeatedly terminated its R worker immediately when the
  full-data joint fit began.
- ZIBNB--Student-t ran stably on one core for more than two hours but did not
  return a fit or checkpoint. It was stopped after exceeding the configured
  two-hour cap inside a long inner calculation.
- ZISICHEL--Student-t was not pursued after these failures. The earlier full-model
  copula screen likewise produced completed results for Gumbel, Gaussian, Frank,
  Clayton, and Joe but no completed Student-t fit.

These are computationally unresolved sensitivity models, not evidence that they
fit better or worse. They do not alter the final margin decision under the selected
Gumbel dependence structure.

## Reproducibility files

- `completed-family-copula-results.csv`: information-criterion comparison.
- `completed-diagnostics-by-time.csv`: 20-seed PIT and worm summaries by visit.
- `checkpoints/`: complete ZISICHEL--Gumbel, ZIBNB--Gumbel, and ZIPIG--Gumbel fits.
- `objects/`: extracted fitted model objects used by the visual report.
- `index.html`: side-by-side PIT, worm, and rootogram review.
- `../run-final-family-copula-screen.R`: checkpointed fitting code.
