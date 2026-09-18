# ACTG visit-1 final polishing screen

Generated 2026-09-10. The accepted shortlist reports and model objects were not
overwritten.

## Decision

Retain shortlist model 3, the ZISICHEL--Gumbel model with `nu ~ visit1`.
It remains the best model by both AIC and BIC. None of the five narrowly targeted
additions improves fit enough to pay for its extra complexity.

The remaining visit-1 worm curvature is modest. Its central fitted curve is mostly
inside the pointwise envelope, and the largest departures are in the sparsely
populated tails. The repeat-seed PIT summaries confirm a residual under-dispersion
signal, but not one for which the tested extensions provide a compelling solution.

## Polishing comparison

The uniform-PIT target SD is 0.2887. Lower CvM and worm RMS values indicate closer
calibration.

| Model | Delta AIC | Delta BIC | Visit-1 PIT SD | Visit-1 CvM | Visit-1 worm RMS |
|---|---:|---:|---:|---:|---:|
| Retained: `nu ~ visit1` | 0.000 | 0.000 | 0.2746 | 0.4928 | 0.0848 |
| Add `visit1` to `tau` | +1.184 | +7.300 | 0.2754 | 0.4912 | 0.0825 |
| Add `visit1 + baseline_cd4:visit1` to `tau` | +3.107 | +15.340 | 0.2754 | 0.4936 | 0.0827 |
| Add `visit1 + baseline_cd4:visit1` to `sigma` | +3.439 | +15.546 | 0.2769 | 0.3967 | 0.0791 |
| Add `I(baseline_cd4^2):visit1` to `mu` | +3.461 | +9.500 | 0.2751 | 0.4923 | 0.0846 |
| Add `visit1` to `sigma` | +7.206 | +13.337 | 0.2760 | 0.4499 | 0.0848 |

All six models converged. Metrics are averages over 20 randomized-PIT seeds.

The `sigma` model with a visit-1 baseline-CD4 interaction produces the largest
numerical diagnostic improvement, reducing CvM by about 19.5% and worm RMS by about
6.7%. However, it is worse by 3.44 AIC and 15.55 BIC points and adds two coefficients.
This is not a good trade-off for the primary model. The one-coefficient `tau ~ . +
visit1` alternative changes the diagnostics only trivially while worsening both
information criteria.

## What has already been ruled out

- Visit-1 treatment, sex, and age effects did not explain the localized residual
  signal in the earlier covariate screen.
- A visit-1-specific quadratic baseline-CD4 mean effect did not improve AIC or BIC.
- Additional visit-1 effects in both remaining ZISICHEL margin parameters (`sigma`
  and `tau`) did not improve AIC or BIC.
- More elaborate `nu` formulas can lower AIC, but worsen BIC materially and were
  visually inferior to the parsimonious `nu ~ visit1` model.

## Interpretation

The residual visit-1 departure is therefore better documented as a small limitation
than removed with another post-selection term. Further pursuit would require a
materially different model class, such as a separately specified first-visit margin
or an additional mixture component. That would change the scientific model rather
than constitute a routine covariate refinement, and would carry a substantial
overfitting risk when motivated by a single diagnostic panel.

## Reproducibility files

- `polish-results.csv`: model formulas, AIC/BIC, and visit-1 diagnostic comparison.
- `all-candidate-diagnostics-by-time.csv`: PIT and worm summaries at every visit.
- `checkpoints/`: the retained model and all five fitted alternatives.
- `aic-best-polish-fit.rds`: a copy of the retained `nu ~ visit1` fit.
- `run-time1-polish-screen.R` in the parent folder: complete fitting script.
