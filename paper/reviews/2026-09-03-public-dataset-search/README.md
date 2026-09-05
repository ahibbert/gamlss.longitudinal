# Public longitudinal data search

Date: 2026-09-03

## Purpose

This audit searches for public data suitable for the main worked application
and/or the simulated workflow example. The preferred data should:

1. contain at least four repeated observations per subject;
2. have enough subjects and strong within-subject dependence for the copula
   component to be visible;
3. favour a supported non-exponential-family GAMLSS margin; and
4. support a non-time covariate in both a non-location marginal parameter and
   the copula parameter.

Information criteria from independent-row GAMLSS fits are screening devices,
not final model comparisons. A candidate is promoted only after the effects
remain defensible in a joint `gamlss.longitudinal` fit.

## Packaging and licensing decision

All three candidates now have compact analysis-ready objects under `data/`,
with reproducible construction in `data-raw/build-public-candidate-data.R` and
installed attribution in `inst/DATA-LICENSES.md`.

- The PLOS workbook is explicitly licensed CC BY 4.0. The packaged derivative
  removes original study identifiers, records all modifications, cites the
  authors/article/dataset DOI, and preserves a link to the license.
- `PatentsRDUS` is distributed as part of pglm 0.2-4 under GPL >= 2. The
  packaged derivative retains that source/license attribution and removes the
  original CUSIP identifier.
- `pbcseq` is distributed as part of survival 3.8-11 under LGPL >= 2, with
  copyright attributed to the Mayo Foundation for Medical Education and
  Research. The packaged derivative retains that notice and replaces patient
  identifiers.

The source-package licenses are compatible with this package's GPL-3 license.
Neither pglm nor survival states a separate, more restrictive data license in
its package metadata. For maximum submission conservatism, the maintainers can
still be asked to confirm that the package license was intended to cover
downstream redistribution of derived data; if either declines, the same
analysis can load the upstream package object at runtime without changing the
model code.

Candidate screens and fits are kept outside the default replication graph:

```text
Rscript paper/run-public-application.R --dataset patents --stage screen
Rscript paper/run-public-application.R --dataset patents --stage fit
Rscript paper/run-public-application.R --dataset patents --stage fit --comparators
Rscript paper/run-public-application.R --dataset patents --stage fit --outputs
Rscript paper/run-public-application.R --dataset patents --stage refine --outputs
```

The `screen` stage follows the vignette workflow directly: `select_margin()`
uses time-specific intercepts, then `select_copula()` screens registered routes
conditional on the retained margin. Both complete selection tables, excluded
routes, and the reason for retaining anything other than the first-ranked
screening result are saved. The `--outputs` command runs the full diagnostic
and publication-output bundle. The `refine` stage runs the covariate-matched
flexible-margin screen and refined patents model; both are deliberately opt-in
because they are more expensive than fitting and saving the base model.
Final application fits compute covariance estimates and validate inference by
default. `--no-vcov` is available only for exploratory timing runs whose saved
fits are not suitable for reporting.

## Candidate ranking

### 1. PatentsRDUS: substantive count candidate

`pglm::PatentsRDUS` contains 346 US firms observed annually from 1970 to 1979
(3,460 observations). The outcome is the number of patent applications in a
year that were eventually granted. It is a balanced panel with no missing
outcomes.

Key audit results:

- exact non-negative integer outcome, range 0--608;
- mean 36.28, variance 5,545.13, and 17.49% zero counts;
- raw adjacent Spearman correlation 0.921;
- residual adjacent Spearman correlation 0.719 after adjustment for year,
  annual R&D expenditure, scientific-sector status, and 1972 capital;
- residual correlation 0.575 below versus 0.874 above median baseline capital;
- Delaporte and NBI were retained by the standard `select_margin()` screen,
  with Delaporte ranked first and NBI second; and
- NBI dispersion depending on scientific-sector status and baseline capital
  strongly preferred to constant dispersion.

Screening specification:

```r
mu.formula    = patents ~ factor(year) + log(rd) + scisect + capital_z
sigma.formula = ~ scisect + capital_z
theta.formula = ~ capital_z
```

where `capital_z` is standardized `log(capital72)`. `sumpat` must not be used as
a predictor because it is derived from the repeated outcome.

Earlier covariate-adjusted independent-margin screening results:

| Model | AIC | BIC |
|---|---:|---:|
| Poisson | 68,643.8 | 68,723.8 |
| NBI, constant sigma | 23,572.6 | 23,658.7 |
| NBI, sigma by scientific sector | 23,502.5 | 23,594.8 |
| NBI, sigma by baseline capital | 23,500.9 | 23,593.1 |
| NBI, sigma by sector and capital | **23,386.4** | **23,484.8** |
| Delaporte, constant sigma and nu | 23,564.1 | 23,656.4 |

After correcting the fitted-score calculation and tightening the RS inner
tolerance, the accepted base joint NBI/Clayton model has log likelihood
-10,417.32, AIC 20,870.65, and BIC 20,981.33. Analytical-Hessian inference is
available; the curvature-standardized fitted score is 0.078.

The standardized vignette screen gives Delaporte AIC 28,270.28 and NBI AIC
28,428.51 after time-specific intercepts. However, the exact joint Delaporte
rectangle fit did not complete an outer iteration inside the routine
application budget. NBI is therefore explicitly retained from the selected
shortlist. Conditional on NBI, `select_copula()` ranks Gumbel before Clayton,
but the full NBI/Gumbel patents fit likewise remained in its initial expensive
optimization phase after 12 minutes. The reproducible patents fit consequently
retains Clayton as the second screened, computationally viable route. These
choices are stored in `distribution-selection.rds`; neither is presented as an
unqualified AIC winner.

A second screen showed that log R&D belongs in the dispersion equation. The
refined joint model

```r
sigma.formula = ~ scientific + capital_z + log_rd
```

has log likelihood -10,374.06, AIC 20,786.12, and BIC 20,902.95. The added
log-R&D dispersion coefficient is -0.234, and inference validation passes with
a standardized fitted score of 0.067. This is the patents model to use if that
application is retained, subject to the documented selection caveat above.

In the covariate-matched independent-margin screen, Delaporte and NBI with
R&D-dependent dispersion are effectively tied by AIC (23,328.96 versus
23,329.15), while BIC prefers NBI (23,433.68 versus 23,439.64). A full joint
Delaporte fit is currently too slow to make part of the routine candidate
workflow and is deferred unless its negligible independent AIC advantage
becomes important.

The packaged derivative retains the upstream GPL >= 2 attribution in
`inst/DATA-LICENSES.md`; downstream redistribution can be replaced by runtime
loading from `pglm` if a final licensing confirmation requires it.

Sources:

- <https://stat.ethz.ch/CRAN/web/packages/pglm/refman/pglm.html>
- Hall, Griliches, and Hausman (1986), *Patents and R&D: Is There a Lag?*

### 2. Mayo PBC sequential data: flexible continuous-margin candidate

`survival::pbcseq` provides public longitudinal laboratory records for 312
randomized patients. Restricting the prothrombin analysis to patients with at
least four observed measurements gives 227 patients, 1,770 observations, a
median of seven visits, and 1,543 adjacent pairs.

Key audit results:

- residual adjacent Spearman correlation 0.611 after time and baseline
  adjustment;
- generalized gamma (GG) decisively preferred to gamma, normal, and log-normal
  margins (BIC 4,674 versus 5,646, 6,292, and 5,408 respectively);
- subject-level five-fold held-out mean log score also favoured GG
  (-1.312 versus -1.607, -1.804, and -1.530);
- baseline histological stage in the GG scale model improved BIC by 11.3; and
- baseline log bilirubin was a stronger scale candidate in a common-case screen
  (estimated sigma ratio 1.166 per SD; BIC improvement 89.0).

This is the strongest distinctively GAMLSS continuous-margin option. Its main
weakness is the irregular and potentially informative visit process. The data
documentation explicitly notes that extra visits and missing laboratory values
may be related to clinical deterioration. Any use must therefore distinguish
actual follow-up day in the marginal model from sequential visit order in the
copula and state the resulting first-order dependence assumption.

The standard workflow selects GG by a wide margin (time-intercept AIC 4,709.71;
the next candidate, log-normal, has AIC 5,335.49), then selects Frank over the
registered Gaussian alternative. The selected joint GG/Frank fit converges and
supports validated analytical-Hessian inference after a tighter continuation
(standardized fitted score 0.039). Its log likelihood is -1,874.48, AIC
3,782.95, and BIC 3,876.09. This is currently the cleanest end-to-end example
of the exact vignette selection and fitting workflow.

Source: <https://stat.ethz.ch/R-manual/R-devel/library/survival/html/pbcseq.html>

### 3. Vietnamese adolescent step counts: clean balanced fallback

The PLOS supporting data contain seven daily pedometer counts for 475
adolescents (3,325 scheduled observations). The article states that values
below 1,000 or above 30,000 steps are invalid and describes weekday/weekend
mean imputation for missed days.

For a strict count analysis, retain all subjects and set only the following
five cells to `NA`:

- two cells outside the stated valid range; and
- three visibly fractional cells, which cannot be unmodified step counts.

This leaves 3,320 observed integer counts. Fit with
`missingness = "segment"`; discarding all five subjects is unnecessary. The
workbook does not identify integer-valued imputations, so these cannot be
removed reliably and must be acknowledged.

`select_margin()` ranks NBI first, ZINBI second with delta AIC 14, and Poisson a
distant third. Conditional on NBI, the registered copula screen selects Gumbel
over Clayton. The selected NBI/Gumbel joint fit converges with cluster-sandwich
inference (standardized fitted score 0.060), log likelihood -31,033.83, AIC
62,093.66, and BIC 62,173.07. However, the best
non-location candidate, `sigma ~ sex`, was not supported after joint fitting
in the earlier analysis (likelihood-ratio p = 0.267; AIC and BIC both favoured
constant dispersion). BMI and PAQ-C dispersion effects were weaker. The data
are therefore useful for count/missingness demonstrations but do not satisfy
the full dual-parameter criterion.

Sources:

- <https://journals.plos.org/globalpublichealth/article?id=10.1371/journal.pgph.0004725>
- <https://figshare.com/articles/dataset/Data_of_the_first_sample_N_35_and_second_sample_N_475_/29258770>

## Other candidates screened

| Candidate | Useful feature | Main reason not promoted |
|---|---|---|
| `JMbayes2::aids` reconstructed CD4 count | 467 subjects; NBI; strong dependence | Heavy dropout; copula covariate separation modest |
| NHANES 2011--12 daily MIMS | 6,507 subjects; seven days | Continuous outcome is nearly symmetric; weak flexible-margin story |
| `Ecdat::Wages` | 595 people x seven years; very high dependence | Identifier/time reconstruction and weak non-location effects |
| `JMbayes2::prothro` | 488 subjects; GG plausible | Irregular schedule and weak treatment effect |
| `geepack::dietox` | Nearly 12 visits; very high dependence | Only 72 subjects |
| `datasets::ChickWeight` | Vivid growth trajectories; high dependence | Only 50 subjects |

## Decision rule

1. Treat `pbcseq` as the strongest end-to-end workflow candidate: both selected
   distributions carry through to a converged joint fit with valid inference.
   Its visit-order and potentially informative follow-up limitations must be
   explained clearly.
2. Retain `PatentsRDUS` as the strongest substantive dual-parameter candidate,
   but disclose that the routine joint fit uses the next viable choices from
   both screens. Do not describe NBI/Clayton as the selection winner.
3. Use the adolescent step data only if a balanced human example is more
   important than satisfying the beyond-mean marginal-parameter criterion.
