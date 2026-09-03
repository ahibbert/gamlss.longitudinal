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
```

The last command runs the full diagnostic and publication-output
bundle and is deliberately opt-in because it is more expensive than fitting
and saving the model.

## Candidate ranking

### 1. PatentsRDUS: primary workflow candidate

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
- NBI strongly preferred to Poisson and constant-parameter Delaporte in the
  independent-margin screen; and
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

Independent-margin screening results:

| Model | AIC | BIC |
|---|---:|---:|
| Poisson | 68,643.8 | 68,723.8 |
| NBI, constant sigma | 23,572.6 | 23,658.7 |
| NBI, sigma by scientific sector | 23,502.5 | 23,594.8 |
| NBI, sigma by baseline capital | 23,500.9 | 23,593.1 |
| NBI, sigma by sector and capital | **23,386.4** | **23,484.8** |
| Delaporte, constant sigma and nu | 23,564.1 | 23,656.4 |

The final joint NBI/Clayton model converged at the prespecified
`outer_tol = 0.001` after a short continuation from the capped screening fit. It had log
likelihood -10,468.64, AIC 20,973.27, and BIC 21,083.96. There were no optimizer
events or model warnings.

Both required components survived joint nested comparisons:

| Joint model | df | log likelihood | AIC | BIC |
|---|---:|---:|---:|---:|
| Full: sigma by sector/capital; theta by capital | 18 | **-10,468.64** | **20,973.27** | **21,083.96** |
| Constant theta | 17 | -10,570.99 | 21,175.97 | 21,280.50 |
| Constant sigma | 16 | -10,535.61 | 21,103.22 | 21,201.60 |

Relative to the full model, constant theta loses 102.35 log-likelihood units
for one parameter; constant sigma loses 66.97 units for two parameters. Thus
the two effects are not merely alternative allocations of the same
heterogeneity.

On the NBI log-sigma scale, scientific-sector status has coefficient -0.594
and one SD of log baseline capital has coefficient -0.358. These correspond to
dispersion ratios of 0.552 and 0.699 respectively. On the Clayton parameter
link scale, the capital coefficient is 0.432. The implied Kendall tau is about
0.426, 0.533, and 0.638 at capital z-scores -1, 0, and +1.

Use the data through `Suggests: pglm` rather than copying it into this package
unless redistribution rights for the underlying study data are confirmed.

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

The first capped joint GG/Gaussian-copula fit did not converge within 15 outer
iterations, so its information criteria were correctly unavailable. A
fixed-margin screen suggested greater dependence at high baseline stage
(ordinal-stage BIC improvement 5.2; subject-block bootstrap interval only just
excluded zero), but a longer, prespecified joint fit is still required before
adoption. The GG shape estimate was also extreme, which may explain some of the
joint optimization difficulty.

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

NBI was overwhelmingly preferred to Poisson and residual dependence was high.
Sex-dependent Clayton association was also supported. However, the best
non-location candidate, `sigma ~ sex`, was not supported after joint fitting
with `theta ~ sex` (likelihood-ratio p = 0.267; AIC and BIC both favoured
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

1. Adopt `PatentsRDUS` as the primary workflow candidate: its joint model has
   converged and both the dispersion and capital-dependent association survive
   nested comparison by wide margins.
2. Retain `pbcseq` as the preferred flexible continuous-margin option, but do
   not put it in the paper until a prespecified joint fit converges and the
   visit-order assumption is explained.
3. Use the adolescent step data only if a balanced human example is more
   important than satisfying the beyond-mean marginal-parameter criterion.
