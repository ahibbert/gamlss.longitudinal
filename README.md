# gamlss.longitudinal

`gamlss.longitudinal` provides a flexible framework for longitudinal GAMLSS
models with copula dependence across repeated measurements. Marginal
distribution parameters can depend on fixed or smooth covariates, and the
copula dependence parameter can also be modelled with covariates.

Supported copula family codes are:

- `"N"`: Gaussian
- `"C"`: Clayton
- `"F"`: Frank
- `"G"`: Gumbel
- `"J"`: Joe
- `"t"`: Student t

## Installation

The core package can be installed from GitHub with CRAN dependencies. Some
optional benchmarking and JSS replication comparisons use `gamlss2`; install it
from the GAMLSS R-universe only when running those opt-in workflows:

```r
install.packages(
  "gamlss2",
  repos = c(
    "https://gamlss-dev.R-universe.dev",
    "https://cloud.R-project.org"
  )
)

devtools::install_github("ahibbert/gamlss.longitudinal")
```

## Article Guide

Start with the article guide when browsing the package site or repository:

- Guide to the articles: [`vignettes/site-guide.Rmd`](vignettes/site-guide.Rmd)

The main worked examples are:

- Standard workflow API, inference, and diagnostics: [`vignettes/standard-workflow.Rmd`](vignettes/standard-workflow.Rmd)
- Native simulation, copula screening, and fitting: [`vignettes/native-simulation-workflow.Rmd`](vignettes/native-simulation-workflow.Rmd)

Focused adoption guides:

- Adoption guide: [`vignettes/replace-gee-glmm.Rmd`](vignettes/replace-gee-glmm.Rmd)
- Adoption decision guide: [`vignettes/adoption-decision-guide.Rmd`](vignettes/adoption-decision-guide.Rmd)
- Inference and uncertainty: [`vignettes/inference-uncertainty.Rmd`](vignettes/inference-uncertainty.Rmd)
- Diagnostics as decisions: [`vignettes/diagnostics-decisions.Rmd`](vignettes/diagnostics-decisions.Rmd)
- Replace a GEE: [`vignettes/replace-a-gee.Rmd`](vignettes/replace-a-gee.Rmd)
- Replace a random-intercept GLMM: [`vignettes/replace-a-random-intercept-glmm.Rmd`](vignettes/replace-a-random-intercept-glmm.Rmd)
- Non-Gaussian longitudinal outcomes: [`vignettes/non-gaussian-longitudinal-outcomes.Rmd`](vignettes/non-gaussian-longitudinal-outcomes.Rmd)
- Missing visits: [`vignettes/missing-visits.Rmd`](vignettes/missing-visits.Rmd)
- Time-varying dependence: [`vignettes/time-varying-dependence.Rmd`](vignettes/time-varying-dependence.Rmd)
- Benchmarking adoption claims: [`vignettes/benchmarking-adoption.Rmd`](vignettes/benchmarking-adoption.Rmd)
- Benchmark investigation: [`inst/benchmarks/gee-glmm-investigation.md`](inst/benchmarks/gee-glmm-investigation.md)

The native simulation vignette demonstrates:

- native longitudinal data simulation,
- exploratory `plotDist()` diagnostics,
- marginal family screening with `gamlss::fitDist()`,
- native copula family screening with `select_copula()`,
- fitting a longitudinal GAMLSS-copula model,
- `summary()`, `plot()`, `plot_terms()`, and `plot_copula_diagnostics()` diagnostics,
- comparison between simulated truth and fitted effects.

To render the two worked examples locally:

```r
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/Quarto/bin/tools")
rmarkdown::render("vignettes/standard-workflow.Rmd")
rmarkdown::render("vignettes/native-simulation-workflow.Rmd")
```

If Pandoc is already available to R, the `Sys.setenv()` line is not needed.

## Standard Workflow

The adoption-facing workflow is:

```r
margin_screen <- screen_margin(dat, response_var = "response")

copula_screen <- select_copula(
  data = dat,
  u_var = "u",
  subject_var = "subject",
  time_var = "time"
)

fit <- fit_longitudinal(
  dataset = dat,
  margin_dist = gamlss.dist::GA(mu.link = "log", sigma.link = "log"),
  time_var = "time",
  subject_var = "subject",
  mu.formula = response ~ treatment + time + s(age_scaled, bs = "ps"),
  sigma.formula = ~ treatment + time,
  theta.formula = ~ time,
  copula_dist = attr(copula_screen, "selected")
)

check_model(fit)
plot(check_model(fit))
copula_time_summary(fit)
wald_test(fit, terms = c("mu.treatment", "sigma.treatment"))
likelihood_compare(fit_reduced, fit)
bootstrap_inference(fit, R = 200, terms = "mu.treatment")
predict(fit, type = "response", se.fit = TRUE, interval = "confidence")
predict(fit, type = "quantile")
predict(fit, type = "probability", q = 10, direction = "above")
simulate(fit, nsim = 10)
marginal_effects(fit, newdata = dat, variable = "treatment", se.fit = TRUE)
effects(fit, newdata = dat, variable = "treatment", se.fit = TRUE)
```

`check_model()` is intended as the first diagnostic pass for applied users. It
summarises convergence, scoring, PIT calibration, tail calibration, residual
dependence, copula summaries, and plain-language diagnostic decisions.

Use [`vignettes/inference-uncertainty.Rmd`](vignettes/inference-uncertainty.Rmd)
when reporting uncertainty. It separates likelihood/Hessian coefficient
intervals, numerical-Hessian fallbacks, delta-method response-mean intervals,
and simulation-based trajectory summaries.
Use [`vignettes/diagnostics-decisions.Rmd`](vignettes/diagnostics-decisions.Rmd)
when deciding whether diagnostic warnings support a primary model, sensitivity
model, or revised specification.

For simulation papers and applied sensitivity analyses, the package also has an
opt-in comparator scaffold:

```r
benchmark_standard_models(
  data = dat,
  formula = response ~ treatment + time + age_scaled,
  subject_var = "subject",
  family = "gaussian",
  comparators = c("gee", "glmm", "gam")
)
```

This fits standard baselines when optional packages are installed and reports
availability, success, runtime, and simple response-scale prediction error.
For repeated simulation evidence, include those comparators explicitly in the
coverage harness:

```r
run_coverage_simulations(
  families = "NO",
  copulas = "N",
  methods = c("rs_separate", "gee", "glmm", "gam"),
  designs = "covariate",
  write_results = FALSE
)
```

Comparator rows include truth-aware mean metrics such as
`benchmark_mean_rmse` when simulations include `true_mu`, plus
distribution-aware `benchmark_q90_mae` and
`benchmark_upper_tail_error_90` diagnostics for Gaussian, Gamma, and Poisson
comparators. Gaussian comparator rows also include residual-calibrated PIT,
tail, and 95% interval diagnostics. Time-varying dependence scenarios include
`benchmark_theta_time_abs_error` for the fitted copula time effect.
Use `summarise_benchmark_results()` to turn repeated simulation rows into
win/tie/loss summaries by metric and estimand:

```r
bench_summary <- summarise_benchmark_results(
  results,
  metrics = c("benchmark_mean_rmse", "benchmark_q90_mae", "benchmark_upper_tail_error_90", "elapsed_sec")
)
bench_summary
```

For a named benchmark plan that is ready to run outside CRAN checks:

```r
adoption_benchmark_scenarios()

bench <- run_adoption_benchmarks(
  scenarios = "gaussian_heteroskedastic",
  reps = 25,
  methods = c("rs_separate", "gee", "glmm", "gam"),
  write_results = TRUE
)
bench$summary$summary
```

For a longer reproducible campaign, run the shipped opt-in benchmark script:

```r
script <- system.file(
  "benchmarks",
  "run-adoption-benchmarks.R",
  package = "gamlss.longitudinal"
)
source(script)
```

Set `GAMLSS_LONGITUDINAL_ADOPTION_REPS`,
`GAMLSS_LONGITUDINAL_ADOPTION_SCENARIOS`,
`GAMLSS_LONGITUDINAL_ADOPTION_METHODS`, and
`GAMLSS_LONGITUDINAL_ADOPTION_OUTPUT_DIR` before sourcing the script to control
the benchmark campaign and output location. The script writes CSV/RDS outputs
and an `adoption_benchmark_report.md` file. If results are already in memory,
use `write_benchmark_report(bench)` to create the Markdown evidence report.
Read the scenario-level report sections when making adoption claims; the
aggregate table is only an orientation because each scenario targets a
different estimand.

## Minimal Fit Skeleton

The core model call looks like this:

```r
fit <- gamlss.longitudinal(
  dataset = dat,
  margin_dist = gamlss.dist::GA(mu.link = "log", sigma.link = "log"),
  copula_dist = "C",
  time_var = "time",
  subject_var = "subject",
  mu.formula = response ~ treatment + time_scaled + s(age_scaled, bs = "ps", k = 8),
  sigma.formula = ~ treatment + time_scaled + s(age_scaled, bs = "ps", k = 5),
  theta.formula = ~ time_scaled + s(age_scaled, bs = "ps", k = 8)
)

summary(fit)
plot_terms(fit, data = dat)
plot_copula_diagnostics(fit, data = dat)
```

See the vignette for a complete reproducible dataset, distribution screening,
copula selection, and diagnostic workflow.

## rOpenSci Standards Notes

This package is being developed against the rOpenSci statistical software
standards for general statistical software, regression software, and
probability-distribution software. The package is in an initially stable state
of development, with active subsequent development expected as the longitudinal
workflow, diagnostics, and benchmark evidence mature.

For a standard-by-standard reviewer crosswalk, see
[`inst/standards/ropensci-srr-compliance.md`](inst/standards/ropensci-srr-compliance.md).
The same compliance claims are also encoded as `srr` roxygen tags in
[`R/srr-stats-standards.R`](R/srr-stats-standards.R), so the package can be
checked with the rOpenSci [`srr`](https://ropensci-review-tools.github.io/srr/)
tools.

### Primary References

The marginal modelling framework follows GAMLSS:

- Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive models for
  location, scale and shape. *Journal of the Royal Statistical Society: Series
  C*, 54, 507-554. <https://doi.org/10.1111/j.1467-9876.2005.00510.x>
- Stasinopoulos, D. M. and Rigby, R. A. (2007). Generalized additive models for
  location scale and shape (GAMLSS) in R. *Journal of Statistical Software*,
  23(7), 1-46. <https://doi.org/10.18637/jss.v023.i07>
- Stasinopoulos, D. M., Rigby, R. A., Heller, G. Z., Voudouris, V., and De
  Bastiani, F. (2017). *Flexible Regression and Smoothing: Using GAMLSS in R*.
  Chapman and Hall/CRC.

The dependence layer uses bivariate copulas and Kendall-tau parameterisations:

- Nelsen, R. B. (2006). *An Introduction to Copulas*, 2nd edition. Springer.
  <https://doi.org/10.1007/0-387-28678-0>
- Brechmann, E. C. and Schepsmeier, U. (2013). Modeling dependence with C- and
  D-vine copulas: The R package CDVine. *Journal of Statistical Software*,
  52(3), 1-27. <https://doi.org/10.18637/jss.v052.i03>

### Prior Art

`gamlss.longitudinal` is not a replacement for `gamlss`; it builds on
`gamlss.dist` family objects and uses `gamlss` as the familiar marginal-model
reference point. It adds a longitudinal first-order copula dependence layer,
workflow helpers, diagnostics, and opt-in benchmarking scaffolds for comparing
against common GEE, GLMM, and GAM baselines.

Comparable R implementations include:

- `gamlss` and `gamlss2` for distributional regression with GAMLSS margins;
- `VineCopula` for bivariate and vine copula densities, CDFs, h-functions, and
  parameter conversions;
- `geepack`, `lme4`, `glmmTMB`, and `mgcv` for standard longitudinal,
  mixed-model, and smooth mean-model baselines.

### Input and Pre-processing Policy

Model input must be a long-format table with one row per observed subject-time
combination, a response in the left-hand side of `mu.formula`, a subject column,
and a time column. The package converts tibbles and other data-frame-like inputs
to a plain data frame for fitting.

`time_var` is used in two ways. Internally, it becomes numeric `time` for margin
ordering and adjacent-pair construction. For formulas it is preserved as
`time_covariate`, so numeric time, factor time, and factor interactions can be
modelled without accidentally treating a categorical visit as a continuous
trend. Numeric-like character time is converted with `as.numeric()`; non-numeric
character time should be converted to a factor by the user before fitting.
Ordered factors are treated with treatment contrasts.

Structurally missing subject-time rows are expanded to explicit rows with
missing responses. Fitting stops if any margin has no observed responses, or if
any adjacent copula pair has no complete response pairs. Predictor missingness
and non-finite predictor values are handled by the model-matrix path; users
should inspect `model.frame(fit, type = "expanded")`,
`model.frame(fit, type = "observed")`, and `fit$convergence` when auditing a
fit. The package does not promise to preserve submitted row names after grid
expansion; case identity is represented by subject, time, and response columns.

### Algorithms and Accessors

Fitted objects have class `gamlss.longitudinal` and expose standard regression
accessors: `coef()`, `confint()`, `formula()`, `terms()`, `nobs()`,
`model.frame()`, `fitted()`, `residuals()`, `vcov()`, `logLik()`, `summary()`,
`predict()`, `simulate()`, and `plot()`. Convergence metadata is stored in
`fit$convergence`.

Most distribution calculations delegate to `gamlss.dist` family functions or
native copula functions. The native copula backend is tested against
`VineCopula`; users can set
`options(gamlss.longitudinal.copula_backend = "vinecopula")` to delegate copula
operations where the optional package is installed. Analytical Hessian paths are
used when available, with documented fallbacks to numerical finite-difference
paths for difficult discrete or near-boundary cases.

## Reference

Motivation for the approach and its performance compared to alternative methods
in the bivariate case is described in the work-in-progress paper:

Sareff-Hibbert, A. *A comparison between copula-based, mixed model, and
estimating equation methods for regression of bivariate correlated data*.
<https://arxiv.org/abs/2410.11892>
