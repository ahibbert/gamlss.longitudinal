# gamlss.longitudinal

`gamlss.longitudinal` fits longitudinal GAMLSS models with flexible marginal
distributions and first-order copula dependence across repeated measurements.
Use it when a longitudinal analysis needs more than a marginal mean: changing
scale, skewness, tails, quantiles, probabilities, or within-subject dependence.

## Install

```r
remotes::install_github("ahibbert/gamlss.longitudinal")
```

Optional comparator and benchmark workflows use additional packages such as
`geepack`, `lme4`, `mgcv`, and `gamlss2`.

## Start Here

The package site is the main entry point:

- [Article guide](https://ahibbert.github.io/gamlss.longitudinal/articles/site-guide.html)
- [Standard workflow](https://ahibbert.github.io/gamlss.longitudinal/articles/standard-workflow.html)
- [Detailed worked example](https://ahibbert.github.io/gamlss.longitudinal/articles/native-simulation-workflow.html)

Support articles:

- [Inference and uncertainty](https://ahibbert.github.io/gamlss.longitudinal/articles/inference-uncertainty.html)
- [Diagnostics as decisions](https://ahibbert.github.io/gamlss.longitudinal/articles/diagnostics-decisions.html)
- [Simulator usage](https://ahibbert.github.io/gamlss.longitudinal/articles/simulator-usage.html)
- [Benchmarking adoption claims](https://ahibbert.github.io/gamlss.longitudinal/articles/benchmarking-adoption.html)

## Workflow Map

```r
margin_selection <- select_margin(dat, response_var = "response")
margin_dist <- best_fit_family(margin_selection)

copula_selection <- select_copula(
  data = dat,
  response_var = "response",
  margin_dist = margin_dist,
  subject_var = "subject",
  time_var = "time"
)
copula_dist <- best_fit_family(copula_selection)

joint_selection <- select_joint_distribution(
  data = dat,
  response_var = "response",
  subject_var = "subject",
  time_var = "time",
  margin_families = head(margin_selection$family, 10),
  copula_families = c("N", "C", "F", "G", "J", "t")
)

missingness <- check_missingness(
  dat,
  response_var = "response",
  time_var = "time",
  subject_var = "subject"
)

fit <- gamlss_longitudinal(
  dataset = dat,
  margin_dist = margin_dist,
  copula_dist = copula_dist,
  time_var = "time",
  subject_var = "subject",
  mu.formula = response ~ treatment + time_scaled + s(age_scaled, bs = "ps"),
  sigma.formula = ~ treatment + time_scaled,
  theta.formula = ~ time_scaled
)

summary(fit)
check_model(fit, include_vcov = TRUE, dependence_cor_cutoff = 0.25)
plot(fit)
plot_terms(fit, data = dat)
plot_copula_diagnostics(fit, data = dat)

confint(fit)
wald_test(fit, terms = "mu.treatment")
likelihood_compare(reduced = fit_reduced, main = fit)
bootstrap_inference(fit, R = 200, terms = "mu.treatment")

predict(fit, type = "response", se.fit = TRUE, interval = "confidence")
predict(fit, type = "quantile", probs = c(0.1, 0.5, 0.9))
predict(fit, type = "probability", q = 10, direction = "above")
marginal_effects(fit, newdata = dat, variable = "treatment", se.fit = TRUE)
simulate(fit, nsim = 10)
```

Core copula family codes are `"N"` Gaussian, `"C"` Clayton, `"F"` Frank,
`"G"` Gumbel, `"J"` Joe, and `"t"` Student t.

## Simulation And Benchmarks

Use the simulation helpers for reproducible examples and method checks:

```r
dat <- simulate_longitudinal_dataset(
  n = 100,
  times = seq(0, 1, length.out = 5),
  margin_dist = GA(mu.link = "log", sigma.link = "log"),
  copula_dist = "C"
)
```

Use opt-in benchmark helpers when you need comparator evidence:

```r
benchmark_standard_models(
  data = dat,
  formula = response ~ treatment + time_scaled + age_scaled,
  subject_var = "subject",
  family = "gaussian",
  comparators = c("gee", "glmm", "gam"),
  fit = fit
)
```

For repeated benchmark campaigns, see
[Benchmarking adoption claims](https://ahibbert.github.io/gamlss.longitudinal/articles/benchmarking-adoption.html).

## References

The marginal modelling framework follows GAMLSS, especially Rigby and
Stasinopoulos (2005), Stasinopoulos and Rigby (2007), and Stasinopoulos et al.
(2017). The dependence layer uses bivariate copulas and Kendall-tau
parameterisations, following standard copula references such as Nelsen (2006).

The motivating bivariate-methods comparison is:

Sareff-Hibbert, A. *A comparison between copula-based, mixed model, and
estimating equation methods for regression of bivariate correlated data*.
<https://arxiv.org/abs/2410.11892>
