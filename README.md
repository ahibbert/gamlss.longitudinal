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

## Worked Example

The main worked example is now maintained as a vignette:

- Online article: <https://ahibbert.github.io/gamlss.longitudinal/articles/native-simulation-workflow.html>
- Source: [`vignettes/native-simulation-workflow.Rmd`](vignettes/native-simulation-workflow.Rmd)

The vignette demonstrates:

- native longitudinal data simulation,
- exploratory `plotDist()` diagnostics,
- marginal family screening with `gamlss::fitDist()`,
- native copula family screening with `select_copula()`,
- fitting a longitudinal GAMLSS-copula model,
- `summary()`, `plot()`, `plot_terms()`, and `plot_copula_diagnostics()` diagnostics,
- comparison between simulated truth and fitted effects.

To render the vignette locally:

```r
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/Quarto/bin/tools")
rmarkdown::render("vignettes/native-simulation-workflow.Rmd")
```

If Pandoc is already available to R, the `Sys.setenv()` line is not needed.

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

## Reference

Motivation for the approach and its performance compared to alternative methods
in the bivariate case is described in the work-in-progress paper:

Sareff-Hibbert, A. *A comparison between copula-based, mixed model, and
estimating equation methods for regression of bivariate correlated data*.
<https://arxiv.org/abs/2410.11892>
