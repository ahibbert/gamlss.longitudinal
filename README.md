# gamlss.longitudinal: a flexible framework for fitting longitudinal models with gamlss margins and copula dependence

`gamlss.longitudinal` fits longitudinal GAMLSS models with flexible marginal distributions and first-order copula dependence across repeated measurements. Any parameter of the marginal distribution or copula can depend on covariates, so strong time or covariate effects on distributinal shape or correlation can be accurately captured. The model provides a great alternative to standard models such as the GEE when a longitudinal analysis needs more than a marginal mean and exponential family distributions, and more control over the dependene structure: changing scale, skewness, tails, quantiles, probabilities, or within-subject dependence.

## Install

```r
remotes::install_github("ahibbert/gamlss.longitudinal")
library(gamlss.longitudinal)
```

## Start Here

For a full overview of available support articles / vignettes: [Article guide](https://ahibbert.github.io/gamlss.longitudinal/articles/site-guide.html)

After installation, you can run a fast end-to-end smoke test to ensure key components have been loaded correctly with:

```r
source(system.file("smoke-tests", "new-user-smoke.R", package = "gamlss.longitudinal"))
```

Two standard workflow examples to get started:

- [Standard workflow](https://ahibbert.github.io/gamlss.longitudinal/articles/standard-workflow.html)
- [Detailed worked example](https://ahibbert.github.io/gamlss.longitudinal/articles/native-simulation-workflow.html)

Additional details on the implemented methods:

- [Inference and uncertainty](https://ahibbert.github.io/gamlss.longitudinal/articles/inference-uncertainty.html)
- [Diagnostics as decisions](https://ahibbert.github.io/gamlss.longitudinal/articles/diagnostics-decisions.html)
- [Simulator usage](https://ahibbert.github.io/gamlss.longitudinal/articles/simulator-usage.html)

## References

The marginal modelling framework follows GAMLSS, Stasinopoulis et al (2024). The dependence layer uses bivariate copulas, following standard copula references, e.g. Nelsen (2006).

The motivating bivariate-methods comparison is:
Sareff-Hibbert, A. *A comparison between copula-based, mixed model, and
estimating equation methods for regression of bivariate correlated data*.
<https://arxiv.org/abs/2410.11892>

The full paper covering the package is a work in progress, link to be added as it becomes available.
