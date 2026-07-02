# gamlss.longitudinal: a flexible framework for fitting longitudinal models with gamlss margins and copula dependence

`gamlss.longitudinal` fits longitudinal GAMLSS models with flexible marginal distributions and first-order copula dependence across repeated measurements. Any parameter of the marginal distribution or copula can depend on covariates, so strong time or covariate effects on distributional shape or correlation can be accurately captured. The model provides a great alternative to standard models such as the GEE when a longitudinal analysis needs more than a marginal mean and exponential family distributions, and more control over the dependence structure: changing scale, skewness, tails, quantiles, probabilities, or within-subject dependence.

The package is in an early CRAN/JSS hardening stage. The exported modelling,
prediction, simulation, diagnostics, and benchmark helper APIs are intended to
be reviewable, but remaining review TODOs are tracked explicitly in the
standards crosswalk.

## Install

```r
remotes::install_github("ahibbert/gamlss.longitudinal")
library(gamlss.longitudinal)
```

## Start here

For a full overview of available support articles / vignettes: [Article guide](https://gamlsslongitudinal.aydins-workbench.com/articles/site-guide.html)

After installation, you can run a fast end-to-end smoke test to ensure key components have been loaded correctly with:

```r
source(system.file("smoke-tests", "new-user-smoke.R", package = "gamlss.longitudinal"))
```

Two workflow examples to get started:

- [Minimal workflow](https://gamlsslongitudinal.aydins-workbench.com/articles/standard-workflow.html)
- [Detailed worked example](https://gamlsslongitudinal.aydins-workbench.com/articles/native-simulation-workflow.html)

Reference guides for specific parts of the workflow:

- [Inference](https://gamlsslongitudinal.aydins-workbench.com/articles/inference-uncertainty.html)
- [Diagnostics](https://gamlsslongitudinal.aydins-workbench.com/articles/diagnostics-decisions.html)
- [Simulator usage](https://gamlsslongitudinal.aydins-workbench.com/articles/simulator-usage.html)

## Review and reproducibility

Reviewer-facing package, statistical-software, and paper-replication guidance is
collected in [REVIEW.md](REVIEW.md). The rOpenSci statistical standards
crosswalk is available at `inst/standards/ropensci-srr-compliance.md`, and the
paper replication entry point is documented in `paper/README.md`.

## References

The marginal modelling framework follows GAMLSS, as described by Stasinopoulos and colleagues (2024). The dependence layer uses bivariate copulas, following standard copula references, e.g. Nelsen (2006).

The motivating bivariate-methods comparison is:
Sareff-Hibbert, A. *A comparison between copula-based, mixed model, and
estimating equation methods for regression of bivariate correlated data*.
<https://arxiv.org/abs/2410.11892>

The full paper covering the package is a work in progress, link to be added as it becomes available.
