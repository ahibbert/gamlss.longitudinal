# gamlss.longitudinal: a flexible framework for fitting longitudinal models with gamlss margins and copula dependence

`gamlss.longitudinal` fits longitudinal GAMLSS models with flexible marginal distributions and first-order copula dependence across repeated measurements. Any parameter of the marginal distribution or copula can depend on covariates, so strong time or covariate effects on distributional shape or correlation can be accurately captured. The model provides a great alternative to standard models such as the GEE when a longitudinal analysis needs more than a marginal mean and exponential family distributions, and more control over the dependence structure: changing scale, skewness, tails, quantiles, probabilities, or within-subject dependence.

Version 0.3.1 is a maintenance update for packaging,
documentation, and checks. This work does not change statistical methods,
defaults, data, or fitted-model recipes from the frozen 0.3.0 review reference.
That version-only reference is commit
`2092e1b`, based on package source `bfb6557a45f4230f306e82a6de7a07375537f6f3`.
This maintenance label does not identify the source of every historical paper
result, certify statistical equivalence, or establish a successful end-to-end
`full` paper run.

The installed guides describe modelling, prediction, simulation, diagnostics
and inference limitations. The standards crosswalk records review evidence and
remaining work separately from paper replication.

## Install

```r
remotes::install_github("ahibbert/gamlss.longitudinal")
library(gamlss.longitudinal)
```

## Start here

For a full overview of available support articles / vignettes: [Article guide](https://gamlsslongitudinal.aydins-workbench.com/articles/site-guide.html)

After installation, run the package new-user smoke test to check its basic
workflow. This is separate from the paper workflow's `smoke`, `paper` and `full`
profiles:

```r
source(system.file("smoke-tests", "new-user-smoke.R", package = "gamlss.longitudinal"))
```

Version 0.3.1 includes six compact installed vignettes. The detailed worked
example and JSS reviewer tour remain repository and website articles, outside
the CRAN source tarball.

Two workflow examples to get started:

- [Minimal workflow](https://gamlsslongitudinal.aydins-workbench.com/articles/standard-workflow.html)
- [Detailed worked example](https://gamlsslongitudinal.aydins-workbench.com/articles/native-simulation-workflow.html)

Reference guides for specific parts of the workflow:

- [Inference](https://gamlsslongitudinal.aydins-workbench.com/articles/inference-uncertainty.html)
- [Diagnostics](https://gamlsslongitudinal.aydins-workbench.com/articles/diagnostics-decisions.html)
- [Simulator usage](https://gamlsslongitudinal.aydins-workbench.com/articles/simulator-usage.html)

## Review and reproducibility

Package and statistical-software review starts in [REVIEW.md](REVIEW.md).
The [standards crosswalk](inst/standards/ropensci-srr-compliance.md) records
standards evidence and remaining work; [method traceability](inst/standards/method-traceability.csv)
maps package source files to tests.

Paper replication is maintained separately. Its repository and reviewer guide
are being prepared for public release; a public download link will be added
when available. Reviewers with a supplied source bundle should use its exact
software and paper revisions, for either clones or archives; do not mix branch
heads.

That guide defines workspace assembly, the `smoke`, `paper` and `full` paper
profiles, outputs and release checks. Package tests and the new-user smoke test
do not establish that those paper profiles reproduce the manuscript. See also
[online JSS reviewer tour](https://gamlsslongitudinal.aydins-workbench.com/articles/jss-start-here.html) for the package tour.

## References

The marginal modelling framework follows GAMLSS, as described by Stasinopoulos and colleagues (2024). The dependence layer uses bivariate copulas, following standard copula references, e.g. Nelsen (2006).

The motivating bivariate-methods comparison is:
Sareff-Hibbert, A. and Heller, G. Z. *A comparison between copula-based, mixed model, and
estimating equation methods for regression of bivariate correlated data*.
<https://arxiv.org/abs/2410.11892>

The full paper covering the package is a work in progress, link to be added as it becomes available.
