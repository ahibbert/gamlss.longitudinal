# JSS Manuscript Blueprint

Working title:

> gamlss.longitudinal: Longitudinal GAMLSS Models With Copula Dependence in R

This file is the decision-complete blueprint for rewriting the current draft as
a Journal of Statistical Software article. It is not a substitute for the final
JSS LaTeX source. Use it as the section-level contract while the manuscript,
examples, and replication workflow are brought into sync.

## Abstract

Longitudinal analyses often require flexible marginal distributions and
interpretable within-subject dependence, especially when scale, skewness, tails,
or association vary with time or covariates. The R package
`gamlss.longitudinal` implements first-order copula dependence for GAMLSS
margins in a formula-based workflow. Any supported marginal parameter and
copula parameter can depend on fixed or smooth covariates, allowing users to
model location, scale, shape, and adjacent-time dependence in one fitted object.
The package provides model fitting, distribution and copula screening,
inference, diagnostics, prediction, simulation, and reviewer-facing replication
workflows. A reproducible worked example and targeted validation studies
illustrate parameter recovery, diagnostic checks, optimizer tradeoffs, and the
relationship to existing R tools for GAMLSS, GEE/GLMM, bivariate copula
regression, and vine copulas.

## 1. Introduction

Purpose: establish the statistical problem and position the package as the
contribution.

Required content:

- Longitudinal data often need more than a conditional mean and constant
  variance model.
- GAMLSS handles flexible marginal distributions but does not, by default,
  provide a simple marginal-regression workflow for within-subject copula
  dependence.
- GEE and GLMM workflows handle dependence but usually target a narrower
  marginal distribution set or change interpretation through random effects.
- Vine-copula tools are highly flexible but can become difficult to specify,
  optimize, and interpret as ordinary longitudinal regression models.
- `gamlss.longitudinal` fills a software gap: formula-based longitudinal GAMLSS
  margins with first-order copula dependence, diagnostics, prediction,
  simulation, and reproducible validation workflows.

Main text to avoid:

- Do not overstate the package as universally simple.
- Do not claim arbitrary vine flexibility.
- Do not center the introduction on the simulation study.

## 2. Statistical Model

Purpose: define enough of the model for users and statistical reviewers without
turning the article into a derivation-heavy methods paper.

Required subsections:

### 2.1 Longitudinal GAMLSS margins

- Define subject index, time index, response, and covariates.
- State that each marginal parameter can have its own formula and link.
- Mention fixed and smooth terms, with smooth details deferred to software
  implementation and references.

### 2.2 First-order copula dependence

- Define adjacent-time copula terms on probability-scale margins.
- State the main continuous-response likelihood.
- Explain first-order truncation: the package models adjacent dependence
  directly; higher-order dependence is implied only through the fitted
  first-order structure and should be checked diagnostically.

### 2.3 Discrete responses and incomplete panels

- Describe rectangle probabilities for discrete margins.
- State the missing-response policy: observed rows contribute marginal terms;
  observed adjacent pairs contribute copula terms; the package does not perform
  imputation.

### 2.4 Inference and uncertainty

- Summarize model-based covariance, available fallbacks, and optional robust or
  simulation-based evidence where used in validation.
- Move Hessian derivations to supplementary material.

Critical wording:

> For continuous margins and the specified first-order copula construction, the
> fitted objective is the likelihood for the model defined by these marginal and
> adjacent-pair dependence assumptions. When the data-generating process has
> material residual higher-order dependence beyond this structure, the fitted
> model should be treated as a simplified first-order longitudinal dependence
> model and checked with residual and tail-dependence diagnostics.

## 3. Software Design and User Interface

Purpose: make this the major JSS software section.

Required subsections:

### 3.1 Package design

- Main fitting function: `gamlss_longitudinal()`.
- Formula-per-parameter interface.
- Supported family objects through `gamlss.dist`.
- Copula backend and family/link metadata.
- Stable fitted-object class and stored metadata.

### 3.2 Minimal fit

Include a short runnable example using the primary reproducible example data:

```r
fit <- gamlss_longitudinal(
  mu.formula = response ~ time + treatment + s(baseline),
  sigma.formula = ~ time + treatment,
  theta.formula = ~ time + treatment,
  dataset = dat,
  time_var = "time",
  subject_var = "id",
  margin_dist = GA(),
  copula_dist = "t"
)
summary(fit)
```

Replace formula and family names with the final primary example.

### 3.3 Selection helpers

- `select_margin()`
- `select_copula()`
- `select_joint_distribution()`
- Explain that selection helpers are workflow aids, not automatic guarantees of
  a final scientific model.

### 3.4 Standard methods

Document and demonstrate:

- `summary()`, `coef()`, `vcov()`, `confint()`
- `predict()`, `simulate()`
- `plot()`, `plot_terms()`, `plot_margin_fit()`,
  `plot_copula_diagnostics()`
- `logLik()`, `formula()`, `terms()`, `nobs()`, `fitted()`, `residuals()`,
  `model.frame()`

### 3.5 Diagnostics and decision workflow

- Marginal diagnostics: PIT, QQ, worm plot, rootogram where applicable.
- Dependence diagnostics: fitted versus observed rank association, Rosenblatt
  residual checks, tail co-occurrence, conditional tail exceedance, residual
  lag dependence.
- Missingness diagnostics: clarify descriptive checks and no imputation.

### 3.6 Optimizer controls

- Explain separate RS as exploratory.
- Explain joint RS/CG as final-fit candidates.
- Report convergence metadata and runtime tradeoffs.

## 4. Reproducible Worked Workflow

Purpose: show the full package workflow from data to interpretation.

Decision: the primary example must use public, package-shipped, or fully
simulated data. Private-data examples may be secondary only.

Required steps:

1. Load or generate data.
2. Explore margins and dependence.
3. Screen margins.
4. Screen copulas.
5. Fit final model.
6. Summarize coefficients and uncertainty.
7. Run marginal and copula diagnostics.
8. Predict distributional quantities.
9. Simulate from the fitted model.
10. State what a practitioner should conclude and what they should check next.

Target main-paper outputs:

- One compact table of candidate families/copulas.
- One final-model coefficient or effect table.
- One diagnostics figure.
- One prediction or simulation summary.

## 5. Relationship to Existing R Software

Purpose: satisfy the JSS expectation that novel software is compared with
similar implementations.

Main comparison table columns:

- Package/software
- Marginal distribution flexibility
- Longitudinal dependence structure
- Covariate-dependent dependence
- Full or staged estimation
- Smooth terms
- Discrete response support
- Diagnostics/prediction/simulation
- Main limitation relative to `gamlss.longitudinal`

Rows:

- `gamlss`
- `gamlss2`
- `geepack`
- `lme4`
- `GJRM`
- `VineCopula`
- `rvinecopulib`, if included
- `gamlss.longitudinal`

Keep the text balanced: each comparator should have a real strength, not only a
limitation.

## 6. Validation and Performance

Purpose: replace broad simulation dominance with targeted evidence.

Main-paper validation set:

- Continuous correctly specified case: recovery, coverage, runtime.
- Discrete correctly specified case: recovery, coverage, runtime.
- Joint versus separate optimization: guidance table.
- First-order misspecification: concise limitation and robust-SE guidance.
- Copula misspecification: tail-shape consequences.

Move to supplement:

- Full grid tables by family/correlation/design.
- Smooth recovery plots not needed for the main narrative.
- Detailed optimizer derivations.
- Full missingness sensitivity until module 04 is final.

## 7. Reproducibility

Purpose: make the paper, package, and replication files one reviewable object.

Required content:

- `paper/replicate.R` is the paper-facing entry point.
- `GAMLSS_LONGITUDINAL_JSS_PROFILE=smoke` runs a quick reviewer check.
- `GAMLSS_LONGITUDINAL_JSS_PROFILE=expanded` regenerates the full paper
  outputs.
- `paper/manifest.csv` maps result IDs to output files and analysis states.
- `results/jss-replication/<profile>/logs/session_info.txt` and
  `output_hashes.csv` record software state and output hashes.
- Private data modules must be skipped or clearly marked unless data access is
  supplied.

## 8. Discussion

Required content:

- The package is most useful when marginal distributional shape or adjacent
  dependence is scientifically important.
- Main limitations: first-order dependence assumption, runtime, model-selection
  burden, and sensitivity to major dependence misspecification.
- Future work: higher-order residual dependence handling, faster derivatives,
  broader public examples, and expanded platform/runtime validation.

## Main Paper Target Outputs

Use the audit file to finalize the exact set, but target no more than:

- 5 to 7 tables.
- 5 to 7 figures.
- One primary worked example.
- One compact validation section.
- Supplementary material for large grids and derivations.

## Submission Bundle Checklist

- JSS-style PDF.
- JSS LaTeX source and references.
- R package source.
- Replication scripts and manifest.
- Data or accepted download instructions.
- Expected smoke and expanded runtime notes.
- Session information and output hashes.
- License and metadata confirmation.
