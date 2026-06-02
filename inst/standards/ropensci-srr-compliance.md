# rOpenSci Statistical Software Standards Compliance Crosswalk

This file documents package compliance with the rOpenSci statistical software
standards for the categories relevant to `gamlss.longitudinal`:

- 6.1 General standards (`G`)
- 6.5 Regression and supervised learning standards (`RE`)
- 6.9 Probability distribution standards (`PD`)

The machine-readable `srr` tags live in `R/srr-stats-standards.R`. This
document is the reviewer-facing companion: it lists every standard ID currently
tracked by `srr`, whether the package is ready for that item, and where a
reviewer should look for evidence or what remains to be done.

Status meanings:

- Met: current code, documentation, or tests provide reviewable evidence.
- Partial/TODO: partly addressed, but additional documentation, tests, or
  release work are needed before claiming full compliance.

## General Standards

| ID | Status | Evidence / change needed |
|---|---|---|
| G1.0 | Met | README lists published GAMLSS, copula, and VineCopula references. |
| G1.1 | Partial/TODO | README positions the package relative to prior art; revisit the final novelty/improvement claim with the manuscript. |
| G1.2 | Met | README includes a lifecycle statement. |
| G1.3 | Partial/TODO | Key terminology is spread across README and help files; add a compact glossary. |
| G1.4 | Met | User-facing functions are roxygen2 documented and Rd files are generated. |
| G1.4a | Partial/TODO | Many internals are documented, but older dense helpers need consistent `@noRd` roxygen blocks. |
| G1.5 | Met | Reproduction and benchmark scaffolds are under `inst/jss-replication` and `inst/benchmarks`. |
| G1.6 | Met | Comparator helpers support GEE, GLMM, GAM, and GAMLSS baselines when optional dependencies are present. |
| G2.0 | Met | Length/scalar assertions exist for core controls and fitting inputs. |
| G2.0a | Partial/TODO | README documents core expectations; all Rd argument entries should be audited for scalar/vector length expectations. |
| G2.1 | Met | Type checks cover formulas, data, family objects/names, subject/time columns, and numeric controls. |
| G2.1a | Partial/TODO | Secondary type documentation should be completed across all help pages. |
| G2.2 | Met | Univariate control arguments are restricted to length-one values where applicable. |
| G2.3 | Partial/TODO | Character controls are constrained; centralize case-sensitivity policy. |
| G2.3a | Met | `match.arg()` or equivalent is used for constrained choices. |
| G2.3b | Partial/TODO | String values are deterministic; document which arguments are case-sensitive. |
| G2.4 | Partial/TODO | Fitting paths perform conversions; broaden explicit conversion tests. |
| G2.4a | Partial/TODO | Audit integer coercions for indices and control inputs. |
| G2.4b | Met | Numeric conversion is used for time/response/numeric covariate paths. |
| G2.4c | Met | Character conversion is used for identifiers and report labels. |
| G2.4d | Partial/TODO | Factor policy is documented; add broader factor conversion tests. |
| G2.4e | Met | Formula/model-matrix machinery handles factor-to-design conversion. |
| G2.5 | Met | README documents ordered/unordered factor behavior and formula/model-matrix routines enforce standard factor handling. |
| G2.6 | Partial/TODO | Main workflow is tabular; one-dimensional helper inputs need class pre-processing audit. |
| G2.7 | Met | Data-frame-like inputs are normalized to base data frames for fitting. |
| G2.8 | Met | Early pre-processing normalizes data, time, subject, formulas, and families. |
| G2.9 | Partial/TODO | Lossy conversions are partly documented; add consistent warnings or explicit policy. |
| G2.10 | Met | Named column extraction avoids default drop behavior. |
| G2.11 | Partial/TODO | Non-standard column classes such as `units` are not yet documented or tested. |
| G2.12 | Partial/TODO | List-column behavior should be explicitly rejected or tested. |
| G2.13 | Met | Required columns, empty margins, and empty adjacent pairs are rejected. |
| G2.14 | Partial/TODO | Missing responses and structural gaps are handled; no user-selectable missingness policy yet. |
| G2.14a | Partial/TODO | Distinct NA/NaN/Inf behavior needs clearer docs and tests. |
| G2.14b | Partial/TODO | Structural gaps are expanded; observed-case retention needs more tests. |
| G2.14c | Partial/TODO | Imputation is not provided; document that explicitly in all workflows. |
| G2.15 | Met | Non-finite likelihood/optimisation states are checked and stored/warned. |
| G2.16 | Partial/TODO | Undefined-value behavior is mostly inherited from model matrices/families; add explicit policy. |
| G3.0 | Partial/TODO | Tolerance checks are being adopted; complete numeric equality audit remains. |
| G3.1 | Met | Analytical and numerical Hessian paths are explicit with guarded fallbacks. |
| G3.1a | Met | Finite-difference controls and convergence diagnostics are accessible. |
| G4.0 | Met | Report output is Markdown and report paths are normalized to `.md`. |
| G5.0 | Partial/TODO | Known-truth fixture exists; add broader canonical/external benchmark datasets if appropriate. |
| G5.1 | Met | testthat and CI cover main workflow, diagnostics, prediction, simulation, and copula parity. |
| G5.2 | Partial/TODO | Many validation branches are tested; not every user-reachable stop/warning/message branch is mapped. |
| G5.2a | Partial/TODO | Error/warning tests should be mapped to validation branches. |
| G5.2b | Partial/TODO | Message/diagnostic branch coverage should be mapped. |
| G5.3 | Met | Numerical tests compare against known or parity expectations. |
| G5.4 | Met | Benchmark and replication directories support performance claims. |
| G5.4a | Met | Comparator packages are optional and opt-in. |
| G5.4b | Met | Benchmark reports are reproducible Markdown summaries. |
| G5.4c | Partial/TODO | Final manuscript outputs should be frozen and cross-referenced. |
| G5.5 | Met | Stochastic tests use fixed seeds where applicable. |
| G5.6 | Met | Simulation tests exercise model fitting and prediction recovery. |
| G5.6a | Met | Tests include Gaussian, positive continuous, and count workflows. |
| G5.6b | Partial/TODO | Expand multi-seed recovery evidence for review. |
| G5.7 | Met | Helper fixtures and known-truth CSV data support reproducible tests. |
| G5.8 | Met | Edge cases include missing margins, invalid copulas, and degenerate paths. |
| G5.8a | Partial/TODO | Add explicit zero-length/empty-data tests. |
| G5.8b | Met | Single-subject/pair edge cases are represented. |
| G5.8c | Met | Invalid and out-of-domain numeric inputs are tested. |
| G5.8d | Met | Invalid class/type inputs are tested for key functions. |
| G5.9 | Partial/TODO | Extended tests exist but need clearer separation from runtime-safe tests. |
| G5.9a | Partial/TODO | Long-running simulations should be guarded by an environment flag and documented. |
| G5.9b | Partial/TODO | Benchmark tests should remain opt-in and linked here. |
| G5.10 | Partial/TODO | Add stress tests for larger panels and harder optimisation cases. |
| G5.11 | Partial/TODO | Document platform/parallel stability if parallel benchmarks are used. |
| G5.11a | Partial/TODO | Note platform-specific numerical tolerances after CI matrix expansion. |
| G5.12 | Partial/TODO | Checks pass with vignette building disabled; full vignette-safe check remains a release gate. |

## Regression And Supervised Learning Standards

| ID | Status | Evidence / change needed |
|---|---|---|
| RE1.0 | Met | Primary model interface uses formulas. |
| RE1.1 | Met | `gamlss.longitudinal()` and `fit_longitudinal()` document formula inputs. |
| RE1.2 | Met | Formula inputs are converted through model-frame/model-matrix machinery. |
| RE1.3 | Partial/TODO | Subject/time metadata is retained; submitted row names are not guaranteed after grid expansion. |
| RE1.3a | Met | Accessor tests cover observed and expanded model-frame reconstruction. |
| RE1.4 | Met | Subject and time variables are explicit arguments. |
| RE2.0 | Met | Pre-processing policy is documented and implemented before model construction. |
| RE2.1 | Partial/TODO | Missingness handling is documented; response/predictor missingness policy is not user-selectable. |
| RE2.2 | Partial/TODO | Non-finite response/predictor behavior needs finer-grained docs/tests. |
| RE2.3 | Partial/TODO | Transformations are user-specified in formulas; default transformation/centering policy should be explicit. |
| RE2.4 | Met | Model-matrix construction and rank checks support standard predictor encoding. |
| RE2.4a | Met | Factor predictors use formula/model-matrix conversion. |
| RE2.4b | Partial/TODO | Rank diagnostics exist; aliased-design/response collinearity messaging should be expanded. |
| RE3.0 | Met | Optimisation controls and convergence state are exposed. |
| RE3.1 | Met | Convergence warnings are emitted and stored. |
| RE3.2 | Met | Hessian and variance-covariance fallbacks are represented in controls/object fields. |
| RE3.3 | Met | Diagnostics expose fit state, dependence estimates, and residual checks. |
| RE4.0 | Met | Returned objects have class `gamlss.longitudinal` and S3 methods. |
| RE4.1 | Met | Stored fields and accessors allow inspection without refitting. |
| RE4.2 | Met | `coef()` is implemented and tested. |
| RE4.3 | Met | `vcov()` is implemented and tested. |
| RE4.4 | Met | `confint()` is implemented and tested. |
| RE4.5 | Met | `summary()` is implemented. |
| RE4.6 | Met | `predict()` is implemented for response/distribution summaries. |
| RE4.7 | Met | `simulate()` is implemented and tested. |
| RE4.8 | Met | `logLik()` is implemented. |
| RE4.9 | Met | `formula()` is implemented. |
| RE4.10 | Met | `terms()` is implemented. |
| RE4.11 | Met | `nobs()` is implemented. |
| RE4.12 | Met | `fitted()` is implemented. |
| RE4.13 | Met | `residuals()` is implemented. |
| RE4.14 | Partial/TODO | Interval predictions exist; forecast interval behavior needs standards-specific tests. |
| RE4.15 | Partial/TODO | Prediction uncertainty exists for supported summaries; horizon/new-panel behavior needs clearer docs. |
| RE4.16 | Partial/TODO | Newdata handling exists; unseen factor levels/groups need explicit validation tests. |
| RE4.17 | Met | `model.frame()` exposes observed and expanded data views. |
| RE4.18 | Met | Convergence metadata is accessible from fitted objects. |
| RE5.0 | Met | Plot methods provide diagnostic and fitted-value workflows. |
| RE6.0 | Met | Prediction supports supplied new data. |
| RE6.1 | Met | Tests cover newdata prediction behavior. |
| RE6.2 | Met | Simulation helpers support longitudinal panels. |
| RE6.3 | Partial/TODO | Future-panel forecasting should be documented separately from ordinary newdata prediction. |
| RE7.0 | Partial/TODO | Regression tests are substantial but not yet mapped one-to-one to all regression standards. |
| RE7.0a | Partial/TODO | Add explicit noiseless response and noiseless predictor tests. |
| RE7.1 | Partial/TODO | Add full method-contract tests for every exported accessor. |
| RE7.1a | Partial/TODO | Subject/time retention is tested; row-name loss should be tested or documented explicitly. |
| RE7.2 | Met | Tests cover parameter recovery and prediction on simulated data. |
| RE7.3 | Met | Tests include categorical, smooth, and time-varying components. |
| RE7.4 | Partial/TODO | Add opt-in benchmark/regression-comparator tests. |

## Probability Distribution Standards

| ID | Status | Evidence / change needed |
|---|---|---|
| PD1.0 | Met | README cites distributional regression and copula references. |
| PD2.0 | Met | Marginal distributions are delegated to `gamlss.dist` family objects and that interface is documented. |
| PD3.0 | Met | Density/CDF/quantile/simulation operations use family objects and copula backends. |
| PD3.1 | Met | Marginal distribution functions follow `gamlss.dist` conventions. |
| PD3.2 | Met | Copula density/CDF/h-function parity is tested against `VineCopula` when available. |
| PD3.3 | Met | Quantile prediction paths are implemented for fitted marginal distributions. |
| PD3.4 | Partial/TODO | Numerical integration and t-copula approximations need a concise methods note. |
| PD3.5 | Partial/TODO | Count workflows are supported through families; discrete summation/truncation policy needs documentation. |
| PD3.5a | Partial/TODO | Add count-distribution edge-case and tail-behavior tests. |
| PD4.0 | Met | Copula backend has numerical parity tests against established implementations. |
| PD4.1 | Met | Simulations validate representative marginal/dependence workflows. |
| PD4.2 | Met | Distribution and copula parameter conversions are tested. |
| PD4.3 | Partial/TODO | Add tests comparing analytic and numerical alternatives where both are available. |
| PD4.4 | Met | Stochastic distribution tests use fixed seeds and tolerances. |
