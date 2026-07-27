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
- Not applicable: intentionally out of package scope, with rationale documented.

Latest local evidence, 2026-06-14 on Windows 11 x64 with R 4.4.1:
`devtools::document()`, `devtools::test(reporter = "summary")`, opt-in extended
`testthat::test_dir()`, `R CMD build .`, `R CMD check --as-cran --no-manual`,
`urlchecker::url_check()`, `spelling::spell_check_package()`, the new-user smoke
test, and the default JSS smoke replication all passed. The CRAN-style check had
0 ERRORs, 0 WARNINGs, and 3 NOTEs: new submission plus optional `gamlss2`
availability via `Additional_repositories`, local time verification, and
unavailable local pandoc for checking `README.md`/`NEWS.md`.

## General Standards

| ID | Status | Evidence / change needed |
|---|---|---|
| G1.0 | Met | README lists published GAMLSS, copula, and VineCopula references. |
| G1.1 | Met | README and CONTRIBUTING position the package relative to prior art. |
| G1.2 | Met | README includes a lifecycle statement. |
| G1.3 | Met | CONTRIBUTING defines the main statistical and workflow terms. |
| G1.4 | Met | User-facing functions are roxygen2 documented and Rd files are generated. |
| G1.4a | Met | Dense numerical/backend internals in `R/likelihood-*.R`, `R/model-fit.R`, `R/model-matrix.R`, `R/model-vcov.R`, `R/hessian-*.R`, and `R/copula-backend-*.R` have concise `@noRd` roxygen context. |
| G1.5 | Met | Reproduction scaffolds are under `inst/jss-replication`; the correlation-misspecification benchmark runner is under `paper/R/08-simulation-sensitivity-correlation-misspecification/standard-model-benchmarking`. |
| G1.6 | Met | Comparator helpers support GEE, GLMM, GAM, and GAMLSS baselines when optional dependencies are present. |
| G2.0 | Met | Length/scalar assertions exist for core controls and fitting inputs. |
| G2.0a | Met | CONTRIBUTING documents scalar/vector expectations for the main APIs. |
| G2.1 | Met | Type checks cover formulas, data, family objects/names, subject/time columns, and numeric controls. |
| G2.1a | Met | CONTRIBUTING documents supported data types for response, predictors, subject, and time. |
| G2.2 | Met | Univariate control arguments are restricted to length-one values where applicable. |
| G2.3 | Met | Character controls are constrained and case-sensitivity policy is centralized in CONTRIBUTING. |
| G2.3a | Met | `match.arg()` or equivalent is used for constrained choices. |
| G2.3b | Met | CONTRIBUTING states string-valued option arguments are case-sensitive unless documented otherwise. |
| G2.4 | Met | Time, subject, character, factor, and integer-control conversions are documented and tested. |
| G2.4a | Met | Integer control/index coercions are explicit and covered by validation or workflow tests. |
| G2.4b | Met | Numeric conversion is used for time/response/numeric covariate paths. |
| G2.4c | Met | Character conversion is used for identifiers and report labels. |
| G2.4d | Met | Factor and ordered-factor handling is documented and tested. |
| G2.4e | Met | Formula/model-matrix machinery handles factor-to-design conversion. |
| G2.5 | Met | README documents ordered/unordered factor behavior and formula/model-matrix routines enforce standard factor handling. |
| G2.6 | Partial/TODO | Main workflow is tabular; one-dimensional helper inputs need class pre-processing audit. |
| G2.7 | Met | Data-frame-like inputs are normalized to base data frames for fitting. |
| G2.8 | Met | Early pre-processing normalizes data, time, subject, formulas, and families. |
| G2.9 | Met | Character time conversion and character predictor treatment warn or are explicitly documented. |
| G2.10 | Met | Named column extraction avoids default drop behavior. |
| G2.11 | Met | CONTRIBUTING documents non-standard column-class policy; `test-p2-srr-input-policy.R` checks that custom predictor classes are rejected. |
| G2.12 | Met | List-columns are explicitly rejected and tested. |
| G2.13 | Met | Required columns, empty margins, and empty adjacent pairs are rejected. |
| G2.14 | Met | Missing-response, structural-missingness, and predictor-missingness policies are documented and tested. |
| G2.14a | Met | Response NaN/Inf and predictor NA/NaN/Inf behavior is explicitly checked and tested. |
| G2.14b | Met | Structural missing rows are expanded and observed/expanded model frames are tested. |
| G2.14c | Not applicable | The package intentionally does not perform statistical imputation; CONTRIBUTING documents this policy. |
| G2.15 | Met | Non-finite likelihood/optimisation states are checked and stored/warned. |
| G2.16 | Met | Response and predictor undefined-value policies are explicit and tested before fitting. |
| G3.0 | Partial/TODO | Tolerance checks are being adopted; complete numeric equality audit remains. |
| G3.1 | Met | Analytical and numerical Hessian paths are explicit with guarded fallbacks. |
| G3.1a | Met | Finite-difference controls and convergence diagnostics are accessible. |
| G4.0 | Met | Report output is Markdown and report paths are normalized to `.md`. |
| G5.0 | Partial/TODO | Known-truth fixture exists; add broader canonical/external benchmark datasets if appropriate. |
| G5.1 | Met | testthat and CI cover main workflow, diagnostics, prediction, simulation, and copula parity. |
| G5.2 | Partial/TODO | Many validation branches are tested; not every user-reachable stop/warning/message branch is mapped. |
| G5.2a | Met | The validation coverage map below links major user-facing errors and warnings to reviewer-oriented tests. |
| G5.2b | Partial/TODO | Message/diagnostic branch coverage should be mapped. |
| G5.3 | Met | Numerical tests compare against known or parity expectations. |
| G5.4 | Met | Benchmark and replication directories support performance claims. |
| G5.4a | Met | Comparator packages are optional and opt-in. |
| G5.4b | Met | Benchmark reports are reproducible Markdown summaries. |
| G5.4c | Partial/TODO | Final manuscript outputs should be frozen and cross-referenced. |
| G5.5 | Met | Stochastic tests use fixed seeds where applicable. |
| G5.6 | Met | Simulation tests exercise model fitting and prediction recovery. |
| G5.6a | Met | Tests include Gaussian, positive continuous, and count workflows. |
| G5.6b | Met | Opt-in extended tests include multi-seed recovery checks. |
| G5.7 | Met | Helper fixtures and known-truth CSV data support reproducible tests. |
| G5.8 | Met | Edge cases include missing margins, invalid copulas, and degenerate paths. |
| G5.8a | Met | Empty-data behavior is explicitly tested. |
| G5.8b | Met | Single-subject/pair edge cases are represented. |
| G5.8c | Met | Invalid and out-of-domain numeric inputs are tested. |
| G5.8d | Met | Invalid class/type inputs are tested for key functions. |
| G5.9 | Met | Extended stochastic/recovery tests are separated from routine tests. |
| G5.9a | Met | Long-running recovery tests are guarded by `GAMLSS_LONGITUDINAL_EXTENDED_TESTS`. |
| G5.9b | Met | Benchmark/stress checks are opt-in and linked from CONTRIBUTING. |
| G5.10 | Met | Opt-in stress tests cover dependence-strength scenarios. |
| G5.11 | Partial/TODO | Document platform/parallel stability if parallel benchmarks are used. |
| G5.11a | Partial/TODO | Note platform-specific numerical tolerances after CI matrix expansion. |
| G5.12 | Met | CONTRIBUTING documents extended-test conditions, runtime expectations, skips, and artifacts. |

## Regression And Supervised Learning Standards

| ID | Status | Evidence / change needed |
|---|---|---|
| RE1.0 | Met | Primary model interface uses formulas. |
| RE1.1 | Met | `gamlss_longitudinal()` documents formula inputs. |
| RE1.2 | Met | Formula inputs are converted through model-frame/model-matrix machinery. |
| RE1.3 | Met | Subject/time metadata is retained and row-name loss after expansion is documented and tested. |
| RE1.3a | Met | Accessor tests cover observed and expanded model-frame reconstruction. |
| RE1.4 | Met | Subject and time variables are explicit arguments. |
| RE2.0 | Met | Pre-processing policy is documented and implemented before model construction. |
| RE2.1 | Partial/TODO | Missingness handling is documented; response/predictor missingness policy is not user-selectable. |
| RE2.2 | Met | CONTRIBUTING documents non-finite response/predictor policy; `test-p2-srr-input-policy.R` checks response NaN/Inf and predictor NA/NaN/Inf. |
| RE2.3 | Met | CONTRIBUTING states that transformations are user-specified through formulas and no automatic centering/scaling is applied. |
| RE2.4 | Met | Model-matrix construction and rank checks support standard predictor encoding. |
| RE2.4a | Met | Factor predictors use formula/model-matrix conversion. |
| RE2.4b | Met | Rank-deficient/noiseless predictor cases warn and constant-response starting behavior is tested. |
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
| RE4.14 | Met | `test-p2-srr-regression-edge-cases.R` checks confidence intervals for future/new-subject panels. |
| RE4.15 | Met | `predict.gamlss.longitudinal()` documents uncertainty and new-panel behavior; standards tests cover new-panel intervals. |
| RE4.16 | Met | `test-p2-srr-regression-edge-cases.R` and `test-p2-srr-error-map.R` check unseen factor-level rejection. |
| RE4.17 | Met | `model.frame()` exposes observed and expanded data views. |
| RE4.18 | Met | Convergence metadata is accessible from fitted objects. |
| RE5.0 | Met | Plot methods provide diagnostic and fitted-value workflows. |
| RE6.0 | Met | Prediction supports supplied new data. |
| RE6.1 | Met | Tests cover newdata prediction behavior. |
| RE6.2 | Met | Simulation helpers support longitudinal panels. |
| RE6.3 | Met | CONTRIBUTING and `predict.gamlss.longitudinal()` document future panels as supplied-covariate `newdata` prediction, not time-series forecasting. |
| RE7.0 | Met | Tests cover exact/noiseless predictor relationships. |
| RE7.0a | Met | Rank-deficient exact predictor input is detected with a warning. |
| RE7.1 | Met | Tests cover constant-response/noiseless response behavior and accessor contracts. |
| RE7.1a | Met | Row/case retention and row-name reset behavior are explicitly tested. |
| RE7.2 | Met | Tests cover parameter recovery and prediction on simulated data. |
| RE7.3 | Met | Tests include categorical, smooth, and time-varying components. |
| RE7.4 | Met | Opt-in extended tests exercise benchmark/recovery checks under an environment flag. |

## Probability Distribution Standards

| ID | Status | Evidence / change needed |
|---|---|---|
| PD1.0 | Met | README cites distributional regression and copula references. |
| PD2.0 | Met | Marginal distributions are delegated to `gamlss.dist` family objects and that interface is documented. |
| PD3.0 | Met | Density/CDF/quantile/simulation operations use family objects and copula backends. |
| PD3.1 | Met | Marginal distribution functions follow `gamlss.dist` conventions. |
| PD3.2 | Met | Copula density/CDF/h-function parity is tested against `VineCopula` when available. |
| PD3.3 | Met | Quantile prediction paths are implemented for fitted marginal distributions. |
| PD3.4 | Met | CONTRIBUTING documents t-copula integration stability assumptions and tests cover finite CDF output. |
| PD3.5 | Met | CONTRIBUTING documents count-family probability handling through finite `gamlss.dist` p/q/d calls. |
| PD3.5a | Met | Count-tail p/q/d consistency is tested for representative Poisson probabilities. |
| PD4.0 | Met | Copula backend has numerical parity tests against established implementations. |
| PD4.1 | Met | Simulations validate representative marginal/dependence workflows. |
| PD4.2 | Met | Distribution and copula parameter conversions are tested. |
| PD4.3 | Met | Tests compare analytical Gaussian copula derivatives with finite-difference alternatives. |
| PD4.4 | Met | Stochastic distribution tests use fixed seeds and tolerances. |

## Validation Coverage Map

| Validation branch | Evidence |
|---|---|
| Empty datasets and missing required columns | `tests/testthat/test-p2-srr-error-map.R` |
| Duplicate subject/time rows and impossible missingness patterns | `tests/testthat/test-p2-srr-error-map.R` |
| List-columns, matrix/data-frame columns, and non-standard predictor classes | `tests/testthat/test-p2-srr-input-policy.R` |
| Response `NaN`/`Inf` and predictor `NA`/`NaN`/`Inf` | `tests/testthat/test-p2-srr-input-policy.R` |
| Character time conversion and character predictor warnings | `tests/testthat/test-p2-srr-input-policy.R` and `tests/testthat/test-p2-srr-error-map.R` |
| Unseen factor levels in `newdata` | `tests/testthat/test-p2-srr-regression-edge-cases.R` and `tests/testthat/test-p2-srr-error-map.R` |
| Future/new-subject prediction panels and interval outputs | `tests/testthat/test-p2-srr-regression-edge-cases.R` |
