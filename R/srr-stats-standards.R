#' srr_stats
#'
#' Reviewer-facing compliance notes for the rOpenSci statistical software
#' standards. These tags are intentionally conservative: standards marked
#' `@srrstats` have code, documentation, or tests that a reviewer can inspect;
#' standards marked `@srrstatsTODO` are partially addressed but still need
#' stronger evidence before submission.
#'
#' The human-readable crosswalk is maintained in
#' `inst/standards/ropensci-srr-compliance.md`.
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstats {G1.0} Met: README lists primary published GAMLSS, copula, and VineCopula references.
#' @srrstats {G1.1} Met: README and CONTRIBUTING document prior art and current package positioning.
#' @srrstats {G1.2} Met: README includes a lifecycle statement for current and anticipated package development.
#' @srrstats {G1.3} Met: CONTRIBUTING defines the main statistical and workflow terms used by the package.
#' @srrstats {G1.4} Met: exported user-facing functions are documented with roxygen2 and regenerated into Rd files.
#' @srrstats {G1.4a} Met: dense numerical/backend internals have concise `@noRd` roxygen context.
#' @srrstats {G1.5} Met: replication and benchmark scaffolds live under `inst/jss-replication` and `inst/benchmarks`.
#' @srrstats {G1.6} Met: opt-in benchmark helpers compare against GEE, GLMM, GAM, and GAMLSS baselines when optional packages are installed.
#' @srrstats {G2.0} Met: scalar controls, family names, copula choices, data columns, and formula inputs are asserted before fitting.
#' @srrstats {G2.0a} Met: CONTRIBUTING documents scalar/vector expectations for the main fitting and workflow APIs.
#' @srrstats {G2.1} Met: type checks cover formulas, data, family objects/names, subject/time columns, and numeric controls.
#' @srrstats {G2.1a} Met: CONTRIBUTING documents supported data types for response, predictors, subject, and time columns.
#' @srrstats {G2.2} Met: univariate control arguments are validated as length-one where applicable.
#' @srrstats {G2.3} Met: character controls are constrained and the case-sensitivity policy is centralized in CONTRIBUTING.
#' @srrstats {G2.3a} Met: copula families, comparison methods, and mode arguments are restricted with `match.arg()` or equivalent.
#' @srrstats {G2.3b} Met: CONTRIBUTING states that string-valued option arguments are case-sensitive unless documented otherwise.
#' @srrstats {G2.4} Met: conversion behavior for time, subject, characters, factors, and integer controls is documented and tested.
#' @srrstats {G2.4a} Met: integer control/index coercions are explicit and covered by validation or existing workflow tests.
#' @srrstats {G2.4b} Met: time, response, and numeric covariate paths use explicit numeric conversion where needed.
#' @srrstats {G2.4c} Met: subject/time identifiers and reporting labels use explicit character conversion.
#' @srrstats {G2.4d} Met: factor and ordered-factor handling is documented and covered by input-policy tests.
#' @srrstats {G2.4e} Met: model-matrix construction handles factor-to-design-matrix conversion through R's formula machinery.
#' @srrstats {G2.5} Met: README documents ordered/unordered factor policy and formula/model-matrix routines enforce standard factor handling.
#' @srrstatsTODO {G2.6} Partial: primary interface is tabular; one-dimensional helper inputs need an audit for class-preserving pre-processing.
#' @srrstats {G2.7} Met: fitting accepts data-frame-like inputs and normalizes them to a base data frame.
#' @srrstats {G2.8} Met: early pre-processing normalizes data, time, subject, formulas, and marginal family objects before lower-level routines run.
#' @srrstats {G2.9} Met: character time conversion and character-predictor treatment warn or are explicitly documented.
#' @srrstats {G2.10} Met: column extraction uses explicit names and preserves data-frame behavior.
#' @srrstats {G2.11} Met: non-standard column classes are documented and custom predictor classes are rejected in input-policy tests.
#' @srrstats {G2.12} Met: list-columns are explicitly rejected before fitting and covered by input-policy tests.
#' @srrstats {G2.13} Met: fitting rejects missing required columns and margins/pairs with no observed data.
#' @srrstats {G2.14} Met: missing-response, structural-missingness, and predictor-missingness policies are documented and tested.
#' @srrstats {G2.14a} Met: response NaN/Inf and predictor NA/NaN/Inf behavior is explicitly checked and tested.
#' @srrstats {G2.14b} Met: structural missing rows are expanded and observed/expanded model frames are tested.
#' @srrstats {G2.15} Met: non-finite optimisation and likelihood values are checked and routed through convergence/warning paths.
#' @srrstats {G2.16} Met: response and predictor undefined-value policies are explicit and tested before fitting.
#' @srrstatsTODO {G3.0} Partial: tolerance comparisons are being adopted; complete floating-point equality audit remains.
#' @srrstats {G3.1} Met: analytical and numerical Hessian paths are explicit and numerical fallbacks are guarded.
#' @srrstats {G3.1a} Met: finite-difference Hessian controls and convergence diagnostics are exposed in returned objects.
#' @srrstats {G4.0} Met: user-facing report output is Markdown and report paths are normalized to `.md`.
#' @srrstatsTODO {G5.0} Partial: known-truth fixture added; broader comparison to canonical external datasets remains useful.
#' @srrstats {G5.1} Met: testthat tests and CI cover the main workflow, diagnostics, simulation, prediction, and copula parity.
#' @srrstatsTODO {G5.2} Partial: many error/warning paths are tested, but not every user-reachable stop/warning/message branch.
#' @srrstats {G5.2a} Met: the reviewer crosswalk maps major user-facing errors and warnings to validation tests.
#' @srrstatsTODO {G5.2b} Partial: message/diagnostic branch coverage needs the same mapping.
#' @srrstats {G5.3} Met: representative tests compare numerical outputs with known or parity expectations.
#' @srrstats {G5.4} Met: benchmark and replication directories support performance and publication claims.
#' @srrstats {G5.4a} Met: benchmark helpers use optional comparator packages and explicit opt-in execution.
#' @srrstats {G5.4b} Met: benchmark reporting writes reproducible Markdown summaries.
#' @srrstatsTODO {G5.4c} Partial: paper-result replication scripts exist, but final manuscript outputs should be frozen and cross-referenced.
#' @srrstats {G5.5} Met: tests use fixed seeds for stochastic workflows where applicable.
#' @srrstats {G5.6} Met: simulation tests exercise model fitting and prediction recovery behavior.
#' @srrstats {G5.6a} Met: tests include representative Gaussian, positive continuous, and count workflows.
#' @srrstats {G5.6b} Met: opt-in extended tests include multi-seed recovery checks guarded by `GAMLSS_LONGITUDINAL_EXTENDED_TESTS`.
#' @srrstats {G5.7} Met: helper fixtures and known-truth CSV data support reproducible tests.
#' @srrstats {G5.8} Met: tests cover edge cases including missing margins, invalid copula choices, and degenerate paths.
#' @srrstats {G5.8a} Met: empty-data behavior is explicitly tested.
#' @srrstats {G5.8b} Met: single-subject/pair edge cases are represented in workflow tests.
#' @srrstats {G5.8c} Met: tests cover invalid and out-of-domain numeric inputs.
#' @srrstats {G5.8d} Met: invalid class/type inputs are tested for key user-facing functions.
#' @srrstats {G5.9} Met: extended stochastic/recovery tests are separated from routine tests.
#' @srrstats {G5.9a} Met: long-running recovery tests are guarded by `GAMLSS_LONGITUDINAL_EXTENDED_TESTS`.
#' @srrstats {G5.9b} Met: benchmark/stress checks are opt-in and linked from CONTRIBUTING.
#' @srrstats {G5.10} Met: opt-in stress tests cover dependence-strength scenarios.
#' @srrstatsTODO {G5.11} Partial: parallel/platform stability evidence should be documented if parallel benchmarks are used.
#' @srrstatsTODO {G5.11a} Partial: platform-specific numerical tolerances should be noted after CI matrix expansion.
#' @srrstats {G5.12} Met: CONTRIBUTING documents extended-test conditions, runtime expectations, skips, and artifacts.
#'
#' @srrstats {RE1.0} Met: primary interface uses R formulas for marginal and dependence components.
#' @srrstats {RE1.1} Met: `gamlss_longitudinal()` documents formula inputs.
#' @srrstats {RE1.2} Met: formulas are converted through model-frame/model-matrix machinery.
#' @srrstats {RE1.3} Met: subject/time metadata is retained and row-name loss after expansion is documented and tested.
#' @srrstats {RE1.3a} Met: accessor tests cover model-frame and observed/expanded data reconstruction.
#' @srrstats {RE1.4} Met: time and subject variables are explicit arguments rather than inferred hidden state.
#' @srrstats {RE2.0} Met: preprocessing is documented in README and implemented before model construction.
#' @srrstatsTODO {RE2.1} Partial: missingness handling is documented; user-selectable response/predictor missingness policies are not exposed.
#' @srrstats {RE2.2} Met: non-finite response/predictor policy is documented and tested.
#' @srrstats {RE2.3} Met: transformation and centering policy is explicit; transformations are user-specified in formulas.
#' @srrstats {RE2.4} Met: model-matrix construction and rank checks support standard predictor encoding.
#' @srrstats {RE2.4a} Met: factor predictors are handled through formula/model-matrix conversion.
#' @srrstats {RE2.4b} Met: rank-deficient/noiseless predictor cases warn and constant-response starting behavior is tested.
#' @srrstats {RE3.0} Met: optimisation controls and convergence state are exposed.
#' @srrstats {RE3.1} Met: convergence warnings are emitted and stored when optimisation fails.
#' @srrstats {RE3.2} Met: Hessian and variance-covariance fallbacks are documented through controls and object fields.
#' @srrstats {RE3.3} Met: diagnostics expose fit state, dependence estimates, and residual checks.
#' @srrstats {RE4.0} Met: returned objects have class `gamlss.longitudinal` and standard S3 methods.
#' @srrstats {RE4.1} Met: fitted model objects can be inspected without refitting via stored fields and accessors.
#' @srrstats {RE4.2} Met: `coef()` is implemented and tested.
#' @srrstats {RE4.3} Met: `vcov()` is implemented and tested.
#' @srrstats {RE4.4} Met: `confint()` is implemented and tested.
#' @srrstats {RE4.5} Met: `summary()` is implemented for fitted objects.
#' @srrstats {RE4.6} Met: `predict()` is implemented for response/distribution summaries.
#' @srrstats {RE4.7} Met: `simulate()` is implemented and tested.
#' @srrstats {RE4.8} Met: `logLik()` is implemented.
#' @srrstats {RE4.9} Met: `formula()` is implemented.
#' @srrstats {RE4.10} Met: `terms()` is implemented.
#' @srrstats {RE4.11} Met: `nobs()` is implemented.
#' @srrstats {RE4.12} Met: `fitted()` is implemented.
#' @srrstats {RE4.13} Met: `residuals()` is implemented.
#' @srrstats {RE4.14} Met: standards tests cover confidence intervals for future/new-subject panels.
#' @srrstats {RE4.15} Met: prediction uncertainty and new-panel behavior are documented and tested.
#' @srrstats {RE4.16} Met: unseen factor levels in `newdata` are explicitly rejected and tested.
#' @srrstats {RE4.17} Met: `model.frame()` is implemented for observed and expanded data views.
#' @srrstats {RE4.18} Met: convergence metadata is accessible from fitted objects.
#' @srrstats {RE5.0} Met: plotting methods provide diagnostic and fitted-value workflows.
#' @srrstats {RE6.0} Met: prediction supports supplied new data through model-matrix reconstruction.
#' @srrstats {RE6.1} Met: prediction tests cover newdata behavior.
#' @srrstats {RE6.2} Met: simulation helpers support longitudinal panels.
#' @srrstats {RE6.3} Met: future panels are documented as supplied-covariate `newdata` prediction rather than time-series forecasting.
#' @srrstats {RE7.0} Met: tests cover exact/noiseless predictor relationships.
#' @srrstats {RE7.0a} Met: rank-deficient exact predictor input is detected with a warning.
#' @srrstats {RE7.1} Met: tests cover constant-response/noiseless response behavior and accessor contracts.
#' @srrstats {RE7.1a} Met: row/case retention and row-name reset behavior are explicitly tested.
#' @srrstats {RE7.2} Met: tests cover parameter recovery and prediction on simulated data.
#' @srrstats {RE7.3} Met: tests include categorical, smooth, and time-varying model components.
#' @srrstats {RE7.4} Met: opt-in extended tests exercise benchmark/recovery checks under an environment flag.
#'
#' @srrstats {PD1.0} Met: README cites primary distributional regression and copula references.
#' @srrstats {PD2.0} Met: package delegates marginal distributions to `gamlss.dist` family objects and documents that interface.
#' @srrstats {PD3.0} Met: density/CDF/quantile/simulation operations use family objects and copula backends.
#' @srrstats {PD3.1} Met: marginal distribution functions are accessed through `gamlss.dist` conventions.
#' @srrstats {PD3.2} Met: copula density/CDF/h-function parity is tested against `VineCopula` where available.
#' @srrstats {PD3.3} Met: quantile prediction paths are implemented for fitted marginal distributions.
#' @srrstats {PD3.4} Met: CONTRIBUTING documents t-copula integration stability assumptions and tests cover finite CDF output.
#' @srrstats {PD3.5} Met: CONTRIBUTING documents count-family probability handling through finite `gamlss.dist` p/q/d calls.
#' @srrstats {PD3.5a} Met: count-tail p/q/d consistency is tested for representative Poisson probabilities.
#' @srrstats {PD4.0} Met: copula backend has numerical parity tests against established implementations.
#' @srrstats {PD4.1} Met: simulations validate representative marginal/dependence workflows.
#' @srrstats {PD4.2} Met: distribution and copula parameter conversions are tested.
#' @srrstats {PD4.3} Met: tests compare analytical Gaussian copula derivatives with finite-difference alternatives.
#' @srrstats {PD4.4} Met: stochastic distribution tests use fixed seeds and tolerances.
#'
#' @noRd
NULL

#' NA_standards
#'
#' @srrstatsNA {G2.14c} Not applicable: the package intentionally does not perform statistical imputation; CONTRIBUTING documents this policy.
#'
#' @noRd
NULL
