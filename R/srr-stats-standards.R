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
#' @srrstatsTODO {G1.1} Partial: README documents prior art and package positioning; final novelty/improvement claim should be revisited with the manuscript.
#' @srrstats {G1.2} Met: README includes a lifecycle statement for current and anticipated package development.
#' @srrstatsTODO {G1.3} Partial: key terms are documented across README and help files; add a compact glossary before review.
#' @srrstats {G1.4} Met: exported user-facing functions are documented with roxygen2 and regenerated into Rd files.
#' @srrstatsTODO {G1.4a} Partial: many internal helpers have roxygen comments, but older dense internals still need full `@noRd` documentation.
#' @srrstats {G1.5} Met: replication and benchmark scaffolds live under `inst/jss-replication` and `inst/benchmarks`.
#' @srrstats {G1.6} Met: opt-in benchmark helpers compare against GEE, GLMM, GAM, and GAMLSS baselines when optional packages are installed.
#' @srrstats {G2.0} Met: scalar controls, family names, copula choices, data columns, and formula inputs are asserted before fitting.
#' @srrstatsTODO {G2.0a} Partial: README documents core length expectations; Rd argument documentation should be checked for every vector/scalar input.
#' @srrstats {G2.1} Met: type checks cover formulas, data, family objects/names, subject/time columns, and numeric controls.
#' @srrstatsTODO {G2.1a} Partial: README documents common data types; complete secondary type documentation should be added to all relevant help pages.
#' @srrstats {G2.2} Met: univariate control arguments are validated as length-one where applicable.
#' @srrstatsTODO {G2.3} Partial: character controls are constrained; case-sensitivity policy should be centralized.
#' @srrstats {G2.3a} Met: copula families, comparison methods, and mode arguments are restricted with `match.arg()` or equivalent.
#' @srrstatsTODO {G2.3b} Partial: supported string values are deterministic; document which arguments are case-sensitive.
#' @srrstatsTODO {G2.4} Partial: conversions are implemented for fitting paths; not every conversion class has explicit tests.
#' @srrstatsTODO {G2.4a} Partial: integer coercions are limited and should be audited where indices/control values are accepted.
#' @srrstats {G2.4b} Met: time, response, and numeric covariate paths use explicit numeric conversion where needed.
#' @srrstats {G2.4c} Met: subject/time identifiers and reporting labels use explicit character conversion.
#' @srrstatsTODO {G2.4d} Partial: factor handling is documented, but explicit factor conversion tests should be broadened.
#' @srrstats {G2.4e} Met: model-matrix construction handles factor-to-design-matrix conversion through R's formula machinery.
#' @srrstats {G2.5} Met: README documents ordered/unordered factor policy and formula/model-matrix routines enforce standard factor handling.
#' @srrstatsTODO {G2.6} Partial: primary interface is tabular; one-dimensional helper inputs need an audit for class-preserving pre-processing.
#' @srrstats {G2.7} Met: fitting accepts data-frame-like inputs and normalizes them to a base data frame.
#' @srrstats {G2.8} Met: early pre-processing normalizes data, time, subject, formulas, and marginal family objects before lower-level routines run.
#' @srrstatsTODO {G2.9} Partial: some conversions are documented; lossy conversions should consistently warn or be explicitly declared.
#' @srrstats {G2.10} Met: column extraction uses explicit names and preserves data-frame behavior.
#' @srrstatsTODO {G2.11} Partial: non-standard column classes such as `units` are not yet explicitly documented or tested.
#' @srrstatsTODO {G2.12} Partial: list-column behavior is not part of the main workflow and should be explicitly rejected or tested.
#' @srrstats {G2.13} Met: fitting rejects missing required columns and margins/pairs with no observed data.
#' @srrstatsTODO {G2.14} Partial: missing responses and structural gaps are handled; user-selectable missingness policy is not implemented.
#' @srrstatsTODO {G2.14a} Partial: missing-value checks exist; distinct NA/NaN/Inf policies need clearer docs/tests.
#' @srrstatsTODO {G2.14b} Partial: structural missing rows are expanded; dropped/retained observed-case behavior needs more tests.
#' @srrstatsTODO {G2.14c} Partial: imputation is not provided; documentation should explicitly say so for all workflows.
#' @srrstats {G2.15} Met: non-finite optimisation and likelihood values are checked and routed through convergence/warning paths.
#' @srrstatsTODO {G2.16} Partial: undefined value behavior is mostly inherited from model matrices and families; add explicit user policy.
#' @srrstatsTODO {G3.0} Partial: tolerance comparisons are being adopted; complete floating-point equality audit remains.
#' @srrstats {G3.1} Met: analytical and numerical Hessian paths are explicit and numerical fallbacks are guarded.
#' @srrstats {G3.1a} Met: finite-difference Hessian controls and convergence diagnostics are exposed in returned objects.
#' @srrstats {G4.0} Met: user-facing report output is Markdown and report paths are normalized to `.md`.
#' @srrstatsTODO {G5.0} Partial: known-truth fixture added; broader comparison to canonical external datasets remains useful.
#' @srrstats {G5.1} Met: testthat tests and CI cover the main workflow, diagnostics, simulation, prediction, and copula parity.
#' @srrstatsTODO {G5.2} Partial: many error/warning paths are tested, but not every user-reachable stop/warning/message branch.
#' @srrstatsTODO {G5.2a} Partial: error/warning tests should be mapped systematically to validation branches.
#' @srrstatsTODO {G5.2b} Partial: message/diagnostic branch coverage needs the same mapping.
#' @srrstats {G5.3} Met: representative tests compare numerical outputs with known or parity expectations.
#' @srrstats {G5.4} Met: benchmark and replication directories support performance and publication claims.
#' @srrstats {G5.4a} Met: benchmark helpers use optional comparator packages and explicit opt-in execution.
#' @srrstats {G5.4b} Met: benchmark reporting writes reproducible Markdown summaries.
#' @srrstatsTODO {G5.4c} Partial: paper-result replication scripts exist, but final manuscript outputs should be frozen and cross-referenced.
#' @srrstats {G5.5} Met: tests use fixed seeds for stochastic workflows where applicable.
#' @srrstats {G5.6} Met: simulation tests exercise model fitting and prediction recovery behavior.
#' @srrstats {G5.6a} Met: tests include representative Gaussian, positive continuous, and count workflows.
#' @srrstatsTODO {G5.6b} Partial: multi-seed recovery evidence exists in simulation scaffolds but should be expanded for review.
#' @srrstats {G5.7} Met: helper fixtures and known-truth CSV data support reproducible tests.
#' @srrstats {G5.8} Met: tests cover edge cases including missing margins, invalid copula choices, and degenerate paths.
#' @srrstatsTODO {G5.8a} Partial: zero-length and empty-data behavior should be explicitly tested.
#' @srrstats {G5.8b} Met: single-subject/pair edge cases are represented in workflow tests.
#' @srrstats {G5.8c} Met: tests cover invalid and out-of-domain numeric inputs.
#' @srrstats {G5.8d} Met: invalid class/type inputs are tested for key user-facing functions.
#' @srrstatsTODO {G5.9} Partial: extended tests exist but should be clearly separable from CRAN/runtime-safe tests.
#' @srrstatsTODO {G5.9a} Partial: long-running simulations should be guarded by an environment flag and documented.
#' @srrstatsTODO {G5.9b} Partial: benchmark tests should remain opt-in and linked to the compliance crosswalk.
#' @srrstatsTODO {G5.10} Partial: stress tests for larger panels and harder optimisation cases should be added.
#' @srrstatsTODO {G5.11} Partial: parallel/platform stability evidence should be documented if parallel benchmarks are used.
#' @srrstatsTODO {G5.11a} Partial: platform-specific numerical tolerances should be noted after CI matrix expansion.
#' @srrstatsTODO {G5.12} Partial: package checks pass with vignette building disabled; full vignette-safe R CMD check remains a release gate.
#'
#' @srrstats {RE1.0} Met: primary interface uses R formulas for marginal and dependence components.
#' @srrstats {RE1.1} Met: `gamlss.longitudinal()` and `fit_longitudinal()` document formula inputs.
#' @srrstats {RE1.2} Met: formulas are converted through model-frame/model-matrix machinery.
#' @srrstatsTODO {RE1.3} Partial: subject/time metadata is retained; submitted row names are not guaranteed after grid expansion.
#' @srrstats {RE1.3a} Met: accessor tests cover model-frame and observed/expanded data reconstruction.
#' @srrstats {RE1.4} Met: time and subject variables are explicit arguments rather than inferred hidden state.
#' @srrstats {RE2.0} Met: preprocessing is documented in README and implemented before model construction.
#' @srrstatsTODO {RE2.1} Partial: missingness handling is documented; user-selectable response/predictor missingness policies are not exposed.
#' @srrstatsTODO {RE2.2} Partial: non-finite response/predictor behavior needs finer-grained documentation and tests.
#' @srrstatsTODO {RE2.3} Partial: transformations are user-specified in formulas; default transformation/centering policy should be explicit.
#' @srrstats {RE2.4} Met: model-matrix construction and rank checks support standard predictor encoding.
#' @srrstats {RE2.4a} Met: factor predictors are handled through formula/model-matrix conversion.
#' @srrstatsTODO {RE2.4b} Partial: predictor-rank diagnostics exist; response collinearity and aliased-design messaging should be expanded.
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
#' @srrstatsTODO {RE4.14} Partial: interval predictions exist; forecasting interval behavior needs explicit standards tests.
#' @srrstatsTODO {RE4.15} Partial: prediction uncertainty is available for supported summaries; horizon and new-panel behavior need clearer docs.
#' @srrstatsTODO {RE4.16} Partial: newdata handling is implemented; unseen factor levels/groups need explicit validation tests.
#' @srrstats {RE4.17} Met: `model.frame()` is implemented for observed and expanded data views.
#' @srrstats {RE4.18} Met: convergence metadata is accessible from fitted objects.
#' @srrstats {RE5.0} Met: plotting methods provide diagnostic and fitted-value workflows.
#' @srrstats {RE6.0} Met: prediction supports supplied new data through model-matrix reconstruction.
#' @srrstats {RE6.1} Met: prediction tests cover newdata behavior.
#' @srrstats {RE6.2} Met: simulation helpers support longitudinal panels.
#' @srrstatsTODO {RE6.3} Partial: forecast-style future panels should be documented separately from ordinary newdata prediction.
#' @srrstatsTODO {RE7.0} Partial: regression tests are substantial, but not yet mapped one-to-one against every regression standard.
#' @srrstatsTODO {RE7.0a} Partial: add explicit noiseless response and noiseless predictor tests.
#' @srrstatsTODO {RE7.1} Partial: accessor tests exist; add full method-contract tests for every exported accessor.
#' @srrstatsTODO {RE7.1a} Partial: row/case retention tests exist for subject/time; row-name loss should be explicitly tested.
#' @srrstats {RE7.2} Met: tests cover parameter recovery and prediction on simulated data.
#' @srrstats {RE7.3} Met: tests include categorical, smooth, and time-varying model components.
#' @srrstatsTODO {RE7.4} Partial: add formal benchmark/regression-comparator tests under an opt-in flag.
#'
#' @srrstats {PD1.0} Met: README cites primary distributional regression and copula references.
#' @srrstats {PD2.0} Met: package delegates marginal distributions to `gamlss.dist` family objects and documents that interface.
#' @srrstats {PD3.0} Met: density/CDF/quantile/simulation operations use family objects and copula backends.
#' @srrstats {PD3.1} Met: marginal distribution functions are accessed through `gamlss.dist` conventions.
#' @srrstats {PD3.2} Met: copula density/CDF/h-function parity is tested against `VineCopula` where available.
#' @srrstats {PD3.3} Met: quantile prediction paths are implemented for fitted marginal distributions.
#' @srrstatsTODO {PD3.4} Partial: numerical integration and t-copula approximations need a concise methods note.
#' @srrstatsTODO {PD3.5} Partial: discrete/count workflows are supported through families; summation/truncation policy needs explicit documentation.
#' @srrstatsTODO {PD3.5a} Partial: add tests for count-distribution edge cases and tail behaviour.
#' @srrstats {PD4.0} Met: copula backend has numerical parity tests against established implementations.
#' @srrstats {PD4.1} Met: simulations validate representative marginal/dependence workflows.
#' @srrstats {PD4.2} Met: distribution and copula parameter conversions are tested.
#' @srrstatsTODO {PD4.3} Partial: add tests comparing analytic and numerical alternatives where both are available.
#' @srrstats {PD4.4} Met: stochastic distribution tests use fixed seeds and tolerances.
#'
#' @noRd
NULL
