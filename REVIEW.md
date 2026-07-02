# Reviewer Guide

This file is designed to give an overview of the `gamlss.longitudinal` codebase for a reviewer or contributor.

We cover where all the various bits of code live, what they do, and the naming scheme for files (Code map), the implemented testing framework for the package including replication for claimed results and mapping to rOpenSci standards (Testing:*), traceability for key functions used in statistical results (Method traceability), and the most recent result of tests (Testing output).

## Code map

All core functionality is contained in the `R/` folders and naming conventions are structured to guide review. Naming follows the convention of tiered naming to group common files, generally the prefix identifies the area the code relates to (e.g. optimizer-) and the suffixes relate the the specific purpose of the code / review question of interest (-reporting.R).

##### High level orchestration
`model-fit.R`: exports the main fitting function `gamlss_longitudinal()` and kicks off the high level routines in `model-fit-entrypoint.R` which structures the overall workflow into three parts:
1. **Preparing for fit**: Preprocessing and validation of input data and formulas, model matrix creation and setting starting values for covariates and optimizer controls. This is done by calling workflow_fn = .gl_prepare_fit_workflow (in `model-fit-workflow.R`)
2. **Optimization / Fit**: Core optimization loop which calls optimizer_fn = .gl_run_prepared_fit_optimizer (in `model-fit-optimizer.R` and `model-optimizer-dispatch*.R`),
3. **Finalization**: Finalization of the fit including collecting all final outputs and model information, as well as kicking off variance/covariance calculations if requested. The call is to finalize_fn = .gl_finalize_prepared_fit (in `model-fit-finalize.R`)

We now outline the files / functions involved in each step and their broad function for reference. 

##### 1. Preparing for fit

**Input validation** and transformation to feed into optimiser:
- `model-fit-setup.R`, `model-preprocess.R`, `model-controls-fit.R`, `model-controls-step.R`,  `model-controls-budget.R`, `model-time.R`, `model-formulas.R`,  `model-missing-panels.R`, `model-column-policy.R`, and `model-data-shape.R`: these handle checks on all the various inputs (and related errors) and transforming these inputs if needed to different formats. Model control inputs for CG / RS, data structure and time / subject identifier provided, missingness checks like 100% missing timepoints or no complete pairs in ajacent timepoints and more. These map to many of the srr standards with regards to input validation and beyond. 

**Model matrix construction** from inputted formulas and dataset, and 
- `model-matrix.R`, `model-matrix-parameter-data.R`, `model-matrix-fixed.R`, `model-matrix-smooth.R`, `model-matrix-data.R`, `model-matrix-formulas.R`,  `model-matrix-bundle.R`, and `model-eta.R`: fixed/smooth design matrix construction from formulas.

**Starting values** initialisation for covariates:
- `model-starting-values.R`, `model-starting-copula.R`, `model-starting-margin.R`, `model-starting-moments.R`, `model-starting-state.R`, `model-parameter-transforms.R` (link to non-link transform), and `model-warm-start.R` (and helpers): initial values for margin and copula covariates, parameter transform helpers, and warm start method.

##### 2. Optimization / fit

Optimizer workflow:
- Model outer loops are run by `optimizer-rs-runner.R` and `optimizer-cg-runner.R` which call all the key helpers for iterations.

The core loop for RS calls:
- .gl_run_rs_parameter_iterations (`optimizer-rs-parameter-loop.R`) which controls the loop for the inner iterations and checks stopping criteria
  - This calls .gl_run_rs_inner_iteration_step (`optimizer-rs-inner-iteration.R`) to run the inner GLM loop for each parameter 
    - This calls .gl_run_rs_inner_backfitting_step (`optimizer-rs-backfitting.R`) to actually run the backfitting 
      - This calls .gl-rs-backfitting-iteration (`optimizer-rs-backfitting-iteration.R`) to do the actual calculation which is where the maths is implemented for the backfitting (**KEY REVIEW POINT**)

- Other functions of interest may be:
  - .gl_should_continue_rs_outer_loop (optimizer-rs-loop-stopping.R) to control outer stopping behaviour
  - .gl_should_continue_rs_inner_loop (optimizer-rs-loop-stopping.R) controls inner stopping behaviour
  - .gl_update_rs_outer_iteration_state (optimizer-rs-loop-stopping.R) to update all the parameters for the outer loop setup like step size, counters, inputs to stopping criteria test
  - .gl_build_rs_optimizer_result (optimizer-rs-result.R) which collects all the outputs from the final finished iterations
  - .gl_evaluate_rs_iteration_likelihood_state (optimizer-rs-likelihood.R) calls functions to calculate current likelihood, eta, derivatives
  - .gl_evaluate_rs_parameter_score_state (optimizer-rs-likelihood.R)

CG:
- .gl_ensure_cg_hessian_available (optimizer-cg-hessian-availability.R)
- .gl_build_cg_runtime_helpers (optimizer-cg-runtime.R)
- .gl_initialize_cg_optimizer_state (optimizer-cg-state-start.R)
- .gl_evaluate_cg_iteration_start (optimizer-cg-state-start.R)
- .gl_prepare_cg_curvature_line_search_state (optimizer-cg-curvature-search.R)
- .gl_apply_cg_line_search_diagnostics_state (optimizer-cg-line-search.R)
- .gl_apply_cg_stop_request_state (optimizer-cg-stop-request.R)
- .gl_finalize_cg_optimizer_state (optimizer-cg-finalization.R)
- .gl_report_cg_optimizer_start (optimizer-cg-reporting.R)
- .gl_build_cg_optimizer_result (optimizer-cg-result.R)
- .gl_prepare_cg_runtime_state (optimizer-cg-runtime-state.R)
- .gl_run_cg_outer_iteration_step (optimizer-cg-outer-iteration.R)
- .gl_run_cg_outer_loop (optimizer-cg-outer-loop.R)




- `optimizer-state.R`, `optimizer-stop-criteria.R`, `optimizer-rs-loop-stopping.R`, and `optimizer-rs-reporting.R`: shared optimizer state and initial values setting (including automated stopping criteria selection if not specified), output printing and implementation of stopping rules.
- `optimizer-rs-*.R`: RS design caches, likelihood contexts, score assembly,
       backfitting, smoothing-parameter updates, acceptance/backtracking,
       diagnostics, and RS outer-loop orchestration. The main score path is split
       between `optimizer-rs-likelihood.R`,
       `optimizer-rs-likelihood-context.R`, `optimizer-rs-score-dispatch.R`,
       `optimizer-rs-score-copula.R`, and `optimizer-rs-score-margin.R`; the
       backfitting update path is split between `optimizer-rs-backfitting.R`,
       `optimizer-rs-backfitting-proposal.R`, `optimizer-rs-acceptance.R`, and
       `optimizer-rs-acceptance-policy.R`.
- `optimizer-cg-*.R`: CG model helpers, curvature/penalty construction,
       lambda updates, trust-region candidate generation, line search,
       step acceptance, best-loglik tracking, stopping, finite-difference
       fallbacks, and finalization.

Likelihood-based calculations (likelihood, hessian, vcov):
- `likelihood-evaluation.R`, `likelihood-evaluation-margin.R`,
       `likelihood-evaluation-copula.R`, and
       `likelihood-evaluation-copula-update.R`: joint likelihood entry point,
       margin CDF/density/derivative evaluation, shared pair-level copula likelihood
       evaluation, and copula-only likelihood refreshes used by optimizer backfitting.
- `likelihood-pair-cache.R`, `likelihood-margin-cache.R`,
       `likelihood-copula-cdf.R`,
       `likelihood-copula-rectangle-bounds.R`,
       `likelihood-copula-rectangle-parameters.R`,
       `likelihood-copula-rectangle-scores.R`,
       `likelihood-copula-parameter-derivatives.R`,
       `likelihood-copula-margin-derivatives.R`,
       `likelihood-copula-margin-indexed.R`,
       `likelihood-copula-gradient-check.R`, `likelihood-numerical-cdf.R`,
       `likelihood-numerical-se-*.R`, `likelihood-gradient.R`,
       `likelihood-gradient-margin-scores.R`, and
       `likelihood-scores.R`: CDF/cache handling, discrete rectangle CDF-bound
       derivatives, copula-parameter rectangle derivatives, discrete score assembly,
       copula derivatives, first-derivative endpoint attribution, indexed
       copula-to-margin score accumulation, optional gradient diagnostics, numerical
       derivative and SE reference paths, analytical margin-score construction,
       analytical gradients, and RS score assembly.
- `hessian-*.R`: semi-analytical Hessian setup, margin derivatives, copula
       contribution assembly, covariate Hessian assembly, theta/zeta copula block
       assembly, and analytical Hessian entry points.
- `model-vcov.R`, `model-vcov-cache.R`, `model-vcov-hessian.R`,
       `model-vcov-setup.R`, `model-vcov-preflight.R`, `model-vcov-primary.R`,
       `model-vcov-solve.R`, `model-vcov-result.R`, and
       `model-vcov-smooth.R`: fitted-object covariance extraction, cache handling,
       vcov setup state, method preflight/fallback logic, Hessian solving, result
       shaping, and smooth-term covariance blocks.

##### Finalization of fit

- `model-fit-finalize.R`, `model-fit-result.R`,  `model-fit-convergence.R`, `model-fit-reporting.R`, and  `model-fit-object.R`: build the final fitted object based on all the outputs of the fit, plus kickoff calculation for variance / covariance and any other final calculations to include post-fit.



##### User-facing methods and review utilities:

- `model-accessors.R`, `model-fitted-residuals.R`, `model-summary*.R`,
  `model-check*.R`, `model-effects.R`, and `model-copula-summary.R`:
  fitted-object accessors, residuals/fitted values, summaries, model checks,
  effects, and copula time summaries.
- `model-predict*.R`, `model-confint.R`, `model-inference-*.R`, and
  `model-simulate*.R`: prediction, prediction value/evaluation helpers,
  coefficient intervals, bootstrap/likelihood inference, and fitted-model
  simulation. Bootstrap validation and result-summary helpers are isolated in
  `model-inference-bootstrap-helpers.R`. Fitted-model simulation `newdata`
  time-grid and dependence parameter alignment rules are isolated in
  `model-simulate-newdata-helpers.R`.
- `model-newdata.R`, `model-newdata-policy.R`, and
  `model-newdata-matrix.R`: newdata entry-point checks, original-to-internal
  column translation, default prediction columns, fitted-factor level
  alignment, and model-matrix column alignment.
- `model-vcov*.R`: active `vcov()` setup, method preflights,
  analytical/numerical computation, Hessian solving, and result assembly. The
  historical derivative-matrix fallback is retained in `model-vcov-legacy.R`
  outside the active documented method path.
- `diagnostics-*.R`, `missingness_diagnostics.R`,
  `missingness-diagnostics-*.R`, and `plot-*.R`: diagnostic data, scoring,
  split-family checks, missingness validation/workflow/model/assessment helpers,
  term plots, distribution plots,
  copula diagnostics, and fit overlays. Copula diagnostics are split into
  reviewer-sized controls, fitted-data, pair-data, grouping, transform, density,
  Rosenblatt metrics, tail metrics, surface metrics, overlay/cut panels,
  Rosenblatt panels, dependence/tail panels, and result helpers in
  `plot-copula-diagnostics-*.R`;
  `plot.gamlss.longitudinal()` dashboard and quantile-panel helpers are
  separated in `plot-method-helpers.R` and `plot-method-quantiles.R`;
  `plot_terms()` dashboard counting and rendering helpers are separated in
  `plot-terms-dashboard.R`; fixed-term plotting helpers are separated in
  `plot-terms-fixed-helpers.R`, `plot-terms-fixed-factor-groups.R`,
  `plot-terms-fixed-groups.R`, `plot-terms-fixed-specs.R`,
  `plot-terms-fixed-data.R`, and `plot-terms-fixed-plots.R`; smooth-term
  x-axis recovery, uncertainty, and plot-data helpers are separated in
  `plot-terms-smooth-helpers.R`; `plot_dist()` data preparation, overlay
  helpers, and panel construction are separated in `plot-dist-data.R`,
  `plot-dist-overlays.R`, and `plot-dist-panels.R`;
  `plot_margin_fit()` fitted/raw density preparation and panel construction
  are separated in `plot-margin-fit-data.R`; `plot_copula_fit()` pair/spec
  preparation, density-grid construction, and overlay assembly are separated in
  `plot-copula-fit-data.R`.
- `copula-backend-*.R`, `copula-family-links.R`, `copula_selection.R`,
  `copula-selection-*.R`, `joint_distribution_selection.R`, and
  `joint-distribution-selection-*.R`: copula backend dispatch, family metadata,
  copula selection, and joint distribution selection.
- `coverage-*.R`, `benchmark-*.R`, `adoption-*.R`, and `coverage_report.R`:
  opt-in coverage simulations, comparator metrics, adoption helpers, and
    reviewer reports. Coverage support is split into case grids/catalogs,
    fitters including longitudinal start/attempt/selection/result-row/dispatcher stages,
    GAMLSS/GAMLSS2 comparator baseline fitters, condition taxonomy, truth
    estimands, eta estimates, recovery summaries, benchmark
  setup/distributions/truth metrics/gamlss metrics, and report generation.
  Benchmark comparator validation/object assembly, coefficient extraction and
  reshaping, report inputs, section assembly, table formatting, interpretation
  object construction, interpretation display formatting, and printing helpers
  are separated in `benchmark-comparator-workflow.R`,
  `benchmark-comparator-coefficient-tables.R`,
  `benchmark-comparator-print.R`, and `benchmark-report-*.R`.
  Adoption margin screening input policy, `fitDist()` orchestration, time-intercept
  rescoring, and print/accessor helpers are separated across `adoption-margin-*.R`.
  Coverage report input/workflow preparation, LaTeX document assembly, LaTeX
  escaping, numeric formatting, table rendering, and report table families are
  isolated in `coverage-report-*.R`.
- `simulation-*.R`: current simulation helpers used by tests and JSS
  replication scaffolding.
- `legacy-simulation-create-longitudinal-dataset.R`,
  `legacy-simulation-load-dataset.R`, `legacy-simulation-se.R`,
  `legacy-optim-outer.R`, and `legacy-fit-nocov.R`: retained historical
  helpers kept out of the main reviewer path.

The machine-readable map in `inst/standards/method-traceability.csv` is the
authoritative source-to-test crosswalk. It is checked by
`tests/testthat/test-p2-method-traceability.R`, including wildcard file-family
coverage for every `R/*.R` file.

## Testing

We provide three lenses for comfirming the quality of the codebase and `gamlss.longitudinal` as a regression solution:
1. Basic build tests against all core functionality as part of `devtools::test(reporter = "summary")`
2. Extended tests against rOpenSci statistical software standards and mappings
3. Replication code for resutls presented to demonstrate package performance

### Testing: CRAN

The CRAN-facing objective is a clean source package with fast routine examples,
tests, and vignettes. Generated pkgdown, result, check, manuscript, and local
development artefacts are excluded from the source tarball.

Recommended commands:

```r
devtools::test(reporter = "summary")
rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"), error_on = "warning")
urlchecker::url_check()
spelling::spell_check_package()
```

Source package hygiene can be checked with:

```sh
R CMD build .
tar -tzf gamlss.longitudinal_*.tar.gz
```

The tarball should not include local state (`.RData`, `.Rhistory`,
`.Rprofile`), generated sites or results (`docs/`, `results/`, `examples/`),
paper sources (`paper/`), check directories, `Rplots.pdf`, or root benchmark
metric outputs. `cran-comments.md` records the latest local and remote check
results immediately before submission.

Optional comparator packages are handled as optional dependencies. `gamlss2` is
not a hard dependency; methods that require it are opt-in and should skip or
report unavailability when the package is absent.

### Testing: rOpenSci standards

The package is prepared against rOpenSci statistical software standards for
general statistical software, regression and supervised learning, and
probability distributions.

Primary evidence:

- `CONTRIBUTING.md`: terminology, input policy, missingness policy, numerical
  assumptions, and extended-test instructions.
- `R/srr-stats-standards.R`: machine-readable `srr` tags.
- `inst/standards/ropensci-srr-compliance.md`: human-readable standards
  crosswalk, including validation/test evidence and remaining TODO items.
- `tests/testthat/`: deterministic unit and integration tests for validation,
  fitting, prediction, simulation, diagnostics, and copula parity.

Routine tests are intended to be CRAN- and CI-friendly. Extended stochastic and
stress tests are opt-in:

```sh
GAMLSS_LONGITUDINAL_EXTENDED_TESTS=true Rscript -e "pkgload::load_all(); testthat::test_dir('tests/testthat')"
```

Current known standards limitations are intentionally left as `Partial/TODO` in
the crosswalk. They mostly concern broader external benchmark datasets, final
manuscript-output freezing, and platform-specific numerical tolerance evidence
after the full CI matrix has run.

### Testing: Replication of results (Paper)

The paper workflow is split into a CRAN-safe installed smoke entry point and a
repository-only full replication workflow.

Recommended smoke workflow:

```r
source(system.file("smoke-tests", "new-user-smoke.R", package = "gamlss.longitudinal"))
source("paper/replicate.R")
```

The default `paper/replicate.R` profile is smoke-sized. The expanded profile is
not CRAN-safe and may take a long time:

```r
Sys.setenv(GAMLSS_LONGITUDINAL_JSS_PROFILE = "expanded")
source("paper/replicate.R")
```

`paper/manifest.csv` is the source claim-to-artifact map. Generated tables,
figures, logs, session information, and output hashes are written under
`results/jss-replication/<profile>/`, which is excluded from source builds.

The JSS paper workflow is organised into six numbered modules covering the two
simulation studies, joint-versus-separate optimisation, missingness/dropout
sensitivity, and the LIPID/RAND applications. Current module artifacts are
clearly labelled as stubs while the paper analyses are being completed.
Private application data are never committed; secure/local paths are supplied
with `GAMLSS_LONGITUDINAL_LIPID_DATA` and `GAMLSS_LONGITUDINAL_RAND_DATA`.

## Method traceability

Reviewer entry points:

- Fitting workflow: `R/model-fit.R` is the current source of truth for
  `gamlss_longitudinal()`. Top-level fitting control, step-size control, and
  elapsed-budget checks live in `R/model-controls-fit.R`,
  `R/model-controls-step.R`, and `R/model-controls-budget.R`;
  the pre-optimizer workflow is sequenced by
  `R/model-fit-workflow.R::.gl_prepare_fit_workflow()`; optimizer setup state is assembled by
  `R/model-fit-setup.R::.gl_initialize_fit_optimizer_context()`;
  prepared workflows are bridged into the optimizer through
  `R/model-fit-optimizer.R::.gl_run_prepared_fit_optimizer()`, and the
  selected optimizer method is routed through
  `R/model-optimizer-dispatch.R::.gl_run_fit_optimizer()`, with branch-specific
  runner wiring in `R/model-optimizer-dispatch-cg.R` and
  `R/model-optimizer-dispatch-rs.R`.
- Data policy and model matrices: submitted data validation, internal
  time/subject/response naming, formula translation, and structural missingness
  expansion live in `R/model-preprocess.R::.gl_prepare_fit_data()`, with
  input-column normalization in
  `R/model-preprocess-input.R::.gl_normalize_fit_input_columns()`; copula-link
  resolution and formula matrix construction enter through
  `R/model-matrix-bundle.R::.gl_build_model_matrix_bundle()`, with design-only
  proxy-data handling in `R/model-matrix-data.R`, parameter-specific design row
  selection in `R/model-matrix-parameter-data.R`, fixed/smooth construction in
  `R/model-matrix-fixed.R` and `R/model-matrix-smooth.R`, and
  formula/column-name normalization in `R/model-matrix-formulas.R`. Input policy evidence is
  in `tests/testthat/test-p0-model-preprocess.R` and
  `tests/testthat/test-p2-srr-input-policy.R`.
- Starting values: fixed-effect starts, smooth coefficient starts, smoothing
  penalty starts, and warm-start smooth carryover enter through
  `R/model-starting-values.R::.gl_build_initial_parameter_state()`.
- Warm starts: the optional separate-RS stabilisation pass used before a joint
  RS fit is isolated in `R/model-warm-start.R::.gl_run_joint_warm_start()`.
- Fitted result assembly: convergence metadata, information criteria,
  reviewer-facing fit summaries, fit-time vcov caching, optimizer traces, and
  final `gamlss.longitudinal` object construction enter through the
  prepared-fit finalization bridge
  `R/model-fit-finalize.R::.gl_finalize_prepared_fit()` and the
  post-optimizer workflow
  `R/model-fit-result.R::.gl_finalize_fit_workflow()`, which delegates to
  `R/model-fit-convergence.R::.gl_build_convergence_info()` and
  `R/model-fit-object.R::.gl_finalize_fit_object()`;
  information criteria live in
  `R/model-fit-convergence.R::.gl_fit_information_criteria()`, and fit-time
  reporting/vcov caching lives in `R/model-fit-reporting.R`.
- Optimizer setup: RS outer/inner loop orchestration now enters through
  `R/optimizer-rs-runner.R::.gl_run_rs_optimizer()`. RS fixed/smooth design-cache construction, smooth penalty
  metadata, eta/score length validation, RS beta-start construction,
  score-path dispatch for discrete, copula, and margin parameters, theta/zeta
  pair-score aggregation, margin score assembly including the optional
  `dlcopdpar` contribution, one-block RS backfitting proposals, GAIC-based
  smoothing-parameter updates, and RS backtracking/trace acceptance enter
  through
  `R/optimizer-rs-design-cache.R::.gl_build_rs_design_cache()`,
  `R/optimizer-rs-eta.R::.gl_validate_rs_eta_lengths()`,
  `R/optimizer-rs-backfitting-inputs.R::.gl_prepare_rs_score_inputs()`,
  `R/optimizer-rs-backfitting-inputs.R::.gl_rs_beta_start()`,
  `R/optimizer-rs-backfitting-inputs.R::.gl_rs_backfitting_inputs()`,
  `R/optimizer-rs-likelihood-context.R::.gl_rs_likelihood_context()`,
  `R/optimizer-rs-likelihood-context.R::.gl_rs_discrete_scores()`,
  `R/optimizer-rs-likelihood-context.R::.gl_rs_copula_derivative_context()`,
  `R/optimizer-rs-score-dispatch.R::.gl_rs_parameter_score()`,
  `R/optimizer-rs-score-copula.R::.gl_rs_copula_parameter_score()`,
  `R/optimizer-rs-score-margin.R::.gl_rs_margin_parameter_score()`,
  `R/optimizer-rs-backfitting-iteration.R::.gl_rs_backfitting_iteration()`,
  `R/optimizer-rs-backfitting-runner-factory.R::.gl_build_rs_backfitting_runner()`,
  `R/optimizer-rs-backfitting-proposal.R::.gl_rs_prepare_backfitting_proposal()`,
  `R/optimizer-rs-lambda.R::.gl_rs_update_smoothing_parameters()`, and
  `R/optimizer-rs-acceptance-policy.R::.gl_rs_accept_backfitting_step()`;
  accepted-step state application enters through
  `R/optimizer-rs-acceptance-policy.R::.gl_apply_rs_acceptance_state()`;
  RS acceptance and backtracking diagnostics enter through
  `R/optimizer-rs-acceptance-policy.R::.gl_report_rs_acceptance()`;
  RS optimisation progress plotting is isolated in
  `R/optimizer-rs-diagnostics.R::.gl_plot_rs_progress()`;
  shared optimizer counters, histories, diagnostics, and likelihood caches enter
  through `R/optimizer-state.R::.gl_initialize_optimizer_state()`; optimizer
  likelihood and parameter history appends enter through
  `R/optimizer-state.R::.gl_append_optimizer_history()`; automatic and
  explicit stopping criteria enter through
  `R/optimizer-stop-criteria.R::.gl_resolve_stop_criteria()`, and negative
  outer-log-likelihood streak handling enters through
  `R/optimizer-rs-loop-stopping.R::.gl_update_outer_negative_streak()`. Shared
  outer-iteration console reporting enters through
  `R/optimizer-rs-reporting.R::.gl_print_outer_iteration_summary()` and
  `R/optimizer-rs-reporting.R::.gl_print_outer_convergence()`; RS outer-iteration
  step-size and likelihood-change state enter through
  `R/optimizer-rs-loop-stopping.R::.gl_update_rs_outer_iteration_state()`. CG outer-loop orchestration now enters through
  `R/optimizer-cg-runner.R::.gl_run_cg_optimizer()`. CG structural helpers for augmented design matrices, beta unpacking,
  smooth penalties, penalized objectives, lambda candidates/scoring, step
  limiting and trust-radius updates (`optimizer-cg-trust-region.R`),
  line-search evaluation (`optimizer-cg-line-search.R`), step acceptance
  (`optimizer-cg-acceptance.R`), best-loglik tracking (`optimizer-cg-best.R`),
  stopping checks, stop-reason selection, and first-lambda-update convergence
  delay live in `R/optimizer-cg-*.R`; the CG
  analytical Hessian availability guard lives in
  `R/optimizer-cg-hessian-availability.R::.gl_ensure_cg_hessian_available()`;
  finite-difference gradient/Hessian checks, observed-Hessian
  analytical/fallback orchestration, smooth EDF calculations, and step/lambda
  trace row construction are isolated in the corresponding CG concern files.
- Likelihood and uncertainty: likelihood evaluation, copula derivatives, margin
  score contributions, Hessian assembly, and vcov paths live in
  `R/likelihood-*.R`, `R/hessian-*.R`, and `R/model-vcov*.R`; parity and
  numerical tests are in the `p0` Hessian/copula/rectangle tests. The exported
  `vcov()` method starts by assembling object fields, optional parameter
  overrides, and eta values in
  `R/model-vcov-setup.R::.gl_prepare_vcov_evaluation()`.
- User-facing review methods: prediction, inference, simulation, diagnostics,
  and checks live in the corresponding `R/model-*.R`, `R/diagnostics-*.R`, and
  `R/plot-*.R` files and are covered by `p1` workflow tests. Prediction
  threshold/density evaluation values are isolated in
  `R/model-predict-values.R::.gl_prediction_eval_values()`, while newdata
  translation, default columns, fitted-factor level alignment, and design
  column alignment are isolated in `R/model-newdata-policy.R` and
  `R/model-newdata-matrix.R`. Fitted-model simulation row/time-grid and
  dependence-parameter alignment policies are isolated in
  `R/model-simulate-newdata-helpers.R`.

The fitted object should be inspected through accessors first: `coef()`,
`vcov()`, `summary()`, `predict()`, `simulate()`, `model.frame()`, `fitted()`,
`residuals()`, `formula()`, `terms()`, `nobs()`, `logLik()`, and
`check_model()`. Internal list fields remain available for method review but
are not the preferred user interface.

The machine-readable traceability table in
`inst/standards/method-traceability.csv` maps major method areas to source
files, tests, paper modules, and refactor status. It is checked by
`tests/testthat/test-p2-method-traceability.R`, which verifies that referenced
source files, exact R entry point functions, wildcard file families, and test
files exist. The same test also verifies that every `R/*.R` source file is
covered by the traceability map, so newly added source files must be classified
for reviewers.

## Testing output (most recent)

The latest local architecture and submission-readiness checkpoint was run on
2026-06-14 on Windows 11 x64 with R 4.4.1 after the reviewer-oriented file
restructuring:

- `devtools::document()`: pass. The local roxygen2 version is 7.3.2 while the
  package was previously documented with 7.3.3; this is reported as a version
  note, not a documentation failure.
- NAMESPACE comparison against the pre-documentation snapshot: unchanged.
- Duplicate top-level R function definition scan: pass.
- `R/common_functions.R` removal check: pass; the file is absent.
- Legacy reachability audit: `loadDataset()`, `create_longitudinal_dataset()`,
  `bvt_norm_true_SE_B0_Bt()`, `optim_outer()`, and `fit_jointreg_nocov()` are
  not exported and are only referenced inside the legacy helper group. They are
  safe to remove from package internals once historical scripts have been
  checked, but are retained for now.
- `devtools::test(reporter = "summary")`: pass, with 2 expected opt-in
  extended-test skips and known coverage convergence warnings.
- Extended `testthat::test_dir()` with
  `GAMLSS_LONGITUDINAL_EXTENDED_TESTS=true`: pass, with 4 CRAN skips and 1
  expected max-iteration warning.
- `R.exe CMD build .`: pass with vignettes built.
- Source-tarball spot-check for generated/local exclusions: pass.
- `R CMD check --as-cran --no-manual gamlss.longitudinal_0.1.0.tar.gz`: 0
  ERRORs, 0 WARNINGs, 3 NOTEs. Remaining NOTEs are new submission plus optional
  `gamlss2` availability via `Additional_repositories`, local time
  verification, and unavailable local pandoc for checking `README.md`/`NEWS.md`.
- `urlchecker::url_check()`: pass
- `spelling::spell_check_package()`: pass.
- New-user smoke test: pass.
- Default JSS smoke replication, `source("paper/replicate.R")`: pass; existing
  smoke targets reported up to date.
