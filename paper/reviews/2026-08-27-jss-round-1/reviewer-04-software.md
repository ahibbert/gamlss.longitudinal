# Independent JSS software/API review

## Recommendation

**Reject in the present form; encourage resubmission after a major software and reproducibility revision.** The package has a substantial statistical implementation, a broad test suite, and unusually good internal traceability. However, the frozen manuscript is visibly an internal draft rather than a reviewable JSS article, its central reproducibility claim is false for the private-data application, and several interfaces advertised as standard R methods do not satisfy their generic contracts. I reviewed the frozen manuscript at source SHA `68c3bad26626ce7c267bd330364cfda8df7a6b76`; the current package checkout is a later commit (`7343d06c986b81e2dbbbf5c2bc2732c3fe539c41`). I did not inspect prior JSS audits or review reports.

Read-only verification included the frozen `main.tex` and PDF build log; package source, namespace, Rd pages, vignettes, tests, README and `REVIEW.md`; the shipped smoke test; and targeted API probes. The new-user smoke test completed. A fitted smoke object nevertheless showed `class(logLik(fit)) == "numeric"`, length three, with no `df` or `nobs`; `AIC(fit)` returned `numeric(0)`; and `vcov(fit)` was a list rather than a matrix.

## Five strengths

1. **Serious statistical scope in one coherent fit.** The public fitter supports parameter-specific formulas for up to four marginal and two copula parameters (`R/model-fit.R:1-24, 152-205`), including smooths and covariate-dependent dependence.
2. **Clear and explicit copula scope.** Six implemented families are validated centrally—Gaussian, Clayton, Frank, Gumbel, Joe, and Student-t (`R/copula-backend-family-map.R:4-15`)—and the paper clearly states the first-order dependence limitation (`main.tex:534-551, 1624-1628`).
3. **Rich post-fit functionality.** Prediction covers means, medians, parameters, quantiles, CDFs, densities, and threshold probabilities (`R/model-predict.R:9-21, 97-164`); `simulate()` preserves the modeled serial copula on supplied panels (`R/model-simulate.R:1-23, 78-109`).
4. **Good input and fitted-object traceability.** Duplicate subject/time rows, unsupported columns, missing predictors, and missing pair support receive explicit checks (`R/model-preprocess.R:54-99`; `R/model-column-policy.R:22-124`), while fitted objects retain submitted and expanded data, formulas, matrices, likelihoods, traces, and convergence metadata (`R/model-fit-object-fields.R:28-74`).
5. **Unusually strong internal review infrastructure.** `REVIEW.md:31-227` maps fitting, likelihood/Hessian, prediction/simulation, diagnostics, and replication to source and tests; the repository includes fixed seeds, tolerance policies, a manifest, a new-user smoke script, and extensive `testthat` coverage.

## BLOCKING issues

### B1. The frozen manuscript is an annotated internal draft, not a submission artifact

- **Exact location:** frozen `main.tex:35-36, 99-112, 143, 146-155, 158-159, 277, 281, 327-348, 580, 745-755, 774, 853, 1031, 1186, 1222, 1610, 1616, 1637`; frozen `main.log:883-901, 944, 998, 1121, 1212, 1231`.
- **Evidence:** `todonotes` is enabled; the source contains 28 live `\todo`, `\jssrevisionnote`, or red-placeholder occurrences, including one in the abstract, a large author planning box after the abstract, “approximately X hours,” and an instruction to report the computational environment. The PDF build log records undefined citations (`czado2015`, `Beareseo2015`, `lambert_copula`, `topmodels`, `simes`, `marra2025?`) and reports a 47-page output.
- **Consequence:** Reviewers cannot distinguish scientific content from author notes, key claims lack resolvable support, and the artifact is not in submission-ready JSS form.
- **Remedy:** Remove all author/reviewer annotations and placeholders; resolve every citation/reference; reduce or move derivations to supplementary material; rebuild in the intended JSS template; visually inspect every page.
- **Completion test:** a source scan returns zero live TODO/revision/placeholder commands; the final LaTeX log contains zero undefined citations/references and zero overfull boxes affecting tables; rendered-page review finds no annotations, clipped tables, or unreadably small figures.
- **Confidence:** High.
- **Classification:** Defect.

### B2. “All results” are not reproducible, and the frozen paper does not identify a reproducible software environment

- **Exact location:** frozen `main.tex:1292-1593, 1595-1616`; `paper/README.md:3-5, 20-22, 32-40, 85-102`; `paper/manifest.csv:22-29, 37-40`; `.Rbuildignore:25`; `DESCRIPTION:1-35`.
- **Evidence:** The paper says, “All results included in this paper can be directly replicated” (`main.tex:1597`) and presents the private LIPID analysis as its clinical example. Repository documentation says LIPID is private and never part of the public graph, and the manifest marks its figures/tables as private static or provenance-only artifacts with no generated path. The paper leaves OS, R/package versions, cores, seed policy, and runtime as a red placeholder. It gives no installation command, package version, license, repository release/tag, or package commit. The full `paper/` workflow is excluded from the source package, while the installed `inst/jss-replication/run-replication.R` is a different lightweight coverage workflow.
- **Consequence:** An independent reviewer cannot reproduce the paper’s principal application or reconstruct the precise package/environment used; the absolute reproducibility claim is materially misleading.
- **Remedy:** Make the main worked application public/package-shipped or fully simulated. Move LIPID to a clearly labeled non-reproducible secondary illustration unless controlled reviewer access is provided. Freeze a release archive/DOI and lockfile, state commit/version/license/install commands and hardware/runtime, distinguish repository-only from installed-package workflows, and narrow all reproducibility claims to what the public graph actually regenerates.
- **Completion test:** on a clean machine, a reviewer can obtain the named release, restore its environment, run one documented command, regenerate every main-text numerical table/figure, and pass manifest/tolerance checks without private files; all non-reproducible secondary assets are explicitly labeled and excluded from “all results” wording.
- **Confidence:** High.
- **Classification:** Defect.

## MAJOR issues

### M1. Advertised “standard methods” violate the `logLik()` and `vcov()` contracts

- **Exact location:** frozen `main.tex:928-953`; `R/model-accessors.R:1-10`; `R/model-fit-object-fields.R:28-60`; `R/model-vcov.R:28-37, 171-180`; `R/model-vcov-result.R:4-24`; `man/logLik.gamlss.longitudinal.Rd:14-16`; `man/vcov.gamlss.longitudinal.Rd:59-61`; `vignettes/jss-start-here.Rmd:89-104`.
- **Evidence:** `logLik.gamlss.longitudinal()` returns all three likelihood components through partial `$calc_lik_out` matching against the stored `calc_lik_out_end`, as a plain numeric vector. It lacks class `logLik`, `df`, and `nobs`. In a smoke fit, `AIC(fit)` therefore returned `numeric(0)`. `vcov.gamlss.longitudinal()` returns a nested list (`vcov`, `se`, method metadata) rather than the covariance matrix expected from `stats::vcov()`. Tests exercise internal structures but do not assert these generic contracts.
- **Consequence:** Base/model-ecosystem tools that consume `logLik()` or `vcov()` can silently fail or misbehave; the manuscript’s compatibility claim is not met.
- **Remedy:** Make `logLik()` return the scalar joint likelihood with class and attributes; expose components through a separate accessor. Make `vcov()` return the fixed-effect matrix by default; expose smooth matrices/SEs/diagnostics through a separate detailed method or an explicit nondefault argument. Add `AIC()`, `BIC()`, and downstream interoperability tests.
- **Completion test:** `inherits(logLik(fit), "logLik")`; `length(logLik(fit)) == 1`; finite `attr(., "df")` and `attr(., "nobs")`; `AIC(fit)` and `BIC(fit)` are finite scalars; `is.matrix(vcov(fit))`; matrix names align exactly with `coef(fit)`; `broom`, `modelsummary`, or equivalent consumers pass smoke tests.
- **Confidence:** High.
- **Classification:** Defect.

### M2. The marginal-family contract is undefined and can silently use an invalid continuous likelihood for unsupported discrete families

- **Exact location:** frozen `main.tex:143, 159, 291-295, 793-795`; `R/model-fit.R:14-19, 35-47`; `man/gamlss_longitudinal.Rd:68-71, 235-252`; `R/package-utils.R:107-134, 172-227`; `R/model-fit-entrypoint.R:70-82`.
- **Evidence:** Documentation implies a generic `gamlss.dist` family interface (“GA(), NO(), PO(), NBI(), etc.”) but provides no authoritative support matrix. `.is_discrete_margin()` explicitly excludes `BB`, `DBI`, `ZABB`, `ZABI`, `ZIBB`, and `ZIBI`; no pre-fit rejection follows, so these discrete families can enter the continuous density-times-copula-density likelihood. The normalizer also silently drops family parameters lacking link metadata and fixes them at inferred defaults.
- **Consequence:** Users can obtain numerically plausible but statistically invalid fits, and “any supported margin” is circular because “supported” is never enumerated by fitting/optimizer/inference capability.
- **Remedy:** Define and enforce a versioned family capability registry covering likelihood type, fitted parameters, simulation, prediction, analytic/numerical vcov, diagnostics, and tested copulas. Reject unsupported families and dropped parameters with an actionable error unless the user explicitly opts into a documented fixed-parameter contract.
- **Completion test:** every installed `gamlss.dist` family is classified supported/partial/unsupported; every unsupported bounded discrete family fails before optimization; representative continuous, count, zero-inflated, and binary families pass likelihood/prediction/simulation/vcov tests; the manuscript and Rd page reproduce the same table.
- **Confidence:** High.
- **Classification:** Defect.

### M3. `check_model()` falsely diagnoses correctly specified discrete margins

- **Exact location:** frozen `main.tex:981`; `R/model-check.R:61-69`; `R/diagnostics-pit.R:1-24`; `R/model-check-thresholds.R:3-9, 28-45`; `R/diagnostics-plot-pithist.R:1-3`; `R/diagnostics-plot-qqrplot.R:1-5`.
- **Evidence:** `check_model()` calls `.gl_pit(object, randomize = FALSE)` and then applies a continuous-uniform KS test. For discrete responses, raw upper-endpoint CDF values are atomized and not Uniform(0,1), even under the true model. A direct known-truth Poisson probe with `n = 500`, `mu = 3` gave KS `p = 1.23e-23` for the non-randomized PIT versus `p = 0.88` for randomized PIT. Individual PIT/QQ/worm methods also default to `randomize = FALSE`, although the QQ axis says “Randomized Quantiles.”
- **Consequence:** The advertised automated diagnostic will systematically warn against valid Poisson/NB/Delaporte/binary fits, undermining a central practitioner workflow.
- **Remedy:** For discrete margins use randomized quantile residuals with explicit seed/RNG handling, or a discrete calibration method; document stochasticity and avoid presenting a single KS cutoff as a definitive pass/fail screen.
- **Completion test:** repeated known-truth discrete simulations produce approximately calibrated diagnostic rejection rates; a fixed seed makes results reproducible; tests cover at least Poisson, negative binomial/Delaporte, and binary cases; labels and defaults agree.
- **Confidence:** High.
- **Classification:** Defect.

### M4. Non-converged CG fits can be returned without a warning

- **Exact location:** frozen `main.tex:853, 981, 1777-1785`; `R/optimizer-cg-stop-request.R:70-96, 152-162`; `R/model-fit-convergence.R:20-47`; `R/model-fit-object.R:36-41`; `R/optimizer-cg-reporting.R:4-27`.
- **Evidence:** CG loop control sets its internal `converged = TRUE` to exit for tolerance, raw-likelihood deterioration, or max stall, while final metadata correctly calls only `stop_reason == "tolerance"` converged. Finalization warns only when `hit_outer_limit` is true. A max-stall or deterioration stop before the iteration limit therefore returns a statistically non-converged object without a warning; max-stall has no verbose stop message either.
- **Consequence:** Users may proceed to inference, diagnostics, prediction, or publication without noticing that optimization did not converge.
- **Remedy:** Separate `stop_requested` from statistical convergence internally and warn for every returned fit with `convergence$converged != TRUE`, including reason and recommended action. Inference methods should warn or require explicit override on non-converged objects.
- **Completion test:** fixtures for max stall, likelihood deterioration, iteration limit, and non-finite criteria each emit a classed warning and retain distinct stop reasons; tolerance convergence is silent; `summary()`, `vcov()`, and `confint()` visibly propagate status.
- **Confidence:** High.
- **Classification:** Defect.

### M5. The paper’s guided example and the producer of its figures are different analyses

- **Exact location:** frozen `main.tex:772-812, 818-849, 859-875, 959-979, 2071-2142`; `paper/R/public-paper-producers.R:21-45`; `paper/manifest.csv:2-4`.
- **Evidence:** The text presents one `n = 500`, four-time-point, covariate-rich BCPE/t example and prints its summary. The registered producer instead simulates `n = 35/80`, three time points, constant marginal/copula parameters; uses only two times for `plot_dist()`; fits intercept-only formulas for four iterations with a 30-second limit; suppresses warnings; and produces the displayed diagnostics figure. The hard-coded console output is not a manifest artifact.
- **Consequence:** Readers cannot execute the displayed code and obtain the displayed figures/output, and a possibly non-converged toy fit is presented as diagnostics from the richer example.
- **Remedy:** Create one public example object and one producer that emits the dataset, fit, printed summary, plots, predictions, simulations, and timings used throughout the section. Do not suppress convergence warnings; fail production unless convergence is established.
- **Completion test:** every code/output/figure block in the software section maps to the same manifest target and seed; generated summary text matches the frozen listing; the producer asserts convergence; the observed dimensions/subjects/times match captions and prose.
- **Confidence:** High.
- **Classification:** Defect.

### M6. User documentation is broad but not executable enough for a JSS software artifact

- **Exact location:** all 53 `man/*.Rd` files; `gamlss.longitudinal.Rcheck/00check.log` line “checking examples ... NONE”; `vignettes/jss-start-here.Rmd:10-15, 73-110`; `vignettes/standard-workflow.Rmd:10-15`; `vignettes/diagnostics-decisions.Rmd:10-15`; `vignettes/inference-uncertainty.Rmd:10-15`; `vignettes/simulator-usage.Rmd:10-15`; `vignettes/site-guide.Rmd:10-15`.
- **Evidence:** No Rd page contains an `\examples{}` section. Five main guidance vignettes and the JSS entry vignette globally set `eval = FALSE`; the JSS fit template uses undefined `dat`/`newdat`. Only the working-install and long detailed vignette are executable. The repository smoke script passes, but it does not exercise `logLik`, the standard `vcov` contract, family screening, discrete diagnostics, or the paper example.
- **Consequence:** Documentation can drift without checks, novice users cannot copy a self-contained minimal example from help pages, and JSS’s executable software-story expectation is not met.
- **Remedy:** Add fast examples to core Rd pages and make the JSS/minimal workflow a self-contained evaluated vignette using package-shipped data or a small simulation. Keep only demonstrably expensive chunks unevaluated, with cached output provenance.
- **Completion test:** `R CMD check` reports examples rather than `NONE`; evaluated minimal/JSS vignettes build in a clean library; every public function shown in the main software story is executed in CI; all referenced objects are defined in-document.
- **Confidence:** High.
- **Classification:** Defect.

### M7. No reviewed package release is tied to the frozen manuscript

- **Exact location:** frozen `main.tex` (no installation/version/license statement; explicit environment placeholder at `1616`); `DESCRIPTION:4, 13-35`; `README.md:5-15, 59`; root `gamlss.longitudinal_0.1.0.tar.gz` metadata (`Packaged: 2026-08-03`); `.Rbuildignore:25, 31`.
- **Evidence:** The checkout declares version 0.1.0 and “early CRAN/JSS hardening,” while the only source tarball predates the August 27 manuscript freeze and omits later material such as `vignettes/jss-start-here.Rmd`. The paper gives neither a tag/DOI nor a GitHub install-at-ref command and does not state GPL-3.
- **Consequence:** Reviewer results depend on an unfrozen working tree, and users cannot retrieve the exact package described by the article.
- **Remedy:** Cut an immutable release from the reviewed commit, archive it, record SHA256/DOI, and cite exact installation and license information in the paper. Generate the submission tarball and replication bundle from that tag.
- **Completion test:** clean-room installation of the archived tarball succeeds; its `DESCRIPTION`, docs, vignettes, tests, and replication manifest match the manuscript; the paper’s version/commit/hash resolve to that artifact.
- **Confidence:** High.
- **Classification:** Defect.

### M8. Computational scaling guidance is insufficient for model and replication decisions

- **Exact location:** frozen `main.tex:561-577, 922, 1008, 1111, 1145-1168, 1239, 1305, 1610, 1626`; `R/model-fit.R:85-139`; `vignettes/native-simulation-workflow.Rmd:406`.
- **Evidence:** The paper says the method is approximately 20 times slower than `gamlss2`, calls runtime the main bottleneck, and uses a 10% clinical subsample, but gives no absolute timings, hardware, memory, scaling curve, or complexity in subjects, time points, coefficients, smooth basis size, family, optimizer, and vcov method. The full-replication runtime remains “X hours.”
- **Consequence:** Practitioners cannot judge feasibility or choose RS/CG, analytical/numerical covariance, screening scope, bootstrap size, or sample size.
- **Remedy:** Add a reproducible benchmark grid and practical table/plot with absolute elapsed time and peak memory, convergence rate, hardware/software, and dominant complexity. Separate fitting from vcov, diagnostics, bootstrap, and full-paper runtime.
- **Completion test:** a public benchmark script regenerates the table/plot; at least `n`, `T`, fixed/smooth dimension, optimizer, family type, and vcov method vary; the paper gives concrete recommendations and a measured end-to-end replication time.
- **Confidence:** High.
- **Classification:** Question.

## MINOR issues

### m1. The manuscript reverses the KS diagnostic inequality

- **Exact location:** frozen `main.tex:981`; `R/model-check-thresholds.R:5, 28-45`.
- **Evidence:** The paper says failure is flagged when the PIT KS “p-value > 0.05”; code flags `p < 0.05`.
- **Consequence:** Readers receive the opposite interpretation of the documented diagnostic.
- **Remedy:** Correct the inequality and describe it as a heuristic screen with multiplicity/sample-size caveats.
- **Completion test:** prose, code, Rd, vignette, and tests all state the same rule.
- **Confidence:** High.
- **Classification:** Defect.

### m2. `response`, `mu`, fitted values, and response means have confusing defaults

- **Exact location:** `R/model-predict.R:9-21, 35-45, 79-93, 97-108`; `R/model-predict-values.R:34-80`; `R/model-fitted-residuals.R:17-45`; `man/predict.gamlss.longitudinal.Rd:28-34, 80-82`; `vignettes/jss-start-here.Rmd:106-110`.
- **Evidence:** `predict()` defaults to `type = "response"`, but that is an alias for the distribution’s `mu`, explicitly not necessarily the response mean. `fitted()` also returns `mu`; response and Pearson residuals subtract or scale by `mu`. For most non-special-cased families, `type = "mean"` is approximated by averaging 199 quantiles.
- **Consequence:** Standard-looking calls can target the wrong estimand, and approximation error is undocumented/untested across supported families.
- **Remedy:** Default to an unambiguous estimand, deprecate the misleading alias, name fitted/residual types explicitly, and document/test mean computation accuracy or mark unavailable families.
- **Completion test:** family-specific tests distinguish `mu` from mean; help and printed output label the estimand; approximation error is bounded against analytic or high-accuracy numerical references.
- **Confidence:** High.
- **Classification:** Suggestion.

### m3. The fitted object is inspectable but not yet a fully conventional model object

- **Exact location:** frozen `main.tex:292-295, 853`; `R/model-fit-object-fields.R:28-74`; `NAMESPACE:7-53`.
- **Evidence:** The object stores useful internals, but no original `call`, complete normalized control bundle, package/schema version, `update()` method, or validation method is stored/registered. Its class is assigned as a single string rather than a documented inheritance chain.
- **Consequence:** Exact refitting, object migration, and reproducible modification are harder than in standard R model workflows, weakening the “single auditable model object” story.
- **Remedy:** Store `call`, normalized controls, data-contract metadata, package/object version, and add `update()` plus an object validator; document which slots are stable public API.
- **Completion test:** `update(fit, ...)` reproduces/refits as expected; saved objects validate after reload; a documented accessor supplies controls without users reaching into internal slots.
- **Confidence:** Medium-high.
- **Classification:** Suggestion.

### m4. Practitioner code should not require an internal triple-colon helper

- **Exact location:** frozen `main.tex:2073-2076, 2130-2139`; `R/copula-family-links.R:5-56`; `NAMESPACE:55-106`.
- **Evidence:** The appendix uses `gamlss.longitudinal:::get_copula_dist("t")` to construct simulation parameters; the helper is deliberately internal and absent from the public namespace.
- **Consequence:** Published code depends on an unstable implementation detail and is not a durable public workflow.
- **Remedy:** Export a small supported copula parameter/link conversion API or express simulation inputs through documented Kendall-tau/df arguments.
- **Completion test:** all manuscript/vignette code runs without `:::`; the public helper has Rd documentation and tests for every supported copula.
- **Confidence:** High.
- **Classification:** Suggestion.

## Proposed software-story sequence for the paper

1. **User problem and scope:** what longitudinal question this solves; first-order, common-family, marginal-interpretation assumptions; when GEE/GLMM/vines remain preferable.
2. **Installable artifact and support contract:** exact version/commit/license/install command, input data contract, supported margin/copula capability matrix, missingness and irregular-time semantics.
3. **One public end-to-end example:** generate or load one package-shipped dataset; explore; screen margin/copula; fit one model; print convergence and timing. Use this same object throughout.
4. **Model object and standard methods:** show `summary`, `coef`, scalar `logLik`, matrix `vcov`, `confint`, `AIC/BIC`, `formula`, `model.frame`, `tidy/glance/augment`, and explicit convergence metadata.
5. **Decisions after fitting:** marginal and dependence diagnostics, with valid discrete/continuous branches and an explanation of what constitutes review versus failure.
6. **What users can predict and simulate:** clearly distinguish `mu`, response mean, quantiles/probabilities, unconditional panel simulation, and the absence of response-conditioned forecasting.
7. **Validation and performance:** concise recovery/calibration evidence tied to advertised families; absolute runtime/memory scaling; RS/CG/vcov/robust-SE decision guidance.
8. **Limitations and alternatives:** first-order dependence, unsupported families, missing-data assumptions, non-convergence, bootstrap cost, private-data constraints.
9. **Exact public reproduction:** one clean-room command, release/lockfile/hardware/runtime, manifest/tolerances; keep private LIPID only as a clearly secondary, non-reproducible illustration.

## Top 10 ordered actions

1. Freeze and archive one package/replication release, then state exact version, commit, hash, license, installation, and environment in the paper.
2. Replace the private LIPID analysis as the main workflow with a public/package-shipped example, or provide accepted reviewer access and narrow reproducibility claims.
3. Remove every draft annotation/placeholder and resolve all citations/references; rebuild and visually QA the JSS-formatted PDF.
4. Repair `logLik()`/`vcov()` contracts and add finite `AIC()`/`BIC()` plus ecosystem interoperability tests.
5. Build and enforce a public family capability registry; hard-error unsupported discrete/bounded families before optimization.
6. Correct discrete diagnostics by using reproducible randomized/discrete calibration methods and fix the manuscript’s KS inequality.
7. Warn on every non-converged return and block or loudly qualify downstream inference.
8. Rebuild the software section from one manifest-tracked, convergence-asserted example that generates all printed output and figures.
9. Add fast Rd examples and evaluated clean-library vignettes for the core fitting, standard methods, diagnostics, prediction, and simulation workflow.
10. Publish reproducible absolute runtime/memory benchmarks and measured smoke/paper/full replication runtimes with practical optimizer/vcov guidance.
