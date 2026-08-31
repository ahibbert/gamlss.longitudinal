# Independent applied-statistics review for *Journal of Statistical Software*

**Manuscript:** *gamlss.longitudinal: Longitudinal GAMLSS Models with Copula Dependence in R*  
**Frozen source reviewed:** `main.pdf` and `main.tex`, source SHA `68c3bad26626ce7c267bd330364cfda8df7a6b76`  
**Review perspective:** whether an applied longitudinal analyst can identify the estimand, determine whether the data and assumptions are compatible with the method, fit and diagnose a model, interpret marginal and copula parameters, make predictions, and report defensible uncertainty. The package source was inspected only where needed to resolve manuscript ambiguity. Existing audits and prior review reports were not used.

## Recommendation

**Reject in the present form, with encouragement to resubmit after major restructuring and statistical clarification.** In editorial terms, this is more than a normal revision: the PDF is visibly an internal annotated draft, and several applied claims cannot presently be interpreted safely. The software addresses a real gap and the public simulated workflow is promising, so a substantially rebuilt paper could become a strong JSS contribution.

Three matters block a scientifically reviewable submission: (1) unresolved draft annotations, placeholders, duplicated text, and missing citations remain in the rendered article; (2) the paper mislabels or under-defines key parameter scales, including an incorrect link for the Student's t copula and an invalid interpretation of its degrees-of-freedom coefficient as a direct test of tail dependence; and (3) the advertised treatment of intermittent missing visits is not the observed-data likelihood of the stated first-order Markov/copula model, because observations separated by a missing intermediate visit are treated as independent rather than integrating over the latent intermediate response.

The LIPID analysis is potentially useful as a secondary illustration, but it is currently too long, too definitive, and insufficiently auditable for a private-data example. Its 10% subsample, post-selection inference, substantial differential-over-time missingness, baseline adjustment, and causal language require resolution. Unless the full-data and missingness analyses can be strengthened, the example should be shortened and framed as a descriptive software illustration, not a scientific re-analysis of treatment effects.

## Five strengths

1. **A useful, identifiable software niche.** The package offers one formula per marginal or dependence parameter while retaining a marginal distributional interpretation. This fills a practical space between mean-focused GEE/GLMM workflows and highly parameterized general vine-copula workflows.

2. **A credible end-to-end public workflow.** Section 3 (pp. 11–17) covers simulation, distribution screening, fitting, standard methods, diagnostics, marginal prediction, and dependent trajectory simulation in one example. That is the right backbone for a JSS article.

3. **The central structural limitation is not hidden.** Sections 2.1–2.3 and 4.3.1 (pp. 6–8 and 21–23) explicitly describe the first-order truncation and compare AR(1)-like and exchangeable dependence. The exchangeable counterexample is especially valuable.

4. **The paper recognizes that both margins and dependence need checking.** The combined PIT/quantile diagnostics, Rosenblatt residual views, fitted-versus-observed rank dependence, tail co-occurrence, simulation, sandwich covariance, and bootstrap options are richer than coefficient-only reporting.

5. **Performance tradeoffs are reported rather than concealed.** Tables 5–8 (pp. 19–21) show recovery, predictive scores, substantial runtime costs, optimizer guidance, and preliminary use/caution advice. This is useful raw material for a more applied decision guide.

## BLOCKING issues

### B1. The rendered manuscript is an internal annotated draft, not a reviewable article

- **Location:** p. 1 abstract and Introduction; pp. 2–3 Section 1.1; pp. 5–7 Section 2.1; p. 9 Section 2.4; p. 11 Section 3; p. 13 Section 3.2.2; p. 17 Section 4; pp. 20–21 Section 4.3; p. 30 Section 6; pp. 31–33 Discussion and Appendix A.
- **Evidence:** The PDF renders red “Next steps” blocks, orange margin comments, “Draft wording” boxes, a duplicated multi-paragraph discussion on p. 5, unresolved `?` citations, “X hours,” “Report operating system, R version, package versions…,” and author questions such as whether optimizer detail should be included. Page 1 alone contains submission-planning notes and a comment inside the abstract. The reference list lacks normal bibliographic fields for several entries.
- **Consequence:** Readers cannot distinguish author claims from internal deliberation, and multiple statements are explicitly provisional. Page count, structure, citation support, and reproducibility claims cannot be evaluated as a finished submission.
- **Concrete remedy:** Produce a clean manuscript build with all revision/todo macros disabled or removed; resolve every placeholder and citation; remove duplicated material; complete the computational-environment and runtime statement; and have a non-author read the rendered PDF page by page.
- **Completion test:** A search of both TeX and extracted PDF text returns zero matches for `Next steps`, `Draft wording`, `todo`, `X hours`, `give some refs`, `Report operating system`, and unresolved `?` citations. A visual page-by-page check shows no colored review annotations and no duplicated paragraphs.
- **Confidence:** High.
- **Classification:** **Defect**.

### B2. The paper does not provide a correct, family-specific estimand map; the LIPID copula link and tail interpretation are wrong

- **Location:** Tables 2–3 and Equation 2, pp. 7–8, Sections 2.2–2.3; Section 3.2.2, pp. 12–13; Table 15, p. 28; LIPID interpretation, pp. 27–30, Section 5.
- **Evidence:** The prose repeatedly calls `mu`, `sigma`, `nu`, and `tau` “mean, variance, skewness and kurtosis,” even though their meanings are family-specific. For BCPE, `mu` is the median, `sigma` is approximately a coefficient of variation only in a restricted regime, and `nu`/`tau` control rather than equal skewness/kurtosis. Table 15 says the Student's t `theta` link is logit; the package implements Fisher-z, `atanh(rho)`, with inverse `tanh(eta)`. Table 15 labels `zeta` “tail correlation,” but it is `log(df - 2)`; tail-dependence strength depends jointly on `rho` and degrees of freedom. Testing a `zeta` coefficient against zero tests `df = 3`, not “no tail dependence” or a Gaussian copula. The manuscript nevertheless says the t parameter is “highly significant” and treats this as evidence for tail correlation. Coefficients in Table 15 are presented only on link scales, with no transformed contrasts in median, scale, Kendall's tau, degrees of freedom, or tail-dependence probability.
- **Consequence:** An applied analyst can report the wrong estimand, direction, or scientific magnitude. In the LIPID example, the treatment coefficient `0.202` for `theta` is not itself a correlation increase, and the `zeta` p-value does not establish tail dependence. This invalidates central scientific interpretations of the example.
- **Concrete remedy:** Add a family-specific “parameter-to-estimand” table for every worked family and copula. For each parameter give its domain, link, natural-scale interpretation, meaningful null, and recommended transformed summary. In examples, report covariate contrasts in response summaries (mean/median/quantiles/scale), `rho`, Kendall's tau, t degrees of freedom, and lower/upper tail-dependence coefficients with uncertainty. Correct Table 15's link. Replace the t-tail significance claim with a comparison of t versus Gaussian fits and uncertainty for scientifically interpretable tail probabilities or tail-dependence coefficients.
- **Completion test:** Every coefficient discussed in Sections 3 and 5 can be traced through a stated link to a named natural-scale estimand. Values reconstructed from the coefficient table agree numerically with package prediction/summary output. No test of `zeta = 0` is described as a test of tail dependence or Gaussianity.
- **Confidence:** High.
- **Classification:** **Defect**.

### B3. Intermittent missing visits are not handled by the likelihood claimed in the paper

- **Location:** Section 3.1, p. 11; Section 3.2.3, p. 15; Table 8, p. 21; Section 4.3.3, pp. 24–25; Section 7, p. 31.
- **Evidence:** The paper says the software fits data with missing responses by using available margins and adjacent observed pairs. It also says that if an intermediate observation is missing, the observations on either side “will not exhibit correlation.” The implementation expands the global subject-by-time grid and includes only complete consecutive time pairs. Under the stated first-order copula/Markov joint model, however, the observed-data likelihood for `Y(t-1)` and `Y(t+1)` when `Y(t)` is missing requires integrating over `Y(t)`; it generally induces a non-independent two-step transition. Dropping both adjacent copula contributions is not that likelihood. Monotone dropout after the last observed response is a different and more defensible case under an ignorable missingness assumption, but the paper pools it with intermittent MAR missingness and calls results “reassuring.”
- **Consequence:** Coefficients, dependence, Hessian-based standard errors, likelihood comparisons, AIC, and simulation claims can be biased or miscalibrated when there are interior gaps. Users are told the package “supports” missing data without knowing that the fitted objective changes the dependence model around gaps.
- **Concrete remedy:** Choose and state one of two scopes. Either (a) implement the proper observed-data likelihood for missing intermediate states, using analytical/numerical integration or an equivalent transition construction, and validate it; or (b) restrict likelihood claims and recommended use to complete or monotone-dropout panels, explicitly label intermittent-gap fitting as a working/composite approximation, disable likelihood-based comparisons that are not valid under that approximation, and give a principled sensitivity procedure. Multiple imputation should not be suggested generically unless the imputation model preserves the multivariate distribution, time ordering, treatment, and dependence estimands and inference is pooled appropriately.
- **Completion test:** A three-visit test with visit 2 intermittently missing compares the package objective to a direct numerical integral of the stated complete-data joint density. They agree within a pre-specified tolerance, or the software and paper issue a prominent restriction/warning and do not call the approximation the joint observed-data likelihood. Separate simulation results are shown for monotone dropout and interior gaps under MCAR, covariate-dependent MAR, and at least one MNAR sensitivity.
- **Confidence:** High.
- **Classification:** **Defect**.

## MAJOR issues

### M1. The applied data contract and pair-level covariate semantics are incomplete

- **Location:** Section 2.2, pp. 6–8; Section 3.1, p. 11; fitting call in Section 3.2.2, p. 12.
- **Evidence:** The article says only that data are longitudinal and that subject/time columns are supplied. It does not consolidate the actual requirements: long format; one unique row per subject-time; non-missing subject and time; observed predictors on submitted rows; a common ordered set of time levels; a single marginal family across visits; and identifiable complete adjacent pairs. It does not explain what “adjacent” means for irregularly spaced visits. The implementation pairs consecutive *global ordered time levels*, regardless of elapsed gap, and evaluates copula formulas using the left member of each pair. Thus a time-varying covariate in `theta.formula` has a specific lagged/left-endpoint meaning that is absent from the paper.
- **Consequence:** Analysts may unknowingly treat visits at 1 month and 24 months like visits one month apart, may misinterpret a time-varying dependence coefficient, or may supply asynchronous subject-specific visit times that create a huge and inappropriate global grid.
- **Concrete remedy:** Add a one-page data contract with required columns/types, uniqueness and missingness rules, balanced/unbalanced panel behavior, global versus subject-specific visit grids, minimum pair support, and exact pair-formula evaluation. Add examples for factor visit, numeric elapsed time, irregular schedules, and a time-varying predictor. State whether elapsed gap must enter `theta.formula`, whether only scheduled panels are supported, and how to construct pair-level covariates when the desired estimand uses both endpoints.
- **Completion test:** For each of four toy datasets—balanced scheduled, monotone dropout, irregular common visits, and asynchronous subject-specific visits—the paper states whether the data are supported, shows the pairs used, and identifies the row from which each copula covariate is taken.
- **Confidence:** High.
- **Classification:** **Question** (requiring explicit scope and documentation).

### M2. The estimand is called “marginal,” but the population and conditioning set are not defined

- **Location:** Abstract and Sections 1.1–1.2, pp. 1–4; Section 2.2, pp. 6–8; Section 3.2.4, p. 15; Section 7, p. 31.
- **Evidence:** “Marginal interpretation” is used mainly to contrast the model with random-effects models. The paper does not define whether coefficients describe the distribution at a scheduled visit conditional on observed covariates, averaged over a covariate distribution, averaged over missingness/visit processes, or standardized to a target population. It also shifts between “mean coefficients,” generic `mu`, and distributional parameters.
- **Consequence:** Users may equate a conditional-on-covariates marginal response distribution with a population-standardized treatment effect or population-average mean effect. This is particularly consequential in LIPID, where baseline covariates, time, dropout, and a 10% subsample all change the population being described.
- **Concrete remedy:** Define the model estimand mathematically as `F_t(y | X_t)` and the pair dependence as a property of adjacent conditional probability transforms. Then distinguish coefficient-scale, covariate-conditional response summaries, and standardized population contrasts. Provide one worked `marginal_effects()`/prediction example that averages over a declared target covariate distribution.
- **Completion test:** The abstract, methods, prediction section, and LIPID section each use the same vocabulary for coefficient, conditional response, and standardized estimands, and every treatment contrast names its target population and time.
- **Confidence:** High.
- **Classification:** **Question**.

### M3. Prediction is marginal and unconditional, but the article invites trajectory/dynamic interpretations

- **Location:** Section 3.2.5, pp. 16–17; Section 1, p. 2; Discussion, p. 31.
- **Evidence:** The paper says the full predictive distribution is available and promotes trajectory simulation. The package contract clarifies an important limitation not stated with equal precision in the article: `predict()` returns marginal row-wise summaries and does **not** condition a future response on that subject's observed history; `simulate()` generates unconditional trajectories on a supplied panel. The shown quantiles are fitted conditional-distribution quantiles, not prediction intervals for a future subject with parameter uncertainty, and the example uses the fitting data as `newdata`.
- **Consequence:** A clinical analyst may use the method for individual dynamic prediction or interpret fitted quantiles as subject-specific forecasts even though the observed history is ignored. Unconditional simulation and conditional forecasting answer different questions.
- **Concrete remedy:** Create a prediction taxonomy: fitted marginal summaries, new-row marginal prediction, standardized population contrast, unconditional new-subject trajectory simulation, and conditional/dynamic prediction given history (currently unsupported). Distinguish confidence intervals for fitted summaries, predictive intervals including outcome variability, and simulation uncertainty. Add a held-out-subject example and an explicit “do not use for dynamic prediction” statement unless conditional prediction is implemented.
- **Completion test:** A code example predicts genuinely held-out subjects and labels every output by estimand and uncertainty source. A test with two subjects sharing covariates but different response histories demonstrates that current `predict()` gives the same future marginal prediction, and the manuscript explains why.
- **Confidence:** High.
- **Classification:** **Suggestion**.

### M4. Inferential guidance does not yet support defensible applied reporting

- **Location:** Section 2.4, pp. 9–11; model summary, pp. 12–13; Sections 3.2.4–3.2.5, pp. 15–17; Tables 5–6, p. 19; Table 15, p. 28; Appendix B, pp. 34–36.
- **Evidence:** Standard errors are defined using `sqrt(abs(diag(-H^{-1})))`. Taking an absolute value can conceal negative variance estimates from an indefinite Hessian rather than declaring inference invalid. The example reports a Hessian condition number of 8309 but does not tell users what is acceptable or show eigenvalue/positive-definiteness diagnostics. Smooth-term variances are computed separately from the fixed-effect covariance and need clearer inferential status. The manuscript demonstrates only 20 parametric bootstrap replicates, which is not adequate for reported percentile intervals, and it gives no rule for failed refits. It does not distinguish model-based, cluster-sandwich, parametric-bootstrap, and cluster-bootstrap targets and assumptions. LIPID intervals are reported after extensive model selection without acknowledging selection uncertainty.
- **Consequence:** Apparently precise p-values and intervals may be reported from unstable curvature, an invalid smooth approximation, too few bootstrap replicates, or a selected model. The existence of several covariance options does not tell an analyst which is primary.
- **Concrete remedy:** Do not mask an indefinite covariance with absolute values; fail or label inferential output unavailable when curvature diagnostics do not pass. Define diagnostics and tolerances. Provide a decision table for model-based, sandwich, cluster bootstrap, and parametric bootstrap inference, including cluster-count conditions, failed-refit rules, and recommended replicate counts. Validate fixed and smooth intervals, transformed dependence summaries, and misspecified-dependence cases. Treat post-selection LIPID intervals as exploratory or use a pre-specified model/validation strategy.
- **Completion test:** Deliberately indefinite/ill-conditioned fits produce no ordinary Wald p-values and a clear diagnostic. Simulation tables report coverage for all advertised interval methods, smooths, and transformed copula summaries with Monte Carlo standard errors. The paper's reporting example uses at least a defensible number of bootstrap replicates or clearly labels the small run as syntax-only.
- **Confidence:** High.
- **Classification:** **Defect**.

### M5. Automated diagnostic “passes” use uncalibrated tests and omit the most important higher-lag check

- **Location:** Section 3.2.5, p. 15; Figure 3, p. 16; Table 8 and Section 4.3.1, pp. 21–23; Figures 15–16, pp. 44–45.
- **Evidence:** `check_model()` declares a marginal failure when a conventional KS p-value for PIT uniformity is below 0.05, but PITs are estimated from the same data and are dependent within subjects. For discrete margins, the implementation uses non-randomized PITs, which are not Uniform(0,1) even under a correct model. The default residual-dependence flag checks only lag 1 at an ad hoc absolute correlation threshold of 0.25. Yet the distinctive model risk is residual dependence at lag 2 and beyond after fitting adjacent copulas. LIPID diagnostics are declared “reasonable” without quantitative criteria, uncertainty bands, or subgroup/sample-size information.
- **Consequence:** Correct discrete models can fail, misspecified higher-order models can pass, and large samples can reject negligible deviations. A green automated result may give false reassurance precisely where the first-order restriction is wrong.
- **Concrete remedy:** Make PIT tests descriptive or calibrate them by subject-level/model-based simulation; use randomized PIT or a discrete PIT method for discrete margins. Show residual dependence by every estimable lag and elapsed-time gap, with simulation envelopes or cluster-aware uncertainty. Label thresholds as review triggers, not pass/fail validity tests. Report quantitative LIPID diagnostic values and pair counts by visit/arm.
- **Completion test:** Type-I flag rates are documented under calibrated continuous and discrete simulations. The exchangeable scenario reliably triggers the higher-lag diagnostic, and an AR(1)-compatible model generally does not. `check_model()` cannot return an unqualified “passed” result when higher-lag diagnostics were unavailable or not assessed.
- **Confidence:** High.
- **Classification:** **Defect**.

### M6. Model-building guidance conflates exploratory screening, joint model selection, and confirmatory inference

- **Location:** Section 3.2.3, pp. 13–14; Section 3.2.4, p. 15; Section 4.2, p. 20; LIPID Tables 12–15, pp. 26–28.
- **Evidence:** Margins and copulas are screened by AIC on the same observations later used for fitting and inference. The text notes that choices may change after adding covariates but gives no stable workflow for revisiting them. The LIPID build begins with “all” covariates for successive parameters, removes non-significant terms, then reports ordinary p-values and 95% intervals from the final selected model. Table 14 calls comparisons sequential likelihood-ratio comparisons even when specifications and effective degrees of freedom change, and the final model is described as reducing non-significant covariates. Hierarchy for treatment-by-time interactions, clinical pre-specification, multiplicity, collinearity, and selection stability are not addressed.
- **Consequence:** The example teaches a high-dimensional stepwise procedure whose standard p-values are anti-conservative and whose scientific narrative can be selected from noise. Copula-family selection can also absorb marginal misspecification.
- **Concrete remedy:** Separate an exploratory workflow from a reporting workflow. Pre-specify scientifically important mean/treatment terms; define candidate marginal/copula families before inspecting results; preserve interaction hierarchy; use joint AIC/cross-validation or an explicit validation set for distributional choices; and treat secondary parameter discovery as exploratory with multiplicity/stability assessment. If stepwise selection remains, state its limitations and do not use unadjusted final-model p-values as confirmatory evidence.
- **Completion test:** A flowchart determines when to refit the margin after dependence terms, when to compare joint families, and which data are used for selection versus assessment. The LIPID model can be reconstructed from a pre-declared candidate set and selection rule, and repeated resampling shows term/family stability.
- **Confidence:** High.
- **Classification:** **Suggestion**.

### M7. The LIPID treatment-effect claims do not match the analysis population or missingness assumptions

- **Location:** Section 5, pp. 25–30; Table 15, p. 28; Figures 8–9, pp. 29–30.
- **Evidence:** The analysis uses a 10% random subsample (902 of 9014 participants) solely “for computational convenience.” Fourteen percent of cholesterol values are missing overall, increasing from 4% to 36% by the last visit. No dropout summary by treatment arm, baseline cholesterol, outcome history, or selected covariates is shown; no inverse-probability, multiple-imputation, joint-dropout, or MNAR sensitivity analysis is performed. Nevertheless the text uses causal language (“treatment has a significant effect,” “effect degrades,” “treatment patients” have increased stability) and implies an intention-to-treat result. Random treatment assignment does not by itself remove bias from post-randomization outcome missingness or outcome-dependent availability. The prose also calls the treatment-by-time interaction “negative” although the reported interaction coefficients are positive; the intended statement appears to be that an initially negative treatment contrast attenuates.
- **Consequence:** The estimand is at best a model-based contrast in the observed outcomes of a random 10% subset, conditional on a selected covariate model. It is not yet a defensible full-trial causal or intention-to-treat treatment effect. Dependence differences are especially sensitive to dropout and changing pair composition.
- **Concrete remedy:** Define the trial estimand and analysis set; report the exact sampling seed/algorithm and arm/baseline balance; analyze the full dataset if computationally feasible or demonstrate repeated-subsample stability; show visit- and arm-specific denominators and adjacent-pair counts; model/weight/impute dropout under a stated MAR strategy and add an MNAR sensitivity analysis. Use “associated with” for non-randomized baseline covariates. Reserve causal treatment language for a properly defined treatment-policy/hypothetical estimand with missing-data assumptions and sensitivity results.
- **Completion test:** The full-data or repeated-subsample treatment contrasts are stable on natural scales; dropout diagnostics and at least two defensible missingness analyses are reported by arm; the treatment-by-time sign is stated correctly; and every causal sentence names the estimand, treatment assignment, and missingness assumption.
- **Confidence:** High.
- **Classification:** **Defect**.

### M8. The LIPID model specification and scientific interpretation are not auditable enough for a private-data example

- **Location:** Section 5, pp. 25–30; Table 14, p. 27; Table 15, p. 28; Figures 6–9, pp. 26–30; reproducibility statement, p. 30; Appendix I, pp. 44–47.
- **Evidence:** The public sanitized LIPID contract specifies only `treatment + time + treatment:time` for `mu`, `treatment + time` for `sigma`, intercepts for `nu/tau/zeta`, and `treatment + time` for `theta`. It does not reproduce the actual Table 14/15 model, which includes a baseline-response smooth and many clinical covariates. Neither the 10% sampling nor data derivation, visit collapsing, covariate coding, missingness handling, fitted coefficients, or selection sequence can be rerun from the public artifact. The manuscript's stated formula includes baseline response while baseline cholesterol is also described as one of the seven longitudinal outcomes; it is unclear whether the baseline row is simultaneously outcome and predictor, which could create deterministic self-adjustment/leakage at baseline. Private figures occupy roughly nine manuscript/appendix pages, yet the example's unique software lesson is not clearly separated from its unverifiable scientific conclusions.
- **Consequence:** Reviewers and readers cannot verify that the static tables/figures arise from the stated model, and the sanitized code describes a materially different analysis. This weakens rather than strengthens the JSS reproducibility story.
- **Concrete remedy:** Either obtain a redistributable/public longitudinal dataset and make it the applied example, or sharply shorten LIPID to a clearly marked secondary illustration. Supply a sanitized but exact data dictionary, derivation and sampling pseudocode, exact formulas/coding/contrasts, model-selection log, software/session metadata, pair counts, and immutable hashes for approved outputs. Explain precisely how baseline cholesterol enters the outcome and covariate definitions. The public simulated example should reproduce the same reporting objects and applied decisions without private data.
- **Completion test:** Starting from an access-controlled LIPID extract, an independent authorized analyst can run the sanitized script unchanged (apart from the input path) and reproduce every published number/hash. The script's formulas match Tables 14–15 exactly. A public reviewer can reproduce an isomorphic analysis and all method demonstrations without LIPID. If that cannot be achieved, LIPID is reduced to a short, explicitly non-reproducible illustration without novel substantive claims.
- **Confidence:** High.
- **Classification:** **Defect**.

### M9. The validation results need denominators, Monte Carlo uncertainty, and a clearer link to applied claims

- **Location:** Section 4.1, pp. 17–19, Tables 4–6; Section 4.2 and Appendix F, pp. 20 and 39–41, Tables 18–21; Section 4.3, pp. 21–25.
- **Evidence:** Table 4 reports 1000 replicates, while Appendix F discusses 100 replications and some prose is ambiguous about retained fits. Main tables give means/SDs but not attempted, converged, retained, or failed fit counts by method. Monte Carlo standard errors are absent. Table 5 excludes seven outer-fence values from one variogram score and puts the raw result in a footnote, creating an asymmetric post hoc summary. Most claims are based on a narrow set of high-signal or correctly specified scenarios; the impact of weak dependence, few subjects, many parameters relative to subjects, and informative dropout is not synthesized into use guidance.
- **Consequence:** Readers cannot distinguish algorithmic failure from statistical performance, assess simulation precision, or know whether favorable averages reflect conditional-on-success results. The guidance is less transportable to small or messy applied studies.
- **Concrete remedy:** Pre-specify the simulation estimands, metrics, failure policy, and robust summaries. Report attempted/converged/retained counts and reasons for failure, Monte Carlo SEs, parameter-to-subject ratios, and hardware/runtime. Keep all finite predictive scores or use the same pre-specified robust summary for every method. Summarize results by applied decision boundaries rather than many appendix tables.
- **Completion test:** Every simulation table has explicit denominators and Monte Carlo uncertainty; no result is silently conditioned on convergence or outlier exclusion; and each use/do-not-use recommendation cites a scenario that tests the relevant failure mode.
- **Confidence:** High.
- **Classification:** **Suggestion**.

### M10. “Use / do not use” guidance needs to cover the full supported domain, not only three sensitivity topics

- **Location:** Tables 7–8, pp. 20–21; Sections 3.1 and 7, pp. 11 and 31.
- **Evidence:** Existing advice covers optimizer choice, exchangeable dependence, copula uncertainty, and missingness. It does not give a consolidated position on irregular/asynchronous visits, interior gaps, very small numbers of subjects, many visit levels or smooth terms, sparse discrete outcomes, zero inflation, negative dependence (not supported by several one-sided copula families), multivariate outcomes, dynamic prediction, time-varying pair covariates, or causal treatment-effect estimation.
- **Consequence:** Analysts may select the method because it is more flexible than GEE/GLMM without recognizing when a simpler model is safer or when a fuller vine/joint missingness/dynamic prediction approach is required.
- **Concrete remedy:** End the methods section with a compact decision table: data pattern/problem, recommended method, `gamlss.longitudinal` status (recommended/caution/not supported), diagnostic, and alternative. Include positive examples and firm exclusions.
- **Completion test:** An analyst can answer “Should I use this package?” for each domain above without consulting source code or inferring from simulations.
- **Confidence:** High.
- **Classification:** **Suggestion**.

## MINOR issues

### N1. Mathematical notation and terminology contain errors that impede implementation

- **Location:** Sections 2.2–2.4, pp. 6–10, Equations 1–4 and the discrete rectangle discussion.
- **Evidence:** Time points are introduced as `t = 1, ..., n` rather than `T` on p. 7; a CDF is called a “cumulative density”; `Y ~ D(mu,sigma,nu,tau)` suppresses subject/time conditioning; pair-specific and shared parameters are interchanged; the discrete rectangle formula switches from `a_it,b_it` to undefined `a1,b1,a2,b2`; the displayed PIT approximation is not readable as an applied recipe.
- **Consequence:** Readers cannot reliably map the mathematical objective to rows, pairs, formulas, and discrete likelihood options.
- **Concrete remedy:** Rebuild notation around subject `i`, scheduled visit `j`, elapsed time, row covariates, and pair covariates. Define continuous and discrete contributions separately and keep every index consistent.
- **Completion test:** A reader can derive the likelihood for one subject with three visits and identify exactly which design row enters each parameter without consulting code.
- **Confidence:** High.
- **Classification:** **Defect**.

### N2. Several figures and tables are too small or insufficiently annotated for applied use

- **Location:** Table 1, p. 2; Figures 3–5, pp. 16, 18, and 24; Tables 10–11, pp. 22–23; Tables 12–15, pp. 26–28; Figures 15–18, pp. 44–47.
- **Evidence:** Dense tables use very small type; Figure 4's labels are difficult to read; the multi-panel diagnostics on pp. 44–47 are full-page but have tiny axes and no in-panel interpretation; Table 14 is hard to parse as a model sequence. Figure 8's caption refers to “true distribution histograms,” although these are empirical observed histograms.
- **Consequence:** The very plots intended to teach diagnosis and interpretation are not usable at normal reading size.
- **Concrete remedy:** Retain only decisive panels in the main text, enlarge labels, add pair/subject counts, and use captions that say what pattern is acceptable or concerning. Move detailed panels to supplementary material with readable dimensions.
- **Completion test:** All text is legible at 100% on an A4/US Letter page and every caption identifies the data, estimand, sample/pair count, and intended diagnostic lesson.
- **Confidence:** High.
- **Classification:** **Defect**.

### N3. Claims about comparator methods are sometimes categorical or inaccurate

- **Location:** Sections 1.1–1.2 and 2.1, pp. 2–6; Table 1, p. 2.
- **Evidence:** The manuscript describes GLMM and GEE as “mean models” unable to handle other distributional features, presents vine workflows as generally requiring manual simulation for adjusted inference, and says “all alternative popular models” are jointly optimized. These statements ignore broader distributional mixed models, marginal models beyond the named packages, robust/two-stage inference procedures, and differences between likelihood and estimating-equation targets. Some claims have unresolved citations.
- **Consequence:** The contribution appears overstated and readers may choose software based on an unfair comparison rather than estimand and data structure.
- **Concrete remedy:** Restrict comparisons to named implementations and default workflows; compare estimands, likelihood status, dependence flexibility, missingness, prediction, and inference without universal claims.
- **Completion test:** Every comparator statement is scoped to a package/version or supported by a citation and could be accepted by a knowledgeable user of that method.
- **Confidence:** Medium-high.
- **Classification:** **Defect**.

### N4. Reporting and copy-editing need a final statistical pass

- **Location:** Throughout; especially pp. 1–7, 15, 25–31, and References pp. 42–47.
- **Evidence:** Examples include `pacakage`, `paramterized`, `oragnized`, `treament`, `dypnoeoa`, “fixed or smooth” as if mutually exclusive, p-values printed as `0.000`, inconsistent Student's t typography, unresolved references, and incomplete bibliographic metadata. Section 6 promises exact reproduction while exempting the private example without a precise static-artifact policy.
- **Consequence:** Terminological slips obscure already difficult distinctions such as linear versus smooth effects, parameter versus response scale, and reproducible versus precomputed outputs.
- **Concrete remedy:** Perform a copy edit after statistical restructuring, standardize notation and p-value reporting, and align reproducibility wording with the actual public/private asset policy.
- **Completion test:** Automated spelling/reference checks pass and a statistical editor verifies every use of mean, median, variance, scale, correlation, Kendall's tau, tail dependence, conditional, marginal, prediction, and causal effect.
- **Confidence:** High.
- **Classification:** **Defect**.

## Ideal applied-reader journey

An effective JSS paper should let the applied reader move through the following sequence without detouring into optimizer derivations:

1. **Start with the problem and the target.** Show one longitudinal dataset where the conditional median/quantiles, scale, and adjacent dependence vary with covariates. Define the response-distribution and dependence estimands on natural scales.
2. **Decide whether the method applies.** Present the data contract and a firm use/caution/do-not-use table, including visit regularity, intermittent gaps, higher-order dependence, discrete responses, and prediction goals.
3. **Understand the model visually.** Use a small diagram: each observation has a GAMLSS margin; each scheduled adjacent pair has one copula whose parameter formula is evaluated at a stated endpoint/pair row; higher-tree copulas are independence copulas.
4. **Prepare and inspect data.** Show long-format validation, subject/visit counts, pair counts, dropout/interior gaps, marginal distributions, dependence by lag/gap, and missingness patterns.
5. **Choose a small candidate set.** Explain how scientific knowledge defines formulas and how marginal/copula screening creates a shortlist rather than a final inferential model.
6. **Fit an interpretable model.** Start simple, preserve hierarchy, show convergence/curvature output, then use the joint optimizer for the reporting fit.
7. **Translate coefficients.** Move immediately from link-scale coefficients to median/mean/quantile/scale contrasts, Kendall's tau, degrees of freedom, and tail probabilities with uncertainty.
8. **Diagnose the exact assumptions.** Check marginal calibration, discrete calibration if relevant, dependence by every lag and elapsed gap, tail behavior, influential subjects, and inference sensitivity.
9. **Predict only what is supported.** Demonstrate held-out marginal prediction and unconditional trajectory simulation; explicitly distinguish these from unavailable dynamic prediction conditional on response history.
10. **Report and reproduce.** Produce a model-specification table, natural-scale results, diagnostic decision record, session/runtime information, and a fully public replication. Only then add a short private LIPID illustration, if retained, to show external realism rather than to carry new scientific claims.

## Top 10 ordered actions

1. **Create a clean, submission-ready build** with all annotations, placeholders, duplicate text, and unresolved citations removed.
2. **Resolve the intermittent-missing-visit likelihood problem** by implementing proper observed-data integration or sharply restricting/documenting the supported missingness scope.
3. **Write the family-specific estimand/link dictionary** and correct the Student's t `theta`/`zeta` interpretations everywhere, especially Table 15.
4. **Define the applied data and pair contract** for scheduled/irregular time, pair construction, pair covariates, predictor missingness, and minimum support.
5. **Rebuild inference safeguards and guidance,** including non-positive-definite Hessians, smooth uncertainty, sandwich/bootstrap choices, failed refits, and transformed-parameter intervals.
6. **Replace automated diagnostic pass/fail claims** with calibrated continuous/discrete diagnostics and explicit lag-2-plus and elapsed-gap checks.
7. **Rewrite the prediction section** around marginal prediction versus unconditional simulation versus unsupported history-conditional prediction, using held-out subjects.
8. **Redesign model building** to distinguish exploratory screening from pre-specified/reporting models and address post-selection uncertainty.
9. **Either rehabilitate or sharply shorten LIPID:** use the full trial or demonstrate subsample stability, analyze dropout by arm, define the causal/descriptive estimand, make exact private-code provenance auditable, and remove unsupported scientific claims.
10. **Condense the article around the applied workflow,** moving optimizer derivations and large diagnostic/simulation panels to supplementary material and ending with a comprehensive use/do-not-use decision table.
