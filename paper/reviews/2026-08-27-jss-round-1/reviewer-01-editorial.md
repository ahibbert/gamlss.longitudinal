# Developmental editorial review

Reviewed the frozen manuscript at source SHA `68c3bad26626ce7c267bd330364cfda8df7a6b76`: `C:/Users/Aydin/AppData/Local/Temp/gamlss-jss-review-round-20260827/paper-repo/main.pdf`.

## Overall developmental recommendation

**Major revision / resubmit for review; do not send to external review in the present form.**

The underlying contribution is potentially well suited to the Journal of Statistical Software: a formula-based R implementation joining flexible GAMLSS margins to covariate-dependent adjacent-time copula dependence, with fitting, inference, diagnostics, prediction, and simulation in one object. I would not recommend a definitive editorial decline on contribution alone.

However, the current manuscript is simultaneously a methods paper, an extensive simulation paper, a software manual, and a clinical analysis. The software contribution is therefore harder to identify and evaluate than it should be. More seriously, the novelty landscape is incomplete, several comparison claims are too categorical, the reproducibility section is unfinished, and the manuscript does not yet establish formal package availability/version/licensing. These are editorial gating issues, not polishing matters.

JSS explicitly requires coverage of implementations with similar scope, encourages empirical comparisons, discourages extensive simulations, and expects exact replication materials that run in roughly one hour or a shorter reviewer profile when the full run is longer ([JSS author requirements](https://www.jstatsoft.org/authors)).

## Five strongest aspects

1. **A plausible and useful software niche.** Shared formulas for flexible marginal parameters and adjacent-time copula parameters offer a more interpretable middle ground between mean-focused longitudinal models and fully parameterized vines.

2. **A genuinely end-to-end fitted object.** The package appears to expose standard methods, diagnostics, prediction, simulation, bootstrap and sandwich inference, and distribution-screening helpers rather than merely wrapping an optimizer.

3. **Unusually candid scope boundaries.** The paper acknowledges first-order dependence, runtime, missingness, and model-misspecification limitations instead of presenting the method as universally preferable.

4. **Attention to user decisions.** Tables 7 and 8 translate optimizer and sensitivity studies into operational guidance, which is more valuable to JSS readers than raw performance tables alone.

5. **Substantial validation effort.** Continuous and discrete margins, fixed and smooth effects, dependence misspecification, missingness, and runtime are all considered. The work needs sharper prioritization, but the validation base is promising.

# BLOCKING issues

## B1. The novelty claim is not supported by a complete current software landscape

- **Type:** Confirmed defect, with a reviewer question about the surviving uniqueness claim.
- **Location:** pp. 2–6, Table 1, §§1.1–1.2, especially the bold claim on p. 3 that “There is currently no direct practitioner-facing software…”.
- **Evidence:** The comparison centers on `geepack`, `lme4`, `gamlss`/`gamlss2`, `GJRM`, and vine packages. It omits the closest recent multivariate distributional-regression work: Kock and Klein’s Gaussian-copula structured additive distributional regression, which allows multivariate GAMLSS-style margins and covariate-dependent dependence ([published description](https://www.tandfonline.com/doi/abs/10.1080/10618600.2024.2434181)). The role of `bamlss`, `gcmr`, multivariate `gamboostLSS`, and direct `copula::fitMvdc()` constructions is also not resolved. The omission is confirmed even if the authors ultimately justify that none offers the same longitudinal shared-formula workflow.
- **Consequence:** A knowledgeable editor or referee may regard the central uniqueness statement as unresearched or overstated. This is the clearest path to editorial rejection.
- **Concrete remedy:** Conduct a dated, explicit landscape review. Compare at least model dimension, supported margins, dependence families, covariate-dependent dependence, shared versus outcome-specific formulas, estimation paradigm, discrete handling, standard methods, diagnostics, and simulation. Rewrite the contribution as a bounded statement, e.g., uniqueness in a **frequentist, shared-formula, first-order longitudinal, joint-likelihood workflow**, if that remains defensible.
- **Completion test:** Every method or package plausibly matching two or more central capabilities is either in the comparison table or explicitly excluded with a documented reason. The abstract and §1.2 contain no unqualified “no software exists” claim.
- **Confidence:** High.

## B2. The paper’s architecture obscures the JSS contribution and includes an extensive simulation paper inside it

- **Type:** Confirmed defect.
- **Location:** Whole manuscript; especially pp. 4–11 (§2), pp. 17–25 (§4), pp. 32–41 (Appendices A–H). The manuscript is 31 main-text pages and 47 rendered pages overall.
- **Evidence:** Seven main-text pages derive the method and optimizer, nine more present simulation and sensitivity studies, and ten appendix pages continue derivations and simulation tables. JSS guidance specifically says extensive simulation studies are discouraged. The software interface begins only on p. 11.
- **Consequence:** The paper’s identity becomes “new statistical method plus several validation studies,” while the package design and practical differentiation receive less narrative weight. Review burden is high and the principal contribution is difficult to summarize.
- **Concrete remedy:** Rebuild the paper around one reproducible software workflow. Retain only the minimal model definition, the likelihood, discrete-data treatment, and user-relevant optimizer distinctions in the main paper. Move derivative/Hessian derivations and most simulation grids to supplementary material or a separate methods paper. Reduce the main manuscript to approximately 24–28 readable pages.
- **Completion test:** A reader reaches executable package use by page 6–7; at least half the main paper concerns interface, workflow, outputs, comparison, and limitations; the full simulation appendix is external to the main article.
- **Confidence:** High.

## B3. Reproducibility is incomplete and internally contradictory

- **Type:** Confirmed defect.
- **Location:** pp. 30–31, §6; also pp. 25–30, §5.
- **Evidence:** The paper says “All results…can be directly replicated,” but the clinical data cannot be redistributed. The runtime is still “X hours,” followed by an author note that it should be near one hour, and p. 31 contains a placeholder instruction to report OS, R version, package versions, cores, seeds, and runtime.
- **Consequence:** This fails a formal JSS submission requirement. The manuscript cannot claim exact reproduction of all displayed results when a major example is unavailable without a documented access path or a redistributable substitute.
- **Concrete remedy:** Provide:
  1. one standalone, commented reviewer script;
  2. a rendered output file including `sessionInfo()`;
  3. a ≤1-hour reviewer profile reproducing all main-paper figures/tables;
  4. a full profile with measured runtime and hardware;
  5. deterministic seeds and platform notes;
  6. either an accessible route to the clinical data, redistributable sufficient derived data, or removal of nonreproducible clinical outputs from the reproducibility claim.
- **Completion test:** A clean-session run regenerates every main-paper artifact and verifies it against a manifest; the manuscript states measured runtimes; no placeholder remains; every exception is explicitly identified.
- **Confidence:** High.

## B4. The submission does not establish the software artifact’s publication state

- **Type:** Confirmed manuscript omission; reviewer question about the package’s actual status.
- **Location:** p. 11, §3.1; pp. 30–31, §6.
- **Evidence:** The text links to a pkgdown site and mentions package source, but gives no release version, release date, persistent repository/DOI, CRAN or r-universe status, license, minimum R version, supported platform information, or citation metadata. The public site currently instructs GitHub installation rather than installation from a standard repository. JSS expects packaged, documented software and a GPL-compatible license ([JSS software preparation](https://www.jstatsoft.org/authors)).
- **Consequence:** An editor cannot determine what immutable software version the article documents or whether the submission meets formal publication requirements.
- **Concrete remedy:** Add a concise “Software availability” subsection/table specifying version/commit, license, repository and archived release DOI, installation command, R/dependency requirements, tested platforms, documentation/support locations, and package citation.
- **Completion test:** A reader can install the exact reviewed release from a persistent source, and the manuscript/package metadata agree on version, license, URL, and citation.
- **Confidence:** High.

# MAJOR issues

## M1. The intended audience and conceptual on-ramp are misaligned

- **Type:** Confirmed developmental defect.
- **Location:** pp. 1–7, Introduction and §§2.1–2.2.
- **Evidence:** The paper promises accessibility to practitioners familiar with `lme4` or `geepack`, but rapidly introduces GAMLSS, PCCs, D-vines, copula tail dependence, marginal versus conditional interpretation, and six distributional predictors. Several concepts are repeated on p. 5 rather than explained progressively.
- **Consequence:** The purported primary audience may not understand the estimand or know when the method is preferable before reaching the fitting interface.
- **Concrete remedy:** Add a one-page conceptual example: “what changes beyond the mean,” “what adjacent dependence means,” and “what a coefficient on θ answers.” State prerequisites and separate the applied-reader route from the methodological detail.
- **Completion test:** A GEE/GLMM user can explain the model’s marginal estimand, dependence assumption, and main limitation after reading the first three pages.
- **Confidence:** High.

## M2. The primary example is too maximalist to teach the software

- **Type:** Confirmed developmental defect.
- **Location:** pp. 12–17, §3.2.
- **Evidence:** The first full fit uses a four-parameter BCPE margin, two-parameter Student copula, six formulas, 24 fixed coefficients, three smooths, a Hessian condition number, screening, diagnostics, prediction, simulation, bootstrap, and sandwich inference. The synthetic covariates have no substantive interpretation.
- **Consequence:** Readers cannot distinguish essential syntax from advanced capability; the example demonstrates breadth but not learnability.
- **Concrete remedy:** Use a staged, reproducible example:
  1. simple margin plus Gaussian copula and two formulas;
  2. add scale and covariate-dependent dependence;
  3. optionally extend to BCPE/t and smooths.
  Prefer a packaged public dataset or a substantively named synthetic longitudinal trial dataset.
- **Completion test:** The first fit is under 12 lines, its coefficients have an immediate interpretation, and each later feature is introduced because a diagnostic or scientific question motivates it.
- **Confidence:** High.

## M3. The software section inventories functions but does not adequately explain design contracts

- **Type:** Confirmed defect.
- **Location:** pp. 11–17, §§3.1–3.2.5.
- **Evidence:** The paper lists many helpers and generic methods, but does not give a compact supported-feature matrix, object anatomy, error/warning behavior, compatibility rules between margins and copulas, or computational scaling.
- **Consequence:** Readers cannot judge package maturity or anticipate whether their data/model is supported.
- **Concrete remedy:** Add:
  - a supported-data/family/copula table;
  - a small fitted-object diagram;
  - convergence and failure semantics;
  - memory/runtime scaling in (n), (T), parameter count, and smooth basis size;
  - a clear distinction between stable public API and internal implementation.
- **Completion test:** A reader can determine support and likely cost without consulting source code or multiple external vignettes.
- **Confidence:** High.

## M4. The empirical comparison does not include the closest software alternatives

- **Type:** Confirmed defect.
- **Location:** Table 1; pp. 17–24, §§4.1–4.3.2, Tables 4–11.
- **Evidence:** The main recovery comparison is against `gamlss2`, which deliberately ignores dependence; the longitudinal benchmark is against GEE under exponential-family margins. Neither comparison tests the principal claimed niche against `GJRM` in (T=2), a staged GAMLSS-plus-vine/conditional-copula workflow, or a multivariate distributional-regression implementation. Official `GJRM` documentation confirms flexible covariate effects and multiple margins for bivariate and some trivariate models ([CRAN documentation](https://search.r-project.org/CRAN/refmans/GJRM/help/gjrm.html)).
- **Consequence:** The manuscript demonstrates that dependence modelling beats independence when dependence exists, but not that this package improves the relevant software frontier.
- **Concrete remedy:** Add a small task-based comparison:
  - (T=2): `GJRM`;
  - (T>2): staged `gamlss2` + `gamCopula`/vine workflow;
  - a current multivariate distributional method where runnable.
  Compare setup complexity, runtime, convergence, prediction/simulation, marginal recovery, dependence recovery, and output usability. If a competitor cannot fit the target model, demonstrate that through a precise capability limitation rather than assertion.
- **Completion test:** At least one reproducible benchmark addresses each nearest-neighbor class, or the manuscript transparently documents why an empirical comparison is infeasible.
- **Confidence:** High.

## M5. Longitudinal indexing and irregular/missing-time semantics are under-specified

- **Type:** Reviewer question with potentially substantive consequences.
- **Location:** p. 12, §3.2.2; pp. 21 and 24–25, Table 8 and §4.3.3; Equation 3 on p. 8.
- **Evidence:** Users supply `time_var` and `subject_var`, but the paper does not state how ties, irregular spacing, unequal visit schedules, absent rows, missing responses, or two observed visits separated by an unobserved scheduled visit are paired. Page 24 says intermittent missingness “splits” trajectories, implying scheduled adjacency rather than observed adjacency, but the data contract is unclear.
- **Consequence:** Two reasonable data encodings may produce different likelihoods. This is central to correctness for real longitudinal data.
- **Concrete remedy:** Define the pair-construction algorithm formally and in code. State whether θ represents dependence per observation step or per unit time and how gaps are treated. Document ordering, duplicate-time policy, unbalanced panels, and missing rows versus `NA` responses.
- **Completion test:** Include a six-row toy panel with a missing middle visit and show exactly which marginal and pair contributions enter the likelihood.
- **Confidence:** High.

## M6. Inferential claims exceed the evidence currently presented

- **Type:** Reviewer question.
- **Location:** pp. 3–4 contribution bullets; pp. 9–11, §§2.4 and variance discussion; pp. 17–20, Tables 5–7; pp. 34–36, Appendix B.
- **Evidence:** The paper repeatedly promises correlation-adjusted standard errors “for all parameters.” Yet p. 36 says smooth-coefficient variances are returned separately using (Q^{-1}\hat\sigma_p^2), rather than clearly deriving them from the full joint covariance. Coverage tables report fixed effects only. Sandwich and bootstrap inference are described but their target, finite-cluster behavior, and smooth-term coverage are not validated.
- **Consequence:** Users may infer that every reported confidence interval has demonstrated joint-likelihood calibration when the evidence supports a narrower statement.
- **Concrete remedy:** Precisely distinguish fixed-coefficient, smooth-function, copula-parameter, numerical-Hessian, sandwich, and bootstrap uncertainty. Add targeted coverage checks for smooth bands and discrete-margin fits, or narrow the claims.
- **Completion test:** Every inferential output named in §3 has an explicit estimator, assumptions, and validation reference; no “all parameters” statement exceeds those results.
- **Confidence:** Medium-high.

## M7. Several methodological comparison statements are inaccurate or too categorical

- **Type:** Confirmed defect.
- **Location:** pp. 2–6, Table 1 and §§1.1–2.1; p. 16, §3.2.5.
- **Evidence:** Examples include:
  - GEE described alongside GLMMs as “jointly optimized,” despite being an estimating-equation/quasi-likelihood approach;
  - vine workflows characterized as inherently two-stage, although joint/full-likelihood constructions are possible;
  - `mgcv` grouped simply with exponential-family mean models despite distributional families;
  - “not possible to preserve correlation when predicting” conflates point predictions with joint predictive simulation;
  - “correlation,” “covariance,” and general copula dependence used interchangeably.
- **Consequence:** These statements invite technical objections and make the comparison appear advocacy-driven.
- **Concrete remedy:** Rewrite comparisons around estimand, objective, and default implementation rather than absolute method capability. Have a GEE/mixed-model and copula specialist audit §§1–2.
- **Completion test:** Each comparative assertion names the package/default/configuration being discussed and avoids universal claims unless mathematically true.
- **Confidence:** High.

## M8. The clinical example is too long, incompletely reproducible, and over-interpreted

- **Type:** Confirmed developmental defect.
- **Location:** pp. 25–30, §5, Tables 12–15 and Figures 6–9.
- **Evidence:** The analysis uses an unexplained 10% subsample “for computational convenience,” conducts extensive sequential covariate selection, and reports a large significance table. Statements about treatment increasing “stability” and unidentified patient characteristics move beyond demonstrating software behavior. The data cannot be redistributed.
- **Consequence:** The example consumes about one-sixth of the main paper, creates inferential distractions, and weakens the reproducible software narrative.
- **Concrete remedy:** Reduce it to one focused scientific question, one model comparison, one effect display, and one diagnostic. Explain the estimand and sampling choice. Prefer the full data if permitted, or treat the example explicitly as nonreproducible supplementary illustration.
- **Completion test:** The clinical case occupies at most 3–4 pages and every central package feature remains demonstrated in the reproducible primary example.
- **Confidence:** High.

## M9. Screening and automated diagnostic rules need stronger qualification

- **Type:** Reviewer question.
- **Location:** pp. 13–15, §§3.2.3 and 3.2.5; `select_margin()`, `select_copula()`, and `check_model()`.
- **Evidence:** Selection relies heavily on AIC after estimated pseudo-observations. `check_model()` uses a KS (p)-value threshold, fixed tail-ratio cutoffs, and a fixed residual-correlation threshold. These rules are sample-size-sensitive and appear to be heuristics, yet “flags failure” suggests a formal adequacy test.
- **Consequence:** Users may treat exploratory screens and diagnostics as automated validation or hypothesis tests.
- **Concrete remedy:** Label these outputs as screening/triage, document their calibration and limitations, avoid “pass/fail” language unless validated, and show how decisions should combine plots, residual checks, and substantive knowledge.
- **Completion test:** Help pages and manuscript use identical qualified language; false-positive behavior is shown across at least two sample sizes or the thresholds are presented only as configurable heuristics.
- **Confidence:** Medium-high.

# MINOR issues

## m1. Table 1’s “maximum formula count” is not a stable or fair comparison metric

- **Type:** Confirmed defect.
- **Location:** pp. 2–3, Table 1.
- **Evidence:** Counts such as (1+T(T-1)/2) for GEE and (PT+QT(T-1)/2) for vine workflows mix parameters, formulas, and possible implementations. Shared covariate models can invalidate these counts.
- **Consequence:** The table gives an impression of exact complexity while comparing unlike interfaces.
- **Concrete remedy:** Replace formula count with “shared formula support,” “pair-specific predictors,” and asymptotic dependence-parameter count under a clearly stated default model.
- **Completion test:** Every table entry is traceable to documented software behavior.
- **Confidence:** High.

## m2. Several central figures and tables are unreadable at normal page size

- **Type:** Confirmed defect.
- **Location:** Table 1, pp. 2–3; Figure 3, p. 16; Tables 10–11, pp. 22–23; Figure 5, p. 24; Table 14, p. 27; Appendix Figures 15–18, pp. 44–47.
- **Evidence:** Multi-panel diagnostics and dense tables use text and labels substantially smaller than surrounding prose; some panels cannot be interpreted without zooming.
- **Consequence:** Even after author notes are removed, readers cannot evaluate the evidence.
- **Concrete remedy:** Split or simplify panels, move detailed tables to supplement, use direct labels, and show only diagnostic panels discussed in the text.
- **Completion test:** All labels and legends remain readable when printed at 100% on A4/Letter without zoom.
- **Confidence:** High.

## m3. Terminology and notation need a controlled vocabulary

- **Type:** Confirmed defect.
- **Location:** Abstract; pp. 2–11, §§1–2; Equation 2 and Tables 2–3.
- **Evidence:** “Fixed” is treated as the opposite of “smooth”; \(\mu\) is sometimes called the mean although it is family-specific location; dependence, correlation, and covariance are interchanged; “CDF” is called “cumulative density”; the discrete rectangle on p. 10 switches from indexed \(a_{it},b_{it}\) to undefined \(a_1,b_1,a_2,b_2\).
- **Consequence:** Interpretation becomes family-dependent and some formulas are harder to verify.
- **Concrete remedy:** Add a short notation/terminology table and distinguish parametric linear terms, smooth terms, location, rank dependence, linear correlation, and covariance.
- **Completion test:** A global search finds consistent usage and every symbol in the discrete likelihood is defined.
- **Confidence:** High.

## m4. The conclusion repeats results rather than consolidating the decision boundary

- **Type:** Optional suggestion.
- **Location:** p. 31, §7.
- **Evidence:** The Discussion restates optimizer and sensitivity results at length after Tables 7–8 already summarize them.
- **Consequence:** The final takeaway is diffuse.
- **Concrete remedy:** Organize the conclusion around: who should use the package, who should not, the main unresolved limitation, and the software roadmap.
- **Completion test:** The Discussion is under two pages and contains no paragraph-length repetition of §4.
- **Confidence:** High.

# Proposed main-paper outline

1. **Introduction and bounded contribution** — 2–3 pages  
   Scientific need, intended user, precise novelty claim, and three nearest software alternatives.

2. **Related software and design choices** — 2–3 pages  
   Capability matrix; explain why shared longitudinal formulas and first-order dependence are deliberate restrictions.

3. **Model and supported data** — 3–4 pages  
   One conceptual figure, linear predictors, joint likelihood, discrete contribution, and assumptions. Move optimizer/Hessian derivations out.

4. **Package architecture and interface** — 3 pages  
   Data contract, supported families/copulas, fitted-object structure, convergence behavior, generics, computational scaling.

5. **Primary reproducible workflow** — 7–8 pages  
   Progressive example: explore → simple fit → diagnose → enrich distribution/dependence → infer → predict/simulate. Use one coherent dataset.

6. **Comparative evaluation** — 4 pages  
   One table of nearest-software comparisons and two focused experiments: correctness/recovery and runtime/robustness. Move simulation grids to supplement.

7. **Applied extension** — 2–3 pages  
   Condensed clinical case or a public-data application, emphasizing one insight unavailable from a mean model.

8. **Limitations, reproducibility, and availability** — 2 pages  
   First-order restriction, irregular/missing visits, inference variants, package version/license/archive, reviewer and full replication profiles.

9. **Discussion** — 1 page  
   Clear use/do-not-use boundary and development roadmap.

Supplementary material should contain optimizer derivations, Hessian formulas, complete simulation grids, detailed clinical selection, and full diagnostic panels.

# Top-10 ordered action list

1. **Rebuild the software/method landscape** and replace the absolute uniqueness claim with a defensible bounded contribution.
2. **Write a one-sentence editorial thesis** naming the user, task, distinctive capability, and intentional first-order limitation.
3. **Choose one progressive, fully reproducible primary example**; begin with a simple fit and add complexity only when motivated.
4. **Re-architect the manuscript around the software workflow** and reduce the main text to roughly 24–28 pages.
5. **Add empirical or task-based comparisons with the nearest copula/distributional software**, not only independence and GEE baselines.
6. **Specify the longitudinal data contract and pair-construction semantics**, including irregular visits, missing middle visits, ties, and gaps.
7. **Audit and narrow inferential claims**, especially smooth-term covariance, discrete likelihood, sandwich inference, and bootstrap coverage.
8. **Finish the exact replication system**, including a measured ≤1-hour reviewer profile and a defensible treatment of unavailable clinical data.
9. **Document the immutable software artifact**: version, archive/DOI, license, installation, dependencies, tested platforms, and citation.
10. **Perform the final technical and presentation pass**: terminology, comparison accuracy, citations, notation, table reduction, and readable figures.
