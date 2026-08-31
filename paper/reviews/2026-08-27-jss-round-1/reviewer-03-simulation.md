# Overall recommendation

**Reject in present form, with encouragement to resubmit after a comprehensive empirical rebuild.** The manuscript’s empirical programme is promising, but several central claims do not correspond to the staged public inputs. Replication counts and designs are inconsistent, the optimizer conclusions contradict the displayed tables, and the alleged dropout study does not simulate dropout. These are evidence-integrity problems rather than editorial details.

Reviewed frozen manuscript at source SHA `68c3bad26626ce7c267bd330364cfda8df7a6b76`. No prior audits/reviews were used and no files were edited during the review.

## Five strengths

1. The validation programme targets important software properties: parameter recovery, prediction, optimizer behaviour, dependence misspecification, missingness, and comparison with standard longitudinal methods.
2. Simulations include known truths for both marginal and copula parameters, including smooth and non-location effects.
3. The authors use held-out log and variogram scores in addition to in-sample likelihood, a sound principle for optimizer comparisons.
4. Public assets include useful per-replicate data for several studies, deterministic seeds, manifests, and provenance metadata.
5. The manuscript recognizes first-order-dependence limitations and sensibly exposes sandwich inference as a diagnostic or robustness option.

# BLOCKING issues

## B1. Main Monte Carlo design does not describe the results presented

- **Location:** PDF p. 18, Section 4.1, Table 4; pp. 18–19, Figure 4 and Tables 5–6.
- **Evidence:** Table 4 declares `n = 250, 1000`, `T = 5, 20`, and **1,000 replicates** for both BCPE/t and negative-binomial/Clayton recovery. The staged BCPE files contain one scenario, `n500_d4`, with `n=500`, `T=4`, and 100 replicates. The staged NBI files likewise contain 100 replicates and an `n500_d4_nbi_signal2` design. Tables 5–6 numerically match the 100-replicate BCPE inputs. The manifest marks the inline design table as `pending_externalization` and manually maintained.
- **Consequence:** Readers cannot determine what experiment generated the paper’s headline recovery, coverage, runtime, and predictive results. Claimed operating characteristics across sample size and number of visits have not been demonstrated by the displayed evidence.
- **Remedy:** Either rerun the complete Table 4 factorial design or rewrite Table 4 and all prose to the actual `n=500, T=4, R=100` experiment. Generate the design table directly from scenario metadata.
- **Completion test:** For every table/figure, an automated check shows identical scenario IDs, `n`, `T`, attempted replicates, successful replicates, and seeds in the manuscript, producer input, and attempt-level output.
- **Type:** **Confirmed defect**
- **Confidence:** **Very high**

## B2. Optimizer study design, replication counts, and conclusions are internally inconsistent

- **Location:** PDF pp. 39–41, Appendix F, Tables 18–21; conclusions on pp. 20 and 31.
- **Evidence:**
  - Page 39 states 100 replications for every example. Public Normal data contain only **10 replicates per case**; Gamma and NBI contain 100.
  - Table 18 presents a common design for all three families. Public Gamma JVS6 is `n=50, T=40, rows=2000`, not `n=100, T=75, rows=7500`. Public NBI JVS1–3/5/6 use `n=20, T=20, rows=400`; JVS4 uses `n=40, T=10`, not the Table 18 designs. NBI JVS6 does not increase sample size at all.
  - The text says JVS3 improves both variogram metrics “for both the median and IQR.” Every JVS3 table has at least one IQR crossing zero; in fact both variogram IQRs cross zero for Normal, Gamma, and NBI.
  - The claimed JVS5 reduction is contradicted by Normal train-LL improvement increasing from 11.2 to 20.4 and Gamma from 8.7 to 9.0.
  - The claimed JVS6 sample-size benefit is contradicted by Gamma train-LL improvement falling from 8.7 to 6.5 and NBI from 3.2 to 1.2.
  - Favorable held-out variogram proportions are modest: for JVS3, 0.60/0.60 in Normal, 0.62/0.60 in Gamma, and 0.62/0.70 in NBI for \(p=0.5/2\).
- **Consequence:** The evidence does not support the stated causal guidance about correlation, visit count, shared effects, or sample size. The central recommendation to prefer joint optimization for prediction is much stronger than the held-out evidence.
- **Remedy:** Replace the selected heterogeneous cases with a declared factorial experiment varying one dimension at a time, using paired seeds and identical designs across families. Report paired effect estimates, uncertainty, sign probabilities, runtime, and failures. Rewrite conclusions from those results.
- **Completion test:** Table 18 is machine-generated from each per-replicate file; every advertised contrast changes only its named factor; narrative claims pass automated direction/sign checks against Tables 19–21.
- **Type:** **Confirmed defect**
- **Confidence:** **Very high**

## B3. The “dropout” experiment does not simulate dropout and omits substantial failures

- **Location:** PDF pp. 24–25, Section 4.3.3; Figures 13–14 on pp. 42–43; Discussion p. 31.
- **Evidence:** The manuscript repeatedly calls `time_mar` “time-based dropout.” In `run_missingness_study.R`, both `mar` and `time_mar` assign missingness independently to each row with `runif(nrow(dat)) < miss_prob`; `time_mar` merely adds time to the missingness logit. It does not enforce monotone post-dropout missingness. Public provenance records only **20 replicates**, with **390 successful and 90 failed fits out of 480** (18.75% failures). The paper shows success-only aggregate curves with no failure counts or uncertainty.
- **Consequence:** The evidence cannot support conclusions about dropout. Success-conditioned curves may be materially optimistic, especially if failure depends on missingness rate, mechanism, or method.
- **Remedy:** Rename the current mechanism “time-dependent intermittent MAR,” add genuine monotone dropout with a subject-level event/hazard, and ideally include an outcome-history-dependent MAR scenario. Increase replication and publish scenario/method-specific attempted, converged, failed, and included counts.
- **Completion test:** Every subject in the dropout scenario has no observed response after its first dropout time; attempt-level outputs reconcile to all summary denominators; failure rates and failure-inclusive sensitivity bounds are reported.
- **Type:** **Confirmed defect**
- **Confidence:** **Very high**

# MAJOR issues

## M1. Monte Carlo uncertainty is largely absent

- **Location:** Tables 5–6, 10–11, 19–21; Figures 4–5 and 13–14.
- **Evidence:** Coverage, RMSE, selection rates, medians, and differences are interpreted without Monte Carlo standard errors or confidence intervals. The correlation benchmark uses 20 replicates per family/scenario; pooled coverage is therefore based on about 60 values, for which MCSE near 0.95 is approximately 0.028. Differences such as 0.950 versus 0.983 are not demonstrably meaningful. Missingness has only 20 attempts per cell and varying successful denominators.
- **Consequence:** Terms such as “highly accurate,” “substantially,” and “better” frequently exceed the precision of the experiments.
- **Remedy:** Predefine precision targets; report MCSE or bootstrap intervals for all summaries and paired intervals for method differences. Show denominators.
- **Completion test:** Each empirical estimate has attempted/successful \(R\) and MC uncertainty; headline comparisons are either supported by intervals or qualified.
- **Type:** **Confirmed defect**
- **Confidence:** **High**

## M2. Copula-selection counts are stale relative to staged results

- **Location:** PDF p. 24, Section 4.3.2, Figure 5 and final paragraph; Discussion p. 31.
- **Evidence:** The manuscript reports Clayton selected in 58/60 cases and, for Gaussian moderate dependence at `n=50, T=4`, selections of Gaussian 1/10, Clayton 3/10, Frank 4/10, and Gumbel 2/10. The public file contains 100 replicates per cell. Recalculation gives Gaussian 38/100, Clayton 12/100, Frank 32/100, Gumbel 14/100, and t 4/100 for that cell. Clayton at high dependence and `n=500` is 100/100.
- **Consequence:** The manuscript is not reporting the staged experiment, and the broad “AIC is generally very reliable” conclusion lacks a current confusion matrix and uncertainty.
- **Remedy:** Recompute all text, heatmaps, and selection rates from the staged data; report selection matrices by generator, dependence, and sample size.
- **Completion test:** Published counts sum to the attempt-level cell sizes and regenerate exactly from `results.csv`.
- **Type:** **Confirmed defect**
- **Confidence:** **Very high**

## M3. Correlation/GEE benchmark overstates differences and uses an intentionally unstable comparator

- **Location:** PDF pp. 21–23, Tables 9–11.
- **Evidence:** Tables average over three response families, masking heterogeneity. The unstructured GEE estimates 190 dependence parameters with only 120 subjects at `T=20`; instability is therefore largely designed in. Nevertheless, the unstable comparator is used rhetorically to strengthen the package’s advantage. In COV-COR, mean RMSE is tied between `gamlss.longitudinal` and AR(1) GEE at 0.062; GEE coverage values of 0.983 cannot be distinguished from 0.950 at the available replication count. The claim of better marginal recovery “than all the GEEs” is too categorical.
- **Consequence:** Comparator fairness and the magnitude of superiority are unclear.
- **Remedy:** Report family-specific results, attempt/failure counts, and uncertainty. Treat high-dimensional unstructured GEE as a stress test, not a primary comparator. Add feasible spline/basis-based time correlation alternatives if available, or explicitly scope the claim to the implemented GEE structures.
- **Completion test:** Main claims hold within families and against feasible comparators with paired uncertainty; unstable methods are labeled and not used as decisive evidence.
- **Type:** **Confirmed defect plus comparator-design question**
- **Confidence:** **High**

## M4. CG optimizer guidance is not empirically demonstrated

- **Location:** PDF p. 20, Table 7 and accompanying text; Discussion p. 31.
- **Evidence:** The manuscript claims CG and joint RS are generally interchangeable, CG is generally slower, and CG is useful for convergence sensitivity. The staged optimizer tables compare only `rs_joint` with `rs_separate`; the public manifest contains no CG comparison table.
- **Consequence:** User-facing recommendations about a principal software option are unsupported.
- **Remedy:** Add a paired RS-joint/CG study covering fixed-only and smooth models, initialization, convergence, objective attained, coefficient/smooth differences, runtime, and failure modes—or clearly label the CG guidance as qualitative.
- **Completion test:** A public per-attempt CG/RS table supports each claim in Table 7.
- **Type:** **Confirmed defect**
- **Confidence:** **High**

## M5. Failure and convergence accounting is incomplete or non-comparable

- **Location:** Throughout Section 4 and Appendices F–H.
- **Evidence:** Main recovery logs report 100/100 package fits converged, but comparator convergence is `NA`, not established. The missingness study loses 90/480 fits. Correlation tables omit attempted/successful counts. Optimizer tables report only paired successes and do not expose warning/stop distributions in the manuscript.
- **Consequence:** Performance may be conditional on being fit successfully, and methods with different failure rates cannot be compared fairly.
- **Remedy:** Publish attempt-level status for every method: error, convergence criterion, gradient/step diagnostics, boundary estimates, warnings, timeout, retry policy, and whether included. Report both conditional and failure-penalized summaries.
- **Completion test:** Attempt totals reconcile exactly across status, estimate, metric, and summary files; no `NA` convergence state is interpreted as success.
- **Type:** **Confirmed defect**
- **Confidence:** **High**

## M6. Parameter-recovery conclusions conceal weak recovery of the t-copula shape parameter

- **Location:** PDF pp. 18–19, Figure 4 and Table 6.
- **Evidence:** The \(\zeta\) terms have RMSEs 0.812–1.331 and coverage 0.98–0.99, while other effects are much more precise. Figure 4 visibly shows very wide \(\zeta\) error ranges. The text nevertheless broadly describes highly accurate recovery.
- **Consequence:** Readers may infer dependable estimation of tail dependence when the study instead suggests weak identification and conservative intervals.
- **Remedy:** Discuss \(\zeta\) separately; report bias and RMSE on both link and natural degrees-of-freedom/tail-dependence scales, boundary rates, interval widths, and sample-size sensitivity.
- **Completion test:** Conclusions explicitly distinguish good location/correlation recovery from weak tail-shape recovery.
- **Type:** **Confirmed defect**
- **Confidence:** **High**

## M7. Clinical claims are exploratory but written as confirmatory treatment effects

- **Location:** PDF pp. 25–30, Section 5, especially Tables 14–15 and text on pp. 27 and 29–30.
- **Evidence:** The analysis uses a 10% subsample of 902 patients, extensive data-driven AIC/LRT reduction, nominal post-selection p-values, 14.3% overall missingness and 36% missingness at the final visit, but no resampling across subsamples, missingness sensitivity, robust/bootstrapped inference, or adjustment for selection. Higher adjacent copula correlation is translated directly into a treatment-induced “increase in stability.”
- **Consequence:** Nominal significance and causal-sounding stability conclusions may reflect subsampling, selection, missingness, or modeling choices.
- **Remedy:** Provide the sampling seed and flow diagram; run the full cohort or repeated-subsample sensitivity; compare model-based, sandwich, and bootstrap intervals; assess dropout by treatment; and describe copula effects as conditional associations unless causal identification is established.
- **Completion test:** Treatment/dependence conclusions remain directionally stable across the full cohort or prespecified resamples and robust inference/missingness analyses.
- **Type:** **Confirmed overclaim**
- **Confidence:** **High**

## M8. Public clinical code does not reproduce the reported model

- **Location:** PDF p. 30, Section 6; clinical Tables 12–15 and Figures 6–9.
- **Evidence:** The paper first states that “all results” are directly replicable, then excludes the private clinical analysis. The public sanitized recipe uses much simpler formulas—e.g. \(\theta \sim\) treatment + time—whereas Table 15 reports treatment, myocardial-infarction history, and dyspnoea effects plus many additional marginal covariates. Clinical assets are static hashes, not regenerated outputs.
- **Consequence:** Even an authorized data holder lacks an exact executable analysis path for the published clinical result.
- **Remedy:** Release immutable, data-free analysis code containing exact preprocessing, sampling, formulas, selection sequence, controls, diagnostics, and table/figure production, plus an input schema and private-data checksum.
- **Completion test:** With authorized data, one command regenerates Tables 12–15 and Figures 6–9 byte-for-byte or within declared numerical tolerances.
- **Type:** **Confirmed defect**
- **Confidence:** **Very high**

## M9. Tail-fit interpretation is not aligned with Table 13

- **Location:** PDF pp. 26 and 29, Table 13 and discussion of the selected t-copula.
- **Evidence:** The empirical lower/upper 5% co-occurrences are 2.63/1.77. The t-copula gives 2.65/2.65, while Gaussian gives 2.31/2.31. The t model closely matches the lower tail but overstates the upper tail more than Gaussian. Claims that it more accurately captures “the observed stronger tail correlation” and opposite-tail behaviour are not supported by the displayed same-tail table.
- **Consequence:** The likelihood advantage is real, but the stated substantive tail explanation is incomplete and potentially misleading.
- **Remedy:** Report lower, upper, and opposite-tail diagnostics separately by treatment and time, with uncertainty; distinguish overall AIC superiority from tail-specific calibration.
- **Completion test:** Every tail claim points to a displayed statistic showing better calibration in that tail.
- **Type:** **Confirmed overstatement**
- **Confidence:** **High**

# MINOR issues

## m1. Data-dependent one-sided trimming of a comparator metric

- **Location:** PDF p. 19, Table 5 footnote.
- **Evidence:** Seven outer-fence `gamlss2` \(p=2\) variogram scores are excluded, giving 139.45 rather than the raw mean 40,881.93; no equivalent prespecified rule is described for both methods.
- **Consequence:** The primary mean is neither fully raw nor symmetrically robust.
- **Remedy:** Make paired median/IQR or a prespecified transformed/trimmed estimator primary; give raw and robust summaries for both methods.
- **Completion test:** The inclusion rule is declared before method labels are inspected and applied identically.
- **Type:** **Confirmed defect**
- **Confidence:** **High**

## m2. “JVS4 only decreases \(T\)” is not literally true

- **Location:** PDF p. 39, Table 18 and following prose.
- **Evidence:** JVS3 to JVS4 changes \(T\) from 75 to 20 and subjects from 40 to 150 while holding total rows fixed.
- **Consequence:** The contrast mixes visit count with number of independent clusters.
- **Remedy:** State that it changes cluster geometry at fixed row count, or construct an experiment separately varying \(T\), subjects, and rows.
- **Completion test:** Factor-isolation language matches the generated design.
- **Type:** **Confirmed defect**
- **Confidence:** **Very high**

## m3. Smooth “50 percent range bands” are insufficiently defined

- **Location:** PDF p. 37, Figure 10; p. 43, Figure 14.
- **Evidence:** It is unclear whether bands are 25th–75th pointwise quantiles, central ranges, or ranges across successful fits. No simultaneous recovery/coverage measure is supplied.
- **Consequence:** Smooth recovery cannot be quantitatively assessed.
- **Remedy:** Define the band and add integrated bias/RMSE plus pointwise and simultaneous coverage where inferential smooth claims are made.
- **Completion test:** Captions state estimand, denominator, quantiles, and failure handling.
- **Type:** **Question/suggestion**
- **Confidence:** **High**

## m4. Clinical sample-size labels are ambiguous

- **Location:** PDF p. 27, Table 14.
- **Evidence:** `n=6314` equals 902×7 scheduled rows despite 14.3% response missingness. BIC also appears to use 6314. The effective observed-response and adjacent-pair counts are not reported.
- **Consequence:** Readers cannot assess likelihood denominators or BIC’s sample-size convention in clustered incomplete data.
- **Remedy:** Report subjects, scheduled rows, observed margins, adjacent pairs, and the BIC \(n\) definition.
- **Completion test:** All fit criteria state their effective likelihood unit and sample-size convention.
- **Type:** **Confirmed ambiguity**
- **Confidence:** **High**

## m5. Reproducibility instructions remain unfinished

- **Location:** PDF pp. 30–31, Section 6.
- **Evidence:** The text calls the deprecated `expanded` profile, retains “X hours,” and contains a reminder to report OS, R/package versions, cores, seeds, and runtime.
- **Consequence:** Submission-level computational expectations are not yet specified.
- **Remedy:** Replace with current `paper`/`full` commands and measured hardware-specific runtimes and session metadata.
- **Completion test:** A clean reviewer run follows the manuscript verbatim and produces the declared outputs and metadata.
- **Type:** **Confirmed defect**
- **Confidence:** **Very high**

# Minimum evidence package required for submission

At minimum, submission should include:

1. A machine-readable protocol listing every scenario, DGP, link-scale truth, `n`, `T`, attempted replicates, seed rule, comparator, tuning, estimand, and metric.
2. Attempt-level status files for every method/scenario/replicate, including errors, warnings, convergence diagnostics, timeouts, boundary estimates, retries, and inclusion decisions.
3. Per-replicate estimates, standard errors, confidence intervals, truths, smooth-grid estimates, and paired predictive scores.
4. Summary tables containing attempted/successful denominators, bias, empirical SD, mean reported SE, RMSE/IRMSE, coverage, runtime, MCSE, and uncertainty intervals.
5. Prespecified paired comparison rules and identical handling of outliers/failures across methods.
6. A true monotone-dropout study plus intermittent MAR, with missing-pattern and retained-pair summaries.
7. Copula-selection confusion matrices by generator, dependence level, and sample size.
8. Family-specific GEE benchmark results and justification for feasible comparator parameterizations.
9. Exact data-free clinical scripts, sampling seed, input checksum/schema, participant/observation flow, model-selection trace, and robust/missingness sensitivity outputs.
10. A one-command manifest that regenerates every public table/figure from attempt-level inputs and verifies scenario metadata and hashes.

# Top-10 ordered action list

1. Reconcile Table 4, Tables 5–6, and all recovery assets to one actual Monte Carlo design.
2. Rebuild the optimizer study with consistent designs, paired seeds, adequate replication, and corrected conclusions.
3. Replace the pseudo-dropout mechanism with genuine monotone dropout and report all 90 current failures by cell.
4. Add attempt-level convergence/failure accounting throughout.
5. Add MCSEs or uncertainty intervals and successful/attempted denominators to every empirical result.
6. Add a real CG-versus-joint-RS comparison or withdraw empirical CG guidance.
7. Regenerate copula-selection text and figures from the current 100-replicate public file.
8. Reframe and strengthen the GEE benchmark with family-specific, feasible, uncertainty-qualified comparisons.
9. Make the clinical analysis exactly executable for authorized users and temper treatment/stability claims pending robust and missingness sensitivity.
10. Finalize the reproducibility section with current profiles, measured runtime, session information, and automated manuscript–asset consistency checks.
