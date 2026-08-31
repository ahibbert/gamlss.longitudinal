# Statistical-methods review

Reviewed the frozen manuscript :codex-file-citation{path="C:\\Users\\Aydin\\AppData\\Local\\Temp\\gamlss-jss-review-round-20260827\\paper-repo\\main.pdf" purpose="source"} and `main.tex` at source SHA `68c3bad26626ce7c267bd330364cfda8df7a6b76`, with selective checks of implementation and staged public assets.

## Recommendation

**Major revision; not publishable in its current form.**

The continuous first-order Markov-copula idea is sound and potentially valuable. However, four issues block acceptance: the implemented likelihood can silently omit invalid contributions; Hessian code can turn invalid curvature into finite standard errors; the reported Monte Carlo design does not match the staged results; and the missingness/dropout claims do not correspond to either the implemented observed-data objective or the actual simulation.

## Five strengths

1. The continuous-data objective is a clear, parsimonious first-tree-truncated D-vine/Markov-copula factorization with distributional regression for both margins and dependence.
2. The package contains a genuine exact discrete–discrete rectangle-mass implementation and targeted derivative/Hessian tests, rather than relying only on midpoint pseudo-observations.
3. The manuscript acknowledges the principal structural limitation—unmodelled conditional dependence above tree 1—and includes useful comparisons with AR(1), exchangeable, and covariate-dependent settings.
4. Subject-cluster sandwich and both parametric and cluster bootstrap machinery are useful additions for sensitivity analysis.
5. The formula interface, prediction/simulation workflow, and separation of marginal and dependence diagnostics could make a technically difficult model accessible to applied users once inferential scope is made precise.

# BLOCKING issues

## B1. Parameter-dependent likelihood contributions can be silently discarded

- **Location:** PDF pp. 8 and 10, Section 2.3, Equation (3) and the discrete-likelihood derivation; TeX lines 536–558 and 699–731.
- **Evidence:** `R/likelihood-evaluation-margin.R:161–162` converts nonfinite or nonpositive observed marginal densities to `NA`; `R/likelihood-evaluation-copula.R:105–107` replaces invalid/nonpositive copula contributions by `1`; `R/likelihood-evaluation.R:64–69` sums only finite terms. Thus an invalid proposal can lose an observation or pair instead of receiving log-likelihood `-Inf`.
- **Consequence:** The optimized function need not be the displayed likelihood, and its effective sample can vary with the parameters. Bias, convergence, likelihood comparisons, AIC/BIC, Hessians, and simulations based on the fitted optimum can all be affected.
- **Concrete remedy:** Evaluate marginal and rectangle contributions on the log scale; define a fixed inclusion mask based only on genuinely missing responses; use stable tail/CDF calculations; reject a parameter vector with `-Inf` whenever an included observed contribution is invalid. Never substitute independence or omit a term because of numerical failure.
- **Completion test:** Add adversarial continuous and count tests with extreme responses and boundary-adjacent parameters demonstrating that the number of included marginal/pair terms remains constant, invalid proposals cannot improve the objective by losing terms, and package log-likelihoods agree with independently calculated values.
- **Confidence:** High.
- **Classification:** Confirmed defect.

## B2. Absolute values mask invalid Hessian-based variance estimates

- **Location:** PDF p. 11, TeX lines 731–740; Appendix B, PDF p. 34 and p. 36, TeX lines 1812–1820 and 1947–1954.
- **Evidence:** The manuscript reports \(\sqrt{|(H^{-1})_{aa}|}\). `R/model-vcov-hessian.R:5–21` computes `vc = -solve(H)` but standard errors as `sqrt(abs(diag(solve(H))))`; diagnostics retain only absolute eigenvalue magnitudes and do not verify that the log-likelihood Hessian is negative definite.
- **Consequence:** A saddle point, nonidentified fit, or incorrect Hessian can be reported with apparently ordinary finite standard errors. This directly undermines the central claim of immediately interpretable adjusted inference.
- **Concrete remedy:** Require symmetry, full numerical rank, acceptable conditioning, and the correct curvature sign before model-based inference. Compute `sqrt(diag(-solve(H)))` without `abs`; return an explicit inference failure when variances are nonpositive. Report signed extreme eigenvalues, rank, and gradient norm.
- **Completion test:** Tests must reject deliberately indefinite and singular Hessians, accept well-conditioned negative-definite Hessians, and show agreement between analytical and numerical Hessians/SEs over representative continuous and discrete fits.
- **Confidence:** High.
- **Classification:** Confirmed defect.

## B3. The published Monte Carlo design does not match the staged evidence

- **Location:** PDF p. 18, Table 4 and the accompanying recovery/coverage claims; TeX lines 1045–1111.
- **Evidence:** Table 4 states \(n=250,1000\), \(T=5,20\), and 1,000 replicates. The staged BCPE/t `fit_run_log.csv` contains one `n500_d4` scenario, 100 replicates per model, and 200 rows total. The staged NBI/Clayton coefficient asset likewise contains 100 unique replicates per model. The displayed coverage values occur in increments consistent with 100 replicates.
- **Consequence:** Readers cannot determine what experiment supports Figures 4 and Tables 5–6. Monte Carlo precision and generalization across \(n,T\) are overstated, and the results are not reproducibly tied to the stated design.
- **Concrete remedy:** Either rerun the declared design or rewrite Table 4 and every associated claim to match the actual scenario. Report attempted, successful, converged, and retained replicates by method/scenario, exclusion rules, and Monte Carlo standard errors.
- **Completion test:** A generated metadata table must agree automatically with the unique `n`, `T`, replicate IDs, methods, and success counts in the public assets; the manuscript tables/figures must be regenerated from those same rows.
- **Confidence:** High.
- **Classification:** Confirmed defect.

## B4. The missingness analysis is not a dropout study and does not use the observed-data likelihood across gaps

- **Location:** PDF pp. 24–25, Section 4.3.3; TeX lines 1280–1290; Discussion p. 31, lines 1624–1630.
- **Evidence:** For \(Y_1,Y_3\) observed and \(Y_2\) missing, the Markov model’s observed-data likelihood integrates over \(Y_2\), inducing a transition between \(Y_1\) and \(Y_3\). The implementation instead removes both adjacent pair factors and treats the two segments as independent. The simulation labelled `time_mar` uses independent Bernoulli deletion with a time-dependent logit (`run_missingness_study.R:401–425`), not monotone dropout. Staged provenance reports only 20 replicates per scenario and 390 successful fits out of 480 attempted, while the manuscript calls the results “quite reassuring” without reporting failures.
- **Consequence:** With intermittent gaps the objective is a composite/pseudo-likelihood, not the stated observed-data joint likelihood. “50 percent time-based dropout” and reassuring MAR/dropout claims are unsupported. The clinical example also requires explicit assessment of intermittent versus monotone missingness.
- **Concrete remedy:** Either implement integration/forward recursion across missing states, or explicitly redefine the gap-case objective as composite likelihood and use Godambe/cluster-bootstrap inference. Run a genuinely monotone dropout design, disclose the missingness model and failures, and separately study intermittent MAR. In the clinical example, report split trajectories and perform robust sensitivity analysis.
- **Completion test:** Verify the likelihood analytically for three- and four-visit missing patterns; show monotonicity of missingness in every simulated dropout trajectory; publish scenario-level attempted/successful counts and Monte Carlo uncertainty.
- **Confidence:** High.
- **Classification:** Confirmed defect.

# MAJOR issues

## M1. The first-order D-vine assumption should be stated as a conditional Markov model

- **Location:** PDF pp. 6–8, Sections 2.1–2.3, especially TeX lines 350, 364–381, and 536–551; PDF p. 21, lines 1224–1226.
- **Evidence:** Tree-1 truncation asserts conditional independence of nonadjacent outcomes given the intervening outcome(s) and covariates. It does not universally imply monotonically decreasing Pearson correlation or Kendall’s \(\tau\), especially with nonstationary families, changing signs, or covariate-dependent transitions. Implementation evaluates pair-copula covariates at the left/earlier row, which is mentioned only in Appendix B.
- **Consequence:** Users may mistake an AR(1) illustration for the defining assumption and may encode time-varying covariates at the wrong endpoint.
- **Concrete remedy:** Write the conditional density \(f(y_{it}\mid y_{i,t-1},X_i)\), state exogeneity/conditioning assumptions, identify the endpoint at which pair covariates are evaluated, and distinguish conditional independence from empirical lag-decay heuristics.
- **Completion test:** The mathematical definition, simulation, prediction, and model-matrix documentation must give the same left/right endpoint mapping and reproduce known nonstationary transition examples.
- **Confidence:** High.
- **Classification:** Confirmed overstatement plus reviewer question.

## M2. Continuous, discrete, and mixed-margin scope is incomplete

- **Location:** PDF p. 7, TeX lines 355–366; PDF p. 10, lines 699–729.
- **Evidence:** The manuscript derives continuous–continuous and discrete–discrete contributions only. In the discrete paragraph, \(\Delta C/(d_1d_2)\) is called a “contribution to the log likelihood,” but its logarithm is required. The implementation accepts one shared marginal family, so continuous–discrete transitions are not supported. For discrete margins, Sklar’s copula is not unique outside the support grid, which is not discussed.
- **Consequence:** “Any supported margin” can be misread as permitting mixed outcome types across visits, and dependence interpretation for counts is stronger than the model identifies.
- **Concrete remedy:** State explicitly that a fit uses one homogeneous marginal family. Either derive and implement continuous–discrete and discrete–continuous contributions using copula partial derivatives, or declare them unsupported. Correct the log-ratio notation and discuss nonuniqueness/parametric identification for discrete copulas.
- **Completion test:** Documentation examples must reject unsupported mixed-family input clearly, or numerical tests must validate all four pair types against independent calculations.
- **Confidence:** High.
- **Classification:** Confirmed defect.

## M3. Parameter links and effective domains are insufficiently specified

- **Location:** PDF p. 7, Equation (2), Table 2, and the paragraph after Equation (2); TeX lines 368–381 and 393–441.
- **Evidence:** The manuscript writes a “link function \(g_p(\eta_p)\),” whereas the parameter is \(p=g_p^{-1}(\eta_p)\). Actual links are \(\exp(\eta)\) for Clayton, identity for Frank, \(1+\exp(\eta)\) for Joe/Gumbel, \(\tanh(\eta)\) for Gaussian/t correlation, and \(2+\exp(\eta)\) for t degrees of freedom. Gumbel is capped at 17 in its inverse link, and Clayton evaluation is capped near 28 in likelihood code, neither matching the displayed domains.
- **Consequence:** Reproduction, coefficient interpretation, gradient/Hessian derivations, and boundary identification are ambiguous; hidden caps create flat or discontinuous effective objectives.
- **Concrete remedy:** Add a link table containing link, inverse link, derivative, implemented numerical bounds, independence point, and warnings. Remove arbitrary caps where stable computation is possible or expose and justify them.
- **Completion test:** Round-trip and derivative tests over the documented domain, with explicit warnings/tests at every numerical cap.
- **Confidence:** High.
- **Classification:** Confirmed defect.

## M4. Discrete and automated diagnostics are not statistically calibrated

- **Location:** PDF p. 15, Section 3.2.5, especially the `check_model()` paragraph; TeX lines 961–981. Clinical diagnostic claim at PDF p. 29, lines 1589–1591.
- **Evidence:** `check_model()` calls `.gl_pit(..., randomize = FALSE)`. For discrete margins this is the upper CDF, not uniform PIT. Copula diagnostics similarly transform the upper CDF and apply continuous h-functions. Even for continuous margins, the KS reference distribution ignores fitted parameters and within-subject dependence. The manuscript also reverses the KS rule, saying failure occurs for \(p>0.05\), while code uses \(p<0.05\).
- **Consequence:** Correct discrete models can be falsely failed, misspecified ones passed, and the advertised automatic PASS/FAIL has no stated size.
- **Concrete remedy:** Use randomized or appropriately nonrandomized discrete residual diagnostics, derive the discrete conditional Rosenblatt transform, and calibrate omnibus thresholds by parametric/cluster bootstrap or label them descriptive screens.
- **Completion test:** Repeated correctly specified continuous and discrete simulations should show calibrated false-flag rates; deliberately misspecified cases should demonstrate power. Fix the p-value inequality.
- **Confidence:** High.
- **Classification:** Confirmed defect.

## M5. Smooth-term and fixed-effect uncertainty are conditional approximations, not full joint inference

- **Location:** Contribution claim on PDF p. 4, TeX line 292; Appendix B.1.5 on PDF p. 36, lines 1957–1985.
- **Evidence:** The analytical/numerical Hessians cover `object$par` fixed coefficients while holding smooth coefficients fixed. Smooth covariance is computed separately as \(Q^{-1}/\bar w\), or with residual variance from the mean curve, ignoring fixed–smooth and cross-parameter/copula covariance and smoothing-parameter uncertainty. Sandwich scores likewise vary fixed coefficients only.
- **Consequence:** The statement that standard errors are supplied for “all parameters” and automatically account for correlation is false for smooth terms; fixed-effect SEs are also conditional on estimated smooths.
- **Concrete remedy:** Assemble and invert the full penalized Hessian for fixed and smooth coefficients, retaining all cross-blocks, and state whether inference is conditional on \(\lambda\). Alternatively label the current smooth output explicitly as a heuristic and remove inferential claims.
- **Completion test:** Full-matrix analytical/numerical agreement and coverage studies for smooth functions and fixed coefficients in models where the same covariate enters margin and copula smooths.
- **Confidence:** High.
- **Classification:** Confirmed defect.

## M6. The CG algorithm appendix does not describe the implemented optimizer

- **Location:** PDF pp. 33–34, Appendix A, TeX lines 1690–1784; optimizer-objective comparison on PDF p. 9, lines 569–575.
- **Evidence:** The manuscript describes solving \((A+\rho I)d=g\), filtering by a quadratic predicted improvement, using an actual/predicted ratio, and updating the radius at 0.25/0.75 thresholds. Code constructs `-solve(H_pen - rho I, g)`, evaluates actual penalized-objective improvement only, applies an absolute Armijo-like threshold, and has no predicted-improvement ratio calculation.
- **Consequence:** The algorithm is not reproducible from the paper, and convergence/optimizer comparisons cannot be interpreted against the stated method.
- **Concrete remedy:** Rewrite the appendix from the implemented equations, including signs, ridge convention, candidate duplication, acceptance rule, objective, smoothing-parameter schedule, raw-likelihood safeguards, and convergence criteria.
- **Completion test:** A deterministic one-iteration test should reproduce every candidate, accepted step, objective change, trust-radius update, and stopping flag from the documented formulas.
- **Confidence:** High.
- **Classification:** Confirmed defect.

## M7. Sandwich and bootstrap claims need narrower estimands and stronger evidence

- **Location:** PDF p. 17, TeX lines 1008–1027; sensitivity guidance pp. 20–21, lines 1204–1245; Discussion p. 31, line 1624.
- **Evidence:** The sandwich estimator covers fixed coefficients only and holds smooth coefficients fixed. The parametric bootstrap simulates from the fitted first-order model, so it is not robust to the dependence misspecification for which the discussion recommends it. The displayed percentile example uses only \(R=20\), and failed refits are omitted from quantile calculations.
- **Consequence:** Users may treat either option as a general robustness correction, including for smooths, copula misspecification, or missingness, when it is not.
- **Concrete remedy:** Distinguish model-based, subject-cluster sandwich, cluster bootstrap, and parametric-bootstrap targets. Recommend materially larger \(R\), report failure rates/effective replicates, and demonstrate coverage under correct and misspecified structures.
- **Completion test:** Coverage tables should separately cover fixed marginal, fixed copula, and smooth estimands, with predefined failure handling and Monte Carlo uncertainty.
- **Confidence:** High.
- **Classification:** Confirmed scope defect.

# MINOR issues

## N1. Comparative claims about GEE, GLMM, and two-stage vine inference are inaccurate

- **Location:** PDF pp. 2–5, particularly TeX lines 275–285 and 327–334.
- **Evidence:** GEE is not jointly optimized by a full likelihood; GLMM coefficients are conditional rather than “adjusted marginal” coefficients; two-stage/IFM copula estimation has established sandwich/Godambe covariance methods and does not inherently require simulation.
- **Consequence:** The contribution is framed against straw-man alternatives.
- **Concrete remedy:** Recast comparisons in terms of estimand, objective, and available covariance correction rather than “joint optimization.”
- **Completion test:** A methods expert should be able to map each comparator row to its objective and estimand without contradiction.
- **Confidence:** High.
- **Classification:** Confirmed defect.

## N2. “Correlation,” “covariance,” and copula parameter are used interchangeably

- **Location:** Throughout Sections 2–5, including TeX lines 355, 366, 793, 1224–1249, and 1591–1593.
- **Evidence:** Clayton, Gumbel, Joe, Frank, and t degrees of freedom parameters are not covariance parameters; only Gaussian/t \(\rho\) is directly a latent correlation.
- **Consequence:** Coefficient interpretation and cross-family comparisons are unclear.
- **Concrete remedy:** Use “dependence parameter” generically and report Kendall’s \(\tau\) or tail-dependence transformations for common interpretation.
- **Completion test:** No Archimedean \(\theta\) is described as correlation/covariance unless accompanied by the exact transformation.
- **Confidence:** High.
- **Classification:** Confirmed terminology defect.

## N3. Copula-selection conclusions are overgeneralized

- **Location:** PDF p. 24, TeX lines 1270–1278; Discussion p. 31, line 1624.
- **Evidence:** “AIC is generally very reliable” is based on a limited Gamma/intercept, \(T=4\), candidate-set experiment; the manuscript itself reports Gaussian selection only 1/10 at \(n=50\), moderate dependence.
- **Consequence:** Readers may underestimate selection uncertainty and post-selection inference effects.
- **Concrete remedy:** Restrict the claim to the simulated high-dependence/high-\(n\) settings and recommend model averaging or sensitivity reporting when AIC differences are small.
- **Completion test:** State scenario-specific selection rates with uncertainty and avoid an unconditional reliability claim.
- **Confidence:** High.
- **Classification:** Confirmed overstatement.

# Claims requiring narrowing or additional evidence

- “Any supported marginal parameter” should mean compatible linked parameters from one shared family per fit.
- “Directly interpretable standard errors for all parameters” should be restricted to regular fixed-coefficient fits with a valid Hessian; current smooth inference is approximate.
- “First-order dependence requires correlation to decay with lag” should be replaced by the conditional Markov-independence assumption.
- “Supports fitting missing observations” must distinguish complete prefixes, intermittent-gap composite likelihood, and an actual observed-data likelihood.
- “Time-based dropout” must be removed unless monotone dropout is simulated.
- “Bootstrap as a robustness check” must distinguish parametric from cluster bootstrap.
- “AIC is generally reliable” must be scenario-specific.
- “Joint models outperform” should identify the exact objective, scenarios, convergence filters, and Monte Carlo uncertainty.
- “All standard tools” should be narrowed where diagnostics/inference are not valid for discrete margins.
- The novelty claim that no direct practitioner-facing software exists needs a documented, current software search rather than an absolute assertion.

# Top-10 ordered action list

1. Make the likelihood contribution mask fixed and reject invalid numerical proposals instead of dropping terms.
2. Remove absolute-value Hessian SEs and enforce negative-definite, full-rank curvature diagnostics.
3. Reconcile Table 4 and all validation results with the actual staged scenario metadata and replicate counts.
4. Redesign/rewrite the missingness analysis: true dropout, failure disclosure, and observed-data or explicitly composite likelihood.
5. Replace smooth and conditional fixed-effect covariance claims with full penalized joint inference or clearly labelled approximations.
6. Formalize the model as a covariate-conditional first-order Markov copula and specify pair-covariate endpoint mapping.
7. Correct and complete the continuous/discrete/mixed likelihood and discrete identifiability discussion.
8. Repair discrete PIT/Rosenblatt diagnostics and calibrate `check_model()` thresholds.
9. Rewrite the optimizer appendix to match the implemented RS/CG objectives and step rules exactly.
10. Add a complete parameter-link/boundary table and narrow comparative, AIC, bootstrap, and novelty claims.
