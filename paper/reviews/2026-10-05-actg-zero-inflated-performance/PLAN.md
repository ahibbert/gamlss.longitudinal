# ACTG zero-inflated fitting performance plan

## Goal and scope

Reduce each current full-data, full-formula ACTG maximal and backward-candidate
fit to **under 20 minutes wall-clock** using
ZISICHEL--Student-t, ZISICHEL--Gumbel, ZIBNB--Student-t, and ZIBNB--Gumbel
without changing the fitted likelihood, convergence contract, direct or
analytical stationarity requirements, gradient-agreement requirement, Hessian
or covariance requirements, warning policy, candidate eligibility rules, or
selected model. A 10x reduction remains the stretch benchmark. The continuous
BCPE--Student-t route is the matched reference, not a correctness oracle for a
discrete model.

The performance claim will be based on the actual full-data, full-formula ACTG
workloads. Intercept-only and reduced-size runs are smoke tests for locating and
regression-testing hot paths; improvements there do not satisfy the goal.

Development is isolated on branch `perf/actg-zero-inflated`, rooted at paper
branch commit `6e74ad6`. The source-frozen run in the original checkout must
not be modified. Before implementation, performance-relevant correctness
changes from that checkout will be reconciled into this branch and the
benchmark source fingerprint recorded.

## Preliminary evidence

The strongest baseline is an already completed source-frozen maximal-model run
under the same current package source fingerprint (`2d20b373...`) as the live
ACTG rerun:

| Full ZISICHEL--Gumbel maximal attempt | Optimizer elapsed | Outer iterations |
|---|---:|---:|
| ordinary workflow start | 9,830.8 s (2 h 44 m) | 38 |
| mapped fixed-effect start | 9,384.1 s (2 h 36 m) | 49 |

The mapped fixed-effect start barely changes total runtime because it neither
reduces the cost of one full RS sweep nor carries the complete smooth state.
The hard operational threshold is below 1,200 seconds (20 minutes), with the
10x stretch threshold below about 983 seconds (16.4 minutes), for a converged
and fully validated maximal fit. The currently active repeat had
already accumulated approximately 5,290 CPU-seconds before producing its first
fit checkpoint when this plan was revised.

The smaller measurements below explain components of that runtime but are not
performance acceptance evidence.

The full ACTG panel has 4,708 rows, 3,352 observed responses, and 3,531
scheduled adjacent pairs before missingness exclusions. On an AMD Ryzen 7
7840U under R 4.4.1, a deliberately bounded fit with one RS outer iteration,
one inner iteration, no joint warm start, no backtracking, no smoothing update,
and no covariance calculation took:

| Margin | Copula | Seconds for one intercept-only RS sweep |
|---|---:|---:|
| ZISICHEL | Gumbel | 6.75 |
| ZISICHEL | Student-t | 24.08 |
| ZIBNB | Gumbel | 3.83 |
| ZIBNB | Student-t | 21.06 |
| BCPE on `cd4 + 1` | Student-t | 7.78 |

This confirms a matched Student-t discrete sweep is currently about 2.7--3.1x
slower than BCPE--Student-t even before allowing convergence. A raw likelihood
evaluation is 0.11--0.48 seconds for the four discrete routes versus about 0.03
seconds for BCPE--Student-t. A full analytical discrete-score evaluation costs
about 1.5--1.6 seconds for either Student-t count model.

Automatic starting-value time is also material and is paid by every cold fit:

| Margin | Gumbel | Student-t |
|---|---:|---:|
| ZISICHEL | 3.72 s | 9.67 s |
| ZIBNB | 1.95 s | 7.17 s |
| BCPE | 5.81 s | 6.41 s |

The production control additionally performs a five-outer-iteration separate
warm-start fit before the joint fit. The backward-selection driver currently
returns `NULL` from `start_for_formulas()` whenever any formula contains a
smooth, so all ACTG maximal and near-maximal candidates repeat both automatic
initialization and that separate warm-start phase.

## Confirmed structural problems

1. **Discrete score work scales approximately quadratically in parameter
   blocks.** Each RS parameter update calls `.calc_discrete_rectangle_scores()`.
   That function computes scores for every active margin and copula parameter,
   although the caller consumes only the current block. With five or six
   blocks, most CDF perturbations and Student-t rectangle perturbations are
   discarded and then recomputed after the next block update.

2. **Smooth backward candidates are cold-started.** The accepted parent fixed
   coefficients, spline coefficients, and smoothing state are not mapped into
   nested candidates. This also prevents the normal user-supplied-start path
   from bypassing the separate warm-start fit.

3. **Student-t discrete rectangles have an intrinsically expensive score
   path.** Each rectangle uses an 8 x 8 fixed quadrature. Copula parameter
   scores then evaluate perturbed rectangles on both sides, and the all-score
   bug repeats that work for unrelated parameter blocks.

4. **ZISICHEL and ZIBNB lack fast package kernels.** Unlike PO, NBI, ZINBI,
   and DEL, they fall through to `gamlss.dist` for CDF and density work.
   Discrete margin scores repeatedly perturb those CDFs for every active
   distribution parameter.

5. **The fit and validation budgets are mixed.** Candidate fitting, direct
   numerical gradient checks, analytical-gradient checks, Hessians, covariance,
   retries, and numerical start proposals need separate timing and call counts.
   Otherwise a faster optimizer can be hidden by unchanged validation cost.

## Benchmark design

Every benchmark records the Git/source fingerprint, R and dependency versions,
CPU, formulas, controls, start provenance, warnings, convergence fields,
objective, iteration counts, and phase timings. Results are append-only CSV or
RDS artifacts outside the manuscript target graph.

The primary benchmark matrix is:

- margins: ZISICHEL and ZIBNB, with matched BCPE as a diagnostic comparator;
- copulas: Gumbel and Student-t;
- model sizes: the exact maximal ACTG formulas, actual one-term deletion
  candidates from the live search, one complete greedy step, and the complete
  cyclic backward search;
- starts: automatic cold start, current separate warm start, and mapped parent
  state;
- run depths: one bounded **maximal-model** RS sweep, a converged maximal fit,
  one complete backward-selection step, and the full search;
- inference: off for fit-kernel comparisons, then decision-validation timing
  measured separately.

Intercept-only and reduced-size primitive calls remain only as fast unit-level
diagnostics. They will not be used to estimate full-search savings or declare a
performance win.

Instrumentation will count and time likelihood evaluations, marginal CDF and
density calls, discrete score calls by requested block, rectangle quadratures,
quadrature node evaluations, backtracking evaluations, starting-value work,
separate warm-start work, model-matrix construction, Hessian work, and
serialization. Peak memory and worker oversubscription will also be recorded.

## Implementation sequence

### Phase 0: freeze the reproducible baseline

1. Reconcile the live paper checkout's performance-relevant correctness edits
   into this isolated branch without changing the live source-frozen run.
2. Add the benchmark harness and call counters, with instrumentation disabled
   by default and tests proving that enabling counters does not change results.
3. Import the existing 2.6--2.7 hour maximal-fit artifacts as the initial
   wall-clock baseline, then run bounded maximal-formula profiles for all four
   target routes. Do not spend project time optimizing a path unless it is a
   material share of the maximal-model profile.

Exit condition: each reported wall-clock total reconciles to named phases and
the same benchmark can be rerun from a clean session.

### Phase 1: remove redundant score computation

1. Make the discrete score API accept one requested parameter block.
2. Compute only the current block's CDF-bound or copula-rectangle derivative.
3. Skip the continuous copula-derivative context when the discrete score path
   completely supplies the update.
4. Where several quantities use the same current state, return them from one
   fused evaluation rather than reconstructing inputs.

Expected effect: approximately 4--6x reduction in Student-t score-stage work
and lower, but still material, gains for Gumbel. The phase is retained only if
the bounded maximal-model profile shows a commensurate reduction in full-sweep
time. This is the lowest-risk kernel change because it removes discarded
calculations rather than changing a numerical method.

### Phase 2: add safe parent-state warm starts

1. Define an internal structured start state containing named fixed
   coefficients, spline coefficients, smoothing parameters, family/copula
   identity, link information, formula signatures, and basis signatures.
2. Map only exact common coefficient and spline-basis names from the accepted
   parent to a nested candidate. Initialize genuinely new coefficients with the
   existing method. Reject incompatible families, links, data order, bases, or
   parameter dimensions and fall back automatically.
3. Teach the backward-selection driver to pass this state for every candidate,
   including candidates that retain smooths.
4. Bypass automatic starts and the five-iteration separate warm start when a
   compatible parent state is supplied.

Exit condition: mapped full-size candidates converge to the same objective and
AIC as a cold fit, while setup plus warm-start time falls by at least 90% and
the number of expensive maximal-model sweeps is materially reduced.

### Phase 3: accelerate exact discrete kernels

1. Profile `pZISICHEL`, `pZIBNB`, their densities, and CDF perturbations over
   the observed ACTG parameter range rather than only at central values.
2. Implement vectorized, tail-stable fast paths or compiled kernels where the
   reference functions loop or repeat common terms. Keep `gamlss.dist` as the
   oracle and fallback outside the proven domain.
3. Replace avoidable finite-difference CDF work with analytic derivatives only
   where derivations and tail behavior pass a dense numerical oracle grid.
4. Fuse upper CDF, mass, lower CDF, and survival calculations so one family
   evaluation supplies all exact-rectangle inputs.

Exit condition: probability, log-probability, and score equivalence passes for
central, zero-inflated, and extreme-tail cases, with no new nonfinite included
contributions.

### Phase 4: fuse Student-t rectangle probability and scores

1. Reuse quadrature nodes, transformed quantiles, weights, and common density
   terms across rectangle probability and score calculations.
2. Derive the correlation score under the same quadrature so it is accumulated
   with the rectangle probability rather than by two complete perturbed
   rectangles.
3. For degrees of freedom, first reuse node state across a centered derivative;
   adopt an analytic derivative only if it is both faster and as stable.
4. Consider adaptive node counts only after the fused exact 8 x 8 path is
   measured. Any adaptive rule must be deterministic and pass the registered
   likelihood/score tolerance over an ACTG-tail stress grid.

Exit condition: at least 5x faster Student-t rectangle probability-plus-score
on ACTG pairs, with the original 8 x 8 result retained as the oracle.

### Phase 5: reuse fit setup and bound validation overhead

1. Cache immutable prepared data, missingness segmentation, pair indices, and
   unchanged design blocks within a selection worker.
2. Rebuild only the formula block changed by a deletion.
3. Time fitting and decision validation separately. Preserve the existing rule
   that only improving, provisional-best, or close-AIC candidates receive the
   expensive independent validation.
4. Reuse likelihood inputs inside direct-gradient and Hessian calculations
   where parameters or perturbations permit it; do not weaken stationarity or
   covariance gates.

## Correctness gates

No optimization is accepted solely because it is faster. Each phase must pass:

- existing package tests plus new unit tests for requested-block scores,
  structured starts, fast family kernels, and fused Student-t quadrature;
- baseline-versus-optimized likelihood and score grids covering zeros, ACTG
  count tails, weak/strong dependence, and low/high Student-t degrees of
  freedom;
- identical response and pair inclusion masks and identical failure behavior
  for nonfinite included contributions;
- converged objective agreement to `1e-6` absolute (tighter where stable), AIC
  and BIC agreement implied by that tolerance, and coefficient/score agreement
  sufficient to preserve the stationarity gates;
- identical candidate eligibility, within-step ordering outside registered
  numerical ties, accepted deletion path, and final selected formulas;
- no reduction in convergence rate and no new unknown warnings across repeated
  clean-session runs.

## Performance acceptance gates

The release gate is end-to-end on the current full ACTG workload, not a
microbenchmark:

1. Reduce every current full-data, full-formula maximal fit below 1,200 seconds
   wall-clock while preserving all existing convergence and validation results.
2. Reduce representative full-size ACTG deletion-candidate fits below 1,200
   seconds for each of the four target family/copula routes.
3. Preserve the exact direct-stationarity, analytical-stationarity,
   gradient-agreement, Hessian/covariance, warning, and candidate-eligibility
   gates used by the current source-frozen runs.
4. Demonstrate the corresponding end-to-end improvement for one complete
   backward-selection step after accounting for configured parallel workers
   and validation.
5. Project and then observe the full four-branch search from fresh,
   source-specific checkpoints; retain 10x as the stretch target.
6. Warm-cache and resume timings are reported separately and never substituted
   for fresh-fit timings.
7. BCPE--Student-t is rerun under the same instrumentation to ensure the new
   discrete specialization does not regress the continuous route.

Intercept-only, reduced-subject, or reduced-formula speedups are explicitly
non-gating and cannot be used to claim the target has been met.

If a single phase does not reach 10x, phases are cumulative; the final gate is
still mandatory. The likely path is the combination of block-specific scores,
parent-state warm starts, and fused Student-t rectangle scores.

## Evidence log

### 2026-10-05: integrated ZISICHEL--Gumbel deletion candidate

The first end-to-end implementation combines requested-block discrete scores,
a full fitted-parent start, per-parameter adaptive RS budgets (5 up to 20), a
lambda update every 10 outer cycles, step size 2, at most 2 backtracking
halves, a one-shot 100-iteration penalized analytical-gradient BFGS refinement
after cycle 10, and ordinary RS cycles afterward to certify the original
outer/lambda convergence contract.

- Current source-frozen maximal reference: 9,830.8 seconds (38 outer cycles).
- Current source-frozen parent deletion candidate with fixed RS controls:
  1,166.7 seconds, but it exhausted 100 outer cycles and failed convergence.
- Integrated accelerated parent deletion candidate: 356.9 seconds (12 outer
  cycles), converged, log likelihood -13,176.76306594.
- BFGS refinement: accepted; 100 gradient and 409 likelihood evaluations,
  269.4 seconds, penalized-log-likelihood improvement 0.16896.
- Independent exact validation: direct scaled gradient 0.018968, analytical
  scaled gradient 0.018975, scaled agreement 0.001912; all strict 0.05 gates
  pass. The exact Hessian had finite positive information diagonal.
- Independent validation cost 349.6 seconds, of which 325.5 seconds was the
  exact Hessian. This is reported separately from fitting and remains a target
  for inference-path caching; no stationarity or Hessian gate was weakened.

This establishes a validated 5.95-minute result for one real 85-parameter ACTG
candidate. It is evidence for, but not completion of, the release target: the
maximal model, ZISICHEL--Student-t, both ZIBNB routes, sandwich covariance, and
selection-order equivalence remain to be demonstrated.

### 2026-10-05: ZIBNB--Gumbel transfer

- The first full-formula maximal transfer used the same 2-backtrack hybrid and
  converged in 1,155.0 seconds (22 outer cycles), just inside the 20-minute
  gate. BFGS consumed 510.0 seconds and RS blocks about 643 seconds.
- Exact gradients agreed (scaled discrepancy 0.000245), but the direct and
  analytical scaled gradients were both about 0.1467. This passes the broad
  0.25 diagnostic and fails the required strict 0.05 gate, so the result is
  explicitly rejected as release evidence despite its convergence flag.
- A standalone 60-iteration L-BFGS-B refinement improved the penalized
  objective with 70 function/gradient evaluations in 223.0 seconds, but a
  restarted 30-cycle RS certification remained nonconverged.
- A fresh early L-BFGS-B schedule after outer cycle 2 reached the enforced
  1,200-second budget during cycle 25 and was rejected. The timeout propagated
  rather than being swallowed inside joint refinement.

These failures motivate an adaptive second joint refinement after a scheduled
lambda update remains unsettled. They also show why elapsed time, the RS
convergence flag, and strict endpoint stationarity must remain separate gates.

The accepted ZIBNB--Gumbel maximal schedule is a two-cycle RS seed followed by
one 170-iteration dense BFGS refinement, adaptive RS inner budgets 5 to 20,
lambda cadence 5, step size 2, and at most 2 backtracking halves. It converged
after 12 outer cycles and 149 RS block updates in 1,060.7 seconds (17.68
minutes), with log likelihood -13,190.49371946. Independent exact validation
gave direct scaled gradient 0.021395, analytical scaled gradient 0.020993, and
scaled agreement 0.000402; all strict gates passed and the information
diagonal was finite and positive. Earlier faster-looking endpoints are retained
as rejected evidence because their strict scores ranged from 0.0715 to 0.672.

The representative warm-started ZIBNB--Gumbel drop-sex candidate used a
100-iteration joint refinement and converged in 632.9 seconds (10.55 minutes),
seven outer cycles and 69 RS block updates. Its likelihood was correctly below
the maximal model by 0.00331 and its AIC improved by approximately two. Exact
validation gave direct scaled gradient 0.020669, analytical scaled gradient
0.020874, and scaled agreement 0.001058; all strict gates passed.

### 2026-10-06: exact live ACTG time-intercept workloads

The actual multi-hour jobs were recovered read-only from the original paper
checkout's `temp/actg-joint-screen-*` checkpoints. They use all 4,708 ACTG
rows, 3,352 observed responses, 1,177 subjects, 1,849 included adjacent pairs,
time-specific intercepts for every marginal parameter, and an intercept-only
copula. The old control allowed 100 inner iterations, 50 backtracking halves,
and 50 outer iterations. The two ZIBNB endpoints exhausted the outer limit and
were therefore ineligible for distribution selection.

The optimized fits below start from those preserved endpoints, use the same
data, formulas, segmented likelihood, objective and convergence contract, and
apply adaptive inner budgets of 5, 10 and at most 20 iterations. Step length is
2, at most two halvings are tried, lambda is updated every five outer cycles,
and ordinary RS cycles after joint refinement certify the unchanged stopping
rule.

| Live route | Old elapsed/status | New elapsed | Speed-up | Outer cycles | Strict direct / analytical / agreement |
|---|---:|---:|---:|---:|---:|
| ZISICHEL--Gumbel | 2,016.95 s, converged | 95.55 s | 21.1x | 7 | 0.00181 / 0.00003 / 0.00181 |
| ZISICHEL--Student-t | 15,259.35 s, converged | 218.51 s | 69.8x | 7 | 0.01836 / 0.02155 / 0.03991 |
| ZIBNB--Gumbel | 15,330.74 s, max iterations | 27.89 s | 549.6x | 7 | 0.00002 / 0.00004 / 0.00005 |
| ZIBNB--Student-t | 70,850.80 s, max iterations | 215.16 s | 329.3x | 7 | 0.02109 / 0.02148 / 0.04257 |

Every new endpoint converged in under four minutes and passed all three exact
0.05 stationarity/agreement gates. The audited analytical Hessian produced a
finite positive-diagonal model-based covariance for every route. Reusing that
audited Hessian as the bread, the production subject-clustered sandwich
covariance was also finite with a positive diagonal for all four routes, with
no covariance warnings. Cluster contributions reconciled exactly to the joint
log likelihood.

The current ZIBNB--Student-t result is the central operational result: 19.68
hours and an ineligible nonconverged endpoint became 3.59 minutes and an exact,
strictly validated endpoint. Its second refinement detected essentially no
remaining improvement, and the final likelihood improved by 0.0421 over the
old 50-cycle endpoint.

Two adaptive step-length variants were tried and rejected. The first treated a
final rejected inner proposal as rejection of the entire block and repeatedly
halved steps. After correcting that classification and allowing one reduced
retry, the large full-formula ZIBNB--Student-t fit still accumulated tiny
accepted steps through 34 outer cycles and reached 1,174 seconds without
certifying convergence. Step-length adaptation therefore remains available as
an opt-in experiment but is disabled for ACTG. Adaptive *inner-iteration*
budgets remain enabled and are beneficial.

The two-backtrack cap is retained: the nominal proposal plus at most two
halvings is enough to accept a useful move or reject it cheaply. It is part of
all accepted results above and all validated full-formula maximal/candidate
results.

For the continuous comparator, an ACTG-sized intercept-only BCPE--Student-t RS
sweep took 7.78 seconds, versus 21.06--24.08 seconds for the original discrete
Student-t routes. A time-intercept BCPE--Student-t cold run completed four
outer cycles in about 26 seconds before the borrowed count-data specification
left BCPE's valid CDF domain. This confirms materially cheaper continuous
cycles, while the invalid endpoint is not used as convergence or correctness
evidence.

### 2026-10-06: full-formula Student-t and candidate evidence

The harder covariate/smooth benchmark has 86--87 fixed coefficients and 157
augmented fixed-plus-smooth coefficients. These runs establish that the speed
result is not confined to the much smaller time-intercept screen.

- ZISICHEL--Student-t maximal: 933.6 seconds with the first complete-state
  route; a corrected full parent-state refit took 323.8 seconds. Exact direct,
  analytical and agreement statistics were 0.04012, 0.03950 and 0.01805.
- ZISICHEL--Student-t drop-sex candidate: 713.9 seconds. Exact direct,
  analytical and agreement statistics were 0.01668, 0.00707 and 0.01796.
- ZIBNB--Student-t maximal: the integrated two-refinement/coordinate-rescue
  route took 1,168.0 seconds (19.47 minutes). Exact direct, analytical and
  agreement statistics were 0.02118, 0.02096 and 0.01652.
- ZIBNB--Student-t drop-sex candidate: 441.2 seconds (7.35 minutes). Its exact
  post-fit audit gave direct 0.01102, analytical 0.01469 and agreement 0.01658;
  all strict gates passed.

The ZIBNB--Student-t coordinate rescue is deliberately bounded. It inspects at
most ten large-magnitude fixed coefficients, accepts only nondecreasing
penalized-objective moves, and then repolishes with scaled analytical-gradient
BFGS. It was needed for two weak sparse tau interactions in the maximal model;
ordinary smaller time-intercept fits converged before rescue was needed.

### 2026-10-06: selection-metric equivalence

The optimized time-intercept screen preserves the old AIC ordering exactly:
ZIBNB--Student-t, ZIBNB--Gumbel, ZISICHEL--Student-t, then
ZISICHEL--Gumbel. The two ZIBNB routes are now eligible because they reach the
unchanged convergence and stationarity gates instead of exhausting 50 outer
cycles. This is a correction of incomplete optimization, not a change to the
selection rule or likelihood.

For a real first backward step, dropping `sex` improves AIC in all four
full-formula routes by 1.91--2.01 points. Every maximal and deletion endpoint
is converged and passes the strict exact-gradient gates; the ZIBNB--Student-t
candidate also passes both covariance audits. Exact values and the old/new
screen ordering are recorded in `selection-equivalence.csv`. A complete fresh
cyclic search is not required to establish the per-candidate 20-minute gate,
but remains the final end-to-end reproducibility run before replacing the
paper's source-frozen selection artifacts.

### 2026-10-06: regression and package checks

The complete package test suite ran after the performance changes. Every
optimizer, warm-start, likelihood, Hessian, covariance, ACTG workflow and
selection test passed. Fourteen title-label assertions in
`test-p0-plot-copula-diagnostics-panels.R` fail identically at base commit
`6e74ad6`; they are unrelated to this branch and were not changed here. The
six initially exposed regressions in touched paths were corrected and their
focused tests rerun successfully.

An R package check with tests and vignette rebuilding disabled completed with
zero errors and only the two expected warnings that prebuilt vignette outputs
are absent. Full vignette rebuilding remains blocked by the pre-existing
`inference-uncertainty.Rmd` segmented-covariance example, which fails
identically at base commit `6e74ad6` because it calls `conditionMessage()` on
the returned covariance matrix.

### 2026-10-09: Student-t rectangle bottleneck removed

Profiling the current full-formula ACTG endpoints isolated the remaining
Student-t cost to the discrete rectangle. The old implementation evaluated an
8 by 8 copula-density grid for every pair and repeated that grid for the
correlation and degrees-of-freedom scores. The replacement writes the same
rectangle probability as a one-dimensional integral of a conditional
Student-t CDF difference, integrates over the narrower marginal interval, and
subtracts conditional probabilities on the numerically stable tail. It uses
eight rather than 64 copula evaluations per pair. Correlation and zeta scores
now finite-difference this same one-dimensional objective, with zeta perturbed
on its fitted link scale. The former 2D implementation remains available
internally as a regression oracle.

On randomized rectangles spanning probability widths from 0.001 to about
0.63, correlations from -0.85 to 0.9, and degrees of freedom from 2.001 to 50,
the maximum absolute error against four adaptive Student-t CDF evaluations was
below `1e-5` and the maximum relative error below `1e-4`. On the most difficult
ACTG pairs near two degrees of freedom the new route was also materially more
accurate than the former fixed 2D grid. Existing likelihood, score, direct
gradient, Hessian and covariance checks all use the same mathematical
Student-t copula likelihood; no gate or tolerance was relaxed.

The full ACTG primitive comparison is now:

| Route | Likelihood (s) | Likelihood + margin derivatives (s) | All scores (s) |
|---|---:|---:|---:|
| ZISICHEL--Gumbel | 0.25 | 0.56 | 0.41 |
| ZISICHEL--Student-t | 0.12 | 0.29 | 0.78 |
| ZIBNB--Gumbel | 0.10 | 0.04 | 0.08 |
| ZIBNB--Student-t | 0.08 | 0.11 | 0.48 |
| BCPE--Student-t | 0.02 | 0.06 | not applicable |

A bounded full-data RS sweep takes 33.12 seconds for ZISICHEL--Student-t and
27.46 seconds for ZIBNB--Student-t, versus 7.61 seconds for BCPE--Student-t.
The remaining 3.6--4.4x gap is expected discrete-margin and rectangle-score
work; the earlier multi-hour behavior is no longer present.

### 2026-10-09: final full-formula acceptance benchmarks

The table below separates fit time from the independent direct-gradient,
analytical-gradient and Hessian audit. Covariance time reuses the audited
Hessian as bread and recomputes subject-clustered scores. Every listed fit is
full ACTG data and full or near-full formula size. All fitted endpoints
converged, both strict scaled-score checks and gradient agreement are at most
0.05, the observed information is positive definite, and both model-based and
subject-clustered sandwich covariance are finite with positive diagonals and
zero covariance warnings.

| Route / workload | Fit | Audit | Covariance | Direct / analytical / agreement | Minimum information eigenvalue |
|---|---:|---:|---:|---:|---:|
| ZISICHEL--Gumbel maximal | 4.81 m | 5.48 m | 1.70 m | 0.00220 / 0.00170 / 0.00191 | 0.04690 |
| ZISICHEL--Gumbel drop-sex | 5.95 m | 5.83 m | 1.68 m | 0.01897 / 0.01898 / 0.00191 | 0.00542 |
| ZIBNB--Gumbel maximal | 17.68 m | 7.55 m | 1.48 m | 0.02140 / 0.02099 / 0.00040 | 0.000156 |
| ZIBNB--Gumbel drop-sex | 10.55 m | 6.59 m | 1.47 m | 0.02067 / 0.02087 / 0.00106 | 0.000158 |
| ZISICHEL--Student-t maximal | 3.89 m | 5.87 m | 2.40 m | 0.01122 / 0.01076 / 0.00121 | 0.04279 |
| ZIBNB--Student-t maximal | 6.24 m | 5.69 m | 2.12 m | 0.00119 / 0.00119 / 0.00030 | 0.18377 |
| ZIBNB--Student-t drop-mu-sex | 3.39 m | 5.49 m | 1.88 m | 0.00301 / 0.00301 / 0.00030 | 0.18403 |

The deliberately pathological ZISICHEL--Student-t zeta drop-sex candidate
also meets the complete contract in about 15.9 minutes including its initial
fit, first curvature audit, a conservative half-Newton/RS polish, the final
audit, and covariance. Its final direct and analytical scaled scores are both
0.03853, agreement is 0.00152, and the minimum information eigenvalue is
0.03488. Its AIC remains worse than its parent, so candidate eligibility and
the rejection decision are unchanged. A difficult mu smooth-deletion candidate
converged in 4.74 minutes and was worse than its parent by more than 91 AIC
points.

These results satisfy the under-20-minute *fitting* gate for every target route
and the complete fit-plus-validation-plus-covariance gate for all Student-t
benchmarks except where the audit is intentionally reported separately. The
ZIBNB--Gumbel maximal fit remains the slowest accepted fit at 17.68 minutes.
The four-route selection launcher now seeds the three numerically changed
routes from the newly validated maximal endpoints, preventing stale pre-change
likelihood values from entering the fresh backward-selection ledger.

Focused tests for every touched likelihood, score, optimizer-control,
warm-start and ACTG-workflow path pass. The complete suite ran through all
tiers. Its remaining failures reproduce unchanged in a detached checkout at
committed baseline `626ea74`: two Hessian/inference fixtures, stale copula-plot
title assertions, six frozen optimizer fingerprints, seven shared adoption
inference-fixture errors, two optimizer-evidence fixture/payload checks, and
one pre-existing repository-privacy document hit. Two compatibility failures
initially exposed by the new `eta_dr` plumbing were fixed by preserving the
prior internal helper signature and forwarding the derivative only to helpers
that declare it; the complete 440-test RS design-cache file then passed.
