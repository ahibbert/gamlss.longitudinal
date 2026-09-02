# Optimizer guidance evidence checklist

Status: empirical CG performance guidance is withdrawn from the Phase 2
optimizer benchmark.

- The registered empirical methods are exactly `rs_separate` and `rs_joint`.
- The v7 design uses Gaussian-copula routes throughout. JVS01--JVS03 use the
  Normal margin and JVS04 changes only the family to ZIP, retaining the same
  panel, predictors, signal strengths, time shape, and Gaussian copula.
- JVS04 registers a Poisson-component mean of 4 and a structural-zero
  probability of 0.20 on the natural parameter scale. Its optimizer-scale
  intercepts are therefore `log(4)` and `qlogis(0.20)`; it must not inherit the
  Normal scale intercept `log(0.75)`.
- Every registered margin/copula route is synchronously capability-checked
  before a run lock is acquired, workers are started, or checkpoints are read
  or written.
- The checkpoint namespace and schema are `paired-one-factor-v7` and 7. The 291
  v6 checkpoints were generated under the rejected NBI/Gaussian route and are
  preserved in the `rejected-unsupported-route` archive; they are not eligible
  for v7 resumption or publication evidence.
- Do not claim that CG and joint RS are interchangeable.
- Do not claim that CG is generally faster or slower.
- Do not recommend CG as an empirical convergence-sensitivity check.
- Any future CG performance statement requires a separately registered,
  precision-driven paired benchmark with attempts, retained fits, failure
  reasons, runtime, objective, recovery differences, MCSEs, and intervals.
- Deterministic implementation evidence is narrower than performance evidence.
  Candidate construction and line-search/acceptance sequencing are exercised by
  `tests/testthat/test-p0-optimizer-cg-helpers.R`, especially the curvature and
  line-search sequencing test and the candidate-step/acceptance helper tests.
  Those tests may support code/document agreement for one iteration; they do not
  support comparative speed, robustness, or interchangeability claims.

## Protected-manuscript blockers

The following edits must be made by the protected-file integrator before the
optimizer wording audit can pass:

- Around line 520, replace the sentence beginning “The RS method has the
  advantage” with: “The RS method uses first derivatives and updates included
  parameter blocks sequentially.” No production timing comparison is available
  yet.
- Around line 522, replace the two sentences beginning “The second method”
  with: “The second method is Cole and Green (CG), which uses a full Hessian to
  update all included parameter blocks within a trust region.”
- Around line 1109, delete the editorial suggestion to use separate RS “when
  speed matters” and joint RS/CG for final fits.
- Around lines 1113–1115, replace the performance and workflow recommendation
  paragraphs with: “The registered comparison evaluates RS with separately
  optimized margin and copula blocks against RS optimization of the joint
  likelihood. CG is not included in this empirical comparison. We report
  convergence, runtime, likelihood, predictive-score, and recovery differences
  with Monte Carlo uncertainty; optimizer guidance is limited to the resulting
  registered contrasts.” This replacement is valid only after production
  evidence is available; until then, state that results are pending. Remove the
  claims that joint methods are substantially slower and that separate RS should
  later be refined to a joint model.
- Around line 1684, replace the paragraph claiming “substantial difficulties”,
  fits “stuck in local optimum locations”, and improved CG “performance” with:
  “The CG implementation constrains each simultaneous Hessian-based update using
  the configured trust-region and line-search controls described below.” This
  states implementation mechanics without unsupported empirical robustness or
  local-optima claims.
- In the conclusion (currently around line 1565 after earlier edits), remove the
  unfrozen claim that the default jointly optimized model “generally
  outperforms” separate optimization for particular correlation, panel, or
  shared-covariate settings. While production is pending, replace it with: “A
  registered paired comparison of separate and joint RS optimization is in
  progress; any directional conclusion will be inserted only from frozen
  attempt-level evidence with Monte Carlo uncertainty.”
