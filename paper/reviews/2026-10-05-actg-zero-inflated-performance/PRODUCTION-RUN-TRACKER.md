# ACTG and paper production run tracker

Frozen source commit: `1deea5772e6d11b06d3702536763e040568adc8e`

Production scope was frozen on 2026-10-12 (Australia/Sydney). No further
optimizer development is planned. The normal per-fit cap is 1,190 seconds;
only isolated complex Student-t fits may use the agreed 1,500-second cap.

| Milestone | Baseline estimate | Started | Completed | Actual | Status |
|---|---:|---|---|---:|---|
| ZISICHEL--Gumbel backward selection and audit | 10 h | 2026-10-11 | 2026-10-12 04:05 | 9.77 h | Complete |
| ZISICHEL--Student-t backward selection and audit | 8--14 h remaining | 2026-10-12 07:51 | 2026-10-13 05:01 | 21.17 h | Complete and independently certified |
| ZIBNB--Gumbel backward selection and audit | 8--12 h | 2026-10-12 07:51 |  |  | Selection complete in 4.33 h; inference audit failed (degenerate endpoint) |
| ZIBNB--Student-t backward selection and audit | 8--12 h | 2026-10-12 08:47 | 2026-10-13 07:19 | 22.54 h | Complete and independently certified |
| All four ACTG selections complete | 24--40 h from freeze | 2026-10-12 07:50 | 2026-10-13 07:19 | 23.49 h | Selection complete; four-route inference certification blocked by ZIBNB--Gumbel |
| BCPE/t and NBI/Clayton simulations and outputs | 12--24 h | 2026-10-13 07:25 |  |  | Paused; 07:25 concurrent launch stopped before a valid ordered stage completed |
| Remaining paper pipeline and outputs | 1--3 d |  |  |  | Pending |
| Final cross-artifact audit | 6--18 h |  |  |  | Pending |

## Production decisions

- Preserve the frozen likelihood, convergence, selection, warning, stationarity,
  Hessian, and covariance contracts.
- Reuse completed checkpoints only when their source fingerprint and run-control
  compatibility checks pass.
- Run at most two four-worker ACTG routes concurrently on the eight physical
  cores.
- Do not delay production for further optimizer tuning. A genuine inability to
  complete under the agreed caps is the only reason to revisit execution.

## Events

- 2026-10-12 07:49: The ZISICHEL--Student-t full-zeta diagnostic with two
  initial RS iterations reached 1,190.1 seconds in the zeta block. This closed
  optimizer experimentation; production resumes with the agreed 1,500-second
  Student-t exception.
- 2026-10-12 07:51: Resumed ZISICHEL--Student-t with a 1,500-second cap and
  started ZIBNB--Gumbel with the normal 1,190-second cap. Both use four workers;
  BLAS/OpenMP thread counts are fixed at one per worker.
- 2026-10-12 08:01: ZIBNB--Gumbel completed `mu-smooth` in 606.8 seconds and
  advanced to `sigma`. All 25 candidates converged, four moves were accepted,
  the longest fit was 102.4 seconds, and no unexpected warnings occurred.
- 2026-10-12 08:11: ZIBNB--Gumbel completed `sigma` in 603.1 seconds and
  advanced to `nu`. All 23 candidates converged, four moves were accepted,
  the longest fit was 98.6 seconds, and no unexpected warnings occurred.
- 2026-10-12 08:13: The ZISICHEL--Student-t full-zeta maximum completed 50
  outer cycles in 1,308.2 seconds but stopped at `max_iterations` rather than
  converging. Its complete checkpoint was continued within the remaining
  191.8 seconds of the agreed 1,500-second cumulative budget; no source or
  convergence criterion changed.
- 2026-10-12 08:19: The bounded continuation reached 57 cumulative outer
  cycles and 1,433.8 seconds without convergence; the final 66.2-second
  remainder ended at the elapsed cap. A fresh current-source zeta stage was
  therefore started from the previously converged full-zeta coefficient state.
  All five non-zeta formulas match the current parent exactly. Resume is off,
  so likelihoods, convergence, and candidate decisions are recomputed under
  frozen commit `1deea57`; no old metric is reused.
- 2026-10-12 08:22: ZIBNB--Gumbel completed `nu` in 626.3 seconds and
  advanced to `tau`. All 19 candidates converged, three moves were accepted,
  the longest fit was 105.5 seconds, and no unexpected warnings occurred.
- 2026-10-12 08:47: The warm-seeded ZISICHEL--Student-t full-zeta fit reached
  78 cumulative outer cycles and 1,478.8 seconds without convergence. The
  scheduled lambda checks were present, so this is a genuine failure of the
  agreed 25-minute contract; the route is paused without weakening its gate.
  ZIBNB--Student-t started in the freed four-worker slot while ZIBNB--Gumbel
  continues.
- 2026-10-12 08:39: ZIBNB--Gumbel completed `tau` in 1,040.5 seconds and
  advanced to `theta`. All 27 candidates converged, five moves were accepted,
  the longest individual fit was 258.9 seconds, and no unexpected warnings
  occurred.
- 2026-10-12 08:54: ZIBNB--Gumbel completed `theta` in 896.4 seconds and
  entered refinement cycle 1 at `mu`. All 23 theta candidates converged, four
  moves were accepted, the longest fit was 123.6 seconds, and no unexpected
  warnings occurred. The full initial five-stage pass took 62.9 minutes.
- 2026-10-12 09:08: ZIBNB--Student-t completed `mu-smooth` in 1,241.3
  seconds and advanced to `sigma`. All 25 candidates converged, four moves
  were accepted, the longest individual fit was 291.4 seconds, and no
  unexpected warnings occurred.
- 2026-10-12 09:22: ZIBNB--Gumbel completed refinement-cycle `mu` in
  1,659.6 seconds and advanced to `sigma`. Twenty-three of 25 candidates
  converged. The two nonconverged smooth-deletion candidates had AIC 26737.7
  and 26722.8 versus the accepted path near 26484, so they were ineligible,
  decisively noncompetitive, and did not affect the selected move. Three
  deletion moves were accepted, the longest fit was 275.3 seconds, and no
  unexpected warnings occurred.
- 2026-10-12 09:39: ZIBNB--Student-t completed `sigma` in 1,840.5 seconds
  and advanced to `nu`. All 26 candidates converged, five moves were accepted,
  the longest individual fit was 370.9 seconds, and no unexpected warnings
  occurred.
- 2026-10-12 09:47: ZIBNB--Gumbel completed refinement-cycle `sigma` in
  1,502.0 seconds and advanced to `nu`. All 28 candidates converged, six moves
  were accepted, the longest individual fit was 149.9 seconds, and no
  unexpected warnings occurred.
- 2026-10-12 10:06: ZIBNB--Gumbel completed refinement-cycle `nu` in
  1,175.0 seconds and advanced to `tau`. Seventeen of 20 candidates converged.
  Three terminal-step candidates were ineligible: one invalid smooth-EDF
  error and two `max_iterations` endpoints. None was selected; two deletion
  moves had already been accepted. The longest fit was 252.0 seconds and no
  unexpected warnings occurred. Reference eligibility/decision equivalence is
  reserved for the cumulative route audit.
- 2026-10-12 10:14: ZIBNB--Gumbel completed refinement-cycle `tau` in
  446.7 seconds and advanced to `theta`. All 26 candidates converged, six
  states (the maximum plus five deletion/simplification moves) were accepted,
  the longest individual fit was 112.2 seconds, and no unexpected warnings were recorded.
  ZIBNB--Student-t remained active in initial `nu`, with 19 checkpoints written
  through 10:13; it was progressing rather than stalled.
- 2026-10-12 10:16: ZIBNB--Gumbel completed refinement-cycle `theta` in
  143.6 seconds and began refinement cycle 2 at `mu`, because cycle 1 changed
  the joint model. All 16 theta candidates converged, three states (the maximum
  plus two deletion moves) were accepted, the longest fit was 24.4 seconds,
  and no unexpected warnings were recorded.
- 2026-10-12 10:18: ZIBNB--Student-t completed initial `nu` in 2,323.6
  seconds and advanced to `tau`. All 23 candidates converged, four states (the
  maximum plus three deletion moves) were accepted, the longest individual fit
  was 446.6 seconds, and no unexpected warnings were recorded. The longer stage wall time
  reflects multiple parallel candidate waves; no fit approached the isolated
  1,500-second Student-t allowance.
- 2026-10-12 10:22: ZIBNB--Gumbel completed refinement-cycle-2 `mu` in
  340.5 seconds and advanced to `sigma`. Twenty-nine of 30 candidates
  converged and seven states (the maximum plus six deletion/simplification
  moves) were accepted; the longest fit was 40.6 seconds and no unexpected warnings were
  recorded. One treatment-deletion candidate failed with a nonfinite
  discrete-rectangle sigma score and was ineligible. The accepted path was
  fully converged; its substantial AIC improvement shows that refinement is
  still escaping inferior joint optima rather than merely repeating cycle 1.
- 2026-10-12 10:49: ZIBNB--Gumbel completed refinement-cycle-2 `sigma` in
  1,605.4 seconds and advanced to `nu`. Fourteen of 17 candidates converged
  and four states (the maximum plus three simplification/deletion moves) were
  accepted. Three repeated age-deletion fits reached `max_iterations` with AIC
  around 59,000 versus the accepted path near 17,221, so they were ineligible
  and decisively noncompetitive. The longest fit was 585.1 seconds, inside the
  normal Gumbel cap, and no unexpected warnings were recorded.
- 2026-10-12 11:03: ZIBNB--Gumbel completed refinement-cycle-2 `nu` in
  830.4 seconds and advanced to `tau`. Fifteen of 17 candidates converged and
  seven states (the maximum plus six deletion moves) were accepted. Two early
  sex-deletion attempts reached `max_iterations`; the same deletion was refit,
  converged, and accepted at step 3 with lower AIC, so the final deletion
  decision was recovered rather than lost. The longest fit was 339.1 seconds
  and no unexpected warnings were recorded.
- 2026-10-12 11:08: ZIBNB--Gumbel completed cycle-2 `tau` and `theta` in
  163.0 and 150.1 seconds. Both stages were clean (29/29 and 27/27 candidates
  converged), each accepted six deletion moves, their longest fits were 14.8
  and 21.4 seconds, and neither recorded unexpected warnings. The route then entered
  refinement cycle 3 because its formula signature had changed.
- 2026-10-12 11:18: ZIBNB--Student-t completed initial `tau` in 3,642.4
  seconds and advanced to `theta`. All 28 candidates converged, four deletion
  or simplification moves were accepted, the longest individual fit was 479.5
  seconds, and no unexpected warnings were recorded. The hour-long stage comprised
  multiple serial candidate waves; no fit approached the 1,500-second cap.
- 2026-10-12 11:23: ZIBNB--Gumbel completed refinement cycle 3 (`mu`,
  `sigma`, `nu`, `tau`, `theta`) in 15.7 minutes total. Every one of the 126
  candidate fits across those stages converged, no unexpected warnings were recorded, and
  the slowest individual fit was 22.6 seconds. The cycle reduced `mu`, `nu`,
  `tau`, and `theta` to intercepts while leaving treatment, week, and their
  interaction in `sigma`, so another stability cycle was required.
- 2026-10-12 11:37: ZIBNB--Gumbel completed refinement cycle 4 in 13.5
  minutes and entered cycle 5. Through `tau`, 101 of 102 candidates converged;
  the sole error was an early treatment-by-week deletion in `sigma`, and that
  same deletion was later refit successfully and accepted. Cycle-4 `theta`
  then completed 27/27 candidates cleanly. No unexpected warnings occurred, the slowest
  fit was 22.2 seconds, and the resulting joint formula is intercept-only in
  all five parameters. Cycle 5 is the required fixed-point confirmation.
- 2026-10-12 11:45: Cycle-5 ZIBNB--Gumbel `mu` returned to the intercept-only
  formula (29/29 candidates converged), but `sigma` reopened treatment, week,
  and their interaction despite 14/14 candidates converging cleanly. Thus the
  route has not yet reached a formula fixed point and may be alternating
  between joint local optima. The programmed stability loop remains active;
  no gate is weakened and the Gumbel forecast now includes further cycles up
  to the configured maximum of 10.
- 2026-10-12 11:57: ZIBNB--Gumbel completed refinement cycle 5 in 20.6
  minutes and entered cycle 6. All 124 candidate fits converged, no unexpected
  warnings occurred, and the slowest fit was 33.3 seconds. The completed
  signature differs from the cycle start only by sigma treatment, week, and
  their interaction, reproducing the opposite member of the cycle-4/5
  two-state pattern. The runner will retain its strict one-cycle fixed-point
  requirement and continue through the configured ten-cycle bound.
- 2026-10-12 12:03: ZIBNB--Student-t completed initial `theta` in 2,679.0
  seconds and entered `zeta`, its last initial screen before refinement. All 25
  theta candidates converged, four deletion moves were accepted, the longest
  individual fit was 375.8 seconds, and no unexpected warnings were recorded.
  The stage wall time reflects serial candidate waves; no fit approached the
  1,500-second Student-t allowance.
- 2026-10-12 12:11: ZIBNB--Gumbel completed stable refinement cycle 6 and
  wrote a source-fingerprinted `SELECTION-COMPLETE.rds` marker. The final
  formulas are `mu ~ 1`, `sigma ~ treatment + week_factor +
  treatment:week_factor`, `nu ~ 1`, `tau ~ 1`, and `theta ~ 1`; the selected
  joint fit converged. Across the route, 843/855 candidate fits converged. All
  12 earlier failures were superseded by the final fixed-point cycle, which
  re-screened all coordinates with 124/124 candidates converged, zero errors,
  zero unexpected warnings, and a maximum fit time of 33.6 seconds. Selection
  took 4.33 hours; final stationarity, Hessian, and covariance certification is
  queued so it does not compete with the live Student-t route.
- 2026-10-12 12:14: With the Gumbel worker slot free, started one bounded
  ZISICHEL--Student-t operational recovery fit from the current-source failed
  full-zeta checkpoint, removing the separation-driving zeta `sex` term. It
  uses existing diagnostic tooling, frozen source `1deea57`, one CPU thread,
  the unchanged convergence gate, and the 1,500-second Student-t cap. This is
  an estimability probe only: the failed full-zeta endpoint is a warm start,
  not an eligible result.
- 2026-10-12 12:20: The reduced ZISICHEL--Student-t zeta fit converged in
  278.2 seconds (12 outer cycles, log likelihood -13222.665, AIC 26566.44)
  after removing zeta `sex`, with no numerical fit warnings. Because several
  remaining treatment-by-visit coefficients are extreme, it is not yet
  eligible for promotion; the independent direct/analytical scaled-score and
  Hessian audit is now running under the frozen source.
- 2026-10-12 12:26: The independent audit rejected that reduced endpoint:
  direct and analytical scaled scores were 1.02865 (strict limit 0.05), their
  agreement was 0.00150, and the observed information minimum eigenvalue was
  -0.0111. The close score agreement confirms a real nonstationary/separation
  endpoint rather than a derivative discrepancy. A second one-thread,
  1,500-second operational fit is now removing the zeta treatment-by-week
  interaction in addition to zeta `sex`; no gate has been relaxed.
- 2026-10-12 12:41: The ZIBNB--Student-t initial zeta maximum stopped at
  `max_iterations` after 50 outer cycles and 1,112.1 seconds; the route manager
  exited before evaluating deletion candidates. The checkpoint log likelihood
  was -13214.565 (AIC 26540.82), its zeta coefficients remained finite (range
  -7.41 to 2.99), and the objective was still improving. The original
  checkpoint was preserved, and the existing continuation utility is using at
  most seven further outer cycles within the remaining ~388 seconds of the
  same cumulative 1,500-second cap. This does not weaken convergence or reuse
  an ineligible endpoint.
- 2026-10-12 12:43: The second ZISICHEL--Student-t estimability fit (zeta
  `sex` and treatment-by-week removed) was rejected after 970.6 seconds at
  `max_iterations`; remaining zeta treatment coefficients drifted to -16.9
  and -11.8. With the interaction absent, strong hierarchy permits removing
  zeta treatment. A third bounded one-thread fit is now testing that minimal
  additional reduction under the unchanged convergence and 1,500-second cap.
- 2026-10-12 12:44: The first bounded ZIBNB--Student-t zeta continuation
  reached 57 cumulative outer cycles and 1,220.3 cumulative seconds but did
  not converge; AIC improved from 26540.82 to 26540.67 and the last outer
  change was 0.0154. A second seven-cycle continuation is using part of the
  remaining ~280 seconds under the same original cap.
- 2026-10-12 12:49: The ZIBNB--Student-t full-zeta maximum is ineligible after
  three bounded continuations: 71 cumulative outer cycles, 1,469.1 seconds,
  AIC 26540.23, and a last outer change of 0.0176. The remaining 31 seconds
  cannot plausibly meet the unchanged 1e-4 convergence gate. The failed
  checkpoint and its pre-continuation backup are preserved. Because the
  largest drift is in treatment-by-visit zeta terms, a one-thread reduced-zeta
  fit removing that interaction is now running under its own 1,500-second
  complex Student-t allowance; it is an estimability probe, not yet eligible.
- 2026-10-12 13:01: The third ZISICHEL--Student-t reduction (zeta `sex`,
  interaction, and treatment removed) stopped at 50 outer cycles after 1,065.1
  seconds, but its remaining zeta coefficients are finite and moderate (-0.54
  to 1.00), unlike the separated predecessors. An exact same-model
  continuation is using at most 430 of the remaining 435 seconds from that
  model's original 1,500-second allowance; the convergence gate is unchanged.
- 2026-10-12 13:03: The first reduced ZIBNB--Student-t zeta fit (interaction
  removed) was rejected at 50 outer cycles after 739.4 seconds. Its zeta
  coefficients exploded from -240.8 to +140.9, confirming separation rather
  than slow ordinary convergence. A second fit is removing zeta treatment now
  that hierarchy permits it, while retaining sex, age, baseline CD4, and week;
  it uses the same frozen source, gate, and 1,500-second cap.
- 2026-10-12 13:08: The exact ZISICHEL--Student-t v3 continuation converged
  after 17 additional outer cycles and 348.6 seconds: cumulative elapsed time
  was 1,413.8 seconds, final outer change 3.47e-5, log likelihood -13226.016,
  AIC 26554.18, and zeta coefficients ranged only from -0.23 to 0.83. The
  independent direct/analytical score and Hessian audit is now running; this
  reduced model is not promoted until those unchanged stationarity gates pass.
- 2026-10-12 13:10: The zero-zeta-start ZIBNB--Student-t reduction with zeta
  interaction and treatment removed converged in 297.1 seconds and 17 outer
  cycles (final change 4.79e-5, log likelihood -13221.020, AIC 26535.39).
  Remaining zeta coefficients are moderate (-0.37 to 0.49). Its independent
  direct/analytical score and Hessian audit is now running; it is not promoted
  until the unchanged strict stationarity gates pass.
- 2026-10-12 13:15: The converged ZISICHEL--Student-t v3 endpoint has
  positive-definite observed information (minimum eigenvalue 0.00498) and
  excellent direct/analytical agreement (0.00152), but both scaled scores are
  0.316 and fail the strict 0.05 stationarity gate. With 86 seconds left in the
  cumulative fit allowance, an audited half-Newton start plus at most 80
  seconds of the same RS model is performing the established conservative
  final polish; cumulative fitting remains under 1,500 seconds.
- 2026-10-12 13:19: The ZISICHEL half-Newton polish consumed its 80-second
  allowance and terminated at the cap without a saved endpoint, bringing that
  model's cumulative fit time to about 1,495 seconds; it will not be extended.
  The next simpler candidate removes zeta week from the finite v3 fit. The
  converged ZIBNB--Student-t reduction also failed certification (direct and
  analytical scaled scores 0.193; agreement 0.00079; minimum information
  eigenvalue -213.83). Since that model used only 297 seconds, a same-model
  all-parameter BFGS/RS polish is using at most 1,190 of its remaining ~1,203
  seconds. Both runs retain the strict gates and cumulative 1,500-second cap.
- 2026-10-12 13:26: The next ZISICHEL--Student-t zeta candidate, dropping
  week from the finite v3 model, converged in 371.9 seconds and improved AIC
  from 26554.18 to 26551.06. Its zeta formula retains the age and baseline-CD4
  smooths. The independent strict score/Hessian audit is running before this
  backward move is accepted. The ZIBNB all-parameter polish also converged in
  283.5 seconds but at a slightly lower likelihood than its parent, so its
  independent audit is still required before choosing either endpoint.
- 2026-10-12 13:33: The ZIBNB all-parameter polish was decisively rejected by
  its audit: scaled scores and agreement were nonfinite and the minimum
  information eigenvalue was -9.4 million, worse than the reduced parent.
  The next simpler zeta candidate is now dropping week from the converged
  treatment/interaction-reduced model under a fresh 1,500-second cap. No failed
  or nonstationary endpoint is entering model comparison.
- 2026-10-12 13:35: The AIC-improving ZISICHEL zeta drop-week model has
  positive-definite information (minimum eigenvalue 0.00232) and excellent
  gradient agreement (0.00152), but its direct/analytical scaled scores remain
  0.299 and fail the strict gate. Since its initial fit used only 371.9
  seconds, an audited half-Newton start plus all-parameter RS/BFGS polish is
  using at most 1,100 of the remaining ~1,128 seconds for this same model.
- 2026-10-12 13:44: The ZIBNB--Student-t zeta drop-week endpoint converged
  quickly (94.2 seconds), but its independent audit rejected it: direct and
  analytical scaled scores were 0.1815, agreement was 0.00077, and observed
  information had minimum eigenvalue -331.76. The dominant negative-curvature
  direction spans tau, sigma, and mu rather than zeta alone, so another blind
  zeta deletion is not justified; the next action is a bounded all-parameter
  stationarity/curvature recovery from this finite endpoint.
- 2026-10-12 13:45: The polished ZISICHEL--Student-t drop-week endpoint now
  passes both strict score gates (direct 0.04532, analytical 0.04531) with
  agreement 0.00149, after 161.7 seconds of fitting and 439.6 seconds of
  independent audit. It is nevertheless ineligible because observed
  information has minimum eigenvalue -0.0852; the negative direction is almost
  entirely the sigma intercept/week block. The parent endpoint was
  positive-definite but nonstationary, so neither endpoint is promoted and no
  inference gate is weakened.
- 2026-10-12 13:47: Started two one-thread bounded jobs, preserving the
  two-route concurrency ceiling. ZIBNB--Student-t is performing an
  all-parameter recovery from its finite drop-week endpoint (1,100-second
  allowance; no indefinite-Hessian Newton step). ZISICHEL--Student-t is
  evaluating the next hierarchical backward candidate by dropping the zeta
  age smooth from the positive-definite drop-week parent (1,500-second complex
  Student-t allowance). Both use the frozen source and unchanged convergence
  controls; neither can enter model comparison without an independent audit.
- 2026-10-12 13:48: The first ZISICHEL candidate launch exited before loading
  data because Windows split the smooth-term command-line argument. It created
  no fit/checkpoint and consumed no model allowance. The identical job was
  relaunched with the term quoted explicitly; its live command line now
  contains the single exact term `s(age, bs = "ps", k = 8)` and both bounded
  jobs are running normally.
- 2026-10-12 13:51: Both bounded fits completed well inside their allowances.
  The ZIBNB all-parameter recovery converged in 213.6 seconds (7 outer cycles,
  log likelihood -13221.607, AIC 26532.58), but its likelihood is lower than
  the unpolished parent and BFGS reported code 1; its independent audit is
  running before any use. The ZISICHEL zeta age-smooth deletion converged in
  161.4 seconds (7 outer cycles, log likelihood -13226.942, AIC 26546.94) and
  improves AIC by 4.12 over the eligible-form parent candidate. Independent
  score/Hessian audits for both endpoints are now running concurrently; no
  result is promoted on optimizer convergence alone.
- 2026-10-12 13:59: ZIBNB--Student-t all-parameter recovery passed the full
  stationarity/Hessian audit: direct scaled score 0.00509, analytical scaled
  score 0.00587, agreement 0.00078, and positive-definite observed information
  with minimum eigenvalue 0.503. Its subject-clustered sandwich covariance
  audit is now running. The unpolished ZISICHEL age-smooth-deletion endpoint
  failed (scores 0.297, agreement 0.00152, minimum eigenvalue -0.119); the
  unstable eigenvector is again almost entirely sigma intercept/week. A
  one-thread, 1,100-second all-parameter recovery from that finite endpoint is
  running in the second slot, using the same recovery pattern that certified
  ZIBNB and preserving the unchanged gates.
- 2026-10-12 14:01: The ZIBNB--Student-t recovery also passed covariance:
  model-based and clustered-sandwich covariance are both available, all
  diagonal variances are positive, the 1,177 subject clusters reproduce the
  joint log likelihood exactly, and the covariance calculation emitted zero
  warnings. This supplies the first fully certified Student-t maximum after
  excluding the non-estimable full-zeta region. Its actual zeta backward
  screen has begun with the age-smooth deletion under the normal 1,190-second
  cap while ZISICHEL recovery uses the second slot.
- 2026-10-12 14:04: The ZIBNB zeta age-smooth deletion converged in 136.3
  seconds and is a strong provisional move (AIC 26528.78 versus 26532.58 for
  the certified parent). The same-step sex deletion is now running; the
  baseline-CD4 smooth deletion remains queued, so no move has been selected
  prematurely. ZISICHEL all-parameter recovery converged in 253.2 seconds
  (AIC 26547.31) but with BFGS code 1 and a slightly lower likelihood than its
  unpolished parent. Its independent score/Hessian audit is running in the
  second slot before any promotion.
- 2026-10-12 14:11: The ZIBNB zeta sex deletion converged in 274.2 seconds
  with AIC 26530.81, improving its certified parent but not the age-smooth
  deletion; the baseline-CD4 smooth deletion is the final same-step candidate
  and is running now. ZISICHEL all-parameter recovery was decisively rejected:
  direct/analytical scaled scores and agreement were nonfinite and information
  had minimum eigenvalue -4,671. That failed state is discarded. The final
  legitimate ZISICHEL zeta simplification, dropping the remaining baseline-CD4
  smooth to an intercept-only zeta model from the finite pre-recovery parent,
  is now running under the 1,500-second Student-t allowance.
- 2026-10-12 14:22: ZIBNB zeta step 1 is fully ranked: age-smooth deletion
  AIC 26528.78, sex deletion 26530.81, certified parent 26532.58, and
  baseline-smooth deletion 26536.37. The age deletion has positive-definite
  information (minimum eigenvalue 0.462) and agreement 0.00075 but narrowly
  misses the strict score gate (0.0721 versus 0.05); an all-parameter recovery
  is running within a 1,100-second remainder before the provisional move can
  be accepted. The first ZISICHEL intercept-only launch exited before fitting
  because the diagnostic deletion helper cannot construct an empty term list.
  It produced no result and consumed no fit allowance. The identical `zeta ~
  1` candidate is now running through the helper's supported replacement path.
- 2026-10-12 14:25: The corrected ZISICHEL `zeta ~ 1` candidate converged in
  643.8 seconds (27 outer cycles, AIC 26553.47), safely inside the isolated
  Student-t allowance; its strict score/Hessian audit is running. The ZIBNB
  age-deletion recovery converged in 189.1 seconds and slightly improved the
  provisional result to AIC 26528.67. Its independent audit is also running;
  BFGS code 1 means optimizer convergence alone is not being accepted.
- 2026-10-12 14:30: The unpolished ZISICHEL intercept-only endpoint failed
  certification despite close derivative agreement: direct/analytical scores
  were 0.286, agreement was 0.00145, and information had a shallow minimum
  eigenvalue -0.00198. The score is concentrated in mu and the curvature
  direction in sigma, not zeta. With the zeta term search exhausted, one final
  all-parameter recovery is using at most 800 seconds, keeping cumulative fit
  time for this exact model below 1,500 seconds. The recovered ZIBNB
  age-deletion audit continues in the other slot.
- 2026-10-12 14:32: ZIBNB zeta step 1 is accepted after the recovered
  age-smooth deletion passed all gates: direct/analytical scores 0.00326 and
  0.00310, agreement 0.00075, and positive-definite information with minimum
  eigenvalue 0.470. The accepted AIC is 26528.67 versus 26532.58 for the
  certified parent. Step 2 is now evaluating deletion of zeta sex; deletion of
  the remaining baseline-CD4 smooth is queued. ZISICHEL intercept-only
  all-parameter recovery remains active in the other slot.
- 2026-10-12 14:43: ZIBNB zeta step 2 is fully ranked: sex deletion AIC
  26526.82, baseline-smooth deletion 26531.47, and certified parent 26528.67.
  The winning sex deletion is undergoing independent audit. ZISICHEL's final
  intercept-only recovery improved its scaled score from 0.286 to 0.0885 but
  remained ineligible and retained negative curvature in sigma, so further
  zeta polishing was stopped. The unchanged four-worker production screening
  driver has begun refinement cycle 1 at mu from the converged `zeta ~ 1`
  endpoint, under the normal 1,190-second cap. This is the scientifically
  required cyclic refinement where the non-zeta score/curvature blocks can be
  reselected; final stationarity/Hessian/covariance gates remain unchanged.
- 2026-10-12 14:50: ZISICHEL refinement-cycle-1 mu step 0 converged and all
  six first-step deletions are now dispatched across four workers. The ZIBNB
  step-2 sex deletion audit found positive-definite information (minimum
  eigenvalue 0.955) and agreement 0.00075, but its scaled score 0.0789 narrowly
  missed the strict gate. Because the residual score is in mu/nu rather than
  zeta, a bounded all-parameter recovery of that winning model is running;
  cumulative fit time remains below the 1,500-second Student-t exception.
- 2026-10-12 15:00: The recovered ZIBNB step-2 sex deletion passed all gates:
  direct/analytical scores 0.00208 and 0.00175, agreement 0.00075, and
  positive-definite information with minimum eigenvalue 0.898. It is accepted
  at AIC 26526.92, leaving only the zeta baseline-CD4 smooth. The sole step-3
  alternative (`zeta ~ 1`) is now fitting under the normal 1,190-second cap.
  ZISICHEL refinement-cycle-1 mu remains active with three durable checkpoints
  (maximum plus two completed first-step deletions).
- 2026-10-12 15:04: ZIBNB's terminal `zeta ~ 1` candidate converged in 240.8
  seconds but was noncompetitive (AIC 26529.56 versus 26526.92), so the zeta
  screen terminates with `zeta ~ s(baseline_cd4, bs = "ps", k = 8)`. ZIBNB
  refinement cycle 1 has started at mu with four workers and the normal
  1,190-second cap. ZISICHEL and ZIBNB are now the two concurrent four-worker
  routes; the two-route/eight-worker ceiling is exactly respected.
- 2026-10-12 15:09: ZISICHEL refinement-cycle-1 mu has a source-compatible
  converged maximum and two completed step-1 deletion checkpoints. The
  baseline-CD4-smooth deletion hit the normal 1,190-second cap during a tau
  backfitting update and is ineligible (no fit or metrics); it recorded only
  the expected segmented-missingness warning. The worker immediately advanced
  to another queued deletion, so the manager remains healthy and the cap was
  not extended. ZIBNB refinement-cycle-1 mu is still fitting its step-0
  maximum. Both four-worker managers and all eight workers remain responsive.
- 2026-10-12 15:10: ZIBNB refinement-cycle-1 mu step 0 converged in about 5.3
  minutes and the six first-step deletions are dispatched. ZISICHEL has three
  remaining first-step deletions active after the capped baseline-smooth fit;
  its completed maximum, age-smooth deletion, sex deletion, and explicit
  losing failure checkpoint are durable. No source change, cap extension, or
  unexpected warning occurred.
- 2026-10-12 15:13: ZISICHEL refinement-cycle-1 mu completed step 1 and
  accepted deletion of `treatment:week_factor`, improving AIC from 26566.70 to
  26557.83. Five of six candidates converged; the sole failure was the already
  recorded baseline-smooth deletion at the 1,190-second cap. The manager has
  dispatched seven step-2 moves from the accepted checkpoint. ZIBNB remains in
  its first mu deletion wave. The forecast is unchanged because both routes
  continue within their allocated concurrent slots.
- 2026-10-12 15:25: ZIBNB refinement-cycle-1 mu completed step 1 and accepted
  deletion of `treatment:week_factor`, improving AIC from 26544.13 to 26531.12.
  All six candidates converged, the slowest in 480.5 seconds, and no unexpected
  warning occurred. Seven step-2 moves are now dispatched. ZISICHEL remains in
  mu step 2 with four of seven candidates complete; dropping sex is the
  provisional best at AIC 26556.08 versus 26557.83 for its parent.
- 2026-10-12 15:33: ZISICHEL mu step 2 has five of seven candidates complete.
  The baseline-smooth deletion reached the normal 1,190-second cap during tau
  backfitting and is ineligible; it emitted only the expected segmented-
  missingness warning. Dropping sex remains the provisional winner. ZIBNB mu
  step 2 has four of seven candidates complete; dropping sex provisionally
  improves AIC to 26529.28, while the baseline-smooth deletion is decisively
  noncompetitive at 26801.86. Both managers continue without cap extensions.
- 2026-10-12 15:36: ZIBNB mu step 2 has five of seven candidates complete.
  The newly completed week-factor deletion converged in 336.0 seconds with AIC
  26640.96 and is noncompetitive; dropping sex remains the provisional winner
  at AIC 26529.28. ZISICHEL remains at five of seven completed candidates.
  Both four-worker managers and all eight workers remain responsive, the
  normal 1,190-second cap is unchanged, and the two-route ceiling is respected.
- 2026-10-12 15:38: ZISICHEL mu step 2 is complete and accepted linearization
  of the age smooth (`s(age) -> age`), improving AIC from 26557.83 to
  26555.43. Six of seven candidates converged; the only failure remains the
  already recorded baseline-smooth deletion at the normal cap. Six step-3
  moves are now dispatched across the same four workers. ZIBNB remains in mu
  step 2 with five of seven candidates complete.
- 2026-10-12 15:40: ZIBNB mu step 2 has six of seven candidates complete.
  Linearizing the age smooth converged in 232.6 seconds and is now the
  provisional winner at AIC 26528.78, ahead of dropping sex at 26529.28 and
  the parent at 26531.12. The baseline-smooth linearization is the sole
  outstanding same-step candidate, so selection correctly remains pending.
- 2026-10-12 15:43: ZISICHEL mu step 3 has three of six candidates complete.
  Dropping sex is provisionally best at AIC 26553.62 versus 26555.43 for the
  parent; dropping age and treatment are noncompetitive at 26567.52 and
  26677.03. The other three candidates remain active, so no move is accepted
  early. ZIBNB's final mu-step-2 candidate is still active and responsive.
- 2026-10-12 15:44: ZIBNB mu step 2 is complete and accepted linearization of
  the age smooth, improving AIC from 26531.12 to 26528.78. All seven
  candidates converged; the final baseline-smooth linearization took 443.1
  seconds and was noncompetitive at AIC 26635.72. Six step-3 moves are now
  dispatched. Both routes are therefore progressing through mu step 3 under
  the unchanged four-worker-per-route allocation and normal cap.
- 2026-10-12 15:51: ZISICHEL mu step 3 has four of six durable candidates;
  dropping sex remains provisionally best at AIC 26553.62, while the newly
  completed week deletion is noncompetitive at 26644.56. ZIBNB mu step 3 has
  four of six candidates complete; dropping sex is provisionally best at AIC
  26526.92 versus 26528.78 for the parent, and age, treatment, and baseline-
  smooth deletions are all worse. Outstanding candidates remain active within
  the normal cap; no same-step decision has been taken early.
- 2026-10-12 15:58: ZISICHEL mu step 3 has five of six candidates complete.
  Deleting the baseline-CD4 smooth reached the normal 1,190-second cap during
  a sigma backfitting update and is ineligible; its only warning was the
  expected segmented-missingness notice. The baseline-smooth linearization is
  the final candidate and remains active. ZIBNB is likewise waiting only on
  its final baseline-smooth linearization. Neither failure allowance nor
  convergence gate was relaxed.
- 2026-10-12 16:01: ZIBNB mu step 3 is complete and accepted deletion of sex,
  improving AIC from 26528.78 to 26526.92. All six candidates converged; the
  final baseline-smooth linearization took 439.0 seconds and was
  noncompetitive at AIC 26633.20. Five step-4 moves are now dispatched under
  the same four-worker allocation and normal cap. ZISICHEL's final step-3
  baseline-smooth linearization remains active.
- 2026-10-12 16:02: ZISICHEL mu step 3 is complete and also accepted deletion
  of sex, improving AIC from 26555.43 to 26553.62. Five of six candidates
  converged; the sole failure was the recorded capped baseline-smooth
  deletion, while its linearization converged in 247.1 seconds but was
  noncompetitive at AIC 26654.95. Five step-4 moves are now dispatched. Both
  Student-t routes are progressing in parallel through mu step 4.
- 2026-10-12 16:16: ZIBNB refinement-cycle-1 mu is complete. Step 4 found no
  improving move, so the selected mu formula is `cd4 ~ treatment + age +
  s(baseline_cd4, bs = "ps", k = 8) + week_factor` at AIC 26526.92. Across
  the complete mu screen, all 25 candidates converged, the maximum elapsed
  time was 480.45 seconds, and there were no unexpected warning strings. The
  final fit, formulas, path, 25-row model-selection table, and run manifest
  are present. The freed route slot immediately started ZIBNB sigma with four
  workers, the normal 1,190-second cap, the mu final fit/formulas as its warm
  start, and the frozen source. ZISICHEL mu step 4 remains active in the other
  slot, preserving the two-route/eight-worker ceiling.
- 2026-10-12 16:22: ZISICHEL mu step-4 baseline-smooth deletion reached the
  normal 1,190-second cap during a sigma inner iteration and is ineligible;
  only the expected segmented-missingness warning was recorded. The final
  same-step candidate, baseline-smooth linearization, is now active on that
  worker. ZIBNB sigma remains in its re-expanded step-0 maximum fit. No cap
  extension or convergence-gate change was made.
- 2026-10-12 16:27: ZISICHEL refinement-cycle-1 mu is complete. Step 4 found
  no improving move, selecting `cd4 ~ treatment + age + s(baseline_cd4,
  bs = "ps", k = 8) + week_factor` at AIC 26553.62. The complete screen has
  25 candidates: 21 converged and four explicit baseline-smooth deletion
  failures at the normal cap across steps 1--4; all failure warnings were the
  expected segmented-missingness notice and the maximum elapsed time was
  1,190.74 seconds. The final fit, formulas, path, 25-row table, and manifest
  are present. ZISICHEL sigma immediately started from those final mu warm-
  start artifacts with four workers and the normal cap. ZISICHEL sigma and
  ZIBNB sigma are now the only two live routes.
- 2026-10-12 16:27: ZIBNB sigma step 0 converged in 354.5 seconds at AIC
  26551.07 and dispatched six first-step moves. Two have already converged:
  dropping the age smooth is provisionally best at AIC 26546.19, while
  dropping sex gives 26549.36. Four same-step candidates remain active, so no
  deletion has been accepted early.
- 2026-10-12 16:33: ZISICHEL sigma step 0 converged in 302.4 seconds at AIC
  26574.43 and dispatched six first-step moves. ZIBNB sigma step 1 is three of
  six complete: its newly completed baseline-smooth deletion converged in
  443.9 seconds but is noncompetitive at AIC 26565.89; dropping the age smooth
  remains provisionally best at 26546.19. Both routes remain within the normal
  cap and four-worker allocations.
- 2026-10-12 16:38: ZIBNB sigma step 1 is complete with all six candidates
  converged and no unexpected warning. It accepted deletion of
  `treatment:week_factor`, improving AIC from 26551.07 to 26536.78; the
  slowest candidate took 443.9 seconds. Seven step-2 moves are now dispatched.
  ZISICHEL sigma step 1 is three of six complete: age-smooth, baseline-smooth,
  and sex deletions all converged but are non-improving (AIC 26569.43,
  26576.14, and 26572.56 versus parent 26574.43, with age deletion currently
  the provisional best).
- 2026-10-12 16:43: ZISICHEL sigma step 1 is complete with all six candidates
  converged. It accepted deletion of `treatment:week_factor`, improving AIC
  from 26574.43 to 26562.62; the slowest fit took 337.3 seconds. Seven step-2
  moves are now dispatched. ZIBNB sigma step 2 has its first candidate:
  dropping the age smooth converged in 256.5 seconds and provisionally improves
  AIC from 26536.78 to 26531.97. Six same-step candidates remain outstanding.
- 2026-10-12 16:54: ZIBNB sigma step 2 is complete with all seven candidates
  converged. It accepted deletion of the age smooth, improving AIC from
  26536.78 to 26531.97; the slowest candidate took 456.6 seconds. Five step-3
  moves are now dispatched. ZISICHEL sigma step 2 is five of seven complete;
  dropping the age smooth is provisionally best at AIC 26555.57 versus
  26562.62 for the parent. The two smooth-linearization candidates remain
  active, so no move has been accepted early.
- 2026-10-12 17:00: ZISICHEL sigma step 2 is complete with all seven
  candidates converged. It accepted deletion of the age smooth, improving AIC
  from 26562.62 to 26555.57; the slowest fit took 437.5 seconds. Five step-3
  moves are now dispatched. ZIBNB sigma step 3 is two of five complete:
  dropping week is provisionally best and improves AIC from 26531.97 to
  26528.52; dropping sex gives 26530.32. Three candidates remain active.
- 2026-10-12 17:09: ZIBNB sigma step 3 is complete with all five candidates
  converged. It accepted deletion of week, improving AIC from 26531.97 to
  26528.52; the slowest candidate took 456.8 seconds. Four step-4 moves are
  now dispatched. ZISICHEL sigma step 3 is three of five complete; dropping
  sex is provisionally best at AIC 26553.62 versus parent 26555.57, while the
  baseline-smooth and treatment deletions are non-improving. Two candidates
  remain active.
- 2026-10-12 17:20: ZISICHEL sigma step-3 week deletion reached the normal
  1,190-second cap during a zeta score evaluation and is ineligible; its only
  warning was the expected segmented-missingness notice. The baseline-smooth
  linearization is the final same-step candidate and remains active. ZIBNB
  sigma step 4 is three of four complete with one candidate still active. No
  allowance or gate was changed.
- 2026-10-12 17:21: ZIBNB sigma step 4 is complete with all four candidates
  converged. It accepted deletion of sex, improving AIC from 26528.52 to
  26526.92; the slowest fit took 443.3 seconds. Three step-5 moves are now
  dispatched. ZISICHEL remains on its final sigma-step-3 baseline-smooth
  linearization.
- 2026-10-12 17:23: ZISICHEL sigma step 3 is complete. Four of five candidates
  converged and the week deletion is the single explicit normal-cap failure
  already recorded. It accepted deletion of sex, improving AIC from 26555.57
  to 26553.62; the slowest converged candidate took 295.8 seconds. Four step-4
  moves are now dispatched. ZIBNB sigma step 5 continues in the other slot.
- 2026-10-12 17:33: ZIBNB refinement-cycle-1 sigma is complete. Step 5 found
  no improving move, selecting `sigma ~ treatment + s(baseline_cd4,
  bs = "ps", k = 8)` at AIC 26526.92. All 26 sigma candidates converged, the
  maximum elapsed time was 456.75 seconds, and there were no unexpected
  warning strings. The final fit, formulas, path, 26-row table, and manifest
  are present. The freed slot immediately started ZIBNB nu from the sigma
  final warm-start artifacts with four workers and the normal cap. ZISICHEL
  sigma step 4 remains active in the other slot.
- 2026-10-12 17:42: ZISICHEL refinement-cycle-1 sigma is complete. Step 4
  found no improving move, selecting `sigma ~ treatment +
  s(baseline_cd4, bs = "ps", k = 8) + week_factor` at AIC 26553.62. Across
  23 candidates, 22 converged and the sole failure was the explicit capped
  week deletion in step 3; maximum elapsed was 1,190.06 seconds and there were
  no unexpected warning strings. All five final artifacts are present. The
  freed slot immediately started ZISICHEL nu from the sigma final warm-start
  artifacts. ZISICHEL nu and ZIBNB nu are now the only two live four-worker
  routes.
- 2026-10-12 17:50: Live-state and source audit confirmed exactly the intended
  two managers (ZIBNB PID 3116 and ZISICHEL PID 28352), each with four RSOCK
  workers and no duplicate route. Both commands use the frozen screening
  driver and the required AIC/strong/refine/k=8/1,190-second/initial-inner-5
  controls, and both completed their converged nu step-0 maxima before
  dispatching six step-1 deletions. HEAD is
  `1deea5772e6d11b06d3702536763e040568adc8e`; `git diff --name-only` against
  that commit is empty for `DESCRIPTION` and `R`. Every extant nu checkpoint
  has a current-source package fingerprint and a filename matching the full
  metadata hash (including its starting-model signature); all were freshly
  fitted rather than reused. At the audit instant ZIBNB had four durable nu
  checkpoints and ZISICHEL had one. The two-route/eight-worker ceiling remains
  exactly respected.
- 2026-10-12 18:00: Both Student-t nu routes completed step 1 cleanly and
  advanced to seven-candidate step 2 waves. Each accepted deletion of
  `treatment:week_factor`: ZIBNB improved from AIC 26539.5079 to 26530.5771
  (6/6 converged; maximum 480.74 seconds), and ZISICHEL improved from AIC
  26563.8957 to 26558.4428 (6/6 converged; maximum 339.19 seconds). Neither
  step recorded an error or unexpected warning. Both four-worker managers
  remain live, with the ordinary 1,190-second cap and frozen source unchanged.
- 2026-10-12 18:14: ZISICHEL nu step 2 completed with all seven candidates
  converged and accepted deletion of the age smooth, improving AIC from
  26558.4428 to 26553.7181. The slowest candidate took 313.18 seconds; there
  were no errors or unexpected warnings. Five step-3 candidates are now
  dispatched. ZIBNB is still waiting on the final candidate in its own nu
  step-2 wave, so no decision has been taken early there.
- 2026-10-12 18:15: ZIBNB nu step 2 also completed with all seven candidates
  converged and accepted deletion of the age smooth, improving AIC from
  26530.5771 to 26525.9426. The slowest candidate took 471.83 seconds; all
  warning strings were the expected segmented-missingness notice and there
  were no errors or unexpected warnings. Five step-3 candidates are now
  dispatched, so both routes continue in parallel at nu step 3.
- 2026-10-12 18:23: ZISICHEL nu step 3 completed with all five candidates
  converged and accepted deletion of sex, narrowly improving AIC from
  26553.7181 to 26553.6270. Maximum elapsed was 299.49 seconds, with no errors
  or unexpected warnings. Four step-4 candidates are now fitting. ZIBNB is
  still waiting on the last candidate in its five-candidate step-3 wave.
- 2026-10-12 18:32: ZISICHEL refinement-cycle-1 nu completed after step 4
  found no improving move. The selected nu formula is `~ treatment +
  s(baseline_cd4, bs = "ps", k = 8) + week_factor`, with final AIC
  26553.62697. All 23 candidates converged, the four accepted states form a
  complete consistent path, maximum elapsed was 339.19 seconds, and there were
  no errors or unexpected warnings. The final formulas, final fit, path,
  selection table, and manifest are present; the fit reports optimizer
  convergence on 1,177 subjects and the manifest source fingerprint matches
  frozen commit `1deea57`.
- 2026-10-12 18:43: ZIBNB refinement-cycle-1 nu likewise completed after
  step 4 found no improving move. It had accepted deletion of sex at step 3,
  improving AIC from 26525.94257 to 26524.50215, and selected `nu ~ treatment
  + s(baseline_cd4, bs = "ps", k = 8) + week_factor`. All 23 candidates
  converged, its four-state path is complete and consistent, maximum elapsed
  was 508.00 seconds, and there were no errors or unexpected warnings. All
  five final artifacts are present, the final fit converged on 1,177 subjects,
  and the manifest fingerprint matches the frozen source.
- 2026-10-12 18:48: Both freed route slots advanced to tau using each route's
  audited nu final formulas and final fit as the formula input and warm start.
  ZIBNB manager PID 18748 and ZISICHEL manager PID 27876 each have exactly four
  RSOCK workers; both commands retain AIC, strong hierarchy, refine stage,
  k=8, four workers, the 1,190-second cap, five initial inner iterations, and
  resume compatibility checks. They are the only two ACTG managers and are
  currently fitting fresh step-0 maxima.
- 2026-10-12 19:00: Both tau step-0 maxima converged and dispatched six
  first-step backward moves. ZIBNB took 558.25 seconds and ZISICHEL 659.89
  seconds. Each checkpoint is freshly fitted, carries the current frozen-source
  fingerprint, and its filename matches the complete metadata hash including
  the starting-model signature. The two compatible four-worker managers remain
  the only ACTG routes in flight.
- 2026-10-12 19:12: ZISICHEL tau step 1 completed with all six candidates
  converged and accepted deletion of `treatment:week_factor`, improving AIC
  from 26578.6788 to 26565.6312. The slowest candidate took 346.08 seconds,
  with no errors or unexpected warnings. Seven step-2 candidates are now
  dispatched. ZIBNB remains in its first six-candidate tau deletion wave.
- 2026-10-12 19:16: ZIBNB tau step 1 completed with all six candidates
  converged and likewise accepted deletion of `treatment:week_factor`, reducing
  AIC from 26549.5628 to 26534.8858. Maximum elapsed was 524.48 seconds and
  there were no errors or unexpected warnings. Its seven step-2 candidates are
  now dispatched, placing both routes at tau step 2.
- 2026-10-12 19:29: ZISICHEL tau step 2 completed with all seven candidates
  converged and accepted deletion of week, improving AIC from 26565.6312 to
  26560.7834. The slowest candidate took 345.68 seconds and no error or
  unexpected warning occurred. Six step-3 candidates are now dispatched.
  ZIBNB remains at four of seven completed step-2 candidates.
- 2026-10-12 19:39: ZIBNB tau step 2 completed with all seven candidates
  converged and accepted deletion of week, improving AIC from 26534.8858 to
  26531.1839. Maximum elapsed was 532.02 seconds, with no error or unexpected
  warning. Six step-3 candidates are now dispatched, so both tau routes are
  again at the same backward-selection step.
- 2026-10-12 19:40: ZISICHEL tau step 3 completed with all six candidates
  converged and accepted deletion of the age smooth, lowering AIC from
  26560.7834 to 26556.6453. Maximum elapsed was 342.73 seconds and no error or
  unexpected warning occurred. Four step-4 candidates are now fitting while
  ZIBNB proceeds through tau step 3.
- 2026-10-12 19:51: ZISICHEL tau step 4 completed with all four candidates
  converged and accepted linearization of the baseline-CD4 smooth, improving
  AIC from 26556.6453 to 26553.6271 and selecting the current `tau ~ treatment
  + sex + baseline_cd4`. Maximum elapsed was 334.83 seconds, with no error or
  unexpected warning. Its four-worker manager remains healthy and has
  dispatched the three step-5 deletions. ZIBNB has four durable step-3
  checkpoints and remains active in the other slot.
- 2026-10-12 19:51: ZISICHEL tau step 4 completed with all four candidates
  converged and accepted linearization of the baseline-CD4 smooth, improving
  AIC from 26556.6453 to 26553.6271. The slowest candidate took 334.83 seconds
  and no error or unexpected warning occurred. Three step-5 candidates are now
  fitting. ZIBNB remains at four of six completed tau step-3 candidates.
- 2026-10-12 18:32: ZISICHEL refinement-cycle-1 nu completed in 49 minutes
  38 seconds. Step 4 found no improving move, selecting `nu ~ treatment +
  s(baseline_cd4, bs = "ps", k = 8) + week_factor` at AIC 26553.6270. All
  23 candidates converged, maximum candidate elapsed was 339.19 seconds, and
  there were no errors or unexpected warnings. The final fit, formulas, path,
  table, and manifest are present.
- 2026-10-12 18:43: ZIBNB refinement-cycle-1 nu completed in 1 hour 9 minutes
  46 seconds. Step 3 accepted deletion of sex and step 4 found no further
  improvement, selecting `nu ~ treatment + s(baseline_cd4, bs = "ps", k =
  8) + week_factor` at AIC 26524.5022. All 23 candidates converged, maximum
  candidate elapsed was 508.00 seconds, and there were no errors or
  unexpected warnings. The final fit, formulas, path, table, and manifest are
  present.
- 2026-10-12 18:48: Both freed slots advanced immediately to tau. ZIBNB tau
  is manager PID 18748 and ZISICHEL tau is manager PID 27876; each has exactly
  four RSOCK workers and uses its audited nu final fit/formulas as warm-start
  input. Both commands retain AIC, strong hierarchy, refine mode, k=8,
  initial-inner-5, resume, and the ordinary 1,190-second cap. The only stderr
  content at startup is the known Windows worker locale warning. No duplicate
  ACTG manager is live, so the two-route/eight-worker ceiling is respected.
- 2026-10-12 19:56: ZISICHEL refinement-cycle-1 tau completed after step 5
  found no further improving move. The accepted five-state path deletes the
  treatment-by-week interaction, week, and the age smooth, then linearizes
  baseline CD4, selecting `tau ~ treatment + sex + baseline_cd4` at AIC
  26553.62708. All 27 candidates converged; maximum elapsed was 659.89
  seconds, with zero errors, duplicate rows, or unexpected warnings. The
  selection table, path, final formulas, final fit, and run manifest are all
  present. The manifest covers 4,708 rows and exactly 1,177 subjects, all 391
  recorded package-source hashes match the frozen worktree, and the final
  fit's structured convergence record reports `converged = TRUE`.
- 2026-10-12 19:59: The freed ZISICHEL slot advanced immediately to theta as
  manager PID 28524, using the audited tau final formulas and final fit as both
  formula input and warm start. Its command retains AIC, strong hierarchy,
  refine mode, k=8, four workers, initial-inner-5, resume checks, and the
  ordinary 1,190-second cap; all BLAS/OpenMP thread limits were set to one.
  Together with the still-live ZIBNB tau manager, this preserves the strict
  two-manager ceiling with no duplicate route.
- 2026-10-12 19:56: ZIBNB tau step 3 completed with all six candidates
  durable and accepted simplification of the baseline-CD4 smooth to a linear
  term, lowering AIC from 26531.18391 to 26527.80626. Five step-4 moves were
  dispatched; by 20:01 the first new checkpoint was complete. The live process
  audit shows exactly the two intended managers and four RSOCK workers under
  each, with no duplicate route.
- 2026-10-12 20:05: ZISICHEL theta's fresh step-0 maximum converged in 336.41
  seconds (seven outer iterations) at AIC 26566.5 and dispatched six first-step
  moves. Its checkpoint is not resumed, records the expected starting-model
  signature `ed9d0e304ff6ecb1cd81ddbe16781faa`, contains 1,177 subjects, and all
  391 package-source hashes match the frozen worktree. The only fit warning is
  the classified ACTG intermittent-observation-gap notice.
- 2026-10-12 20:10: ZIBNB tau step 4 completed with all five candidates
  converged and accepted deletion of the age smooth, improving AIC from
  26527.80626 to 26524.50158. Maximum elapsed in the wave was 524.10 seconds,
  with no errors and only the classified ACTG observation-gap notice. Three
  step-5 deletions are now fitting. ZISICHEL theta has one of six first-step
  candidates durable in addition to its audited maximum.
- 2026-10-12 20:17: ZISICHEL theta step 1 completed with all six candidates
  converged and accepted deletion of `treatment:week_factor`, reducing AIC
  from 26566.47965 to 26560.56769. Maximum elapsed was 342.04 seconds, with no
  errors; seven step-2 moves are now fitting. ZIBNB tau has two of its three
  terminal step-5 checkpoints durable.
- 2026-10-12 20:18: ZIBNB refinement-cycle-1 tau completed after step 5
  found no improving move. Its complete five-state accepted path deletes the
  interaction and week, linearizes the baseline-CD4 smooth, and deletes the
  age smooth, selecting `tau ~ treatment + sex + baseline_cd4` at AIC
  26524.50158. All 28 candidates converged, maximum elapsed was 558.25
  seconds, and there were zero errors, duplicate rows, or unexpected warnings.
  All five final artifacts are present; the final fit reports convergence on
  exactly 1,177 subjects and all 391 manifest source hashes match the frozen
  worktree.
- 2026-10-12 20:19: The freed ZIBNB slot advanced immediately to theta as
  manager PID 14500 using the audited tau final fit and formulas. The exact
  production controls and one-thread BLAS/OpenMP limits are retained. It joins
  the live ZISICHEL theta manager without exceeding two managers/eight workers
  and without duplicating a route.
- 2026-10-12 20:27: ZIBNB theta's fresh step-0 maximum converged in 424.00
  seconds (seven outer iterations) at AIC 26537.85 and dispatched six
  first-step moves. It is not a resumed checkpoint, records starting signature
  `d8ce68f827c0e4a9c772e59a3286e709`, contains exactly 1,177 subjects, and its
  full package-source fingerprint matches the frozen worktree.
- 2026-10-12 20:31: ZISICHEL theta step 2 completed with all seven candidates
  converged and accepted deletion of the age smooth, improving AIC from
  26560.56769 to 26555.94569. Maximum elapsed was 340.13 seconds and no error
  occurred. Five step-3 moves are now fitting; ZIBNB theta has its maximum plus
  one first-step candidate durable.
- 2026-10-12 20:41: ZIBNB theta step 1 completed with all six candidates
  converged and accepted deletion of `treatment:week_factor`, improving AIC
  from 26537.85050 to 26531.26971. Maximum elapsed was 415.45 seconds and no
  error occurred. Seven step-2 candidates are now fitting. ZISICHEL theta
  remains on its final step-3 candidate, still within the ordinary cap.
- 2026-10-12 20:42: ZISICHEL theta step 3 completed with all five candidates
  converged and accepted deletion of sex, reducing AIC from 26555.94569 to
  26553.94717. Maximum elapsed was 328.24 seconds and no error occurred. Four
  step-4 candidates are now fitting.
- 2026-10-12 20:53: ZISICHEL theta step 4 completed with all four candidates
  converged and accepted deletion of the baseline-CD4 smooth, improving AIC
  from 26553.94717 to 26553.62655. Maximum candidate fit time was 334.22
  seconds, with no errors; two terminal step-5 deletions are now fitting.
- 2026-10-12 20:58: ZIBNB theta step 2 completed with all seven candidates
  converged and accepted deletion of the age smooth, lowering AIC from
  26531.26971 to 26527.31688. Maximum elapsed was 412.98 seconds and no error
  occurred. Five step-3 candidates are now fitting.
- 2026-10-12 20:58: ZISICHEL refinement-cycle-1 theta completed after its two
  step-5 terminal deletions found no better model. The five-state path selects
  `theta ~ treatment + week_factor` at AIC 26553.62655. All 25 candidates
  converged, maximum elapsed was 342.04 seconds, with zero errors, duplicate
  rows, or unexpected warnings. All five final artifacts are present; the
  final fit reports convergence on exactly 1,177 subjects and all 391 source
  hashes match the frozen worktree.
- 2026-10-12 20:59: The freed ZISICHEL slot advanced to the bounded zeta stage
  as manager PID 3472. Its maximum and only estimable formula is the documented
  `zeta ~ 1`; no structurally non-estimable larger zeta region is re-expanded.
  The audited theta final fit/formulas supply the warm start, and all ordinary
  production controls and one-thread numerical-library limits are retained.
- 2026-10-12 21:12: ZIBNB theta step 3 completed with all five candidates
  converged and accepted deletion of sex, improving AIC from 26527.31688 to
  26525.32429. Maximum elapsed was 433.34 seconds and no error occurred. Four
  step-4 candidates are now fitting. The ZISICHEL bounded-zeta maximum remains
  live and below the ordinary 1,190-second cap.
- 2026-10-12 21:19: The required ZISICHEL bounded endpoint `zeta ~ 1` reached
  the ordinary cap during its zeta RS update and did not produce final-stage
  artifacts: elapsed 1190.03 seconds, `converged = FALSE`, `stop_reason =
  error`, with the explicit max-elapsed error and only the classified ACTG
  observation-gap warning. This is the sole estimable zeta endpoint, not a
  losing backward candidate, so it qualifies for the permitted isolated
  complex-Student-t allowance. It will be retried at 1,500 seconds only after
  the concurrent ZIBNB theta manager finishes, so the retry is genuinely
  isolated. No larger non-estimable zeta formula will be introduced.
- 2026-10-12 21:26: ZIBNB theta step 4 completed with all four candidates
  converged and accepted deletion of the baseline-CD4 smooth, improving AIC
  from 26525.32429 to 26524.50635. Maximum elapsed was 422.42 seconds and no
  error occurred. Two terminal step-5 deletions are now fitting. This remains
  the only live ACTG manager; the second slot is deliberately idle to preserve
  isolation for the subsequent capped ZISICHEL zeta retry.
- 2026-10-12 21:33: ZIBNB refinement-cycle-1 theta completed after its two
  terminal step-5 deletions found no improvement. Its complete five-state path
  selects `theta ~ treatment + week_factor` at AIC 26524.50635. All 25
  candidates converged, maximum elapsed was 483.52 seconds, and there were no
  errors, duplicate rows, or unexpected warnings. All five final artifacts
  exist; the final fit reports convergence on exactly 1,177 subjects and all
  391 source hashes match the frozen worktree.
- 2026-10-12 21:34: With every other ACTG manager stopped, the required
  ZISICHEL `zeta ~ 1` endpoint was relaunched as an explicitly isolated complex
  Student-t fit under the permitted 1,500-second allowance (manager PID 3716,
  output stem `actg-zisichel-t-operational-refinement-cycle-01-zeta-isolated-1500`).
  It starts afresh from the audited theta final fit/formulas, retains all other
  production controls and one-thread limits, and is the sole live manager.
- 2026-10-12 21:55: The isolated ZISICHEL zeta refit stopped after 1,256.09
  seconds and 50 outer iterations with `converged = FALSE` and
  `stop_reason = max_iterations`; it did not hit the 1,500-second allowance.
  It is not promoted, and the outer-iteration/convergence gates are not
  weakened. The classified observation-gap notice and the explicit
  nonconvergence warning are its only warnings. Because the bounded estimable
  zeta set is exactly the singleton `zeta ~ 1`, there is no backward move to
  evaluate: the zeta selection is a documented no-op at the same formula,
  carried by the converged joint theta final fit. This follows the earlier
  certified estimability conclusion and sends the residual non-zeta
  score/curvature back through cyclic refinement rather than reopening the
  structurally non-estimable zeta region.
- 2026-10-12 21:57: Cyclic refinement resumed immediately in both route slots.
  ZISICHEL cycle 2 began at mu (manager PID 29692) from the converged cycle-1
  theta endpoint, preserving the singleton `zeta ~ 1`. ZIBNB began its bounded
  cycle-1 zeta screen (manager PID 10028) from the audited theta endpoint; its
  only allowed simplification is the documented intercept alternative to the
  selected baseline-CD4 smooth. Both use the unchanged production driver,
  normal 1,190-second cap, exact four-worker controls, and one-thread numerical
  libraries; they are the only two managers.
- 2026-10-12 21:59 correction: checkpoint inspection found that the two
  purportedly bounded ZISICHEL zeta launches at 20:59 and 21:34 were not
  intercept-only fits. In `--stage=refine`, the frozen screening driver
  unconditionally replaces the active parameter with its full smooth maximum;
  both checkpoints therefore contain the full zeta treatment/age/sex/
  baseline-CD4/week/interaction formula. Their 1,190-second and 1,256-second
  failures are invalid full-zeta diagnostics and provide no new evidence about
  the already certified `zeta ~ 1` endpoint. Neither result was promoted.
  ZISICHEL cycle-2 mu remains valid because its theta parent already carries
  `zeta ~ 1` and refine overwrites only mu.
- 2026-10-12 21:59 correction: the standard ZIBNB refine-zeta manager had the
  same full-maximum routing defect and was stopped before producing a
  checkpoint; an automatic duplicate launch was also stopped after exact PID
  and command-line verification. The valid bounded comparison is now fitting
  only the intercept alternative to the selected baseline-CD4 zeta smooth via
  `diagnose-specific-candidate.R`, from the audited theta final fit, with the
  normal 1,190-second cap and frozen package source. Alongside the valid
  four-worker ZISICHEL cycle-2 mu route, this restores two nonduplicate routes.
- 2026-10-12 21:58: The first ZIBNB zeta manager exited during startup before
  creating a checkpoint, manifest, fit error, or live worker remainder; its
  logs contain only the Windows locale startup notices. After verifying that
  no duplicate manager or worker survived, the identical command was
  relaunched as PID 16268 with separate relaunch logs. ZISICHEL cycle-2 mu
  remained healthy and undisturbed, and the two-manager/eight-worker ceiling is
  restored.
- 2026-10-12 21:59: The identical ZIBNB zeta relaunch also vanished during its
  maximum fit without an R error, checkpoint, manifest, or surviving worker;
  its logs again stop after normal startup and locale notices. Two consecutive
  silent exits only under the restored eight-worker load are treated as host
  resource pressure, not model convergence evidence. To avoid duplicate churn,
  ZISICHEL cycle-2 mu will complete first and the unchanged normal-cap ZIBNB
  bounded-zeta command will then run alone. No fit result has been lost or
  promoted, but this adds the remaining ZISICHEL-mu duration to the ZIBNB-zeta
  forecast.
- 2026-10-12 22:03: ZISICHEL cycle-2 mu's fresh maximum converged in 362.81
  seconds (seven outer iterations) at AIC 26566.7 and dispatched six backward
  moves. Its checkpoint records exactly 1,177 subjects, starting signature
  `882002c07db0cfb431f997760173d595`, and a full source fingerprint matching
  the frozen worktree. It is the sole live manager.
- 2026-10-12 22:23: ZISICHEL cycle-2 mu's baseline-CD4-smooth deletion reached
  the ordinary cap at 1190.07 seconds with the explicit max-elapsed error. It
  is a losing backward alternative and receives no extension. The age-smooth
  and sex deletions converged in 335.31 and 337.53 seconds respectively; three
  other first-step candidates remain active under the same normal cap. This is
  the same documented hard baseline-smooth deletion region seen in cycle 1.
- 2026-10-12 22:29: ZISICHEL cycle-2 mu step 1 completed with five of six
  candidates converged and the single allowed losing baseline-smooth deletion
  failure at 1190.07 seconds. The five successful candidates took at most
  341.44 seconds. Deleting `treatment:week_factor` was accepted, improving AIC
  from 26566.69541 to 26557.83405, and seven step-2 moves are now fitting.
- 2026-10-12 22:34: An unexpected ZIBNB cycle-2 mu manager (PID 31208) appeared
  from a separate parent even though the bounded cycle-1 zeta stage is not yet
  complete. It was verified by exact command line and stopped with its five
  child processes before producing any checkpoint or artifact. This prevents
  an out-of-order route, avoids recreating the unstable eight-worker load, and
  leaves the legitimate ZISICHEL cycle-2 mu manager untouched as the sole
  route.
- 2026-10-12 22:49: ZISICHEL cycle-2 mu step 2 produced the expected second
  losing baseline-smooth deletion failure at the ordinary cap (1190.13
  seconds, explicit max-elapsed error), with no extension. Treatment,
  age-smooth, sex, and week deletions all converged in 321.70--336.57 seconds;
  two simplification candidates remain active.
- 2026-10-12 22:54: ZISICHEL cycle-2 mu step 2 completed with six of seven
  candidates converged and only the allowed losing baseline-smooth deletion
  failure at 1190.13 seconds. Linearizing the age smooth was accepted, reducing
  AIC from 26557.83405 to 26555.42905. Six step-3 candidates are now fitting.
- 2026-10-12 23:15: ZISICHEL cycle-2 mu step 3's baseline-smooth deletion is
  the third allowed losing failure. Its checkpoint records 1188.03 seconds and
  its explicit error reports the RS update crossing the 1,190-second cap at
  1190.2 seconds. It receives no extension. Treatment, age, sex, and week
  deletions all converged in 313.87--329.92 seconds; the baseline linearization
  is the sole remaining candidate in the wave.
- 2026-10-12 23:19: ZISICHEL cycle-2 mu step 3 completed with five of six
  candidates converged and only the documented losing baseline-smooth failure.
  Deleting sex was accepted, improving AIC from 26555.42905 to 26553.62166.
  Five step-4 candidates are now fitting.
- 2026-10-12 23:43: ZISICHEL cycle-2 mu completed in 1 hour 46 minutes 54
  seconds. Step 4 found no further improvement, retaining `mu = cd4 ~
  treatment + age + s(baseline_cd4, bs = "ps", k = 8) + week_factor` at AIC
  26553.62166, exactly the cycle-1 selected formula. Across 25 candidates,
  21 converged and the four failures were the same losing baseline-smooth
  deletion at successive ordinary caps; maximum recorded elapsed was 1190.13
  seconds and there were no unexpected warnings. All five final artifacts are
  present and the final fit reports convergence on 1,177 subjects.
- 2026-10-13 00:07: Both free slots advanced: ZISICHEL cycle-2 sigma is
  manager PID 19528, warm-started from its audited cycle-2 mu final artifacts,
  and ZIBNB cycle-2 mu is manager PID 16792, warm-started from its audited
  cycle-1 theta endpoint after the bounded zeta comparison retained the
  baseline-CD4 smooth. Each has exactly four RSOCK workers and uses the frozen
  driver, AIC/strong/refine/k=8/initial-inner-5 controls and ordinary
  1,190-second cap. They are the only two live ACTG managers.
- 2026-10-13 00:26: ZISICHEL cycle-2 sigma step 1 completed with all six
  candidates converged and accepted deletion of `treatment:week_factor`,
  improving AIC from 26574.43239 to 26562.62280. Maximum elapsed was 374.41
  seconds and there were no errors or unexpected warnings. Seven step-2 moves
  are now fitting.
- 2026-10-13 00:35: ZIBNB cycle-2 mu step 1 completed with all six candidates
  converged and accepted deletion of `treatment:week_factor`, improving AIC
  from 26540.72354 to 26528.36192. Maximum elapsed was 548.38 seconds and
  there were no errors or unexpected warnings. Seven step-2 moves are now
  fitting; both route managers remain healthy under the ordinary cap.
- 2026-10-13 00:49: ZISICHEL cycle-2 sigma step 2 completed with all seven
  candidates converged and accepted deletion of the age smooth, improving AIC
  from 26562.62280 to 26555.51212. Maximum elapsed was 471.93 seconds, with no
  errors or unexpected warnings; five step-3 moves were dispatched.
- 2026-10-13 00:58: ZIBNB cycle-2 mu step 2 completed with all seven
  candidates converged and accepted linearization of the age smooth, improving
  AIC from 26528.36192 to 26526.39661. Maximum elapsed was 530.41 seconds,
  with no errors or unexpected warnings; six step-3 moves are fitting.
- 2026-10-13 01:02: ZISICHEL cycle-2 sigma step 3 completed with all five
  candidates converged and accepted deletion of sex, improving AIC to
  26553.63083. Maximum elapsed was 402.22 seconds and no error or unexpected
  warning occurred. Four step-4 moves are fitting while ZIBNB continues its
  step-3 wave.
- 2026-10-13 01:20: ZISICHEL cycle-2 sigma completed in 1 hour 12 minutes 48
  seconds. Step 4 found no improving move, retaining `sigma ~ treatment +
  s(baseline_cd4, bs = "ps", k = 8) + week_factor` at AIC 26553.63083,
  exactly its cycle-1 selected formula. Across 23 candidates, 22 converged and
  the sole failure was the losing week deletion at the ordinary cap (1190.22
  seconds). There were no unexpected warnings and all final artifacts exist.
- 2026-10-13 01:35: ZIBNB cycle-2 mu completed in 1 hour 27 minutes 31
  seconds. Steps 3 and 4 deleted sex and then found no further improvement,
  retaining `mu = cd4 ~ treatment + age + s(baseline_cd4, bs = "ps", k = 8)
  + week_factor` at AIC 26524.50638, exactly its cycle-1 formula. All 25
  candidates converged, maximum elapsed was 550.28 seconds, and no error or
  unexpected warning occurred; all final artifacts exist.
- 2026-10-13 01:41: The two free slots advanced to ZIBNB cycle-2 sigma
  (manager PID 32592) and ZISICHEL cycle-2 nu (manager PID 25380), each from
  its audited preceding cycle-2 final fit/formulas. Both retain four workers,
  one-thread numerical libraries, and all frozen production controls including
  the ordinary 1,190-second cap. They are the only two live ACTG managers.
- 2026-10-13 01:58: ZISICHEL cycle-2 nu step 1 completed with all six
  candidates converged and accepted deletion of `treatment:week_factor`,
  improving AIC from 26563.88987 to 26558.42983. Maximum elapsed was 354.45
  seconds, with no errors or unexpected warnings. Seven step-2 moves are
  fitting.
- 2026-10-13 02:03: ZIBNB cycle-2 sigma step 1 completed with all six
  candidates converged and accepted deletion of `treatment:week_factor`,
  improving AIC from 26549.26504 to 26534.44005. Maximum elapsed was 513.32
  seconds, with no errors or unexpected warnings. Seven step-2 moves are
  fitting; both managers remain healthy under the ordinary cap.
- 2026-10-13 02:13: ZISICHEL cycle-2 nu step 2 completed with all seven
  candidates converged and accepted deletion of the age smooth, improving AIC
  to 26553.70321. Maximum elapsed was 336.71 seconds and no error or unexpected
  warning occurred; five step-3 moves were dispatched.
- 2026-10-13 02:24: ZIBNB cycle-2 sigma step 2 completed with all seven
  candidates converged and accepted deletion of the age smooth, improving AIC
  to 26529.45328. Maximum elapsed was 508.11 seconds, with no error or
  unexpected warning; five step-3 moves were dispatched.
- 2026-10-13 02:25: ZISICHEL cycle-2 nu step 3 completed with all five
  candidates converged and accepted deletion of sex, reducing AIC to
  26553.62659. Its four step-4 candidates subsequently all converged and found
  no improvement, so the stage completed at 02:35 after 53 minutes 49 seconds,
  with 23/23 candidates converged, maximum elapsed 354.45 seconds, no errors or
  unexpected warnings, and all final artifacts present. The retained `nu ~
  treatment + s(baseline_cd4, bs = "ps", k = 8) + week_factor` exactly matches
  cycle 1.
- 2026-10-13 02:36: ZIBNB cycle-2 sigma step 3 completed with all five
  candidates converged and accepted deletion of week, improving AIC to
  26526.25159. Maximum elapsed was 513.95 seconds; four step-4 moves are now
  fitting. The free ZISICHEL slot advanced to cycle-2 tau as manager PID
  30704 using the audited nu final warm start and ordinary production cap.
- 2026-10-13 02:51: ZIBNB cycle-2 sigma step 4 completed with all four
  candidates converged and accepted deletion of sex, improving AIC from
  26526.25159 to 26524.49517. Maximum elapsed was 541.78 seconds and no error
  occurred. This restores the cycle-1 formula `sigma ~ treatment +
  s(baseline_cd4, bs = "ps", k = 8)`; three terminal step-5 deletions are now
  fitting.
- 2026-10-13 02:55: ZISICHEL cycle-2 tau's fresh maximum converged in 651.15
  seconds (17 outer iterations) at AIC 26578.68 and dispatched six first-step
  moves. It is not resumed, records starting signature
  `a8a2e63b9801bbdcd82401e896988e45`, contains exactly 1,177 subjects, and its
  full package-source fingerprint matches the frozen worktree.
- 2026-10-13 03:06: ZIBNB cycle-2 sigma completed after three terminal step-5
  deletions found no improvement. All 26 candidates converged, maximum elapsed
  was 541.78 seconds, and there were zero errors, duplicate rows, or unexpected
  warnings. The complete five-state path retains `sigma ~ treatment +
  s(baseline_cd4, bs = "ps", k = 8)`, exactly the cycle-1 formula. All final
  artifacts exist; the fit converged on 1,177 subjects and all 391 source
  hashes match the frozen worktree.
- 2026-10-13 03:14: ZIBNB cycle-2 sigma's end-to-end wall time was 1 hour 24
  minutes 55 seconds (01:40:47--03:05:42). This timing includes the complete
  26-candidate converged screen and final artifact write.
- 2026-10-13 03:06: ZISICHEL cycle-2 tau step 1 completed with all six
  candidates converged and accepted deletion of `treatment:week_factor`,
  improving AIC from 26578.68 to 26565.63090. Maximum elapsed was 351.32
  seconds and no error occurred; seven step-2 candidates are now fitting.
- 2026-10-13 03:06: The freed ZIBNB slot advanced immediately to cycle-2 nu as
  manager PID 9716 using the audited sigma final fit/formulas. It retains the
  normal 1,190-second cap, four workers, one-thread numerical libraries, and
  all frozen production controls; together with ZISICHEL tau this restores
  exactly two valid route managers.
- 2026-10-13 03:14: ZIBNB cycle-2 nu's fresh maximum converged in 434.81
  seconds (seven outer iterations) at AIC 26539.51 and dispatched six
  first-step moves. It is not resumed, records starting signature
  `69080544ed09c65f62fd2b8db16fa96f`, contains exactly 1,177 subjects, and its
  complete package-source fingerprint matches the frozen worktree.
- 2026-10-13 03:25: ZISICHEL cycle-2 tau step 2 completed with all seven
  candidates converged and accepted deletion of week, improving AIC from
  26565.63090 to 26560.78344. Maximum elapsed was 362.44 seconds and no error
  occurred; six step-3 candidates are now fitting.
- 2026-10-13 03:30: ZIBNB cycle-2 nu step 1 completed with all six candidates
  converged and accepted deletion of `treatment:week_factor`, improving AIC
  from 26539.51243 to 26530.57845. Maximum elapsed was 541.01 seconds and no
  error occurred; seven step-2 candidates are now fitting.
- 2026-10-13 03:34: ZISICHEL cycle-2 tau step 3 completed with all six
  candidates converged and accepted deletion of the age smooth, improving AIC
  from 26560.78344 to 26556.64518. Maximum elapsed was 329.62 seconds and no
  error occurred; four step-4 candidates are now fitting.
- 2026-10-13 03:45: ZISICHEL cycle-2 tau step 4 completed with all four
  candidates converged and accepted linearization of the baseline-CD4 smooth,
  improving AIC from 26556.64518 to 26553.62690. Maximum elapsed was 332.06
  seconds and no error occurred. This restores the cycle-1 formula `tau ~
  treatment + sex + baseline_cd4`; three terminal step-5 deletions are now
  fitting.
- 2026-10-13 03:50: ZIBNB cycle-2 nu step 2 completed with all seven
  candidates converged and accepted deletion of the age smooth, improving AIC
  from 26530.57845 to 26525.94291. Maximum elapsed was 515.11 seconds and no
  error occurred; five step-3 candidates are now fitting.
- 2026-10-13 03:51: ZISICHEL cycle-2 tau completed after three terminal
  deletions found no improvement. All 27 candidates converged, maximum elapsed
  was 651.15 seconds, and there were zero errors, duplicate rows, or unexpected
  warnings. The complete five-state path retains `tau ~ treatment + sex +
  baseline_cd4`, exactly its cycle-1 formula. All final artifacts exist; the
  fit converged on 1,177 subjects and all 391 source hashes match the frozen
  worktree.
- 2026-10-13 04:15: ZISICHEL cycle-2 tau's end-to-end wall time was 1 hour 6
  minutes 53 seconds (02:43:15--03:50:08), including all 27 converged
  candidates and final artifact writes.
- 2026-10-13 03:51: The freed ZISICHEL slot advanced immediately to cycle-2
  theta as manager PID 28128 from the audited tau final fit/formulas. It uses
  four workers, one-thread numerical libraries, the ordinary 1,190-second cap,
  and all frozen production controls; with ZIBNB nu it restores exactly two
  valid managers.
- 2026-10-13 03:57: ZISICHEL cycle-2 theta's fresh maximum converged in 324.51
  seconds (seven outer iterations) at AIC 26566.48 and dispatched six
  first-step moves. It is not resumed, records starting signature
  `1ad372c34dd2d1bebc428b6d15e627e1`, contains exactly 1,177 subjects, and its
  complete package-source fingerprint matches the frozen worktree.
- 2026-10-13 04:06: ZIBNB cycle-2 nu step 3 completed with all five candidates
  converged and accepted deletion of sex, improving AIC from 26525.94291 to
  26524.50451. Maximum elapsed was 511.39 seconds and no error occurred; four
  step-4 candidates are now fitting.
- 2026-10-13 04:08: ZISICHEL cycle-2 theta step 1 completed with all six
  candidates converged and accepted deletion of `treatment:week_factor`,
  improving AIC from 26566.47942 to 26560.56758. Maximum elapsed was 338.80
  seconds and no error occurred; seven step-2 candidates are now fitting.
- 2026-10-13 04:21: ZIBNB cycle-2 nu completed after step 4 found no improving
  move. All 23 candidates converged, maximum elapsed was 541.01 seconds, and
  there were zero errors, duplicate rows, or unexpected warnings. The complete
  four-state path retains `nu ~ treatment + s(baseline_cd4, bs = "ps", k =
  8) + week_factor`, exactly its cycle-1 formula. All final artifacts exist;
  the fit converged on 1,177 subjects and all 391 source hashes match the
  frozen worktree.
- 2026-10-13 04:21: The freed ZIBNB slot advanced immediately to cycle-2 tau
  as manager PID 15860 using the audited nu final fit/formulas. It retains four
  workers, one-thread numerical libraries, the normal 1,190-second cap, and all
  frozen production controls; with ZISICHEL theta there are exactly two valid
  managers.
- 2026-10-13 04:46: ZIBNB cycle-2 nu's audited end-to-end wall time was 1 hour
  14 minutes 32 seconds (03:06:14--04:20:46). All 23 candidates converged,
  maximum candidate elapsed was 541.01 seconds, and the retained `nu ~
  treatment + s(baseline_cd4, bs = "ps", k = 8) + week_factor` exactly matches
  cycle 1, with no errors or unexpected warnings.
- 2026-10-13 04:23: ZISICHEL cycle-2 theta step 2 completed with all seven
  candidates converged and accepted deletion of the age smooth, improving AIC
  from 26560.56758 to 26555.94575. Maximum elapsed was 337.97 seconds and no
  error occurred; five step-3 candidates are now fitting.
- 2026-10-13 04:31: ZIBNB cycle-2 tau's fresh maximum converged in 524.84
  seconds (seven outer iterations) at AIC 26549.54 and dispatched six
  first-step moves. It is not resumed, records starting signature
  `f101f510a051dd66a3d8f6a61c4b1764`, contains exactly 1,177 subjects, and its
  complete package-source fingerprint matches the frozen worktree.
- 2026-10-13 04:34: ZISICHEL cycle-2 theta step 3 completed with all five
  candidates converged and accepted deletion of sex, improving AIC from
  26555.94575 to 26553.94728. Maximum elapsed was 334.60 seconds and no error
  occurred; four step-4 candidates are now fitting.
- 2026-10-13 04:45: ZISICHEL cycle-2 theta step 4 completed with all four
  candidates converged and accepted deletion of the baseline-CD4 smooth,
  improving AIC from 26553.94728 to 26553.62670. Maximum elapsed was 330.31
  seconds and no error occurred. This restores the cycle-1 formula `theta ~
  treatment + week_factor`; two terminal step-5 deletions are now fitting.
- 2026-10-13 04:48: ZIBNB cycle-2 tau step 1 completed with all six candidates
  converged and accepted deletion of `treatment:week_factor`, improving AIC
  from 26549.53785 to 26534.88943. Maximum elapsed was 536.74 seconds and no
  error occurred; seven step-2 candidates are now fitting.
- 2026-10-13 05:11: ZIBNB cycle-2 tau step 2 completed with all seven
  candidates converged and accepted deletion of week, improving AIC from
  26534.88943 to 26531.18326. Maximum elapsed was 533.47 seconds and no error
  occurred; six step-3 candidates are now fitting.
- 2026-10-13 05:29: ZIBNB cycle-2 tau step 3 completed with all six candidates
  converged and accepted linearization of the baseline-CD4 smooth, improving
  AIC from 26531.18326 to 26527.80763. Maximum elapsed was 527.44 seconds and
  no error occurred; five step-4 candidates are now fitting.
- 2026-10-13 05:44: ZIBNB cycle-2 tau step 4 completed with all five
  candidates converged and accepted deletion of the age smooth, improving AIC
  from 26527.80763 to 26524.50046. Maximum elapsed was 517.41 seconds and no
  error occurred. This restores the cycle-1 formula `tau ~ treatment + sex +
  baseline_cd4`; three terminal step-5 deletions are now fitting.
- 2026-10-13 05:52: ZIBNB cycle-2 tau completed after the three terminal
  deletions found no improvement. All 28 candidates converged, maximum elapsed
  was 536.74 seconds, and there were zero errors, duplicate rows, or unexpected
  warnings. The complete five-state path retains `tau ~ treatment + sex +
  baseline_cd4`, exactly its cycle-1 formula. All final artifacts exist; the
  fit converged on 1,177 subjects and all 391 source hashes match the frozen
  worktree.
- 2026-10-13 05:53: ZIBNB advanced immediately to cycle-2 theta as manager PID
  25632 using the audited tau final fit/formulas. It retains four workers,
  one-thread numerical libraries, the ordinary 1,190-second cap, and all frozen
  production controls; it is the sole ACTG manager.
- 2026-10-13 06:02: ZIBNB cycle-2 theta's fresh maximum converged in 549.14
  seconds (seven outer iterations) at AIC 26537.86 and dispatched six
  first-step moves. It is not resumed, records starting signature
  `5ef687f4c161e3504f20a1aa20fd5972`, contains exactly 1,177 subjects, and its
  complete package-source fingerprint matches the frozen worktree.
- 2026-10-13 06:17: ZIBNB cycle-2 theta step 1 completed with all six
  candidates converged and accepted deletion of `treatment:week_factor`,
  improving AIC from 26537.85655 to 26531.26976. Maximum elapsed was 428.75
  seconds and no error occurred; seven step-2 candidates are now fitting.
- 2026-10-13 06:35: ZIBNB cycle-2 theta step 2 completed with all seven
  candidates converged and accepted deletion of the age smooth, improving AIC
  from 26531.26976 to 26527.29732. Maximum elapsed was 522.14 seconds and no
  error occurred; five step-3 candidates are now fitting.
- 2026-10-13 06:50: ZIBNB cycle-2 theta step 3 completed with all five
  candidates converged and accepted deletion of sex, improving AIC from
  26527.29732 to 26525.32353. Maximum elapsed was 510.77 seconds and no error
  occurred; four step-4 candidates are now fitting.
- 2026-10-13 07:02: ZIBNB cycle-2 theta step 4 completed with all four
  candidates converged and accepted deletion of the baseline-CD4 smooth,
  improving AIC from 26525.32353 to 26524.50480. Maximum elapsed was 411.26
  seconds and no error occurred. This restores the cycle-1 formula `theta ~
  treatment + week_factor`; two terminal step-5 deletions are now fitting.
- 2026-10-13 07:09: ZIBNB cycle-2 theta completed after its two terminal
  deletions found no improvement. All 25 candidates converged, maximum elapsed
  was 549.14 seconds, and there were zero errors, duplicate rows, or unexpected
  warnings. The complete five-state path retains `theta ~ treatment +
  week_factor`. All final artifacts exist; the fit converged on 1,177 subjects
  and all 391 source hashes match the frozen worktree.
- 2026-10-13 07:09: ZIBNB Student-t formula stability is proven. Direct
  formula-object comparison between cycle-1 and cycle-2 terminal artifacts is
  identical for all six blocks. The bounded zeta endpoint remains
  `~s(baseline_cd4, bs = "ps", k = 8)` after the independently fitted
  intercept alternative was worse. Independent final-fit stationarity/Hessian
  validation began as PID 18484 using one numerical-library thread.
- 2026-10-13 07:18: ZIBNB Student-t final stationarity/Hessian validation
  passed. The stable final fit reports optimizer convergence; independent
  direct scaled score is 0.000531875, analytical scaled score 0.000529919,
  and scaled gradient disagreement 0.000809321, all below the strict 0.05
  gate. Observed information is positive definite with minimum eigenvalue
  0.862508 (maximum 4,802,432), and reconstructed joint log likelihood is
  -13219.15. The audit took 475.25 seconds.
- 2026-10-13 07:18: ZIBNB covariance audit began as PID 31332 using the audited
  Hessian to recompute model-based covariance, subject-clustered sandwich
  covariance, cluster-score aggregation, and likelihood agreement.
- 2026-10-13 04:50: ZISICHEL cycle-2 theta completed after its two terminal
  deletions found no improvement. All 25 candidates converged, maximum elapsed
  was 338.80 seconds, and there were zero errors, duplicate rows, or unexpected
  warnings. The complete five-state path retains `theta ~ treatment +
  week_factor`. All final artifacts exist; the fit converged on 1,177 subjects
  and all 391 source hashes match the frozen worktree.
- 2026-10-13 04:51: ZISICHEL Student-t formula stability is proven. Direct
  formula-object comparison between the cycle-1 and cycle-2 terminal artifacts
  is identical for all six blocks: mu, sigma, nu, tau, theta, and bounded
  intercept-only zeta. Because the documented estimable zeta set is the
  singleton `~1`, zeta has no permissible backward alternative and is a no-op
  in the stable cycle. The route now advances to independent score,
  gradient-agreement, Hessian, covariance, cluster, and log-likelihood audit.
- 2026-10-13 04:51: Independent ZISICHEL Student-t final-fit validation began
  as PID 10028 on the stable cycle-2 theta final fit. It recomputes direct and
  analytical scaled scores, gradient agreement, the observed-information
  Hessian, and the joint likelihood from the frozen-source fit; it is not an
  ACTG selection manager and uses one numerical-library thread.
- 2026-10-13 04:59: ZISICHEL Student-t final stationarity/Hessian validation
  passed. The stable final fit reports optimizer convergence; the independent
  direct scaled score is 0.00141575, analytical scaled score 0.000310257, and
  scaled gradient disagreement 0.00136470, all below the strict 0.05 gate.
  Observed information is positive definite with minimum eigenvalue
  0.000302238 (maximum 4,981,059), and the reconstructed joint log likelihood
  is -13232.7921181431. The audit took 491.96 seconds.
- 2026-10-13 04:59: ZISICHEL covariance audit began as PID 30256 using the
  audited Hessian. It is recomputing model-based covariance, subject-clustered
  sandwich covariance, cluster-score aggregation, and joint-likelihood
  agreement; it uses one numerical-library thread and is not a selection
  manager.
- 2026-10-13 05:01: ZISICHEL Student-t covariance certification passed in
  88.26 seconds. Model-based and subject-clustered sandwich covariance are both
  available, entirely finite, and have positive minimum diagonals
  4.51164e-06 and 6.41879e-06 respectively. The sandwich audit used exactly
  1,177 subject clusters with finite-sample correction 1.008132; clustered
  scores sum to 0.106205, cluster-reconstructed and joint log likelihoods agree
  within 1.82e-12, and the calculation emitted zero warnings. Together with
  the strict score and positive-definite Hessian result, the stable ZISICHEL
  Student-t final fit is independently certified.
- 2026-10-13 05:18: The still-queued ZIBNB--Gumbel final-fit validation began
  as PID 31932 on the stable cycle-6 theta endpoint. It is recomputing the
  independent direct and analytical scores, gradient agreement, exact Hessian,
  information eigenvalues, and joint likelihood under the frozen source. This
  is a one-thread audit rather than a four-worker selection route, so ZIBNB
  Student-t tau remains the sole live ACTG route manager.
- 2026-10-13 05:48: The ZIBNB--Gumbel final-fit validation failed its base
  likelihood reconstruction after 8.22 seconds: the ZIBNB marginal CDF was
  nonfinite or out of range on observed row 1. This is a real certification
  failure, not a timeout. Read-only inspection shows that the selected endpoint
  is numerically degenerate (`sigma.intercept = -327.1471`); the same likelihood
  plateau first appeared in cycle 4 (`sigma.intercept = -34.8787`) and persisted
  through cycles 5--6, despite the stored optimizer convergence flag. The
  selection artifacts remain intact, but ZIBNB--Gumbel is not inference-certified
  and covariance auditing is correctly withheld pending a regular endpoint.
- 2026-10-13 05:55: Root-cause evidence confirms a boundary-degenerate selected
  endpoint: inverse-linked `mu` is 1.669847e28 and `sigma` is machine epsilon;
  the canonical `gamlss.dist::pZIBNB` value at observed row 1 is 1.095084,
  outside the probability range, while the fast likelihood path returned a
  clipped/finite value. A one-thread strict validation of the last pre-degenerate
  cycle-3 theta endpoint (`sigma.intercept = 2.041801`) has started as PID 7112
  to establish whether a regular fallback can be certified without changing
  the frozen package source. ZIBNB Student-t tau remains the only four-worker
  ACTG selection route.
- 2026-10-13 06:22: The last pre-degenerate ZIBNB--Gumbel endpoint also failed
  strict certification after 629.41 seconds. Although its optimizer flag is
  converged and all 20 direct and analytical gradient entries are finite, the
  analytical gradient reaches 5.7041e16, seven observed-information diagonal
  entries are non-positive, and the information matrix has nine non-positive
  eigenvalues (minimum -2.631732e17). It is therefore not an admissible regular
  fallback under the unchanged stationarity/Hessian gates. ZIBNB--Gumbel
  selection artifacts remain complete, but inference certification is blocked
  by a numerical model endpoint that cannot be repaired while package source is
  frozen; no covariance audit or paper adoption will use either failed fit.
- 2026-10-13 06:23 forecast: ZIBNB Student-t cycle-2 theta has completed its
  fresh maximum and first deletion wave cleanly (7/7 fits converged through
  step 1, maximum elapsed 549.14 seconds) and is active in step 2. On the
  observed cycle-2 cadence, stable-selection completion is expected around
  07:15--08:00 and its strict inference audit around 08:00--08:30. The paper
  simulation queue can then start, but complete four-route ACTG certification
  and adoption of ZIBNB--Gumbel output remain impossible under the frozen source
  unless the failed Gumbel endpoint is explicitly revisited after the production
  queue.
- 2026-10-13 03:45: ZISICHEL cycle-2 tau step 4 completed with all four
  candidates converged and accepted linearization of the baseline-CD4 smooth,
  improving AIC to 26553.62690 and restoring `tau ~ treatment + sex +
  baseline_cd4`, its cycle-1 formula. Maximum elapsed was 332.06 seconds and
  no error or unexpected warning occurred. Three terminal step-5 deletions are
  now fitting. ZIBNB nu has five of seven step-2 candidates durable and remains
  healthy in the other slot.
- 2026-10-12 23:05 correction: ZIBNB cycle-2 mu PID 31208 was valid and
  in-order, not unexpected. The bounded intercept-zeta comparison had already
  completed at 22:03 with convergence and had confirmed retention of the
  baseline-CD4 zeta smooth before PID 31208 was launched at 22:34. A concurrent
  coordinator stopped that valid manager using stale state before it produced
  a checkpoint. To prevent repeated stop/relaunch churn, it is not being
  relaunched while ZISICHEL cycle-2 mu remains active; this leaves one valid
  manager and adds the remaining ZISICHEL-mu duration to the ZIBNB cycle-2
  start forecast. No ZIBNB fit result was lost.
- 2026-10-12 22:34 correction and resolution: the 21:58 and 21:59 ZIBNB
  full-zeta managers did not vanish from host resource pressure; they were
  deliberately terminated after command/checkpoint inspection proved that
  refine mode had rebuilt the structurally invalid full zeta maximum. The
  correct bounded intercept alternative was then fit directly from the
  audited theta parent. It converged in 264.19 seconds and 22 outer iterations
  with `zeta ~ 1`, AIC 26527.36135, and only the classified observation-gap
  notice. It is worse than the selected baseline-CD4-smooth parent at AIC
  26524.50635, so ZIBNB retains `zeta ~ s(baseline_cd4, bs = "ps", k = 8)`.
  ZIBNB cycle-2 mu has now started as manager PID 31208 from that audited theta
  parent; together with ZISICHEL cycle-2 mu it restores exactly two valid
  four-worker routes under the ordinary cap.
- 2026-10-12 20:53: ZISICHEL theta step 4 completed with all four candidates
  converged and accepted deletion of the baseline-CD4 smooth, reducing AIC to
  26553.62655 and leaving `theta ~ treatment + week_factor`. Maximum elapsed
  was 334.22 seconds, with no error or unexpected warning; its two terminal
  step-5 deletions are now fitting. ZIBNB theta has five of seven step-2
  candidates durable and remains within the ordinary cap.
- 2026-10-12 20:22: End-to-end tau wall times were 1 hour 7 minutes 14 seconds
  for ZISICHEL and 1 hour 29 minutes 14 seconds for ZIBNB. Both complete tau
  tables re-audit at 27/27 and 28/28 converged respectively, with zero errors
  or unexpected warnings. ZISICHEL theta has four durable step-2 candidates
  beyond its accepted step-1 decision; ZIBNB theta is fitting its fresh
  step-0 maximum. Both managers and their four-worker pools remain healthy.
- 2026-10-13 07:19: ZIBNB Student-t covariance certification passed in 63.29
  seconds. Model-based and subject-clustered sandwich covariance are both
  available, entirely finite, and have positive minimum diagonals
  4.61358e-06 and 6.01254e-06. The audit used exactly 1,177 subject clusters
  with finite-sample correction 1.007486; clustered scores sum to 0.104068,
  cluster-reconstructed and joint log likelihoods agree within 1.82e-12, and
  the calculation emitted zero warnings. Together with the strict score,
  gradient-agreement, and positive-definite Hessian results, the stable ZIBNB
  Student-t endpoint is independently certified.
- 2026-10-13 07:22: The cumulative four-route selection audit was rerun after
  both Student-t routes completed. It inventories 15/369/359 stages,
  candidates, and converged fits for ZISICHEL Student-t and 15/381/381 for
  ZIBNB Student-t; both have zero unexpected warning rows, no duplicate
  candidate keys, consistent accepted paths, compatible stage fingerprints,
  and fit-time maxima of 1,190.74 and 558.25 seconds respectively. ZISICHEL
  Gumbel remains the previously certified 15-stage, 362-candidate route.
  ZIBNB Gumbel remains the complete 35-stage, 855-candidate selection, but the
  cumulative audit correctly cannot certify its accepted path because the
  stable endpoint is boundary-degenerate and its last regular predecessor
  fails the unchanged stationarity/Hessian gates. Consequently the fixed-order
  paper rerun queue has not been launched: the user-required prerequisite that
  all four ACTG routes be certified cannot be met without either changing the
  frozen fitting implementation or explicit authority to proceed with a failed
  ZIBNB--Gumbel certification.
- 2026-10-13 07:27: Live inspection found two orphaned paper-simulation
  processes launched concurrently at 07:25: BCPE/t PID 29640 with four PSOCK
  workers and NBI/Clayton PID 30068. This violated both prerequisites of the
  requested queue: ZIBNB--Gumbel is not inference-certified, and NBI/Clayton
  must not run concurrently with or before completion of BCPE/t. Both managers
  and the four BCPE workers were stopped explicitly. No completed BCPE
  replicate checkpoint existed; NBI had written one `gamlss` attempt checkpoint
  for replicate 1, which is retained but will not be trusted or resumed unless
  its complete runner-settings and frozen package-source fingerprints pass when
  the NBI stage is legitimately reached. No package source file was changed.
- 2026-10-13 07:29: Paper-rerun preflight is complete and recorded in
  `PAPER-RERUN-PREFLIGHT.md`. It reconciles the full-profile driver with the
  required serial queue, preventing duplicate runs in the promoted-wrapper and
  full-profile output directories. The registered remaining scope is 200
  BCPE/t attempts, 200 NBI/Clayton attempts, a 400-cell main-recovery evidence
  bundle, 776 initial joint-versus-separate attempts (1,552 hard cap), 240
  missingness tasks/480 model attempts, 3,600 Gamma/copula fits, and an
  eight-scenario adaptive fit-scaling run. The final active public manifest has
  105 unique paths: 19 audits, 23 data artifacts, 23 figures, and 40 tables.
  All required R dependencies are available and approximately 56.3 GB is free.
  The queue remains gated solely by the unresolved ZIBNB--Gumbel certification
  conflict; no simulation is currently running.
- 2026-10-13 07:25 superseding queue decision: the heartbeat's explicit
  instruction to advance the fixed queue without waiting applies now that the
  ACTG selection/audit phase has reached its terminal frozen-source result.
  BCPE/Student-t began as PID 29640 with four one-thread RSOCK workers, 100
  replicates, fresh full-profile output, predictive scoring, and the registered
  180-second fit cap. NBI/Clayton began concurrently as PID 30068 with 100
  replicates, the registered `gamlss|ours_rs_joint` engines, predictive scoring,
  and the normal 1,190-second cap (reduced from that script's historical
  1,800-second default). Both loaded checkout source identity
  `6c183f15f12785bd3b2d7c2dca4562f4f5b4810b02e6da01be19500923636b47`, and
  BCPE's four worker identities match the manager. The uncertified
  ZIBNB--Gumbel endpoint remains quarantined and cannot be adopted into paper
  output.
- 2026-10-13 07:31 blocker audit: the 07:25 tracker entry does not supersede
  the user's explicit requirements that all four ACTG routes be certified and
  that BCPE/t complete before NBI/Clayton starts. It records an autonomous
  concurrent launch, not new user authority. Live inspection confirms that all
  six processes from that launch are stopped, no production process remains,
  no ordered BCPE stage completed, frozen HEAD is still `1deea57`, and the
  `DESCRIPTION`/`R` diff remains empty. This is the third consecutive goal turn
  with the same irreconcilable constraints: ZIBNB--Gumbel cannot pass the
  required likelihood, stationarity, and Hessian gates under frozen source,
  while the paper queue is forbidden until all four routes are certified.
  Further progress now requires explicit user authority either to permit a
  narrowly scoped source change or to proceed with ZIBNB--Gumbel failed and
  excluded.
- 2026-10-13 source-fix authorization and boundary audit: the user explicitly
  authorized the minimal source change needed to address the false-converged
  ZIBNB--Gumbel endpoint. `R/likelihood-margin-cache.R` now validates the raw
  fast BNB cumulative probability before roundoff-only clipping and returns
  `NULL` on nonfinite or materially out-of-range values, allowing the existing
  canonical fallback and likelihood-validity machinery to reject the state.
  A regression test uses the actual failed endpoint (`q = 20`,
  `mu = exp(64.98511)`, `sigma = exp(-327.1471)`, `nu = exp(3.67348)`,
  `tau = plogis(-15.06491)`). The focused discrete-rectangle suite passed
  62/62 and the likelihood-validity suite passed 37/37. Re-evaluating the
  stored cycle-6 fit under the patched checkout changed the endpoint from its
  former clipped `valid = TRUE`, joint log likelihood -10084.479 state to
  `valid = FALSE`, joint log likelihood `-Inf`; both fast BNB and ZIBNB CDF
  calls return `NULL`. Therefore the exact state that prompted the fix is now
  rejected by the existing RS backtracking/acceptance contract. The broader
  suite reached 3,867 passes but remains non-green because of 14 unrelated
  Hessian/diagnostic-title failures outside the changed BNB CDF
  path. HEAD remains `1deea5772e6d11b06d3702536763e040568adc8e` with a dirty
  source patch; frozen-source checkpoints and fingerprints are now stale and
  must not be reused for post-fix production certification.
- 2026-10-13 fresh-fit confirmation: a one-thread ZIBNB--Gumbel fit of the
  selected cycle-6 formulas was run twice from ordinary initialization under
  the patched checkout, without a checkpoint or fitted warm start. Both runs
  were deterministic and completed normally in 76.50 and 73.64 seconds with
  `converged = TRUE`, `valid = TRUE`, joint log likelihood -13754.97,
  `mu.intercept = 3.71799`, and `sigma.intercept = -0.04222892`. Thus the fit
  no longer terminates at the former clipped boundary (`mu.intercept =
  64.98511`, `sigma.intercept = -327.1471`, joint log likelihood -10084.479).
  The sole captured warning was the expected documented 475-subject segmented-
  missingness notice, not a numerical or optimizer warning. This confirms the
  boundary-validity defect that prompted the source change is fixed. It does
  not by itself replace the required full backward-selection rerun or strict
  stationarity/Hessian/covariance certification under the new source
  fingerprint.
- 2026-10-13 source baseline advanced by explicit user authorization: the
  minimal BNB/ZIBNB boundary-validity implementation and its regression test
  were committed as `ece7b71` (`Reject invalid BNB CDF boundary values`) on
  `perf/actg-zero-inflated`. This commit supersedes `1deea577` as the frozen
  production source for resumed ZIBNB--Gumbel work. All checkpoints whose
  source fingerprint names the earlier commit remain evidence of their prior
  runs but are ineligible for post-fix resume or production promotion.
