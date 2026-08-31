# Adversarial meta-review of the round-1 synthesis

## Verdict

**The synthesis is directionally strong and substantially faithful to the eight
reports, but it is not yet decision-complete. Do not use the current issue
register and next-steps document as the implementation backlog without first
correcting them.**

The central factual findings are supported: all eight reviewers advise against
submission of the frozen draft; the likelihood, curvature, missing-gap,
simulation-metadata, optimizer-study, software-contract, reproducibility, and
submission-form defects are real; and the proposed software-first editorial
thesis is a sensible JSS target. I found no invented headline defect in the
synthesis.

The problems are in normalization and resolution. One high-severity inference
finding has been incorrectly absorbed into the Hessian issue and therefore has
no adequate action or acceptance test. Major package-documentation and runtime
guidance findings are absent. Several proposed remedies do not by themselves
resolve the cited evidence. The dependency graph contains direct and indirect
cycles. Some acceptance tests contradict their own action branches or are too
vague to observe. Finally, the optimizer and comparator programme is larger than
the focused JSS software-paper case requires.

## Confirmed strengths of the synthesis

1. **The core diagnosis is accurate.** The bottom line correctly distinguishes a
   potentially valuable contribution from a presently untrustworthy alignment
   among implementation, evidence, and claims. The counts in the register are
   internally correct: 29 rows, comprising 10 P0, 16 P1, and three P2 items.

2. **The most serious numerical defects are preserved.** `JSS-001` and
   `JSS-002` accurately retain the evidence that parameter-dependent failures can
   remove or neutralize likelihood contributions and that absolute values can
   conceal invalid curvature. Their proposed fail-fast direction is supported by
   the methods review.

3. **The empirical discrepancies are not minimized.** `JSS-004`, `JSS-005`,
   `JSS-015`, `JSS-016`, and `JSS-017` correctly capture the mismatch between
   Table 4 and the staged `n = 500`, `T = 4`, `R = 100` assets; the heterogeneous
   optimizer cases; omitted failure counts and Monte Carlo uncertainty; stale
   copula-selection counts; and unsupported CG guidance. The decision not to
   recreate the oversized Table 4 design is proportionate to a JSS software
   paper.

4. **Several disagreements are resolved sensibly.** It is reasonable not to add
   a general integration engine across intermittent gaps solely for this paper,
   not to require full fixed-smooth joint covariance if claims can be narrowed,
   and not to let private LIPID data carry the public software story. Keeping
   LIPID short and secondary, while making the public example primary, reflects
   the reviewer consensus.

5. **The editorial thesis is appropriately bounded.** The proposed
   frequentist, shared-formula, first-order longitudinal positioning is more
   defensible than the manuscript's absolute novelty language. The software-first
   order, early executable fit, capability table, and use/do-not-use boundary are
   all well supported.

6. **The visual summary is mostly accurate.** The stated totals correspond to
   the audit: ten figures are marked `redesign`, Figure 6 needs local minor
   improvement, and six tables need redesign while Tables 19--21 need relocation.
   The broad instruction to defer visual polishing until the evidence set is
   frozen is sound.

7. **The release direction is correct.** The synthesis accurately recognizes
   that the improved replication machinery exists only in a dirty working state
   and that clean release, provenance, tolerance, resume, and platform evidence
   must be archived before submission.

## Required corrections

### 1. Restore the missing fixed-smooth inference issue

`JSS-002` incorrectly deduplicates three distinct findings. Its title, action,
and acceptance test concern only signed curvature of the fixed-parameter
Hessian, but its reviewer list also cites R01-M6 and R06-M4, and the synthesis
claims that those concerns have been resolved. They have not. R02-M5 establishes
that the fixed Hessian holds smooth coefficients fixed, smooth covariance is
computed separately, cross-block covariance and smoothing-parameter uncertainty
are omitted, and even fixed-effect inference is conditional on estimated
smooths. Removing `abs()` and checking negative definiteness does not repair that
scope problem.

Create a separate P0/P1 inference-scope row, or materially expand `JSS-002`, with
an explicit choice:

- implement and validate the full penalized fixed-smooth covariance, including
  cross-blocks; or
- state that both fixed and smooth inference is conditional on the fitted smooth
  structure, label smooth uncertainty approximate, remove “all parameters” and
  general joint-inference claims, and specify exactly what sandwich and bootstrap
  outputs cover.

The acceptance test must inventory every inferential output and its estimand,
assumptions, fitted coefficient blocks, failure behavior, and validation
evidence. A search must find no broader claim. The current Phase 1 sentence
“Validate model-based fixed-coefficient inference only” is insufficient because
it can still imply unconditional fixed-effect inference.

### 2. Do not imply that robust standard errors solve the intermittent-gap model

The chosen segmented/composite-likelihood scope is defensible, but the synthesis
overstates what “use robust inference” resolves. A sandwich covariance can
address variance under appropriate estimating-equation conditions; it does not
make the segmented objective equal the stated observed-data Markov likelihood,
recover the omitted two-step transition, or remove bias from informative gap
patterns.

Revise `JSS-003` and Phase 1 so the supported primary scope is complete panels or
ignorable monotone dropout/observed prefixes. For interior gaps, either implement
integration or explicitly define a different segmented estimand/objective,
disable ordinary likelihood/AIC interpretations across incompatible objectives,
warn at fit time, and show operating characteristics under the claimed gap
mechanism. Robust covariance may be a sensitivity tool, not the resolution.

The current acceptance phrase “three/four-visit missing-pattern likelihood tests
pass” is ambiguous under the no-integration branch. It must say whether the test
compares against the integrated observed-data likelihood or verifies the
documented segmented objective. The two branches need separate, observable
tests.

### 3. Repair the invalid dependency graph before assigning work

The register is not executable as a dependency graph:

- `JSS-006` depends on `JSS-001:JSS-023 substantive blockers`. This is not a
  defined dependency syntax, includes `JSS-006` itself, and includes `JSS-023`,
  which in turn depends on `JSS-006`.
- There is an indirect cycle
  `JSS-007 -> JSS-026 -> JSS-022 -> JSS-019 -> JSS-007`.

Separate prerequisite work from terminal release gates. A workable order is:
correctness and supported scope; regenerated evidence; manuscript and public
example; final artifact mapping and visual/style work; clean reproducibility
validation; immutable tag/archive. If `JSS-007` continues to combine workflow
implementation with final release, split it into “replication workflow and
wording” and “archive the final accepted state.” Replace the range expression in
`JSS-006` with explicit, acyclic prerequisites or no dependency.

### 4. Split or sharpen incorrectly deduplicated issue rows

Several rows group defects whose remedies and tests differ enough that closure of
one can conceal the others:

- `JSS-009` combines broken `logLik()`/`vcov()` generics, unsupported-family
  likelihood routing, and silent CG nonconvergence. Preserve at least three
  separately checked acceptance clauses. The generic tests should include a
  scalar classed `logLik` with finite `df`/`nobs`, finite scalar `AIC`/`BIC`, and
  a named covariance matrix aligned with `coef`. Family tests must enforce a
  versioned capability registry and reject unsupported bounded/discrete and
  mixed-family paths before optimization. Convergence tests must cover every stop
  reason and block or visibly qualify downstream inference.
- Add R02-M2 to the family-scope evidence. The present normalization omits the
  one-homogeneous-family-per-fit restriction, the unsupported mixed continuous-
  discrete transitions, the missing logarithm in the discrete likelihood prose,
  and the discrete-copula identifiability qualification.
- `JSS-006`, `JSS-022`, and `JSS-023` duplicate template conversion, float order,
  and JSS presentation. Define `JSS-006` as submission cleanliness and resolved
  draft decisions, `JSS-022` as object-level legibility/placement, and `JSS-023`
  as JSS markup/bibliography/build compliance. Do not require the same work three
  times.
- `JSS-007` and `JSS-020` both require creating the immutable tagged release.
  Make `JSS-020` the manuscript metadata/installability check and reserve the
  actual tag/archive gate for one row only.

### 5. Add the omitted major software-paper findings

R04-M6 is not closed by the public worked-example row. The package has no Rd
examples, its principal guidance vignettes are globally unevaluated, and the JSS
template vignette contains undefined objects. Add a P1 documentation action with
observable `R CMD check` and clean-library evaluated-vignette tests. It may share
the same small public dataset as `JSS-013`, but it is a package deliverable rather
than only a manuscript producer.

R04-M8 is also absent. Measured end-to-end replication runtime in Phase 4 does
not tell a practitioner how fit cost changes with subjects, visits, parameter
dimension, smooth basis, family, optimizer, and covariance method. Add a modest
P1 feasibility benchmark or explicitly narrow the manuscript to representative
measured timings and qualitative complexity. Do not require an exhaustive grid.

Also preserve the software review's prediction-default issue somewhere in
`JSS-008`/`JSS-025`: `response`, `mu`, fitted values, and the response mean are
not synonymous for many families, and the current generic mean approximation
needs a documented support/accuracy contract.

### 6. Make alternative action branches have alternative acceptance tests

Several rows allow a narrow documentation remedy in their action but require an
implementation study in acceptance:

- `JSS-010` permits descriptive/configurable diagnostics, yet acceptance always
  requires prespecified false-flag rates and demonstrated power. If automated
  pass/fail is removed, acceptance should instead verify reproducible randomized
  discrete residuals where offered, descriptive labels, no unqualified pass, and
  documented higher-lag unavailability. Calibration/power is required only if a
  formal test remains.
- `JSS-017` permits qualitative CG guidance, yet acceptance always requires a
  public paired study supporting every Table 7 claim. Under the narrow branch,
  remove empirical claims and Table 7 recommendations; retain only the
  deterministic implementation/documentation test.
- `JSS-019` requires an authorized one-command private rerun. Reviewers explicitly
  allowed a fallback if data access or custodian rerun is unavailable: reduce
  LIPID to a short, clearly non-reproducible descriptive illustration with no new
  scientific claim. Its acceptance test must include that branch.
- `JSS-025` should not imply that history-conditioned prediction must be
  implemented. Acceptance can demonstrate that current predictions are
  history-unconditional, document that two subjects with identical covariates
  receive the same result despite different histories, and state dynamic
  prediction as unsupported.

### 7. Restore specific high-severity evidence lost in broad actions

The visual action must explicitly include Figure 2. It is a blocking move because
it interrupts the `margin_screen` output, but Phase 5 names only Figures 13--18
and Tables 19--21. Also include the empty/misowned Appendix C--E headings and the
Appendix G code interruption in the observable float-order check.

For LIPID, either remove the tail-calibration interpretation or require every
lower-, upper-, and opposite-tail claim to point to a displayed statistic. The
current broad instruction to temper language does not resolve R03-M9, where the
t copula matches the lower tail but overstates the upper tail more than Gaussian.

For artifact verification, dimensions, approved hashes, and a visual-review flag
are not enough. Require a deliberate mutation test: changing a plotted series
while preserving a valid PNG and dimensions must fail. Record the plotted-data
and plot-specification hashes, define the canonical byte-identical platform or a
numeric/perceptual tolerance, and distinguish renderer drift from changed
evidence.

### 8. Make acceptance criteria quantitatively observable

Replace terms such as “acceptable conditioning,” “small gradient,” “readable,”
“full test suite,” “within policy,” and “material warnings” with a registered
threshold, named test, artifact, or human approval record. In particular:

- predeclare curvature/rank/gradient tolerances and test their boundary behavior;
- define supported family/coupla combinations and the expected error class;
- state denominator and MC-precision rules per simulation estimand;
- record actual-size print, 100% screen, grayscale, and accessibility approval
  for each visual row;
- require zero unresolved citations/references and zero overfull boxes affecting
  content, while allowing individually reviewed benign underfull/float warnings
  rather than the unqualified and possibly unattainable “warning-clean” gate;
- identify which full run must pass on a canonical platform and which paper/smoke
  profiles must pass on both Windows and Ubuntu.

### 9. Reduce the empirical programme to the minimum that supports the paper

The proposed optimizer benchmark is no longer “compact”: two families, three
methods, five scenarios, and at least 100 attempts per cell imply at least 3,000
fits before retries and uncertainty variants. The mandatory comparator programme
also demands a fully reproducible accuracy/dependence/runtime/usability benchmark
for every near-neighbor class. This conflicts with the synthesis's own
software-first, limited-simulation thesis.

Use a precision-driven design. A base scenario plus the two or three factors that
actually change user decisions is enough; include CG only to support retained CG
guidance. For comparators, a dated capability table plus one or two fair,
tractable task-based benchmarks is sufficient. Document infeasible/non-equivalent
comparisons rather than forcing them. The page target of 24--28 pages should be
an editorial aim, not a pass/fail criterion; the observable criterion is that the
software workflow is central and secondary grids are outside the main article.

### 10. Clarify priority semantics

The ten P0 rows mix statistical correctness, submission cleanliness, release
state, and manuscript interpretation. That is defensible only if P0 means “must
close before submission,” not “execute first.” State this explicitly. Within the
work sequence, `JSS-001`--`JSS-003` and the unsupported-family part of `JSS-009`
precede regenerated evidence; manuscript-template and archive gates come later.
Conversely, the conditional fixed-smooth inference defect must not remain hidden
as a P1-like prose task while JSS formatting is a P0.

## Optional improvements

1. Add columns for `resolution_branch`, `acceptance_artifact`, and
   `closure_evidence` rather than embedding alternatives in prose. This would
   make the register auditable without inflating it into many near-duplicates.

2. Distinguish `scientific blocker`, `software blocker`, `submission blocker`,
   and `release blocker` from priority. This would explain why a JSS template
   defect is submission-blocking but should not be performed before likelihood
   repair.

3. Treat `figure-table-audit.csv` as a visual disposition only. Table 4 is marked
   `keep` visually even though its content is empirically false; the next-steps
   text should state that `keep` never overrides evidence regeneration.

4. Add a traceability matrix from every blocking/major reviewer item to exactly
   one primary issue and, where needed, secondary linked issues. This would have
   exposed the lost fixed-smooth inference and documentation findings immediately.

5. Preserve the useful editorial thesis, but avoid promising that the public
   BCPE/t workflow “reproduces every software step” until the single producer,
   evaluated documentation, and standard-method contracts actually pass.

## Final readiness check

| Check | Status | Reason |
|---|---|---|
| Headline evidence fidelity | **Pass with qualifications** | Core factual findings are supported; no headline hallucination was found. |
| Blocking/major issue coverage | **Fail** | Fixed-smooth inference, executable package documentation, and practical scaling guidance lack adequate normalized actions. |
| Deduplication quality | **Fail** | Several unrelated defects are bundled, while template/release/float work is duplicated. |
| Action resolves evidence | **Fail** | Robust SEs do not repair the intermittent-gap objective; conditional fixed inference remains overstated. |
| Dependency graph | **Fail** | Direct/self and indirect cycles make the register impossible to schedule. |
| Acceptance observability | **Fail** | Several tests are vague or contradict the permitted narrow action branch. |
| Severity/order | **Needs correction** | P0 mixes must-close gates with execution order and omits a distinct inference-scope gate. |
| JSS-proportionate scope | **Fail** | Optimizer and comparator requirements are broader than necessary for the bounded software-paper thesis. |
| Ready to guide implementation | **No** | Correct the register and next-steps controls first; then begin likelihood and scope work. |

After the required corrections above, the synthesis would be a strong basis for
revision. In its present form, it is a good narrative summary but not a safe
closure plan.

## Post-correction verification — 2026-08-27

### Overall result

**Substantive corrections: pass. Bundle safety: fail pending one mechanical CSV
repair.** The corrected synthesis and next-steps document now resolve almost all
of the adversarial findings and are proportionate enough to guide the scientific,
software, evidence, manuscript, and release-candidate phases. However,
`JSS-035`, the terminal immutable-archive gate, has one field fewer than the CSV
header. Its values are shifted left from `consequence` onward, leaving
`dependencies` null. Consequently, the register currently treats the archive
test as an action, reviewer citations as acceptance, `release` as reviewers, `M`
as owner, and the intended dependency list as effort. The final gate therefore
has no valid action/acceptance/ownership/dependency contract. The bundle is not
fully safe to guide implementation until that row is repaired and the dependency
check is rerun.

### Verification results

| Required correction/check | Result | Verification |
|---|---|---|
| Register counts and IDs | **Pass** | The corrected register parses to 35 unique IDs (`JSS-001`--`JSS-035`): 14 P0, 18 P1, and three P2, exactly matching `synthesis.md`. |
| 1. Distinct fixed-smooth inference scope | **Pass** | New `JSS-032` separately records conditional fixed and smooth inference, omitted cross-block/smoothing uncertainty, sandwich/bootstrap block scope, failure handling, targeted conditional coverage, and a search for prohibited all-parameter/unconditional claims. `JSS-002` is now correctly limited to curvature validity. |
| 2. Intermittent-gap scope and robust-inference limitation | **Pass** | `JSS-003`, the synthesis, and Phase 1 consistently restrict ordinary likelihood to complete panels/observed prefixes, mark interior gaps as a segmented composite objective, disable ordinary AIC/BIC and model-based Hessian inference, and state that sandwich/bootstrap sensitivity does not restore the omitted transition or remove informative-gap bias. Its acceptance test now distinguishes the segmented objective from integrated likelihood. |
| 3. Acyclic and schedulable dependencies | **Partial fail** | A graph traversal over all dependency strings that actually parse finds no unknown IDs and no cycles; the earlier direct and indirect cycles are gone. But `JSS-035.dependencies` is null because the row is malformed, so the declared graph omits the terminal dependencies `JSS-006,JSS-007,JSS-020,JSS-023,JSS-026,JSS-033`. Acyclicity and final ordering must be reverified after putting those values in the correct column. |
| 4. Incorrect deduplication | **Pass** | Family routing is separated as `JSS-030`, convergence as `JSS-031`, and inference scope as `JSS-032`; `JSS-009` now covers only standard generic contracts. `JSS-006`, `JSS-022`, and `JSS-023` have distinct cleanliness, visual-order, and JSS-style responsibilities. `JSS-020` prepares availability/installability metadata while `JSS-035` is intended to own the final tag/archive. |
| 5. Omitted documentation and scaling findings | **Pass** | New `JSS-033` requires clean-library Rd examples, evaluated core vignettes, justified unevaluated chunks, and `R CMD check`. New `JSS-034` adds a modest measured scaling benchmark. `JSS-025` now distinguishes `response`, location/`mu`, response mean, family support, approximation, and unsupported dynamic prediction. |
| 6. Branch-consistent actions and acceptance | **Pass** | `JSS-010` removes formal diagnostic pass/fail rather than requiring a full power study; `JSS-017` requires empirical CG evidence only for retained empirical claims; `JSS-019` uses exact data-free code-path tests plus author-run static provenance without claiming public regeneration; and `JSS-025` tests/documentarily excludes dynamic prediction rather than requiring its implementation. |
| 7. Lost high-severity visual, LIPID-tail, and artifact evidence | **Pass** | `JSS-022` and Phase 5 explicitly include Figure 2, Figures 13--18, Tables 19--21, Appendix C--E ownership, and the interrupted Appendix G/reference order. `JSS-019` requires every retained lower/upper/opposite-tail statement to match a displayed statistic. `JSS-026` and Phase 4 require plotted-data/specification hashes and a mutation test that changes values while preserving PNG validity and dimensions. |
| 8. Observable acceptance criteria | **Pass** | Curvature thresholds must be registered and boundary-tested; family and optimizer failures use named/classed conditions; Monte Carlo precision is registered; visual approval records cover 100% screen, actual-size print, grayscale, and color accessibility; build criteria distinguish prohibited errors/overfull content from approved benign warnings; and platform/profile responsibilities are explicit. |
| 9. JSS-proportionate empirical scope | **Pass** | The optimizer study is reduced to a base plus two or three decision-relevant factors, uses MC-precision rather than a fixed 100 attempts, and includes CG only if guidance remains. Comparator work is reduced to a dated capability table plus one or two fair tractable benchmarks. Scaling is intentionally modest, and 24--28 pages is an aim rather than a gate. |
| 10. Priority semantics and execution order | **Pass** | The synthesis explicitly defines P0 as “must close before submission,” while `next-steps.md` controls execution. Correctness precedes regenerated evidence; manuscript/artifact work follows; release-candidate validation precedes style; and the immutable archive is placed last in Phase 6. |
| Explicit final archive ordering | **Pass in next steps; fail in register** | Phase 6 correctly creates the tag/archive/DOI only after Phases 1--5 and prohibits post-tag analysis or asset changes. The intended `JSS-035` dependencies express the same order, but they are currently stored in the `effort` field because of the malformed row. |

### Required final repair

Repair `JSS-035` so all 14 columns align:

- `consequence` should state why the missing immutable archive matters;
- `action` should contain the instruction to create the tag/archive/DOI after all
  prior scientific, manuscript, artifact, and platform gates pass;
- `acceptance` should contain the clean-tag, resolving DOI/checksum, identity-
  agreement, and archived-evidence checks;
- `reviewers` should be `R01-B4;R04-M7;R05-B1/M6`;
- `owner` should be `release`;
- `effort` should be `M`; and
- `dependencies` should be
  `JSS-006,JSS-007,JSS-020,JSS-023,JSS-026,JSS-033`.

Then rerun the CSV schema/count/unknown-dependency/cycle checks. If the corrected
row preserves those intended dependencies, the graph remains acyclic and the
bundle will be safe to guide implementation.

**Final closure — 2026-08-27.** `JSS-035` is now correctly aligned, with a distinct consequence, action, acceptance test, reviewer provenance, owner, effort, and terminal dependencies. Strict revalidation confirms 35 unique rows; no blank action, acceptance, reviewer, owner, or dependency fields; only allowed effort values; no unknown dependency IDs; and no dependency cycles. Together with the substantive passes recorded above, this closes the remaining meta-review defect: **the corrected bundle is now safe to guide implementation**, subject to closing each registered acceptance test in the stated phase order before submission.
