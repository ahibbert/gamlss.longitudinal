# Decision-complete next steps

Work in the order below. Do not begin the manuscript-wide visual or copy-editing
pass until phases 1-3 have frozen the supported model and evidence set.

## Phase 1 - correctness and supported scope

1. **Repair likelihood failure handling (`JSS-001`).** Keep a contribution mask
   determined only by observed/missing data. Compute included marginal and pair
   contributions on stable log scales. Any invalid included contribution returns
   `-Inf`; never replace it with independence or omit it. Add extreme-value and
   boundary tests for continuous and count margins.
2. **Repair Hessian inference (`JSS-002`).** Remove absolute-value SEs. Require
   symmetric, full-rank, negative-definite log-likelihood curvature and a small
   gradient before model-based inference. Register the exact rank, signed-
   eigenvalue, condition-number, and gradient thresholds and test their boundary
   behavior. Return a structured inference failure otherwise.
3. **Lock the missing-data contract (`JSS-003`, `JSS-018`).** Define contiguous
   observed blocks as likelihood units. Ordinary likelihood applies only to
   complete panels and observed prefixes/ignorable monotone dropout. An interior
   gap starts a new block and sets a composite-objective flag; AIC/BIC and model-
   based Hessian inference error with named conditions. Document ordering,
   duplicate times, irregular spacing, scheduled versus observed adjacency, and
   which row supplies pair covariates. Sandwich/bootstrap is sensitivity only.
4. **Lock the family capability registry (`JSS-030`).** Version the supported
   family/copula/pair/diagnostic combinations, enforce one homogeneous family per
   fit, and reject unsupported bounded, discrete, or mixed transitions before
   optimization using named error classes.
5. **Repair standard methods and convergence (`JSS-009`, `JSS-031`).** Correct
   scalar classed `logLik()` with `df`/`nobs`, scalar AIC/BIC where valid, named
   `vcov()` aligned with `coef()`, and structured optimizer stop reasons. Every
   nonconverged fit must be visible and downstream inference qualified or blocked.
6. **Repair diagnostic scope (`JSS-010`).** Implement reproducible randomized
   discrete residuals where supported. Remove formal pass/fail claims and fixed
   thresholds; return descriptive flags and state when higher-lag diagnostics are
   unavailable.
7. **Narrow inference scope (`JSS-032`).** Inventory every model-based, smooth,
   sandwich, parametric-bootstrap, and cluster-bootstrap output. State that fixed
   and smooth results are conditional on fitted smooth structure, label smooth
   covariance approximate, define coefficient blocks and failure handling, and
   remove unconditional/all-parameter claims.

**Phase gate:** independent likelihood values match the package; invalid/singular
fits cannot return ordinary inference; every capability, stop reason, and
inference output has a named contract test; the full test suite passes with zero
unexpected failures and an approved ledger for any expected numerical event.

## Phase 2 - rebuild and reconcile evidence

1. **Main recovery study (`JSS-004`).** Make the actual `n=500, T=4, R=100`
   BCPE/t and NBI/Clayton designs authoritative. Generate the design table from
   attempt metadata. Report attempted, converged, retained, failure reasons, bias,
   RMSE/IRMSE, empirical SD, mean SE, coverage, runtime, predictive metrics, and
   MCSE/intervals. Discuss weak t-copula shape recovery separately.
2. **Optimizer benchmark (`JSS-005`, `JSS-017`).** Retire the heterogeneous JVS
   case interpretation. Use one base design plus only the two or three factors
   that change user decisions, with paired seeds and one factor changed per
   contrast. Compare RS-separate and RS-joint; include CG only if empirical CG
   guidance remains. Choose attempts per cell from a registered MC-precision
   target rather than a fixed 100-replicate rule.
3. **Missingness study.** Rename the current `time_mar` mechanism to
   time-dependent intermittent MAR. Add genuine subject-level monotone dropout
   with no observations after dropout. Report all attempt/failure counts and
   failure-inclusive sensitivity; make monotone dropout the headline analysis.
4. **Copula selection and standard-model benchmarks (`JSS-014`-`JSS-016`).**
   Regenerate selection confusion matrices from current attempt-level data. Show
   family-specific GEE results and treat high-dimensional unstructured GEE as a
   stress test. Build a dated capability table and run only one or two fair,
   tractable nearest-neighbor benchmarks; document non-equivalence/infeasibility.
5. **Automate evidence checks.** Generate scenario tables and data-backed inline
   statements from inputs; test narrative directions, denominators, and scenario
   IDs against the produced artifacts.
6. **Add modest fit-scaling evidence (`JSS-034`).** For one representative
   continuous and discrete model, vary subjects, visits, and smooth basis only.
   Generate median/IQR timings, failures, hardware and version; keep other
   complexity guidance qualitative and bounded by measurements.

**Phase gate:** every empirical claim maps to attempt-level data; scenario and
replication metadata agree automatically; all differences have uncertainty;
failures are visible; no stale counts remain.

## Phase 3 - rewrite the manuscript around the software

1. **Adopt the bounded thesis and landscape (`JSS-011`).** Complete the current
   software review and replace absolute novelty/comparator claims.
2. **Use this main-paper order (`JSS-012`):** introduction/landscape; concise
   conditional Markov-copula model; data contract and package design; progressive
   public workflow; focused comparative evaluation; short LIPID illustration;
   limitations/reproducibility/availability; one-page discussion.
3. **Replace the primary workflow (`JSS-013`).** Use one named public/simulated
   dataset. Begin with an under-12-line simple fit, then add scale, shape, smooth,
   and dependence terms only when a plotted diagnostic or scientific question
   motivates them. The exact displayed code must produce Figures 2-3 and all
   printed output.
4. **Add an estimand and capability table (`JSS-008`, `JSS-018`, `JSS-025`,
   `JSS-027`).** Cover data requirements, one-family-per-fit scope, link/natural
   scales, numerical bounds, pair semantics, prediction conditioning, inference
   status, and use/do-not-use conditions.
5. **Condense LIPID (`JSS-019`).** Limit to 2-4 pages and one scientific question.
   Publish exact data-free formulas and producers for authorized rerun. Remove
   unsupported causal/stability language and present it as a private secondary
   illustration.
6. **Move technical material.** Put CG/RS derivations, Hessian formulas, full
   grids, detailed parameter tables, clinical selection trace, and extended
   diagnostics in supplementary material. Treat 24-28 main pages as an editorial
   aim; the gate is that the public software workflow is the largest component.
7. **Make package documentation executable (`JSS-033`).** Use the same public
   example dataset in core Rd examples and the fitting, diagnostics, prediction,
   and simulation vignettes. Run them in a clean library; justify every expensive
   unevaluated chunk explicitly.

**Phase gate:** executable software appears by page 6-7; the public workflow is
the narrative center; every main claim is supported by final Phase 2 evidence;
LIPID is explicitly secondary; no unsupported inference/likelihood claim remains.

## Phase 4 - release-candidate replication validation

1. Track and commit every workflow file, public-derived input, reviewer guide,
   seed/tolerance registry, CI definition, and producer used by the manuscript.
2. Make `seeds.csv` authoritative and include every executed script/file in the
   targets dependency graph and input hashes.
3. Externalize all public data-backed inline tables and align every manifest row
   with the final TeX label/path.
4. Strengthen figure approval beyond PNG dimensions: hash plotted data and plot
   specifications, define canonical byte or numeric/perceptual policy, record
   human visual approval, and add a mutation test where changed plotted values
   must fail despite valid PNG dimensions.
5. Run cold Windows and Ubuntu `paper` profiles, then `full` from empty
   checkpoints and resume mode on the canonical platform. Wire and pass metric-
   level tolerance comparisons. Smoke and paper must pass on both platforms.
6. Prepare version/commit, license, installation, R/dependency requirements,
   tested platforms, citation, and measured runtimes for the release candidate;
   defer the final tag/archive/DOI until Phase 6.

**Phase gate:** a context-free reviewer can install and reproduce every public
main-paper artifact from the release candidate without intervention; the exact
LIPID exception is stated; full/resume and CI evidence are archived.

## Phase 5 - visual, JSS-style, and submission pass

1. Move to the current JSS LaTeX template and apply `\pkg`, `\proglang`, `\code`,
   code-input/output, caption, heading, and cross-reference conventions.
2. Apply every row of `figure-table-audit.csv`. In particular, move Figure 2
   after its complete code/output; redesign dense Figures 3-5 and 7-12; relocate
   Figures 13-18 and Tables 19-21; repair empty/misowned Appendix C-E headings
   and the interrupted Appendix G code; split Tables 1, 3, 6, 10, 13, and 14.
3. Use accessible palettes plus line/shape redundancy; require caption-sized
   labels. Record per-artifact approval at 100% screen, actual-size print,
   grayscale, and color-accessibility review.
4. Resolve all bibliography keys, remove database debris, and obtain a clean
   BibTeX/LaTeX build with zero unresolved references/citations and zero
   overfull boxes affecting content. Record and approve any remaining benign
   underfull/float warnings individually.
5. Remove all author notes/TODOs and copy-edit using the terminology sheet in the
   writing review. Make captions self-contained and the Discussion a one-page
   use/do-not-use boundary.
6. Repeat a page-by-page visual inspection and a final independent JSS review.

## Phase 6 - immutable archive gate

After Phases 1-5 pass, create the final clean Git tag, source archive, and DOI.
Update the manuscript, package metadata, `CITATION`, manifest, and reviewer guide
to the same immutable version without changing analyses. Re-run hash-only release
checks and archive the Windows/Ubuntu smoke/paper plus canonical full/resume
evidence with the tag.

**Archive gate:** the tag is clean; DOI and archive checksum resolve; every
identity field agrees; archived evidence matches the tag; no post-tag analysis or
asset difference exists.

**Submission gate:** all P0 and P1 register rows are closed with their acceptance
tests; the final PDF is readable in print/grayscale; the manuscript and immutable
release agree; no reviewer needs unpublished context to interpret or reproduce a
claim.
