# JSS developmental review synthesis

## Bottom line

All eight independent reviewers recommend **do not submit in the current form**.
Their common view is not that the contribution lacks value. The package occupies
a plausible JSS niche and already exposes an unusually complete fit-diagnose-
predict-simulate workflow. The problem is that the manuscript currently asks
readers to trust claims that are either ahead of the released artifact, described
by stale simulation metadata, or broader than the implemented likelihood and
inference support.

The next revision should therefore start with correctness and evidence alignment,
not prose or figure polish. After adversarial normalization, the register contains
35 actions: 14 P0 submission blockers, 18 P1 major changes, and three P2
finishing tasks. `P0` means **must close before submission**, not “execute first”;
the ordered phases in `next-steps.md` control execution.

## Strongest consensus positives

1. **Useful software niche.** A shared-formula longitudinal distributional model
   with covariate-dependent adjacent-time copula dependence is practically and
   statistically interesting.
2. **End-to-end software ambition.** Standard methods, diagnostics, prediction,
   simulation, model screening, bootstrap, and sandwich tools make this more than
   an optimizer wrapper.
3. **Good user-guidance instinct.** Optimizer and sensitivity results are framed
   as decisions users must make, even though the current evidence needs repair.
4. **Broad validation foundation.** Continuous and discrete margins, smooths,
   dependence misspecification, missingness, and standard-model benchmarks are
   all represented.
5. **Substantial reproducibility progress.** The current repository has separated
   public/private workflows, isolated profiles, provenance logs, checkpointing,
   and CI definitions; the remaining problem is release and acceptance evidence.

## Key themes

### 1. Correctness gates precede paper revision

The most serious findings concern the implemented objective and inference:
invalid likelihood contributions can be omitted or neutralized, Hessian signs can
be hidden by absolute values, intermittent gaps are treated as segment breaks
rather than integrated observed-data likelihood, discrete diagnostics are not
calibrated, and some standard-method/family/convergence contracts are unsafe.

These are scientific/software P0s because they can change estimates, likelihood
comparisons, standard errors, and the validity of simulations. The package should
adopt fixed missingness-only contribution masks, fail invalid proposals, enforce
curvature validity, publish a versioned family-capability registry, and expose all
convergence failures before any new headline simulation.

For this submission, the recommended missing-data scope is explicit rather than
algorithmically expansive: contiguous observed blocks are valid likelihood units;
an intermittent gap starts a new block and is described as a segmented composite
likelihood. The paper should not claim an integrated observed-data likelihood.
Interior-gap fits must carry a composite-objective flag, ordinary AIC/BIC and
model-based Hessian inference must be disabled, and robust covariance must be
described only as a sensitivity tool that does not restore the omitted transition
or eliminate informative-gap bias. Use true monotone dropout for the headline
missingness study; keep intermittent MAR as a separately qualified sensitivity.

### 2. The manuscript and staged evidence have diverged

Table 4 does not describe the 100-replicate `n=500, T=4` recovery assets. The
optimizer appendix combines Normal, Gamma, and NBI inputs with different sample
sizes, visit counts, and replication counts, then draws factor-isolation
conclusions contradicted by several displayed intervals. Copula-selection prose
uses stale 10/60 counts rather than the current 100-replicate cells. The current
missingness curves omit 90 failed fits.

Do not rerun the oversized design advertised in Table 4. JSS discourages
extensive simulations, so the main recovery study should be rewritten around the
actual focused design, with generated metadata, denominators, MCSEs, and narrower
claims. In contrast, the optimizer guidance needs a new compact paired factorial
benchmark because the current heterogeneous cases cannot support causal guidance.

### 3. Reframe as a software paper

The paper presently combines a methods paper, software manual, large simulation
paper, and clinical analysis. The preferred target is a 24-28 page main article:

1. bounded contribution and current software landscape;
2. concise conditional Markov-copula model and data contract;
3. package architecture and supported feature matrix;
4. one progressive fully reproducible workflow;
5. focused recovery and comparative evidence;
6. a short secondary LIPID illustration;
7. availability, reproducibility, limitations, and decision boundary.

Move optimizer derivations, Hessian details, full grids, and detailed clinical
selection to supplementary material. The first executable fit should appear by
page 6-7 and be simple enough to interpret immediately.

### 4. Narrow claims to validated estimands

The revision needs a family-specific map of link scale, natural parameter,
independence value, Kendall/tail interpretation, and numerical domain. Replace
generic uses of mean, correlation, and covariance with location and dependence
language where appropriate.

For the next submission, do not claim full joint inference for all smooth and
fixed parameters. State that fixed and smooth inference is conditional on the
fitted smooth structure; model-based fixed-coefficient output also requires a
valid Hessian; label current smooth uncertainty approximate; inventory exactly
which blocks sandwich and bootstrap methods cover; and remove every broader
statement. Subject-level bootstrap or sandwich output remains sensitivity
analysis, not a substitute for unimplemented full fixed-smooth covariance.

### 5. Keep LIPID secondary and transparent

Reviewers agreed that LIPID should not carry the reproducible software story.
Retain at most a 2-4 page secondary illustration with one scientific question,
one model comparison, one effect display, and one diagnostic. Publish exact
data-free preprocessing contracts, formulas, controls, selection steps, and
producers for authorized users, but no data-derived intermediates. State the
private exception plainly and avoid causal or treatment-stability language not
supported by the 10% subsample, missingness, and selection procedure.

### 6. Release and visual polish come after evidence repair

The improved replication system is still in a dirty/untracked state, full
tolerance validation is incomplete, executed producer scripts are not all in the
input graph, and several figure hashes or inline-table mappings are unresolved.
Build and test a clean release candidate only after the P0 evidence is
regenerated. Execute cold Windows/Ubuntu paper runs, canonical-platform full plus
resume, and the complete artifact approval workflow. Create the final tag/archive
only after the manuscript, visual, and style gates also pass.

The visual audit found that 11 of 18 figures require redesign or local
improvement and nine tables require redesign or relocation. Figures 13-18 and
Tables 19-21 currently break appendix/reference order. Apply those changes only
after the analysis set is final.

## Recorded disagreements and resolutions

| Question | Reviewer positions | Resolution for the next revision |
|---|---|---|
| Implement a fully integrated likelihood across intermittent gaps? | Methods review preferred integration or an explicit composite likelihood; applied review treated the current claim as blocking. | Do not add a general integration engine for this paper. Limit ordinary likelihood to complete panels/observed prefixes, flag interior gaps as a different segmented composite objective, disable ordinary AIC/BIC and model-based Hessian inference there, and run true monotone dropout separately. Robust covariance is only a sensitivity tool. |
| Implement full fixed-smooth joint covariance now? | Methods review preferred a full penalized Hessian; editorial review allowed narrower claims. | Narrow the submission claim. State that both fixed and smooth inference is conditional on fitted smooth structure, label smooth covariance approximate, inventory sandwich/bootstrap blocks, and remove unconditional/all-parameter claims. |
| Keep or remove LIPID? | Some reviewers preferred a public example only; others saw value in a brief clinical illustration. | Keep LIPID as a short secondary private-data example. The public BCPE/t workflow remains primary and reproduces every software step. |
| Rerun the large Table 4 design? | Reviewers required either rerun or reconciliation. | Reconcile to the focused `n=500, T=4, R=100` design, add MC uncertainty, and remove unsupported n/T generalization. Do not expand the main recovery grid. |
| How much comparator work is necessary? | Editorial review requested nearest-neighbor software; simulation review emphasized fair GEE settings. | Add a focused task-based comparison against one bivariate joint copula tool, one staged longitudinal copula workflow, and feasible GEE baselines. Document infeasible comparisons rather than forcing them. |

## Recommended editorial thesis

`gamlss.longitudinal` provides a frequentist, shared-formula workflow for fitting
flexible longitudinal marginal distributions jointly with covariate-dependent
adjacent-time copula dependence, together with diagnostics, inference,
prediction, and simulation; it is intended for ordered repeated measurements
whose remaining dependence is plausibly first-order and is not a replacement for
exchangeable or general higher-order dependence models.
