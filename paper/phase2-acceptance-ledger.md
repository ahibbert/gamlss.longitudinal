# Simplified Phase 2 acceptance ledger

Status: in progress. Final main-recovery promotion and clean-room verification
remain outstanding.

## Accepted scope

Phase 2 establishes portable, reviewer-facing evidence for five public modules:

1. BCPE/t and NBI/Clayton main recovery;
2. the three-scenario paired RS optimizer comparison;
3. monotone-dropout and intermittent-MAR missingness analyses;
4. Gaussian-versus-Clayton copula misspecification; and
5. fit scaling in its recorded execution environment.

The ordinary reviewer gate is `Rscript paper/replicate.R --profile paper` from
versioned analysis-ready inputs. `smoke` is the fast installation/workflow
check. `full` remains an optional checkpointed author recomputation and is not
a routine reviewer or submission gate.

## Gate ledger

| Gate | Status | Evidence |
|---|---|---|
| Phase 1 prerequisite | Pass | `reviews/2026-08-27-jss-round-1/phase1-acceptance-ledger.md` |
| Main recovery | In progress | BCPE/t has 200 complete method attempts and a no-repeat resume pass; NBI/Clayton production is still running. |
| Optimizer evidence | Pass | 3 scenarios, 97 paired replicates per scenario, 582/582 attempts converged and retained; producer and submission validators pass. |
| Missingness evidence | Pass | 240 scenario/replicate tasks, 480 fits, 399 retained, 81 named failures, 298 expected warnings, and zero unexpected warnings. |
| Copula misspecification | Pass | 2,400 registered paired candidate fits pass the submission schema and publication validator. |
| Fit scaling | Pass | 62 registered attempts pass the submission schema and reconstruction checks; timing claims remain environment-specific. |
| Package suite | Pass | All ordinary P0/P1/P2 tests pass; only two explicitly opt-in extended SRR stress tests are skipped. |
| Reviewer profiles | In progress | Strict smoke is clean; the paper profile awaits final main-recovery promotion. |
| Manuscript asset classification | Pass | Read-only canonical audit classifies 40/40 assets with zero unclassified items. |
| Publisher/compiled-paper isolation | Pending | Run after the final clean paper profile. |
| Independent final audit | Pending | Run after all evidence, tests, and documentation are final. |

## Explicit non-gates

- The exploratory multivariate/VineCopula Module 09 study is deferred.
- Detached signatures, author-held keys, and checkout-specific attestations are
  not required by the reviewer workflow.
- Exact cross-platform PNG byte identity is not required; dimensions,
  nonblank content, hashes for traceability, publisher isolation, and visual QA
  are required.
- A publication-scale `full` rerun is not required of reviewers. The retained
  production studies must instead demonstrate fixed seeds, complete lifecycle
  accounting, and checkpoint resume behavior before their portable CSV bundles
  are promoted.

## Phase 3 boundary

The canonical manuscript wording is not a Phase 2 gate. Its separate audit is
expected to remain red until obsolete empirical prose, pending claim markers,
asset references, citation issues, and the public/private reproducibility text
are revised in Phase 3.
