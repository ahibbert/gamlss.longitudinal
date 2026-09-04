# Simplified Phase 2 acceptance ledger

Status: technical gates complete. The final context-free reviewer audit is the
only remaining close-out check.

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
| Main recovery | Pass | BCPE/t and NBI/Clayton each contain 200/200 successful, converged, publication-retained method attempts (400 total). Exact-settings resume accepted every checkpoint with zero pending work and no hash or timestamp changes. The portable 26-CSV bundle passes all 31 evidence checks. |
| Optimizer evidence | Pass | 3 scenarios, 97 paired replicates per scenario, 582/582 attempts converged and retained; producer and submission validators pass. |
| Missingness evidence | Pass | 240 scenario/replicate tasks, 480 fits, 399 retained, 81 named failures, 298 expected warnings, and zero unexpected warnings. |
| Copula misspecification | Pass | 2,400 registered paired candidate fits pass the submission schema and publication validator. |
| Fit scaling | Pass | 62 registered attempts pass the submission schema and reconstruction checks; timing claims remain environment-specific. |
| Portable public boundary | Pass | Public-input scan finds no checkout, user-home, cloud-storage, or private-data paths. Generic interpreter paths remain as explicit execution provenance. Public profiles have no private-data application dependency, and Module 09 is deferred. |
| Package suite | Pass | The complete ordinary P0/P1/P2 source-tree suite passes after correcting the documentation privacy-boundary regression. The only skips are two explicitly opt-in extended SRR stress tests. |
| Reviewer profiles | Pass | Strict paper completed from a detached clean checkout with all 10 Phase 2 gates true, 83 exact public manifest paths, 17 verified claims, and no target warnings/errors (3 minutes 50.6 seconds pipeline). Strict smoke also passed from a detached clean checkout in 1 minute 4 seconds with no target warnings/errors and nine structured optimizer-nonconvergence events. |
| Manuscript asset classification | Pass | Read-only canonical audit classifies 40/40 assets with zero unclassified items. |
| Publisher/compiled-paper isolation | Pass | Dry-run and apply staged exactly 17 renderable assets, all under `charts/` or `tables/`, with no commit or push. A disposable canonical clone compiled to 47 nonblank pages after three already-documented manuscript-source fixes; the remaining 2 unresolved references and 8 citation occurrences are Phase 3 work. |
| Independent final audit | Repeat pending | The first context-free run passed smoke but exposed CRLF-bound main-recovery hashes. Hashes/byte counts are now newline/BOM canonicalized, hash logs are terminally verified, and a new detached checkout passes the complete paper profile. Repeat the context-free audit against this repair before closure. |

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
