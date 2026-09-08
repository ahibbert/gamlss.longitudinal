# JSS Replication Workflow

This directory contains the public, CRAN-excluded replication workflow for the
JSS paper. The workflow has no private-data application dependency.

## JSS Submission Preparation

The paper is being reshaped as a JSS-first software article. Use these files as
the authoring front door before changing the manuscript:

- `jss-submission-audit.csv`: current manuscript items mapped to code,
  reproduction status, and recommended placement.
- `jss-manuscript-blueprint.md`: section-by-section replacement outline for the
  JSS article.
- `jss-submission-checklist.md`: blocking gates for manuscript, package,
  replication, and submission packaging readiness.
- `phase2-acceptance-ledger.md`: the simplified public-evidence gate and its
  executed acceptance checks.
- `manuscript/jss-skeleton.tex`: JSS-style manuscript source skeleton.

The main-paper example must be public, package-shipped, or fully simulated.

## Quick Start

```text
Rscript paper/replicate.R --profile smoke
Rscript paper/replicate.R --profile paper
Rscript paper/replicate.R --profile full --workers N
```

The old `expanded` name temporarily maps to `paper`. Each profile has a separate
target store and results directory. The bootstrap restores `paper/renv.lock`,
then installs the checked-out package source before the graph starts. Its
isolated library is stored under the platform user cache, keyed by lockfile hash
and R major/minor version, to avoid cloud-sync locking in OneDrive checkouts.
Set `GAMLSS_LONGITUDINAL_JSS_LIBRARY` to choose another isolated location. See
`REVIEWER.md` for the clean-room path.
Current executed acceptance checks and the remaining full/submission gates are
recorded in `verification-status.md`.

Audit and publish only allowlisted public assets to a clean paper clone:

```text
Rscript paper/audit-manuscript.R --paper-repo <clone>
Rscript paper/publish-assets.R --paper-repo <clone> --profile paper --dry-run
Rscript paper/publish-assets.R --paper-repo <clone> --profile paper --apply
```

## Outputs

Generated files are written under `results/jss-replication/<profile>/`:

- `tables/`: CSV summaries for paper tables.
- `figures/`: PNG diagnostics and summary figures.
- `logs/`: current session/run metadata, concise target timings, structured fit
  and target events, actual-input hashes, artifact hashes, figure-reference
  comparisons, portable target-graph vertices/edges, a reviewer summary,
  hashes of the provenance logs themselves, and the enforced full-profile
  metric tolerance audit.
- `manifest.csv`: mapping from paper result IDs to generated artifacts.

`paper/data/public-derived/` is deliberately tracked: it contains the compact,
approved analysis-ready inputs used by the reviewer workflow. In contrast,
`results/` contains rebuildable run outputs and is deliberately ignored. These
directories therefore have different roles and should not be merged.

The source manifest template is `paper/manifest.csv`. Treat this file as the
authoritative map from manuscript result IDs to workflow targets and generated
artifacts. The generated manifest is validated at the end of the workflow, and
`logs/output_hashes.csv` records hashes for generated artifacts (not mutable log
files), while `logs/provenance-hashes.csv` hashes the completed run logs.

## Public analysis modules

The public graph contains the introductory and BCPE/t software-workflow figures
plus five evidence modules:

1. Main recovery: BCPE/t and NBI/Clayton results from tracked replicate-level
   public inputs (`paper`), based on 100 data-generating replicates per study.
2. Joint-versus-separate optimization: a supplementary three-scenario paired
   RS benchmark with 97 retained pairs per scenario.
3. Missingness/dropout: monotone dropout as the headline analysis and
   time-dependent intermittent MAR as a sensitivity analysis.
4. Gamma-margin copula-family misspecification summaries and heatmap.
5. Fit-scaling evidence, with runtimes explicitly limited to the recorded
   execution environment.

`smoke` instead runs small new public simulations and representative workflow
figures. `paper` never reads ignored local results and is the portable Phase 2
acceptance build. `full` remains available as an optional, checkpointed
author-side recomputation; it is not required for an ordinary reviewer pass.
Reviewer-facing validation checks versioned schemas, complete scenario and
replicate keys, lifecycle denominators, finite numeric fields, and reconstruction
of the reported summaries from their analysis-ready inputs. It does not require
an author-held signing key or compare thousands of incidental intermediate
quantities. The deferred multivariate comparison is documented in
`DEFERRED-MULTIVARIATE-STUDY.md` and is not part of the submission gate.
For approved main-recovery results, the `paper` profile validates the recorded
source manifest and internal provenance carried by the bundle. The stricter
`source_policy = "current"` check is reserved for production and approval: it
requires the bundle to match the exact current producer and package source tree.
This distinction keeps historical results verifiable after interface-only
package releases without representing them as newly computed results.
Manuscript claim wording and exact display selectors are tracked in the generated
Phase 3 checklist but do not block the Phase 2 evidence gate; the separate
manuscript audit continues to fail until those authoring tasks are resolved.

The manifest records which artifact belongs to each module, its input bundle,
verification policy, optional reference hash, and publication status. Public inline
tables that still need conversion to generated `\\input{tables/...}` assets are
listed in `overleaf-inline-table-checklist.md`.

absent from every public profile.
