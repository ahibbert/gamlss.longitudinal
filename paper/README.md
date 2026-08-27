# JSS Replication Workflow

This directory contains the public, CRAN-excluded replication workflow for the
JSS paper. LIPID is a separately identified private-data example and is never a
prerequisite of the public target graph.

## JSS Submission Preparation

The paper is being reshaped as a JSS-first software article. Use these files as
the authoring front door before changing the manuscript:

- `jss-submission-audit.csv`: current manuscript items mapped to code,
  reproduction status, and recommended placement.
- `jss-manuscript-blueprint.md`: section-by-section replacement outline for the
  JSS article.
- `jss-submission-checklist.md`: blocking gates for manuscript, package,
  replication, and submission packaging readiness.
- `manuscript/jss-skeleton.tex`: JSS-style manuscript source skeleton.

The main-paper example should be public, package-shipped, or fully simulated.
Private LIPID and RAND data modules must remain secondary unless accepted data
access instructions are supplied for reviewers.

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
Current executed acceptance checks and the remaining full/CI gates are recorded
in `verification-status.md`.

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

The source manifest template is `paper/manifest.csv`. Treat this file as the
authoritative map from manuscript result IDs to workflow targets and generated
artifacts. The generated manifest is validated at the end of the workflow, and
`logs/output_hashes.csv` records hashes for generated artifacts (not mutable log
files), while `logs/provenance-hashes.csv` hashes the completed run logs.

## Public analysis modules

The public graph contains the introductory and BCPE/t software-workflow figures
plus six analysis modules:

1. BCPE/t recovery from tracked replicate-level public inputs (`paper`) or a
   checkpointed 100-replicate rerun (`full`).
2. NBI/Clayton recovery from tracked replicate-level public inputs (`paper`) or
   a checkpointed 100-replicate rerun (`full`).
3. Joint-versus-separate tables rebuilt from tracked per-replicate deltas, with
   the published 10/100/100 replicate designs rerun with checkpoints in `full`.
4. Missingness/dropout figures rebuilt from tracked aggregate simulation inputs,
   with the published 20-replicate checkpointed study rerun in `full`.
5. Gamma-margin copula-family misspecification summaries and heatmap.
6. Correlation-structure misspecification and standard-model benchmarking.

`smoke` instead runs small new public simulations and representative workflow
figures. `paper` never reads ignored local results. `full` uses fixed seeds and
resumable checkpoints before rebuilding the same publication interfaces. Its
final validation compares 2,000+ grouped estimands against the approved public
inputs using `paper/tolerances.csv`; missing comparison rows and values outside
the registered absolute-plus-relative bounds fail the pipeline.

The manifest records which artifact belongs to each module, its input bundle,
verification policy, approved paper hash, and publication status. Public inline
tables that still need conversion to generated `\\input{tables/...}` assets are
listed in `overleaf-inline-table-checklist.md`.

LIPID data, synthetic analogues, fitted objects, pseudo-observations, residuals,
design matrices, and data-derived intermediates must not be committed. Its
sanitized data contract and formulas are in `paper/R/05-application-lipid.R`.
Approved clinical assets are hash-only private publication artifacts and are
never copied by `publish-assets.R`. RAND is likewise absent from every public
profile.
