# JSS Replication Workflow

This directory contains the full, CRAN-excluded replication workflow for the
JSS paper. The installed package keeps only a lightweight smoke entry point in
`inst/jss-replication/`; this directory is the paper-facing workflow used to
regenerate expanded simulations, tables, figures, logs, session information,
and output hashes.

## Quick Start

```r
source("paper/replicate.R")
```

By default the workflow runs the smoke profile. To run the expanded paper
profile:

```r
Sys.setenv(GAMLSS_LONGITUDINAL_JSS_PROFILE = "expanded")
source("paper/replicate.R")
```

The expanded profile is intentionally not CRAN-safe. It may take a long time and
can optionally use `gamlss2` when that package is installed.

## Outputs

Generated files are written under `results/jss-replication/<profile>/`:

- `tables/`: CSV summaries for paper tables.
- `figures/`: PNG diagnostics and summary figures.
- `logs/`: session information, timing, and output hashes.
- `manifest.csv`: mapping from paper result IDs to generated artifacts.

The source manifest template is `paper/manifest.csv`. Treat this file as the
authoritative map from manuscript result IDs to workflow targets and generated
artifacts. The generated manifest is validated at the end of the workflow, and
`logs/output_hashes.csv` records hashes for generated outputs so reviewers can
confirm that tables and figures were regenerated rather than hand-edited.

## Paper Modules

The replication workflow is organised into seven numbered modules:

1. `01-simulation-bcpe-t.R`: continuous BCPE margin with t-copula simulation.
2. `02-simulation-delaporte-clayton.R`: Delaporte margin with Clayton copula
   simulation and parameter recovery.
3. `03-joint-vs-separate-optimization.R`: joint versus separate optimisation.
4. `04-missingness-dropout-sensitivity.R`: missingness and dropout sensitivity.
5. `05-application-lipid.R`: LIPID clinical-trial application.
6. `06-application-rand-doctor-visits.R`: RAND doctor-visits application.
7. `07-gamma-copula-misspecification.R`: Gamma-margin copula
   mis-specification simulation.

Each module writes at least one CSV artifact and one PNG figure artifact. Most
module outputs are currently marked as stubs so the target graph, manifest
validation, and hash logging are testable before final paper analyses are
filled in. Module 03 now writes the current joint-versus-separate optimisation
candidate review and simulation comparison artifacts.

Module 03 expanded runs default to a practical review grid. Heavier follow-up
runs can set `GAMLSS_LONGITUDINAL_JSS_JVS_REPS`,
`GAMLSS_LONGITUDINAL_JSS_JVS_FAMILIES`, or
`GAMLSS_LONGITUDINAL_JSS_JVS_DESIGNS`.

Module 07 expanded runs default to a 10-replicate pilot over the full
generating/fitted copula grid. Set `GAMLSS_LONGITUDINAL_JSS_MISSPEC_STAGE=full`
to continue from completed pilot checkpoints up to the 100-replicate run.

The LIPID and RAND application data are private and must not be committed to
this repository. Future secure/local runs can provide them through
`GAMLSS_LONGITUDINAL_LIPID_DATA` and `GAMLSS_LONGITUDINAL_RAND_DATA`. The
default workflow passes when these variables are unset and records that the
external data were unavailable.
