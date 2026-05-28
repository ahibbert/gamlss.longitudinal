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

The source manifest template is `paper/manifest.csv`; the generated manifest is
validated at the end of the workflow.
