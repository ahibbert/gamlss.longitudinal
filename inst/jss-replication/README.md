# JSS Replication Workflow

This directory is the package-shipped entry point for reproducing paper-facing
simulation and figure outputs. It is intentionally lightweight for CRAN: source
scripts and instructions live here, while generated CSV, RDS, PDF, and image
outputs should be written outside the installed package, usually under
`results/`, which is excluded from the source tarball.

## Recommended Workflow

1. Install the package and its CRAN dependencies.
2. Optionally install `gamlss2` from the GAMLSS R-universe if comparisons using
   method `"gamlss2"` are required. The CRAN-safe smoke workflow uses the
   `gamlss` marginal baseline plus native RS/CG methods.
3. Run `system.file("jss-replication", "run-replication.R",
   package = "gamlss.longitudinal")` from a working directory where `results/`
   can be created.
4. Record the generated `session_info.txt`, result tables, and figures with the
   manuscript artefacts.

The replication script uses fixed seeds and calls opt-in simulation helpers. It
is not run by CRAN checks.
