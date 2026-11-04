# Package coverage simulation example

This package-shipped example runs a small coverage simulation for the NO, GA,
PO and NBI families using the `gamlss`, `rs_separate`, `rs_joint` and `cg`
methods. Generated outputs go to `results/jss-replication/` beneath your working
directory, outside the installed package.

For the paper's `smoke`, `paper` and `full` profiles, use the separate
[paper reviewer guide](https://github.com/ahibbert/working-paper-gamlss-long/blob/main/paper/REVIEWER.md)
and the exact source pair named in its release record. That guide is currently
a release-preparation draft; this package example does not reproduce the full
paper workflow.

## Run the example

1. Install the package and its CRAN dependencies.
2. From a writable working directory, run:

   ```r
   source(system.file("jss-replication", "run-replication.R",
                      package = "gamlss.longitudinal", mustWork = TRUE))
   ```

3. Retain the generated CSV files and `session_info.txt` from
   `results/jss-replication/` when reviewing the results.

The replication script uses fixed seeds and calls opt-in simulation helpers. It
is not run by CRAN checks.
