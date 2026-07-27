# Final Simulation Runners

This directory contains the promoted runner scripts used to reproduce the final
BCPE/t and NBI/Clayton simulation artifacts included in the JSS paper workflow.

The formal paper build does not rerun these simulations by default. Instead,
`paper/R/01-simulation-bcpe-t.R` and
`paper/R/02-simulation-delaporte-clayton.R` copy the frozen final artifacts from
`results/jss-exploratory/` into `results/jss-replication/<profile>/` and record
them in `paper/manifest.csv`.

Use these runners only when the final simulation outputs need to be regenerated:

- `bcpe-t/run_bcpe_t_rs_joint_se_diagnostics.ps1`
- `nbi-clayton/run_nbi_highsignal_se_diagnostics.ps1`

Both runners default to the final 100-replicate settings used for the paper.
