# Final Simulation Runners

This directory contains the promoted runner scripts used to reproduce the final
BCPE/t and NBI/Clayton simulation artifacts included in the JSS paper workflow.

The ordinary `paper` profile does not rerun these simulations. It validates the
compact, versioned replicate-level inputs in
`paper/data/public-derived/main-recovery/` and rebuilds the manuscript tables
and figures under `results/jss-replication/paper/`. The `full` profile uses the
checkpointed production runners when a complete recomputation is requested.

Use these runners only when the final simulation outputs need to be regenerated:

- `bcpe-t/run_bcpe_t_rs_joint_se_diagnostics.ps1`
- `nbi-clayton/run_nbi_highsignal_se_diagnostics.ps1`

Both runners default to the final 100-replicate settings used for the paper.
