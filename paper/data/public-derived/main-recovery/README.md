# Legacy JSS-004 main-recovery reconciliation (not authoritative)

These files reconcile the historical staged public runs to their recorded
design: 500 subjects, 4 scheduled time points, and 100 data-generating
replicates per study. They are retained for workflow testing and traceability;
they are not publication-authoritative Phase 2 evidence.

The historical BCPE runner used a Student's t copula. Despite its directory
name, the historical NBI runner generated and fitted a Gaussian copula. The
attempt ledger and design table therefore label that route NBI/Gaussian (`N`),
not NBI/Clayton. Only fresh post-Phase1 BCPE/t and NBI/Clayton (`C`) production
runs may replace this legacy reconciliation as authoritative evidence.

`attempt_metadata.csv` is the attempt-level source within this explicitly
non-authoritative legacy fixture. `design_table.csv`, all denominators,
failure-reason counts, runtime summaries, and retention counts are generated
from that ledger. Parameter, smooth, prediction, and diagnostic summaries are
joined to publication-retained attempts by retry-aware attempt identity. All
reported Monte Carlo estimates include an MCSE and/or interval appropriate to
the statistic; coverage uses Wilson intervals.

`weak_t_copula_shape_recovery.csv` deliberately separates the Student's
t-copula shape result from general coefficient recovery. It reports link-scale
zeta, implied degrees of freedom, correlation, and lower-tail dependence at the
registered reference profile. The very wide implied-degrees-of-freedom
intervals are evidence of weak shape recovery, not a numerical formatting
error. The hardened producer retains and counts finite, missing, infinite, and
overflow events, with the registered descriptive boundary `zeta >= 7.5`.

The existing staged BCPE comparator log records 100 successful fits but no
convergence indicator. Its attempts are not eligible for ordinary retained-fit
recovery or prediction summaries; `convergence_not_reported = 100` remains
explicit and the attempts remain available to failure-inclusive sensitivity.
The seed values for both older staged logs are reconstructed from the exact
deterministic rules in their runners and labeled as such. A fresh full-profile
run with the updated runners will replace those two legacy metadata fields with
runner-native values. `bundle_status.csv` is the machine-readable publication
gate. A publication bundle additionally requires normalized per-attempt metric
files, paired method differences, `runner_settings_identity.csv`, a complete
SHA-256 `output_manifest.csv`, and `bundle_checkpoint.csv`. Public validation
recomputes authority and summaries rather than trusting status CSVs.

Regenerate with
`paper/scripts/final-simulations/main-recovery/run_main_recovery_evidence.R`.
