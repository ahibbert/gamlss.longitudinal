# Authoritative main-recovery evidence

This directory is the portable, reviewer-facing evidence layer for the two
post-Phase-1 recovery studies:

- BCPE margins with first-order Student's t copula dependence; and
- NBI margins with first-order Clayton copula dependence.

Each study uses 500 subjects, four scheduled observations, and 100 fixed-seed
data-generating replicates. Two methods are attempted for each replicate, for
400 method attempts in total. The attempt ledger records execution,
convergence, retention, named failure reasons, seeds, runtime, and source
identity. All summaries and paired comparisons are reconstructed from the
normalized attempt-level CSVs.

The public layer contains exactly the 26 CSV artifacts registered by
`jss_main_recovery_artifact_names()`. It does not contain source datasets,
fitted objects, raw residuals, pseudo-observations, design matrices, or RDS
checkpoints. `output_manifest.csv` and `bundle_checkpoint.csv` provide portable
artifact integrity records; reviewer validation does not require a detached
signature, author-held key, or the original checkout path.

The expensive fit runners use append-only per-attempt checkpoints outside this
directory. Their no-repeat resume behavior was verified before this bundle was
promoted. The ordinary reviewer command consumes these versioned CSVs:

```text
Rscript paper/replicate.R --profile paper
```
