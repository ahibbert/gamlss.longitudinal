# Reviewer replication guide

Run from the repository root with R 4.4.1 initially recommended:

```text
Rscript paper/replicate.R --profile smoke
Rscript paper/replicate.R --profile paper
Rscript paper/replicate.R --profile full --workers N
```

The first run restores the lockfile and installs this checkout into an isolated project library. Outputs and target stores are separated by profile under `results/jss-replication/` and `paper/_targets/`.

- `smoke` uses only newly simulated public data and is intended to finish in 15 minutes.
- `paper` rebuilds every public manuscript asset from versioned analysis-ready inputs and is intended to finish in 60 minutes.
- `full` reruns the retained public Monte Carlo studies with fixed seeds and
  checkpoint/resume support. The exploratory multivariate benchmark is deferred
  from this submission and is not in the target graph.

The public workflow has no private-data application dependency. The fully simulated BCPE/t workflow reproduces the complete model-building, diagnostic, prediction, and visualization sequence.

Every run writes session information, Git/dirty state, lockfile and actual-input
hashes, seeds, cores, target timings, structured target/fit events, artifact
hashes, figure-reference comparisons, and provenance-log hashes beneath the
profile's `logs/` directory. Platform startup messages emitted before the script
runs (for example a host's invalid locale setting) cannot be captured there and
do not count as target warnings.

Start with `logs/reviewer-summary.md`. The two `target-graph-*.csv` files provide
a portable vertex/edge dependency graph, and `input-hashes.csv` covers the
package R source, workflow definitions/configuration, and only the public data
bundle actually read by the selected profile.
For text inputs, `sha256` is the newline/BOM-canonicalized digest used by the
cross-platform verifier; `raw_sha256` is also supplied for comparison with
ordinary file-hashing tools. Binary inputs use the same raw digest in both.

Reviewer-facing evidence is checked using versioned schemas, complete scenario
and replicate keys, explicit attempt/convergence/retention denominators, finite
numeric fields, and reconstruction of publication summaries. It requires no
author-held signing key or detached attestation. PNG hashes are recorded and
compared with the reference, but a changed hash is informational because
plotting-library and renderer versions can change binary bytes; PNG validation
requires valid dimensions and nonblank content, followed by the publisher
dry-run and visual manuscript check.

The generated profile manifest contains public rows only and uses paths relative
to its result directory. Private static manuscript assets are hash-audited
separately and are never copied, regenerated, or published by these commands.
Do not provide private data to any public command.
