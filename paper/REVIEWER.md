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
- `full` reruns all public Monte Carlo studies with fixed seeds and checkpoint/resume support.

The LIPID clinical example is an explicit exception: its source data cannot be redistributed. Sanitized analysis contracts and formulas are available, while its approved figures/tables are static private-data publication assets. The fully simulated BCPE/t workflow reproduces the same model-building, diagnostic, prediction, and visualization sequence without clinical data.

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

For CSV and TeX publication assets, validation requires the approved canonical
hash. PNG hashes are recorded and compared with the approved reference, but a
changed hash is informational because plotting-library and renderer versions can
change binary bytes; PNG validation requires valid dimensions and nonblank
content, followed by the publisher dry-run and visual manuscript check.

The generated profile manifest contains public rows only and uses paths relative
to its result directory. Do not provide private data to any public command.
