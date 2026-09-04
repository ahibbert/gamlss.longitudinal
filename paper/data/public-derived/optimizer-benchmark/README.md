# Lean optimizer benchmark

This directory contains the reviewer-facing outputs from the supplementary
three-scenario paired comparison of separate and joint RS optimization. Each
scenario contains 97 paired replicates for `rs_separate` and `rs_joint`.

Runtime summaries are descriptive only: the methods were run sequentially, so
timings can be affected by execution order and machine load. The author-side
production workflow is checkpointed, but this public bundle intentionally
contains only analysis-ready CSV and figure outputs. It omits host names,
process identifiers, checkout paths, checkpoint signatures, and fitted model
objects.

The broader multivariate Module 09 benchmark is not part of the submission
evidence gate and is not included here.
