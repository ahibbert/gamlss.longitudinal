# Deferred multivariate benchmark

The exploratory Module 09 multivariate/copula-selection benchmark is not part of
the active JSS submission evidence set. It is not run by the reviewer `smoke`,
`paper`, or `full` profiles, does not support an active manuscript claim, and is
not included in the Phase 2 artifact allowlist.

The implementation and its earlier exploratory inputs are preserved in the
annotated Git tag `archive/pre-submission-exploratory-2026-09-08`; they are not
kept in the submission tree. Recover the former module with:

```text
git restore --source archive/pre-submission-exploratory-2026-09-08 -- paper/R/09-simulation-multivariate-longitudinal
```

Re-activating it requires a separately reviewed design, versioned input schema,
tractable reviewer runtime, manuscript scope decision, and new artifact-contract
and manifest entries.
