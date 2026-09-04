# Archived Phase 2 evidence-signing design

This ceremony is no longer part of the submission or reviewer workflow. The
earlier design required checkout-external Ed25519 attestations for several
public evidence bundles. It was retired because it added operational machinery
without improving a reviewer's ability to reconstruct or assess the reported
results.

The current workflow uses versioned input schemas and direct scientific checks:
complete scenario and replicate grids, explicit attempt/convergence/retention
denominators, finite numeric fields, and reconstruction of reported summaries
from the public analysis-ready inputs. Author-side production checkpoints still
record source and settings fingerprints to prevent incompatible resume, but no
private key, signature, attestation file, or approval environment variable is
required by `smoke`, `paper`, or `full`.

The historical helper remains in the repository only as implementation history
and is not a reviewer prerequisite. The multivariate benchmark that used the
most elaborate form of the ceremony is deferred from this submission; see
`paper/DEFERRED-MULTIVARIATE-STUDY.md`.
