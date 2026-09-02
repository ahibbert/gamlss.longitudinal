# Public Phase 2 evidence approvals

This directory is the release bundle of public, detached evidence approvals.
It contains serialized attestations, Ed25519 signatures, and their SHA-256
manifest. These are public verification inputs, not secret key material.

`paper/replicate.R` verifies every bundled byte against `manifest.csv`, copies
the files to a temporary checkout-external directory, and supplies those paths
to the existing pinned-key validators. This preserves the validators' rule that
they never trust an editable in-checkout path while allowing a fresh reviewer
clone to run without receiving separate files from the authors.

The private signing key must never be added here or anywhere in the repository.
The bundle is complete only after all five independently reviewed production
studies have an attestation and signature and `manifest.csv` records their exact
hashes. External author overrides remain supported study by study, but each
override must provide both the attestation and its detached signature.
