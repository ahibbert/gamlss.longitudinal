# Phase 2 evidence approval ceremony

The public main-recovery, optimizer, missingness, copula-misspecification, and
multivariate-benchmark validators accept only a detached Ed25519 signature made
by the production key whose public half is pinned in all five validators.
Caller-provided hashes, keys, or editable checkout trust files are
not accepted. The current pinned public-key SHA-256 is
`cb73c05cede55bfd56357b1780b90c1bf254413d82765c7a682fc9db3a0d8587`.
The corresponding private key is retained only in an access-restricted, checkout-external
author-controlled location. It must never be copied into the
repository, a cloud-synchronised folder, logs, tickets, or test fixtures. The
publication gate remains closed until an independent reviewer validates and
signs each production bundle.

1. The reviewer runs `Rscript paper/scripts/phase2-evidence-approval.R
   generate-key C:/checkout-external/phase2-approval` using a directory outside
   every checkout. The helper refuses overwrite, applies mode 0600 and attempts
   a Windows ACL limited to the current account. The reviewer must verify the
   ACL and keep the private key outside Git, cloud-synchronised folders, logs,
   tickets, and test fixtures.
2. A different reviewer compares the displayed public-key fingerprint with the
   external public-key file. In a code-only change, replace the 32 pinned public
   key bytes in all five validators and in the signing helper. Review and merge that
   change before any evidence approval. This pin-only change does not alter the
   producer fingerprints or invalidate durable attempt checkpoints.
3. Run the full candidate validation without promotion. The independent
   reviewer then invokes the helper's `sign` command explicitly. It refuses a
   private key that does not match the pin, reruns candidate validation, and
   signs the exact bundle hash, package and producer hashes, complete canonical
   checkpoint/artifact manifest where applicable, exact registered configuration,
   production audit, runtime/source provenance, approval instant, and reviewer name.
4. Set the study-specific `..._ATTESTATION` and
   `..._ATTESTATION_SIGNATURE` environment variables to the checkout-external
   files and run the public validator. Approval must postdate every checkpoint.
   The exact prefixes are `GAMLSS_LONGITUDINAL_JSS_MAIN_RECOVERY`,
   `GAMLSS_LONGITUDINAL_JSS_OPTIMIZER`,
   `GAMLSS_LONGITUDINAL_JSS_MISSINGNESS`, `GAMLSS_LONGITUDINAL_JSS_COPULA`,
   and `GAMLSS_LONGITUDINAL_MVT`. Module 07 validates resumable per-fit
   checkpoints in the registered `ga-nc-v2` namespace and records rejected
   checkpoints in a quarantine ledger; its signature binds the resulting exact
   six-file evidence bundle: the 2,400-row full-R100 attempt CSV, paired-effect
   and selection artifacts, structured warning audit, and immutable
   generation-time execution manifest
   (1,200 independently seeded generated datasets, each fitted twice)
   for the Phase-1-supported Gamma/Gaussian and Gamma/Clayton routes,
   immutable configuration, checked-out package source, producer, clean Git
   state, dependency versions, registered seeds, workers, timestamps, and
   warning totals. Main recovery binds the exact canonical bundle,
   output manifest/checkpoint, runner settings, control-source manifest,
   package/producer identities, n=500/T=4/R=100 configuration, and every
   execution completion instant. Module 09 binds the immutable root snapshot plus its
   artifact/checkpoint manifests and complete 28-check production audit.
5. Promote the validated evidence in a separate data-only review/commit. Never
   combine a public-key rotation and evidence promotion in one commit.

Rotation repeats steps 1–5. The fit runners do not read the approval key, so a
pin rotation does not by itself change their attempt/checkpoint fingerprints.
For main recovery, however, the validator and signing helper are deliberately
part of the control-source/producer identity: rotation therefore requires a
fresh candidate-bundle build from the still-valid normalized attempts before
an authorized reviewer signs it. Every signature made under the old key becomes
unacceptable after the pin change.
