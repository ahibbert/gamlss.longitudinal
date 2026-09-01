# JSS-004 main recovery evidence

Run `run_main_recovery_evidence.R` from the repository root after the BCPE/t
and NBI/Clayton production runners finish. The script reads only public
attempt-level and per-replicate inputs, validates the authoritative
`n = 500`, `T = 4`, `R = 100` contract, and writes the reconciled evidence
bundle under `paper/data/public-derived/main-recovery/`. The bundle is marked
authoritative only when both routes carry fresh post-Phase1 runner metadata;
historical staged inputs are explicitly labeled legacy reconciliation.

The builder is replace-safe and resumable. It records input checksums and skips
regeneration only when the schema, confidence level, all inputs, and all outputs
match the existing checkpoint. Set `JSS_MAIN_RECOVERY_RESUME=false` to force a
rebuild. The two fit runners remain the expensive per-replicate checkpointing
layer.

The attempt ledger is the sole source for design dimensions, family identity,
seed provenance,
attempted/successful/converged/retained denominators, failure reasons, and
runtime. Older staged BCPE/t comparator logs did not record convergence or seed
columns; the evidence bundle exposes convergence as `not_reported` and labels
the deterministic seed reconstruction. A fresh full-profile rerun will replace
both with runner-native metadata. In particular, the historical staged NBI run
is Gaussian even though its directory is named `nbi-clayton`; the adapter does
not relabel it as Clayton.

Ordinary recovery summaries use exactly one publication-retained attempt per
method/replicate cell. Retention requires successful execution, explicit
optimizer convergence, and complete fixed, smooth, predictive, and diagnostic
outputs. Nonconvergence and incomplete post-fit output are named failed
attempts; they remain in the attempt ledger and failure-inclusive sensitivity.
Retries increase the attempt-row count but never change the registered 400
planned cells (two studies by two methods by 100 replicate seeds).

The BCPE `gamlss2` reference is a marginal-only comparator. Its fixed-effect
contract contains only `mu`, `sigma`, `nu`, and `tau`; dependence parameters
`theta` and `zeta` are not applicable and are never represented by fabricated
rows. Paired outputs label such non-common estimands explicitly with a zero
paired denominator instead of manufacturing a comparison.

`gamlss2` does not expose a trustworthy `converged` member or the final RS
tolerance. The registered `registered_gamlss2_rs_v1` rule therefore uses its
actual result API: both `logLik` and deviance, every coefficient, and every
fitted value must be finite; iterations must be a positive integer strictly
below `control$maxit` (default 20 when absent). A cap hit is conservatively
`outer_iteration_cap_reached_or_unverified`. Each attempt records the named
status and basis; execution failures and nonfinite values are ineligible.
Every method writes the exact `raw-convergence-2026-09-01.1` payload into each
attempt: API, registered basis and named status, method-specific indicator,
log-likelihood and deviance, coefficient/fitted-value counts and nonfinite
counts, iteration count/cap, and optimizer limit/stall/deterioration flags.
Public validation requires exactly that field set and recomputes convergence,
eligibility, and generic status from it. Derived booleans and named statuses are
checked for contradiction and are never accepted as authority by themselves.

Student-t shape recovery retains non-finite events. The output reports separate
finite, missing, infinite, and numerical-overflow counts for zeta, degrees of
freedom, and their interval widths. The registered descriptive near-Gaussian
boundary is `zeta >= 7.5`; crossing it is counted and never used as an exclusion
rule.

Every authoritative runner loads the checked-out package source on the parent
process and PSOCK workers. Attempts and settings record the source path, package
version, package-source SHA-256, runner SHA-256, settings signature, and settings
SHA-256. These identities are checkpoint fingerprints, so changing package or
producer source invalidates resume.
The recorded checkout identity is the portable repository-relative `.` rather
than an absolute machine path; an absolute or parent-traversing identity is
rejected. Package version and the sorted `DESCRIPTION`, `NAMESPACE`, and `R/*.R`
byte hash establish identity across clones.
The BCPE settings artifact contains one ordered canonical row covering every
registered DGP, optimizer, tolerance, step/backtracking, CG/lambda, inference,
predictive, retry, checkpoint-outcome, worker/core, thread-library, verbosity,
PSOCK scheduler, and R interpreter control. The interpreter binary SHA-256 and
version are included. Its signature is the exact ordered `field=value`
serialization of that row, so changing `N_CORES` (including 1 to 7) invalidates
checkpoints and prevents incomparable runtimes from being pooled.

The public bundle is protected by `output_manifest.csv` and
`bundle_checkpoint.csv`. Validation verifies every artifact hash and then
recomputes authority, design/cardinality, denominators, lifecycle counts, and
publication summaries from attempt-level and normalized per-attempt files;
`bundle_status.csv` and `evidence_validation.csv` are never trusted as authority.
The normalized attempt payload and both canonical runner-settings rows are
the bundled public inputs, and `input_provenance.csv` hashes those actual files.
Paper-profile use additionally requires checkout-external serialized
attestation bytes and a detached Ed25519 signature made by the shared pinned
Phase 2 production key (fingerprint
`cb73c05cede55bfd56357b1780b90c1bf254413d82765c7a682fc9db3a0d8587`).
There is no caller key override and an editable approval CSV has no authority.
A full-profile run emits only a candidate identity for independent signing.
The bundle also records exact-byte hashes for the evidence builder, both fit
runners, target wrapper, central publication gate, target graph, artifact
allowlist, claims registry, signing helper, and signing ceremony documentation.
Git commit/state and a hash of every tracked worktree file are checked against
the validating checkout; self-reported provenance rows are not authoritative.
Candidate building and public validation additionally require Git state
`clean`; any modified or untracked file is a publication-policy failure.
The signed bytes bind the exact canonical bundle, output manifest, checkpoint,
both runner settings files, control-source manifest, package and producer
identities, exact n=500/T=4/R=100 design, every attempt completion instant,
approval instant, and named approver. Approval must postdate every execution.
Validation is repeated on an immutable staged copy and the installed copy, with
before/after source hashes and the same attestation checked each time.

The hardened code gate and the production-evidence gate are separate. No
production attestation currently exists, and manuscript assets are conditional,
until a fresh n=500, T=4, R=100 run under the current runner/schema fingerprints
has completed and its candidate bundle has been independently validated. Runs
started under earlier fingerprints are not publication evidence.
Full profile execution stops after writing and validating an unapproved
candidate bundle and candidate identity; it writes no manuscript tables or
figures. Only paper-profile validation against the checkout-external pinned-key
attestation/signature may generate those assets or satisfy the central gate.
