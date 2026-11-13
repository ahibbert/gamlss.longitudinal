## Submission candidate

gamlss.longitudinal 0.3.1, first CRAN submission. The checked source archive is
1,164,323 bytes, SHA-256:
70135fb11fd8181e4a470d5d3477048e425fa8ae82a65e7cf15192f4829778b1.
It was built from commit f347e19b2c0a5056d3804e3904fee084834bb9be.

## Test environments and results

Full R CMD check --as-cran completed on 13 November 2026, with tests, examples,
vignette rebuilding, and PDF/HTML manuals enabled:

* Windows 11 x64, R 4.6.2: 0 errors, 0 check warnings, 1 note.
* Ubuntu 24.04.5 x86_64 under WSL2, R 4.6.2: 0 errors, 0 check warnings, 1 note.
* Windows 11 x64, current CRAN R-devel build, R 4.7.0 development revision
  r90881 (runtime source date 2026-11-12; distribution build 2026-11-13):
  0 errors, 0 check warnings, 1 note.

The R-devel check uses the recorded R-release-built dependency binaries; all 32
stated dependency namespaces loaded without warnings/errors before the full
check. Linux verified all 28 recorded non-base dependency versions. No statistical
dependency was patched or pinned to avoid a failure.

Each incoming note combines "New submission" with optional gamlss2 availability
via the declared Additional_repositories entry. gamlss2 is in Enhances and was
not installed for these checks; its use remains conditional.

## Recorded test conditions

Windows release and R-devel each reported 8,888 successful expectations; Linux
reported 8,887. All three reported zero failures, 24 test warning conditions and
18 skips. Actual warning/skip messages, tests, source locations and warning
namespace sequences match across platforms; their full records are retained.

All 24 warning conditions originate in upstream gamlss.dist vector subassignment
in the coverage benchmark's probability/density calls. Static attribution supports
six pZAP, six dZAP and twelve dZINBI conditions; the recorded traces abbreviate
numeric inputs, so this family-level attribution is explicitly source-supported.
These warnings do not validate the affected benchmark values. Upstream report:
https://github.com/gamlss-dev/gamlss.dist/issues/14.

The 18 skips comprise seven existing slow-test CRAN guards, five installed-source
or documentation guards, one separate-paper-source guard, two explicit extended
opt-ins, and three ZISICHEL integration contexts. Source-tree checks retain a
separate verification route; paper and opt-in studies are outside CRAN checking.

The three ZISICHEL contexts retain their original assertions and tolerances. They
skip only after a scalar probability-contract probe recognizes the exact reported
upstream positive-density defect; unexpected contracts or probe warnings fail.
A valid provider automatically reenables the integrations. This unavailable
coverage is not reported as a pass. Independent NBI/ZINBI probability references
replace invalid upstream test oracles; scalar ZAP reference calls preserve the
original inputs and comparison tolerances. No warning suppression is introduced.

## Maintenance scope and limitations

The authenticated locked-review package was originally labelled 0.2.0 at bfb6557
and was subsequently frozen, with only the version changed, as research reference
0.3.0 at 2092e1b. This 0.3.1 maintenance descendant changes packaging, metadata,
documentation and test harnesses. All 400 executable statistical R-source
expressions, seven datasets, defaults, runtime dependency declarations, namespace
and public computational APIs are preserved. The only R-source byte edit is a
separately verified roxygen URL correction.

Six compact guides ship with the package. Detailed articles remain online or in
the repository. The separately revised website example explicitly requests the
existing CG method for its working fit/bootstrap; package defaults and paper
computations are unchanged. Website execution retains its convergence, bootstrap
failure and exploratory uncertainty limitations.

Known gamlss.dist probability/density limitations are described in README, NEWS
and the installed article guide. Native fast paths protect some NBI/ZINBI uses;
this is not a general workaround or qualification of every affected consumer.
The existing Gamma DSS scoring limitation is tracked separately:
https://github.com/ahibbert/gamlss.longitudinal/issues/14.

Original review archives, audits and historical paper-result producer identities
remain unchanged. Same-environment preservation checks have 39 exactly equal
successful cases plus one matching backend refusal retained as a nonpass; retained
fit consumers have ten exactly equal successful cases plus one legacy-newdata
refusal. Preservation and successful package checks do not establish an overall
scientific pass. Paper FULL replication and manuscript/scientific acceptance
remain separate.
