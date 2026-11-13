## Working submission status

Maintenance candidate 0.3.1 is not ready for submission. No package has been
uploaded. These are actual preparation results, not final submission evidence.
The exact qualified tarball requires Aydin's separate approval.

Build 05 was built from commit 1303a5c44db78a8139bf615747951138d9c1f9f1.
Its SHA256 is
6922a7ad5ed0da5609a2ce1674c246bd52e3b7dfb2f82a639b6bbc2100ab0648.
This working notes file is excluded from that tarball by .Rbuildignore.

## Scope and provenance

The locked review package was originally labelled 0.2.0 at commit bfb6557.
It was subsequently frozen as version-only 0.3.0 (2092e1b). This maintenance
descendant preserves executable statistical R expressions, data, defaults and
public computational interfaces. Historical paper results and audit records
remain separate and unchanged. This does not certify paper replication or
scientific acceptance; FULL remains incomplete and producer revisions differ.

## Current CRAN qualification

Windows 11 x64, R 4.6.2 (2026-10-29 ucrt), current CRAN dependencies including
gamlss.dist 6.1-11: the full R CMD check --as-cran of build 03 completed with
1 ERROR, 0 top-level WARNINGs and 3 NOTEs. Tests reported 8,815 passes,
7 failures, 612 warnings and 15 skips. Test warnings are separate from
the top-level warning count and require classification and review.

Four failures involve NBI/ZINBI CDF-oracle comparisons; three involve ZISICHEL
score, Hessian and fitting contexts. Independent fixed-parameter probability
identities confirm defects in the current dependency's NBI CDF and ZISICHEL
density. The separate investigation opened upstream issue
https://github.com/gamlss-dev/gamlss.dist/issues/14.
Aydin chose to keep statistical code frozen. No dependency pin, patch,
workaround, weakened test or tolerance change has been adopted. Await a
verified upstream correction and compatibility/impact assessment before
resuming affected qualification. The older dependency is not universally
certified correct, and two corrected statements alone need not restore all
previous behavior.

Build 05 passed a packaging-only R CMD check --as-cran on Windows R 4.6.2
with 0 ERRORs, 0 WARNINGs and one incoming NOTE for a new submission and
optional gamlss2 availability. This check explicitly skipped tests, examples
and vignette execution. PDF and HTML manuals passed; private HTML Tidy was
available, and duplicate Author/Maintainer fields were removed. This scoped
pass does not replace a full check. Full current Linux and R-devel checks
remain unperformed.

## Preservation and documentation evidence

* All 400 R files preserve parsed executable expressions; all seven data files,
  NAMESPACE and runtime dependency declarations are unchanged. The only R-file
  byte edit is a documented roxygen URL correction.
* In separate R 4.4.1 processes with retained gamlss.dist 6.1-1 dependencies,
  39 successful numerical cases are exactly equal. One matching VineCopula
  Gumbel domain refusal remains a nonpass. Original failed/partial comparisons
  remain retained. This is preservation evidence, not current-CRAN qualification.
* One actual retained BCPE/t fit has 10 successful consumers exactly equal and
  one matching legacy newdata refusal. Its sole cross-version serialization
  byte difference is namespace version 0.3.0 to 0.3.1; input/in-process integrity
  checks pass. No refit, historical-SE certification or whole-paper rerun occurred.
* Focused repaired harness checks passed 1,412 expectations. Build 05 source-tree
  checks passed 1,291 expectations with one separate-paper-source skip. With
  numDeriv and VineCopula absent, 286 native expectations passed and four
  numerical-reference cases skipped before dependency use.
* A normal build rendered all six installed guides. Build 05 installed examples
  passed separately. All six exact delivered inst/doc/*.R scripts passed in
  fresh sequential R 4.6.2 processes; the wrapper observed no unhandled warnings.
  Commented templates and illustrative nonexecuted chunks remain unperformed.
* The repository-only JSS start guide rendered successfully. The detailed native
  simulation article remains unvalidated: an explicit process-local ragg device
  passed the earlier plotting failure, then AIC correctly rejected a model with
  numerical_failure at fit-sigma-compare. Fits, seeds and controls were unchanged.
  This failure is not being attributed to the discrete dependency defects.
* URL and spelling reviews passed after documented valid-name exceptions.
  The exact archive retains dataset notices, with no local paths, executable
  binaries or paper/output trees found in the targeted inspection.

## Remaining work before final submission notes

Review the warning classification, obtain the upstream correction and assess
compatibility/impact, resolve the detailed website example without changing its
protected fitting recipe, and complete exact-candidate full checks on current
Windows and Linux R-release plus R-devel. Refresh affected checks after changes.
Then assemble final notes, complete the evidence pack and obtain exact-package
approval. DSS stays unchanged as the owner's accepted nonblocking known issue
(https://github.com/ahibbert/gamlss.longitudinal/issues/14); it is not a release
hold. Paper review and manuscript approval remain separate open tasks.
