## Working submission status

Maintenance candidate 0.3.1 is not ready for submission. No upload or release
approval has occurred. These are actual preparation results, not final notes.

Exact build 06: commit b97580bf000a5fa118496f4b9415bced2e9c5424.
Tarball SHA256:
4dad1dd89f5f1febe8f6caac838f6c0d4abd731a46313a4408f64c4b7334d56d
(1,160,839 bytes). This notes file is excluded by .Rbuildignore.

## Scope and lineage

The locked review package was originally 0.2.0 at bfb6557. It was subsequently
frozen as version-only 0.3.0 (2092e1b). This maintenance descendant preserves
executable statistical R expressions, data, numerical defaults and public
computational interfaces. Historical paper results, producers and audits remain
separate and unchanged. FULL and overall scientific acceptance remain incomplete.

The only R-source byte edit is a documented roxygen URL correction. Test harness,
metadata and documentation changes are listed separately. The latest additions
retain raw test warnings and report unavailable worked-example comparisons.
The CRAN tarball changes no statistical methods, fitting recipes, seeds or
controls. The excluded website article has the separately authorised CG recipe
exception described below.

## Observed checks

| Source and environment | Actual scope and outcome |
| --- | --- |
| B6, Ubuntu 24.04.5, R 4.6.2 | Packaging-only --as-cran: 0 errors, 0 warnings, 1 incoming note; indexed PDF and HTML manuals passed. |
| B6, Windows, R-devel 2026-11-11 r90876 (4.7.0) | Packaging-only --as-cran: 0 errors, 0 warnings, 1 incoming note; PDF and HTML manuals passed. All 32 declared dependencies loaded; dependency binaries were reused from the recorded R-release library. |
| B5, Windows, R 4.6.2 | Earlier packaging-only --as-cran: 0 errors, 0 warnings, 1 incoming note. |
| B3, Windows, R 4.6.2/current gamlss.dist 6.1-11 | Full --as-cran: 1 error, 0 top-level warnings, 3 notes. Tests: 8,815 passes, 7 failures, 612 warnings, 15 skips. |

Packaging-only checks explicitly skipped tests, examples and vignette execution
and rebuilding. They do not replace full exact-candidate qualification. The
single incoming note concerns a new submission and optional gamlss2 availability
through the additional repository. The earlier redundant-author and missing-Tidy
notes were resolved. Linux initially lacked Inconsolata and its deferred upquote
dependency; existing local TeX resources repaired that isolated environment.
The original warning/failing font proofs are preserved. A GNU cp portability
notice came from the external evidence-copy wrapper, outside the package/check
warning count; exports and archive hashes were verified.

## Numerical qualification remains open

Four full-check failures concern NBI/ZINBI CDF-oracle comparisons; three involve
ZISICHEL score, Hessian and fitting contexts. Independent probability identities
confirm defects in current gamlss.dist. The separate investigation reported these
at https://github.com/gamlss-dev/gamlss.dist/issues/14. CRAN still lists 6.1-11
at the latest review. Aydin chose frozen statistical code; no dependency patch,
pin, workaround, weakened assertion or tolerance change was adopted.

Await an identifiable upstream correction and compatibility/impact assessment
before affected full checks. Some older dependency defects and approximation
differences are separate; two typo corrections do not guarantee every test passes.

The original 612 warning conditions were not retained by testthat's CRAN
reporter. They remain unclassified. B6 now saves that reporter's raw conditions,
test/source references and session information even on failure. Synthetic
CRAN-mode probes preserve the original printed output and nonzero failure exit,
including repeated warnings. This prepares the next justified full run; it does
not recover or classify the original warnings.

## Preservation and focused evidence

* All 400 statistical R files retain identical parsed expressions; all seven
  data files, NAMESPACE and runtime dependency declarations are unchanged.
* Separate original-environment R 4.4.1/gamlss.dist 6.1-1 processes produced
  39 successful exactly equal numerical cases, plus one matching backend refusal
  retained as a nonpass. This is preservation evidence, not current qualification.
* One actual retained BCPE/t fit produced ten successful equal consumers and
  one matching legacy-newdata refusal. The sole cross-version serialization byte
  is namespace version 0.3.0 to 0.3.1. No historical-SE certification or refit.
* Focused repaired harness checks passed 1,412 expectations. B5 source-tree
  checks passed 1,291 with one separate-paper-source skip. Optional-absence
  checks passed 286 native expectations with four numerical-reference skips.
* B6 built all six guides normally. Its six delivered extracted scripts are
  byte-identical to B5's six successfully executed fresh-process scripts. All
  407 statistical source/data files and NAMESPACE are byte-identical between
  B5 and B6. The only generated guide differences are three elapsed-time values.
  Earlier installed examples and source-tree results retain their actual scope.
* B6 URL and spelling checks passed. Archive boundaries, attribution and targeted
  privacy checks are retained; the new raw-condition output is excluded.

## Website example and remaining steps

The unchanged JSS start guide rendered successfully. Aydin authorised an
example-only revision to the excluded detailed native article. Its working model
now requests the existing CG optimizer, and bootstrap refits receive CG explicitly.
The data, formulas, seeds, replicate count and public default controls remain
unchanged. Of 50 protected calls across the guides, 48 remain identical and only
these two optimizer requests differ. Package statistical code remains frozen.
Unavailable/nonconverged comparisons remain visible and no complete-set winner
is claimed.

The revised article completed 75/75 knitting steps in 675.3 seconds within a local
900-second validation budget on Windows R 4.6.2, using process-local ragg output
and null devices. This is not an established website-service timeout or an Ubuntu
pkgdown validation. CG converged and the inference guard passed; diagnostics,
prediction and simulation executed. Of five bootstrap smoke attempts, three
failed convergence and are reported. This is not adequate bootstrap inference.
Final prose was checked against the new diagnostic figures and clarified without
changing executable chunks after the run. Illustrative eval=FALSE templates
remain unexecuted. No website deployment was performed.

Complete full exact-candidate Windows/Linux R-release and R-devel checks after
the dependency assessment, classify retained warning conditions, validate the
actual website build when appropriate, assemble final notes and obtain
exact-package approval.
DSS remains the owner's unchanged, accepted nonblocking issue:
https://github.com/ahibbert/gamlss.longitudinal/issues/14.
Paper review and manuscript approval remain separate and open.
