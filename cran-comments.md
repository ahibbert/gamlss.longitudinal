## Submission status

Maintenance candidate 0.3.1 is not ready for submission. No package has been
uploaded. The exact qualified tarball requires the maintainer's separate approval.

## Scope and provenance

The locked review package was originally labelled 0.2.0 at commit bfb6557.
It was subsequently frozen as version-only 0.3.0 (2092e1b). This maintenance
descendant preserves executable statistical R expressions, data, defaults and
public computational interfaces. Historical paper results and audit records
remain separate and unchanged.

## Actual check evidence

Windows 11 x64, R 4.6.2 (2026-10-29 ucrt), current CRAN dependencies including
gamlss.dist 6.1-11: full R CMD check --as-cran of build 03 completed with
1 ERROR, 0 WARNINGs and 3 NOTEs. Seven test failures comprise four NBI/ZINBI
CDF-oracle comparisons and three ZISICHEL score, Hessian and fitting contexts.
The test summary is 8,815 passes, 7 failures, 612 warnings and 15 skips.
The test warnings are distinct from the top-level check warning count and
remain to be fully classified.

Fit-free independent probability identities confirm defects in the installed
current dependency: pNBI at zero is zero despite positive mass at zero, and
two positive ZISICHEL masses each equal one. The retained dependency 6.1.1
satisfies both identities. Published CRAN source contains the corresponding
statements. No dependency pin, patch, workaround, weakened test or tolerance
change has been adopted. The owner impact decision is pending.

Other stages passed installation, package loading/unloading, code analysis,
S3 consistency, Rd/code consistency, examples, data checks, vignette rebuilding,
PDF manual compilation and filesystem checks. HTML manual validation was
skipped because Tidy was absent; a verified private Tidy installation is now
prepared for the next affected check.

The three notes were new submission/optional gamlss2 availability, a redundant
Author field differing from Authors@R, and missing HTML Tidy. The source metadata
now derives author and maintainer from Authors@R; that change awaits a rebuilt
candidate and affected verification.

## Preservation and focused checks

* All 400 R files preserve parsed executable expressions; all seven data files,
  NAMESPACE and runtime dependency declarations are unchanged.
* The six compact installed guides build; long worked/reviewer articles stay
  outside CRAN. Protected fit, control, simulation and seed calls are unchanged.
* Focused harness fixes passed 1,412 expectations in the retained environment.
* Separate current-R processes produced exact results in 34 successful selected
  cases. Four added segmented NBI optimizer cases failed identically; the
  original 38-case comparison therefore failed overall and is not claimed as
  complete equivalence. Supplemental discrete cases are prepared, not run.
* With numDeriv and VineCopula absent, 286 native expectations pass and four
  numerical-reference cases skip before dependency use.
* URL checks pass. Spelling findings were reviewed as valid names/technical
  terms. Archive inspection retains dataset notices and finds no local path
  leaks, executable binaries or generated paper/output trees.

## Remaining qualification

Resolve the dependency impact decision, classify warnings, complete the focused
comparison, rebuild the final candidate, and complete full current Linux,
Windows and R-devel qualification. Validate the website-only articles separately.
Repeat affected checks after changes. Replace these working notes with the
actual exact-candidate evidence before submission; no historical pass for a
different source or dependency environment substitutes for those checks.
