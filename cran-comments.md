## Working submission status

Maintenance candidate 0.3.1 is not yet ready for submission. A new exact source
archive and full Windows/Linux R-release and R-devel checks are being prepared.
No package approval, upload, confirmation or publication has occurred.

## Scope and lineage

Original locked-review package source bfb6557 was labelled 0.2.0 and subsequently
frozen as version-only 0.3.0 at 2092e1b. Original archives, audits, installation
records and paper producer identities remain preserved. This maintenance release
keeps all 400 statistical R expressions, seven data files, defaults, runtime
dependency declarations and computational APIs unchanged. The only R-source byte
edit is the recorded roxygen URL correction. Paper FULL, scientific acceptance
and manuscript approval remain separate and incomplete.

The owner accepts the reported gamlss.dist defects as an outstanding known issue;
preparation does not wait for an upstream release. The dependency is unmodified.
Known affected distribution outputs are not scientifically qualified by a passing
package check. README, NEWS and the installed article guide explain that limit.
Report: https://github.com/gamlss-dev/gamlss.dist/issues/14.

## Narrow test-reference and dependency-contract handling

Four NBI/ZINBI CDF/log-survival comparisons use independent finite PMF sums and
the defining zero-inflation identity, preserving the existing Poisson-limit
branch, input grid and 1e-12/1e-10 tolerances. Scalar density/log-density checks
are retained. ZAP CDF references use scalar rows to avoid the upstream vector
assignment warning; their parameters and comparison tolerances are unchanged.
No warning suppression was introduced.

Three original ZISICHEL score/Hessian/fit contexts remain intact behind a
precondition that detects only the exact documented density defect. The probe
uses the unmodified dependency at fixed scalar inputs and checks SICHEL mixture,
CDF, zero-mass/log-mass and positive-count density identities. Unexpected
contracts or probe warnings fail normally. A corrected provider automatically
reenables the original assertions. The three unavailable integrations are
reported as skips, not passes. NBI/Delaporte/NO/BCPE controls run separately;
supplemental NBI/Gumbel controls retain coverage of that copula route without
claiming ZISICHEL validity. Tests cover recognized, fixed and unexpected-invalid
provider signatures without patching a namespace.

Focused run 12: 452 passed, zero failures, three exact known-defect skips and two
upstream ZAP reference warnings. After the scalar-reference correction, affected
run 13 passed 223 assertions with zero failures, warnings or skips. These are
focused source-test results against the unchanged B6 statistical namespace;
full rebuilt-archive qualification remains required. Original failed full-check
logs are retained. The new test driver saves complete warning/skip/failure
conditions without weakening assertions or changing failure status.

## Existing evidence and unfinished qualification

The previous B6 archive, SHA256
4dad1dd89f5f1febe8f6caac838f6c0d4abd731a46313a4408f64c4b7334d56d,
is preserved as historical evidence, not the new candidate. Its Linux R 4.6.2
and Windows R-devel packaging-only checks each passed with zero errors, zero
warnings and one incoming note. They skipped tests/examples/vignette execution.
Earlier full Windows B3 had seven failures and 612 unretained warning conditions;
those historical warnings cannot be retrospectively classified from counts.
The next full run will retain actual conditions.

Original-environment preservation checks have 39 successful exactly equal cases
and one matching backend refusal retained as a nonpass. A retained fitted object
has ten successful equal consumers and one legacy-newdata refusal. Preservation
is not independent correctness certification. Six compact installed guides and
examples have scoped execution proof; the new archive will build all six normally.

The excluded website example explicitly requests existing CG for its working
model and bootstrap refits. Other formulas, data, seeds and controls are retained.
Windows/ragg and ordinary Ubuntu/pkgdown runs passed in 675.3 and 372.5 seconds.
Three of five bootstrap smoke refits failed and remain reported. Exploratory
mu-only GG response-mean uncertainty is explicitly distinguished from validated
response-mean inference. Final prose/comments preserve all executed expressions.
No website deployment or whole-paper rerun occurred.

DSS remains unchanged as the owner's separate accepted known issue:
https://github.com/ahibbert/gamlss.longitudinal/issues/14.

Before submission, replace these working notes with exact final archive identity,
full platform results, actual warning/skip explanations and any incoming notes.
Obtain separate approval of that exact qualified package; manuscript approval is
separate. No checks or release gates are waived by accepting the upstream issue.
