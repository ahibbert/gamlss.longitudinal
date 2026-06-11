## R CMD check results

Initial CRAN/JSS hardening baseline for `gamlss.longitudinal` 0.1.0.

Current readiness targets:

* no hard dependency on non-CRAN `gamlss2`;
* curated public API with low-level numerical internals unexported;
* no generated pkgdown, check, result, or manuscript artefacts in the source
  tarball;
* fast CRAN examples/tests, with long JSS simulations moved to opt-in
  replication scripts;
* user-facing diagnostic wrappers (`plot_terms()` and
  `plot_copula_diagnostics()`) with compatibility shims for older plotting
  entry points;
* CRAN-only smoke replication in `inst/jss-replication/` and a full,
  CRAN-excluded JSS workflow in `paper/`.

This file should be updated again with platform-specific win-builder, r-hub,
macOS, Linux, and R-devel results before submission. GitHub Actions contains a
Linux/macOS/Windows/R-devel matrix, but remote CI status should be recorded here
from the submission commit.

Current local baseline, run on 2026-06-11 on Windows 11 x64 with R 4.4.1:

* `devtools::test(reporter = "summary")`: pass, with 2 opt-in extended-test
  skips and expected max-iteration warnings in coverage-simulation stress
  fixtures.
* `GAMLSS_LONGITUDINAL_EXTENDED_TESTS=true testthat::test_dir("tests/testthat")`:
  pass, 0 failures, 1 expected max-iteration warning, 4 CRAN skips, 1099
  passing expectations.
* `R CMD build .`: pass; source tarball excludes local state, check outputs,
  generated sites/results, paper sources, root benchmark artifacts, generated
  plot files, and reviewer-only repository documents.
* `R CMD check --as-cran --no-manual gamlss.longitudinal_0.1.0.tar.gz`: 0
  ERRORs, 0 WARNINGs, 2 NOTEs.
  Remaining NOTEs are new submission plus optional `gamlss2` availability via
  `Additional_repositories`, and local time verification.
* `urlchecker::url_check()`: pass.
* `spelling::spell_check_package()`: pass.
* New-user smoke test,
  `source(system.file("smoke-tests", "new-user-smoke.R", package = "gamlss.longitudinal"))`:
  pass.
* JSS smoke replication, `Rscript -e "source('paper/replicate.R')"`: pass;
  default smoke profile selected and existing up-to-date targets skipped.
