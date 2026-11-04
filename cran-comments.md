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

Current local baseline, run on 2026-06-14 on Windows 11 x64 with R 4.4.1:

* `devtools::document()`: pass. Local roxygen2 7.3.2 reports that the package
  was previously documented with roxygen2 7.3.3; no documentation generation
  failure.
* `devtools::test(reporter = "summary")`: pass, with 2 opt-in extended-test
  skips and expected max-iteration warnings in coverage-simulation stress
  fixtures.
* `GAMLSS_LONGITUDINAL_EXTENDED_TESTS=true testthat::test_dir("tests/testthat")`:
  pass, 0 failures, 1 expected max-iteration warning, and 4 CRAN skips.
* `R CMD build .`: pass with vignettes built. The resulting source tarball
  excludes local state, check outputs, generated sites/results, paper sources,
  root benchmark artifacts, generated plot files, and reviewer-only repository
  documents.
* `R CMD check --as-cran --no-manual gamlss.longitudinal_0.1.0.tar.gz`: 0
  ERRORs, 0 WARNINGs, 3 NOTEs. Remaining NOTEs are new submission plus optional
  `gamlss2` availability via `Additional_repositories`, local time
  verification, and unavailable local pandoc for checking `README.md`/`NEWS.md`.
* `urlchecker::url_check()`: pass.
* `spelling::spell_check_package()`: pass.
* New-user smoke test,
  `source(system.file("smoke-tests", "new-user-smoke.R", package = "gamlss.longitudinal"))`:
  pass.
* JSS smoke replication, `Rscript -e "source('paper/replicate.R')"`: pass;
  default smoke profile selected and existing smoke targets reported up to date.
