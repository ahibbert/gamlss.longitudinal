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
macOS, Linux, and R-devel results before submission.

Current local baseline:

* `devtools::test(reporter = "summary")`: pass.
* `R CMD check --as-cran`: 0 ERRORs, 0 WARNINGs, 4 NOTEs.
  Remaining NOTEs are new submission plus optional `gamlss2` availability via
  `Additional_repositories`, local time verification, README/NEWS Pandoc
  detection in the check subprocess, and skipped HTML manual math rendering
  because optional package `V8` is unavailable.
* `urlchecker::url_check()`: pass.
* `spelling::spell_check_package()`: pass.
* JSS smoke replication, `Rscript paper/replicate.R`: pass; generated tables,
  figures, session information, manifest, and output hashes validate.
