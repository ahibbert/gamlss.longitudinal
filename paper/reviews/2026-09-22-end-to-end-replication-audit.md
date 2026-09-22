# End-to-end replication audit

Status: in progress  
Audit date: 2026-09-22 (Australia/Sydney)

## Benchmark

The manuscript benchmark is the 49-page PDF
`gamlss_longitudinal__Longitudinal_GAMLSS_Models_With_Copula_Dependence_in_R (4).pdf`,
created at 2026-09-22 15:42:10 and identified by SHA-256
`C1D832970B707828E41F92B526760D0A31CF3361A2FA9910D47F0E7A6859E436`.
It contains 18 numbered figures and 25 numbered tables.

The closest available TeX source is the 2026-09-22 14:41 attachment used for
the asset inventory. The canonical `working-paper-gamlss-long` repository is
not yet the source of this PDF: its current head is
`68c3bad26626ce7c267bd330364cfda8df7a6b76` (Initial Overleaf Import,

## Current verdict

The promoted ACTG computational workflow is reproducible, but the latest
manuscript is not yet wired end to end to that workflow. The remaining blockers
are manuscript-source and publication-contract issues, rather than missing
ACTG analysis code.

## Blocking gaps

1. **Canonical source mismatch.** The latest ACTG TeX must become the canonical
   paper repository source before the read-only audit, publisher, and compile
   check can validate the same document as the submitted PDF.
2. **Retired analyses remain in the manuscript.** The latest source references
   eight artifacts deliberately marked `retired`: the joint-versus-separate
   scenario table and three result tables, the correlation-misspecification
   scenario and two result tables, and the main simulation-design table. The
   fastest submission path is to remove those retired appendix analyses and
   their narrative references rather than reactivate their replication scope.
3. **Generated ACTG tables are not used by TeX.** Seven ACTG tables remain
   inline. The workflow generates complete table environments under
   `results/jss-replication/paper/tables/`; each entire inline table block must
   be replaced by the corresponding `\input{tables/actg-....tex}` call.
4. **The reproducibility section is unfinished.** It currently contains author
   comments rather than the `smoke`, `paper`, and optional `full` commands,
   measured runtimes, environment statement, artifact coverage, and public
   ACTG data statement.
5. **The PDF contains unresolved references and visible drafting notes.** In
   particular, two correlation-benchmark references render as `Table ??`, the
   optimizer appendix callout renders as `section ??`, and draft notes remain
   visible in the algorithm appendix.
6. **Reference hashes are not yet frozen.** Twenty-six of the 29 active public
   assets referenced by the latest source have no approved reference hash.
   Output hashes are recorded and recomputed, but the final manuscript-facing
   figure/table baseline still needs author approval and freezing.
7. **Some automated claim tests inspect a stale local TeX copy.** Several tests
   and Phase 2 claim checks still read `paper/manuscript/main.tex`, although the
   repository policy says the external paper repository is the only live TeX
   source. Passing those tests does not currently establish that the latest
   ACTG manuscript obeys the same scope contract. These checks should move to
   a source-independent claim/asset contract, with `audit-manuscript.R` used
   against the actual paper clone.
8. **One selection-path count needs reconciliation.** The manuscript describes
   the ACTG selection history as 51 steps. The promoted ledgers contain 41
   accepted non-initial changes and 18 targeted candidate fits (59 rows when
   these are combined, or 71 rows including stage baselines). The prose should
   state the intended denominator explicitly rather than retain the unsupported
   count of 51.
9. **Four simulation table producers do not reproduce the manuscript tables.**
   The files published as the BCPE/t and NBI/Clayton `fit_characteristics` and
   `fixed_parameter_bias_rmse` tables currently contain lifecycle counts and
   aggregate recovery summaries. The PDF's Tables 11, 16, 17, and 18 instead
   contain fit/predictive metrics and full parameter-level recovery results.
   Applying the publisher would silently replace those four manuscript tables
   with different content under the same filenames. This is a P0 producer bug;
   the formal producer must emit the exact manuscript table schemas from the
   versioned main-recovery evidence bundle.
10. **The current recovery bundle does not support the manuscript's coverage
    columns for `gamlss.longitudinal`.** For the BCPE/t package fits, standard
    errors are available for only 1 of 100 replicates per coefficient in the
    promoted `fixed_by_attempt.csv`; the resulting summary coverage values of
    zero or one are therefore not meaningful Monte Carlo coverage estimates.
    The current aggregate producer nevertheless averages them. The leanest
    defensible fix is to publish bias and RMSE only for this benchmark unless a
    new inference-enabled recovery run is authorized.
11. **Four simulation figure producers do not match the manuscript figures.**
    The formal BCPE/t and NBI/Clayton fixed-effect producers currently draw
    point estimates only, whereas the manuscript figures show a fitted-bias
    range for each coefficient. The formal smooth-effect producers draw
    integrated-RMSE bar charts, whereas the manuscript figures show the fitted
    smooth functions and their between-replicate bands. Publishing the current
    files would therefore replace Figures 7, 13, 14, and 15 with plots answering
    different questions. The producers and captions must be brought back to a
    single agreed display contract before their hashes are approved.

## Audit hardening completed during this review

- Registered `tab:package_workflow` as a non-data static table, giving 45/45
  classified assets for the latest available TeX.
- Changed the manuscript audit to write its CSV report before failing, so a
  failed audit no longer leaves a stale success report.
- Changed the audit to fail if a retired artifact is still referenced. Against
  the latest source it now reports all eight retired references explicitly.
- Corrected the Overleaf checklist: generated ACTG files contain complete table
  environments, so the entire inline table must be replaced rather than nesting
  the generated file inside an existing environment.
- Corrected the BCPE/t generated-table labels to the labels registered in the
  manifest and used by the manuscript. A regression test now enforces the
  label contract.
- Replaced checkout-absolute paths in the Module 07 and central Phase 2 audit
  CSVs with output-relative paths. This removes two avoidable cross-checkout
  hash differences and prevents local directory names leaking into the public
  audit bundle.
- Registered the four ACTG random streams in `paper/seeds.csv` and added tests
  for their presence.

## Clean reviewer run

A fresh local clone at package commit `a5510d650ab0faa5e0d1461618d35eca77f87bca`
was confirmed to contain no target store and no `results` directory before
running:

```text
Rscript paper/replicate.R --profile paper --workers 2
```

The run completed successfully in 39 minutes 6 seconds on Windows R 4.4.1
with two workers. All 27 targets completed from an empty store; the checkout
was clean, the checked-out package source was installed, and the target event
ledger contains zero warnings or errors. The run recorded 530 input hashes,
129 output hashes, and 23/23 passing figure dimension/nonblank checks.

The ACTG warning ledger contains six expected, structured events rather than
leaked target warnings: four segmented-likelihood notices, one analytical to
numerical Hessian fallback, and the final model's documented
`fitted_gradient_too_large` inference caution. The last item is a substantive
interpretation caveat and must remain visible in the manuscript rather than
being described as an execution failure.

Comparison with the preceding local paper run gave identical hashes for 125
of 129 outputs, including every manuscript table and figure. The four expected
nonsemantic differences were two audit CSVs containing checkout-absolute
paths, the run manifest timestamp, and the saved ACTG model object's elapsed
timings. Loading the two model objects showed that only
`convergence$elapsed_sec` and `rs_block_trace$elapsed_sec` differed. These
volatile fields should either be canonicalized or explicitly excluded from
cross-run exact-hash expectations; manuscript-facing outputs are unaffected.

## Publication and compile check

The asset publisher staged 36 allowlisted public assets into a fresh clone: 14
updates to existing simulation/software assets and 22 new ACTG or Phase 2
assets. No source, bibliography, style, commit, or remote state was changed by
the publisher.

For diagnosis only, the closest available ACTG TeX attachment was then placed
in that disposable clone and compiled with the staged assets. It produced a
47-page PDF. The initial compile exposed two BCPE/t table-label mismatches;
these were corrected in the producer and the regenerated document resolved
those table references. The remaining compile issues are six bibliography keys
missing from the stale canonical `sample.bib` and a reference to the retired
`sec:sim_jvs` section. The exact 49-page benchmark PDF is newer than the
available TeX attachment, so a final exact-source compile remains impossible
until that source and bibliography are committed to the canonical paper
repository.

## Close-out order

1. Put the latest ACTG `main.tex` and bibliography in the canonical paper
   repository.
2. Remove or explicitly reactivate the eight retired analyses. Removal is the
   recommended lean submission choice.
3. Replace the seven inline ACTG tables with generated `\input{...}` calls.
4. Decide whether to remove coverage from the recovery tables (recommended
   lean path) or authorize a new inference-enabled simulation run, then correct
   the four main-recovery table producers and four figure producers against an
   agreed manuscript display contract.
5. Complete the reproducibility section and resolve visible draft markers and
   cross-references.
6. Replace stale local-TeX assertions with source-independent scope contracts.
7. Run the manuscript audit, publish allowlisted assets into a clean clone,
   compile, and visually inspect the exact candidate PDF.
8. Approve and freeze hashes for the final 29 public manuscript assets, then
   repeat the clean `paper` profile and publisher/compile check.
