# JSS reproducibility and privacy review

**Review date:** 2026-08-27 (Australia/Sydney)  
**Frozen manuscript:** `main.pdf` / `main.tex` at manuscript repository commit `68c3bad26626ce7c267bd330364cfda8df7a6b76`  
**Frozen file digests:** `main.tex` SHA-256 `b34dfe7e8ed13b6e12a1868d5b9ff330d4e6a317257dafb3d04fcfcfb27648ee`; `main.pdf` SHA-256 `3c4e4cfe08e65cd79dcd8e95581d933267f62716de4165b5bc00cfc45f6609b9`  
**Package working tree inspected:** Git HEAD `7343d06c986b81e2dbbbf5c2bc2732c3fe539c41`, with extensive modified and untracked replication work  
**Scope:** Independent review of the frozen manuscript and current smoke/paper/full workflows, manifest, bootstrap, target graph, public-derived inputs, LIPID boundary, CI, and run logs. Earlier audits and review reports were not consulted.

## Recommendation

**Do not submit the frozen manuscript to JSS in its current form; major reproducibility revision is required.** The repository now contains the core of a credible JSS replication system, but it is not yet an immutable, cloneable release, its full recomputation has not completed, and the frozen manuscript describes an older workflow in ways that are materially false. These are release and documentation blockers rather than evidence that the statistical package is intrinsically irreproducible.

Conditional acceptance from a reproducibility perspective would be reasonable after: (1) all replication files and public inputs are committed in a tagged/archived release; (2) a clean `paper` run and complete `full` run pass on a declared platform; (3) the manuscript states the three profiles and LIPID exception truthfully; (4) provenance includes every executed script; and (5) artifact verification is strengthened and aligned with what “exact” means.

## Five strengths

1. **Strong public/private graph boundary.** `paper/_targets.R:1-83` does not source the LIPID or RAND modules, and generated reviewer manifests are filtered to public, active rows at `paper/R/replication-helpers.R:367-378`. A scan of the 45 staged public-derived files (about 40.5 MB) found simulation metrics and no LIPID/RAND records or obvious clinical identifiers.

2. **Useful three-level reviewer design.** `paper/replicate.R:11-39` implements `smoke`, `paper`, and `full`, with separate stores at `paper/_targets/<profile>` and outputs at `results/jss-replication/<profile>`. Long-running modules use checkpoint/resume patterns.

3. **Good structured run evidence.** `paper/R/replication-helpers.R:718-813` records Git state, R/platform, lock hash, workers, seeds, target timings/events, input/output hashes, graph vertices/edges, reviewer summary, and hashes of the provenance logs.

4. **Appropriate distinction between exact text and stochastic recomputation.** Public TeX/CSV assets use newline/BOM-canonical SHA-256, while the full profile compares recomputed statistical metrics against registered absolute/relative tolerances (`paper/tolerances.csv:2-8`; `paper/R/replication-helpers.R:657-692`).

5. **Promising portability and clean-install design.** The bootstrap restores `paper/renv.lock` and installs the checked-out package source (`paper/bootstrap.R:23-42`). A Windows/Ubuntu, R 4.4.1/current, smoke/paper matrix is defined at `.github/workflows/paper-replication.yaml:15-44`. The local status file reports one successful clean-snapshot Windows smoke run.

## Already-implemented repository fixes versus stale manuscript wording

The following improvements already exist in the working tree but are absent from the frozen manuscript:

- `paper/REVIEWER.md:3-17` gives correct command-line invocations for smoke, paper, and full, target runtimes, and an explicit LIPID exception.
- `paper/LIPID-REPRODUCIBILITY-WORDING.md:1` provides substantially better exception wording than the manuscript.
- `paper/manifest.csv` classifies public, private, and non-data-static assets and assigns verification policies.
- `paper/bootstrap.R`, `paper/seeds.csv`, `paper/tolerances.csv`, the public-derived bundle, full-profile tolerance validation, and the CI workflow are new or revised working-tree components.
- The current paper log shows a successful public-asset build with 16 recorded targets and no target warnings/errors; the current smoke log records nonconvergence as structured events rather than hiding it.

These fixes should not be credited to the frozen manuscript until they are committed, archived, and described accurately in the submitted source.

## BLOCKING issues

### B1. The replication bundle is not distributable from an immutable repository state

- **Exact location:** Package Git HEAD `7343d06c...`; `git status`; key paths `.github/workflows/paper-replication.yaml`, `paper/REVIEWER.md`, `paper/bootstrap.R`, `paper/seeds.csv`, `paper/tolerances.csv`, `paper/R/public-paper-producers.R`, and all `paper/data/public-derived/**`.
- **Evidence:** The listed paths are untracked (`??`), while `paper/replicate.R`, `paper/_targets.R`, `paper/renv.lock`, `paper/manifest.csv`, and many producer/package files are modified. `git ls-files` confirms that the staged public-derived files, bootstrap, reviewer guide, and CI definition are absent from HEAD. A reviewer cloning the public Git URL in `DESCRIPTION:33-35` would not receive the workflow that was reviewed.
- **Consequence:** No commit currently binds package source, replication code, public inputs, manifest hashes, and manuscript commit `68c3bad...`. Local success cannot be reproduced by an external reviewer.
- **Remedy:** Commit the complete public replication bundle, remove unrelated/local-only state, tag it, and create a preservation archive with a DOI. Record both the package release commit and manuscript commit in the manuscript and manifest.
- **Completion test:** From a new machine or disposable VM, clone only the declared tag/archive, verify a published top-level SHA-256 manifest, and complete the documented `paper` workflow without copying any local files.
- **Confidence:** High.
- **Classification:** **BLOCKING — release-state defect; repository fixes exist locally but are not released.**

### B2. The frozen manuscript's replication instructions and scope claims are materially false

- **Exact location:** Frozen `main.tex:1595-1616`; implementation `paper/replicate.R:11-16`, `paper/_targets.R:34-60`; improved guide `paper/REVIEWER.md:3-17`.
- **Evidence:** `main.tex:1597` says “All results” are directly replicable and that the default `source("paper/replicate.R")` generates all tables and figures. The default is `smoke`, which generates only three workflow figures plus coverage/convergence audit outputs. `main.tex:1603-1608` calls `expanded` the “full run,” but `paper/replicate.R:12-15` deprecates `expanded` and maps it to `paper`, not `full`. `main.tex:1610` leaves “X hours” in red and claims exact hash matching; `main.tex:1616` leaves all platform/version/core/seed metadata as a red instruction. The same section simultaneously includes a private clinical example, contradicting “All results.”
- **Consequence:** A JSS reviewer following the paper will run the wrong profile, will not reproduce the claimed asset set, cannot estimate review time, and is misled about the private-data exception and the meaning of hashes.
- **Remedy:** Replace the entire section with the implemented three-profile contract. State that smoke is an execution check, paper is the quick public-asset rebuild from versioned analysis-ready simulation inputs, and full is the long public Monte Carlo recomputation. State the LIPID exception before any universal reproducibility claim, give measured runtimes and commands, and define verification precisely.
- **Completion test:** A context-free reviewer follows only the revised manuscript, runs each named command, and the resulting manifest contains exactly the asset subset promised by the prose.
- **Confidence:** High.
- **Classification:** **BLOCKING — stale manuscript wording; the repository guide already largely fixes it.**

### B3. The required clean public reproduction and full recomputation have not yet been demonstrated end to end

- **Exact location:** `paper/verification-status.md:6-38`; current `results/jss-replication/paper/logs/run-metadata.csv`; current `results/jss-replication/full/`; `.github/workflows/paper-replication.yaml:15-44`.
- **Evidence:** The preserved current paper run reports `git_state=dirty`, `restored=FALSE`, and `source_installed=FALSE`; it therefore is not evidence of the reviewer bootstrap path. `paper/verification-status.md:29-38` explicitly states that the full run, no-repeat resume test, acceptance of the 2,275 metric comparisons, remote CI, and live-manuscript update remain pending. At review time the full store had only `settings`, `manifest`, `public_input_files`, and `public_workflow_figures` completed and 12 worker processes were still active. The status file reports a clean Windows smoke, not a preserved clean paper/full result.
- **Consequence:** The repository has not yet shown that all public claims can be rebuilt from source on at least one clean platform, nor that the long workflow reaches its tolerance gate. A quick derived-input build cannot validate its own upstream Monte Carlo lineage.
- **Remedy:** Complete full from an empty store, then rerun to prove checkpoint resume; run paper from a clean clone and empty library/cache with restore and source installation enabled; execute the CI matrix; archive all logs and artifacts against the release commit. If full exceeds one hour, state its measured duration and make paper the explicit similar-results reviewer path.
- **Completion test:** Archived logs show a clean Git tag, `restored=TRUE`, `source_installed=TRUE`, zero execution errors, a passing full-tolerance audit, a successful no-repeat resume, and a passing paper build on at least one platform; CI corroborates portability.
- **Confidence:** High.
- **Classification:** **BLOCKING — incomplete verification gate, not merely stale wording.**

## MAJOR issues

### M1. The LIPID boundary is privacy-safe, but the public description overstates the available analysis recipe

- **Exact location:** Frozen `main.tex:1294-1593` and `1597-1612`; `paper/REVIEWER.md:17`; `paper/R/05-application-lipid.R:1-65`; private rows in `paper/manifest.csv:22-33`.
- **Evidence:** The public target graph excludes LIPID, and public-derived inputs contain no detected LIPID data. However, the manuscript presents eight LIPID figures, four inline numeric tables, exact sample/missingness summaries, and a many-covariate fitted model. The “sanitized” recipe requires only `subject`, `time`, `response`, and `treatment`, and its formulas use treatment/time only (`05-application-lipid.R:22-30`); it is not the exact model specification shown in manuscript table `tab:covariate-estimates`. The manifest provides only static image hashes and “approved_static_provenance” labels, with no public provenance record for the inline values.
- **Consequence:** Privacy is protected, but the clinical results cannot be reproduced or even specification-checked, and reviewers may wrongly believe exact sanitized formulas/code are available.
- **Remedy:** Publish the exact privacy-safe analysis specification (all transformations, variable definitions, formulas, selection rules, seeds/sample draw, software versions, and output code) without records. State whether controlled data access is possible. Label the current four-column recipe “schematic” if retained. Have a data custodian rerun the exact release code and sign/hash an output provenance report.
- **Completion test:** A custodian rerun from the private source regenerates every LIPID table/figure hash or documented semantic equivalent; the public bundle contains the exact code/specification and signed provenance but no row-level or disclosive intermediate data.
- **Confidence:** High.
- **Classification:** **MAJOR — transparency defect at an otherwise sound privacy boundary.**

### M2. Provenance hashes and the target graph omit scripts that are actually executed

- **Exact location:** `paper/R/replication-helpers.R:740-765`; invocations at `paper/R/01-simulation-bcpe-t.R:23-32`, `paper/R/02-simulation-delaporte-clayton.R:5-12`, `paper/R/04-missingness-dropout-sensitivity.R:7-16`, and `paper/R/public-paper-producers.R:245-269`.
- **Evidence:** `control_inputs` hashes package `R/`, selected `paper/R/` files, and top-level configuration, but not `paper/scripts/final-simulations/**` or nested correlation-benchmark scripts. Those omitted scripts perform the full BCPE/t, NBI/Clayton, missingness, and correlation computations. They are invoked through `sys.source()`/helper calls rather than declared as `format="file"` targets, so changing them need not invalidate the target graph and will not appear in `input-hashes.csv`.
- **Consequence:** Two runs can report identical logged control-input provenance while executing different analysis code; stale cached targets may survive a producer-script change.
- **Remedy:** Declare every executed script/configuration as a file target or source it directly into the target dependency graph, and construct `input-hashes.csv` from actual graph/file dependencies. Include publisher and manuscript-build scripts if they are part of the claimed workflow.
- **Completion test:** Modify one byte in each producer script in a disposable clone; `tar_outdated()` must identify the correct downstream targets and the input-hash log/root release manifest must change.
- **Confidence:** High.
- **Classification:** **MAJOR — repository implementation defect in deterministic provenance.**

### M3. The seed registry is documentary, unused, and contains at least one incorrect value

- **Exact location:** `paper/seeds.csv:1-9`; `paper/R/replication-helpers.R:1-19`; `paper/R/07-gamma-copula-misspecification.R:21-50`; `paper/R/replication-helpers.R:743`.
- **Evidence:** No workflow reads `seeds.csv`; it is only included in the input-hash list. The copula-misspecification row records base seed `20360528` and says “settings seed plus 700000,” but settings seed `20260528` plus `700000` is `20960528`. Staged results contain seeds beginning `20972244`, consistent with the code, not the registry. Several other simulation scripts hard-code their own seed bases.
- **Consequence:** The human-readable seed policy can silently diverge from execution, weakening auditability and making the manuscript's requested seed report unreliable.
- **Remedy:** Make `seeds.csv` executable configuration consumed by every module, or generate it from code and test exact agreement. Remove duplicate hard-coded seed sources.
- **Completion test:** A unit test derives representative per-replicate seeds for all studies from the registry and matches the checkpoint/result seed columns exactly.
- **Confidence:** High.
- **Classification:** **MAJOR — repository implementation/provenance defect.**

### M4. Figure verification does not establish the manuscript's claimed exact match

- **Exact location:** Frozen `main.tex:1610`; `paper/tolerances.csv:3`; `paper/R/replication-helpers.R:411-448` and `694-700`; current `results/jss-replication/paper/logs/figure-verification.csv`.
- **Evidence:** Ten public figures are checked, but five current Windows paper outputs have SHA-256 values different from their approved references (the three software-workflow figures and two missingness figures). Mismatch is informational. `jss_png_nonblank()` checks only PNG signature, file size, and width/height greater than 10; it does not test pixels for blankness or visual/numeric equivalence. The manuscript nevertheless says outputs “match ... exactly using hash comparisons.”
- **Consequence:** A materially wrong but valid PNG can pass. Reviewers cannot distinguish benign renderer-byte drift from changed data, scales, labels, or plotted estimates.
- **Remedy:** Either enforce byte identity on one declared canonical platform, or use semantic figure validation: hash the plotted data/specification, check dimensions and expected panels/labels, and apply a documented perceptual/pixel threshold. Keep visual QA as an additional, not sole, gate. Revise “exact” wording accordingly.
- **Completion test:** Deliberately replace a plotted series while retaining valid PNG dimensions; validation must fail. On the canonical platform, all approved reference checks must pass under the declared policy.
- **Confidence:** High.
- **Classification:** **MAJOR — repository verification defect plus stale manuscript claim.**

### M5. The manifest is not yet a complete, accurate mapping of manuscript outputs

- **Exact location:** `paper/manifest.csv:2-41`; frozen graphic labels at `main.tex:318-320`, `872-874`, `976-978`, `1098-1100`, `1265-1267`, `1994-2014`, `2149-2158`; inline tables at `main.tex:1053-1080` and `2029-2069`.
- **Evidence:** Two public inline tables (`inline_sim_design`, `inline_jvs_design`) are `pending_externalization` with `manual_until_externalized` and disappear from the generated public manifest. The `manuscript_label` values for public figures do not match the actual frozen labels (for example `fig:intro-copulas` versus `fig:example-copulas`, `fig:bcpe-fixed` versus `fig:sim_fixed_recovery`). Private inline LIPID tables have no approved digest. Paths provide partial mapping, but the label/provenance claim is inaccurate.
- **Consequence:** Coverage audits can report complete classification while excluding presented public content or pointing to nonexistent labels. Inline numeric transcription is not mechanically checked.
- **Remedy:** Generate all data-derived tables with `\\input{}`; use actual manuscript labels; give every public result a producer, input set, output, and verification rule; give private results an explicit exception/provenance record.
- **Completion test:** A parser over the final TeX finds every `\\includegraphics`, `\\input`, and data-derived inline table exactly once in the manifest, with matching path/label/access, and no `pending_externalization` rows.
- **Confidence:** High.
- **Classification:** **MAJOR — manifest/release implementation defect.**

### M6. Source/package identity and archival metadata are missing from the manuscript

- **Exact location:** Frozen `main.tex` (no repository/install/version/commit availability statement); `DESCRIPTION:1-35`; `README.md:10-15`.
- **Evidence:** The package metadata identifies version 0.1.0, GPL-3, maintainer, website, and GitHub URL, but the manuscript does not give a source URL, package version, immutable commit/tag, archive DOI, license, or tested installation command. The manuscript commit and package commit are separate and unbound.
- **Consequence:** A future reviewer cannot know which source implements the article or retrieve the exact dependency/analysis snapshot after the branch changes.
- **Remedy:** Add a Software and Data Availability section naming package version, source archive DOI, Git commit/tag, license, supported R/platform, install command, replication entry point, and manuscript-source commit.
- **Completion test:** Every identifier resolves from a blank browser/session, archived contents hash to the declared release manifest, and the article can be mapped one-to-one to the package release.
- **Confidence:** High.
- **Classification:** **MAJOR — stale/absent manuscript metadata with a release-process dependency.**

### M7. Bootstrap isolation and cross-platform portability are not yet fully controlled

- **Exact location:** `paper/bootstrap.R:14-27`; `paper/renv.lock:2-9` and the `Formula`/`gamlss2` entries; `.github/workflows/paper-replication.yaml:20-40`.
- **Evidence:** The bootstrap comment says a restored run sees only the pinned library, but `.libPaths(unique(c(lib, inherited_libs)))` retains inherited libraries. Optional availability can therefore depend on host state; for example `glmmTMB` is queried by the correlation benchmark but is not in the lockfile. The lock has a GitHub SHA for `gamlss2`, which is good, but `Formula` is marked `R-Forge` while only CRAN is declared in the lock's repository list. Remote CI has not run.
- **Consequence:** A host-installed package or unavailable repository source can alter behavior or break restore. Local cache success does not prove a cold install.
- **Remedy:** Start the reviewer subprocess with an explicitly minimal `.libPaths()`, pin all optional packages that affect executed branches, declare every repository/source URL, and document system requirements. Test with empty renv/download caches.
- **Completion test:** Cold Windows and Ubuntu runs use only base/site libraries plus the lock-keyed library, restore without a package cache, and produce the same canonical text assets and passing statistical tolerances.
- **Confidence:** Medium-high.
- **Classification:** **MAJOR — portability risk pending CI/cold-install evidence.**

## MINOR issues

### m1. Smoke is an execution test, not evidentiary model validation

- **Exact location:** `paper/REVIEWER.md:13`; `paper/R/replication-helpers.R:199-256`; current smoke `reviewer-summary.md`.
- **Evidence:** All 24 smoke fits ended in optimizer nonconvergence under deliberately small limits, though there were no execution failures and the guide now calls smoke a short test.
- **Consequence:** Readers could mistake smoke summaries for statistical evidence if the manuscript says it reproduces paper results.
- **Remedy:** Label every smoke output “execution-only; estimates not for scientific interpretation,” and keep it out of manuscript result comparisons.
- **Completion test:** Smoke logs and documentation carry that label, while paper/full acceptance does not rely on smoke convergence.
- **Confidence:** High.
- **Classification:** **MINOR — documentation/interpretation issue; structured logging is already good.**

### m2. The reviewer summary points to a log that does not exist for paper/full

- **Exact location:** `paper/R/replication-helpers.R:774-796`; current `results/jss-replication/paper/logs/`.
- **Evidence:** When `fit-events.csv` is absent, the helper silently uses an empty data frame but still tells the reviewer to inspect `fit-events.csv`. The current paper log directory has no such file.
- **Consequence:** Small loss of trust and a broken audit trail link.
- **Remedy:** Always write an empty, schema-valid fit-events file or mention it only when present.
- **Completion test:** Every reviewer-summary link resolves after each profile.
- **Confidence:** High.
- **Classification:** **MINOR — repository logging defect.**

### m3. `input-hashes.csv` is broader than “only inputs actually read” for paper/full

- **Exact location:** `paper/_targets.R:26-31` and `43-51`; `paper/R/replication-helpers.R:752-765`; `paper/REVIEWER.md:26-32`.
- **Evidence:** All public CSV/TeX/Markdown files are declared as one file target and every non-log public-derived file is hashed, while modules read subsets. Conversely, executed scripts outside the selected control list are omitted (M2).
- **Consequence:** The log is an inventory, not a precise read-set, despite the guide's wording.
- **Remedy:** Record per-module file dependencies/read sets and describe the aggregate log as such.
- **Completion test:** Each logged input has a consumer edge, and every consumer file appears in the log.
- **Confidence:** High.
- **Classification:** **MINOR — provenance precision/documentation issue.**

### m4. The frozen manuscript source is not submission-clean

- **Exact location:** Examples include `main.tex:143-158`, `747-774`, `1031`, `1186`, `1222`, `1610-1616`, `1637`, and `1694`; frozen `main.blg:12-17`.
- **Evidence:** Red TODO/revision/placeholder material remains, including the runtime/environment placeholders. The compile log has undefined citations including `czado2015`, `Beareseo2015`, `lambert_copula`, `topmodels`, `simes`, and `marra2025?`.
- **Consequence:** Even after code fixes, the supplied source is not a stable scholarly record and obscures the final reproducibility contract.
- **Remedy:** Remove editorial macros/notes, resolve references, compile twice from a clean directory, and archive only the final source/PDF.
- **Completion test:** Clean compile with no undefined citations/references and no visible TODO/revision text.
- **Confidence:** High.
- **Classification:** **MINOR for this reproducibility review, but submission-blocking editorial cleanup.**

## Minimum JSS replication bundle

The minimum acceptable bundle should contain:

1. An immutable package source release and manuscript-source release, with tag/commit, DOI, license, and a signed/top-level SHA-256 manifest binding both.
2. Final `main.tex`, `main.pdf`, bibliography, JSS/style files, and every static publication asset required for a clean compile.
3. A single reviewer start page with tested smoke, paper, and full commands; measured wall times; hardware/cores; expected disk/memory; and output locations.
4. `DESCRIPTION`, exact source-install instructions, R version range, system requirements, `renv.lock`, and all declared package repositories/remote SHAs.
5. The complete `targets` graph and every executed producer/publisher/build script as declared file dependencies.
6. `seeds.csv` and `tolerances.csv` as executable configuration, not parallel documentation.
7. All 45 non-clinical public-derived inputs with origin script, source release commit, creation platform/session, design/replicate count, and hashes.
8. A complete manifest mapping every final figure/table/inline result to manuscript path and actual label, producer, inputs, access class, profile, and verification rule.
9. A LIPID exception statement, exact privacy-safe analysis specification, access route or reason none exists, and independently attested static-output provenance; no private records/intermediates in public commands.
10. Archived clean-run artifacts for paper and full from at least one platform, plus CI smoke/paper evidence on Windows and Ubuntu, including target events, fit events, input/output hashes, figure checks, tolerance audit, session info, and checkpoint-resume proof.

## Top-10 actions, in priority order

1. Commit/tag/archive the complete working-tree replication bundle and bind it to manuscript commit `68c3bad...`.
2. Replace frozen `main.tex:1595-1616` with the accurate three-profile and LIPID-exception text; remove all placeholders.
3. Finish the empty-store full run, pass all 2,275 comparisons, and perform the no-repeat resume test.
4. Run and archive a cold, clean paper build with lock restore/source install enabled; then run the CI matrix.
5. Add every executed `paper/scripts/**` file to the target dependency graph and provenance hashes.
6. Make the seed registry authoritative and fix the copula-misspecification seed discrepancy.
7. Strengthen figure validation; either enforce a canonical-platform hash or validate plotted data and perceptual/semantic equivalence.
8. Externalize the two pending public inline tables, correct every manifest label, and verify complete TeX-to-manifest coverage.
9. Publish the exact privacy-safe LIPID model specification and custodian-attested provenance, while explicitly retaining the no-private-data public boundary.
10. Add an immutable Software/Data Availability section and produce a clean manuscript build with resolved citations and no editorial annotations.
