# Reviewer 08 — writing, references, and JSS style

**Review date:** 2026-08-27  
**Recommendation:** **Do not submit in the present form; major revision is required before JSS review.** The paper has a credible software contribution and a useful applied story, but the frozen artifact is administratively non-compliant with JSS, bibliographically broken, and still contains unresolved scientific/editorial decisions. These are pre-submission blockers rather than matters that should be left to copy-editing.

## Scope and review basis

I independently reviewed the frozen 47-page `main.pdf` (31 pages before the appendices), `main.tex`, `main.log`, `main.blg`, `main.bbl`, and `sample.bib` in the supplied review bundle. I visually inspected every PDF page and did not consult prior reviews or internal audits.

Frozen file identifiers:

- `main.pdf`: SHA-256 `3C4E4CFE08E65CD79DCD8E95581D933267F62716DE4165B5BC00CFC45F6609B9`
- `main.tex`: SHA-256 `B34DFE7E8ED13B6E12A1868D5B9FF330D4E6A317257DAFB3D04FCFCFB27648EE`
- `sample.bib`: SHA-256 `4F6E5D0385E9FF14F46D87CE178E380ACA616BC1A49C9D898521DB09D404A02D`

The governing guidance was the current [JSS Style Guide](https://www.jstatsoft.org/style) and [JSS Submission Guide](https://www.jstatsoft.org/guides/submission). In particular, JSS requires its LaTeX style, `\proglang`/`\pkg`/`\code`, title style for the article and BibTeX titles, sentence style for headings/captions, labelled and referenced artifacts, captions below both figures and tables and ending in a period, properly cited software, code within text width without comments in verbatim blocks, and exact reproducibility of displayed output.

Red revision notes were treated as known draft context. The issues below assess the unresolved decisions and claims that the notes expose, not merely their red appearance.

## Five strengths

1. **The software purpose is identifiable.** The package is positioned around a clear practical gap: distributional longitudinal regression with covariate-dependent adjacent-time copula dependence and marginal interpretation.
2. **The paper contains a real user workflow.** Sections 3.2.1–3.2.5 move from simulation through fitting, screening, standard methods, diagnostics, prediction, simulation, and bootstrap inference, which is appropriate for a JSS software article.
3. **Limitations are not wholly hidden.** The first-order dependence restriction, runtime cost, exchangeable-correlation weakness, and missingness caveats are discussed explicitly (especially `main.tex` lines 1204–1214, 1224–1251, 1282–1290, and 1624–1628).
4. **The validation agenda is broad and user-oriented.** Parameter recovery, optimizer choice, dependence misspecification, copula misspecification, and missingness are connected to practical guidance tables rather than presented only as algorithmic benchmarks.
5. **Artifact scaffolding is mostly present.** Every figure/table environment in `main.tex` has a label, the figures render without clipping, and the manuscript attempts to map paper artifacts to replication outputs. The remaining reference and readability defects are therefore repairable rather than structural absences.

## BLOCKING issues

### B1. The manuscript is not in JSS LaTeX style

- **Location/evidence:** `main.tex` lines 1–3 use `\documentclass{article}` and `\usepackage{arxiv}`. The frozen PDF has an “A PREPRINT” header on every page. JSS declarations such as `\Plainauthor`, `\Plaintitle`, and `\Plainkeywords` are absent; `plainnat` is selected at lines 2192–2193.
- **Consequence:** JSS states that its style files are required and that formally non-conforming manuscripts can be returned before review. Pagination, title block, headers, abstract/keywords, code environments, citations, and bibliography cannot be meaningfully certified while the wrong class remains in place.
- **Remedy:** Rebase the source on the current JSS article template (`\documentclass[article]{jss}`), use its front-matter declarations and bibliography setup, remove `arxiv.sty` and redundant/incompatible packages, and compile with pdfLaTeX.
- **Completion test:** A clean build uses the current `jss` class, contains no “A PREPRINT” furniture, has all required JSS front matter, and visually matches the official template.
- **Confidence:** High.
- **Classification:** Formal JSS compliance.

### B2. Citations and the compiled bibliography are broken

- **Location/evidence:** `main.log` reports unresolved citations `czado2015`, `Beareseo2015`, `lambert_copula`, `topmodels`, `simes`, and `marra2025?` (log lines 883–898, 944, 998, and 1121). `main.blg` contains **72 warnings**, including the six missing keys, three unsupported entry types, and repeated missing `journal`/`year` fields. The cause is visible in `sample.bib`: many entries use biblatex fields such as `journaltitle` and `date` while the document runs BibTeX with `plainnat`. Consequently, the PDF reference list omits journals and years for many articles; for example, Aas et al., Czado and Nagler, Liang and Zeger, and Smith et al. appear only with volume/issue/pages. `Rigby2005` compiles as volume 54(3), pages 1–2, with no journal or year.
- **Consequence:** Readers cannot reliably identify sources, claims appear unsupported, package citations do not meet JSS policy, and the build is not publication-ready.
- **Remedy:** Build a small, clean JSS-compatible `.bib` containing only cited records; repair the six missing keys; convert fields to the bibliography engine actually used; verify every compiled record against a primary source/DOI; use title-style titles and JSS markup inside BibTeX where appropriate.
- **Completion test:** BibTeX and LaTeX complete with zero missing-entry, empty-required-field, or undefined-citation warnings; every in-text citation resolves; every compiled article has author, year, article title, journal, volume/issue/pages or article number, and DOI/URL as appropriate.
- **Confidence:** High.
- **Classification:** Bibliographic integrity.

### B3. Material scientific/editorial choices remain unresolved

- **Location/evidence:** The source contains 12 `\todo` calls, nine `\jssrevisionnote` calls, and one `\jssplaceholder`. They concern substantive matters: how to describe covariate effects in the abstract (line 143), introduction of GAMLSS (159), missing support for a GLMM claim (277), a possibly future extension (281), the interpretation of GLMM estimates (327), model notation (338, 344), optimizer detail (580, 853, 1694), software-workflow content (747–755, 774), simulation emphasis (1031, 1186, 1222), and appendix scope (1637). These notes visibly occupy pages 1, 3–6, 9, 11–13, 17, 20–21, and 31–33.
- **Consequence:** The draft does not yet have a stable argument or stable description of the implemented method. Removing the colored boxes alone would conceal unresolved claims.
- **Remedy:** Resolve each decision against the implemented software and evidence, revise the prose/notation, then remove all annotation machinery (`todonotes`, custom revision macros, and colored draft text).
- **Completion test:** Searches for `\\todo`, `\\jssrevisionnote`, `\\jssplaceholder`, draft-only red text, “Next steps,” and placeholder language return zero; an independent read finds no decision deferred to future drafting.
- **Confidence:** High.
- **Classification:** Manuscript completeness.

### B4. The reproducibility claim is stronger than the documented protocol

- **Location/evidence:** The abstract promises a “reproducible worked example” (line 143), and Section 6 says all results can be replicated (1597–1614), but the manuscript leaves runtime as “approximately X hours” and explicitly lacks OS, R version, package versions, cores, seed policy, and runtime (1610, 1616). The primary main-text code does not show package installation/loading or a pinned package version. Appendix code calls the internal API `gamlss.longitudinal:::get_copula_dist()` (2075). The clinical data are not redistributable (1294, 1612), yet figures from that analysis are mixed into the article and bibliography sequence.
- **Consequence:** A reviewer cannot know what environment reproduces the frozen numbers, whether the public workflow uses stable APIs, or which artifacts are exempt because of data restrictions. This conflicts with JSS’s requirement that displayed output be fully and exactly reproducible on at least one platform.
- **Remedy:** Document a tested, version-pinned clean-session command; OS/R/package versions; CPU/cores; seeds and parallel RNG policy; runtime and storage; expected outputs; hash verification; and the exact package commit/release. Replace `:::` use with an exported/stable interface or justify it as package-internal replication code. Explicitly separate reproducible simulated artifacts from non-redistributable clinical artifacts and provide a lawful synthetic or processed substitute where possible.
- **Completion test:** A clean-machine run from the archived submission regenerates every reproducible table/figure and passes the stated hashes; the manuscript reports the tested environment and measured runtime; the manifest clearly marks any legally non-reproducible clinical outputs and their provenance.
- **Confidence:** High.
- **Classification:** Reproducibility and claim accuracy.

## MAJOR issues

### M1. The title and abstract do not yet do the full editorial job

- **Location/evidence:** The title at line 114 capitalizes the prepositions “With” and “in” contrary to JSS title style and does not apply JSS markup to the package or R. A compliant form would be `\pkg{gamlss.longitudinal}: Longitudinal GAMLSS Models with Copula Dependence in \proglang{R}`. The abstract (line 143) uses GAMLSS without expansion, inventories features, contains a live terminology dispute (“fixed or smooth”), reports no quantitative validation result, omits the first-order limitation, and calls the workflow reproducible before the protocol is complete.
- **Consequence:** A reader learns what the package exposes but not the strength, boundary, or cost of the evidence. The title/abstract also fail visible JSS conventions.
- **Remedy:** Expand GAMLSS at first use; state the problem, implemented approach, principal distinction, one or two quantitative validation/runtime results, first-order-dependence boundary, and exact software/reproduction availability. Remove generic feature inventory and unverified adjectives.
- **Completion test:** The title passes JSS title-style rules and markup; the abstract is self-contained, contains no undefined acronym or draft note, includes at least one concrete result and one limitation, and every reproducibility claim is demonstrably true.
- **Confidence:** High.
- **Classification:** Front matter and narrative framing.

### M2. The introduction and theory repeat themselves and delay the software story

- **Location/evidence:** Lines 323–325 and 330–332 repeat nearly the same two paragraphs on multivariate decomposition, Sklar’s theorem, vines, applications, and higher-order dependence. The introduction first summarizes competitors in Table 1, then re-summarizes them in lines 275–285, then gives another contribution summary at 289–306. Main text reaches the software interface only on PDF page 11. The paper is 31 pages before appendices and 47 pages overall.
- **Consequence:** The contribution appears less crisp than it is, readers encounter duplicated literature claims, and the software article reads initially as a methods survey. The extra length also works against JSS’s readability guidance for manuscripts over 30 pages.
- **Remedy:** Delete the duplicate block; reduce competitor discussion to the decision-relevant contrasts in Table 1; introduce GAMLSS, copulas, and first-order PCC once; move derivations and secondary simulations to a supplement; bring the package contract, installation/version, and minimal example earlier.
- **Completion test:** No paragraph-level duplication remains; each background concept is introduced once; a reader reaches a minimal working package example within the first third of the main paper; the main text has a clear problem → gap → software → validation → application sequence.
- **Confidence:** High.
- **Classification:** Narrative coherence and length.

### M3. Several novelty and comparative claims are too absolute or insufficiently sourced

- **Location/evidence:** Line 287 asserts in bold that “There is currently no direct practitioner-facing software...” without a systematic search or bounded date/scope. Lines 275 and 1228 call approaches “the most common” or “most popular” without evidence. Line 277 claims “Many meta-analyses” show GLMMs are often misapplied but supplies no references. Line 327 says the two-stage approach does not provide the “most optimal fit,” conflating an estimation target with a superlative. Line 1243 says the package is “uniquely designed”; 1245 uses “substantially” twice without an effect-size summary; 1620 calls the approach “novel” although the introduction also locates closely related joint and first-order vine methods.
- **Consequence:** The tone invites avoidable challenges to novelty and fairness, especially because several competitor packages are uncited and the comparison scope is not reproducibly defined.
- **Remedy:** Date- and scope-bound the software search; replace absolutes with precise capability comparisons; cite every comparative factual claim; distinguish joint maximum penalized likelihood from two-stage estimation without “optimal”; pair qualitative adjectives with effect sizes and uncertainty.
- **Completion test:** Every market/novelty/prevalence claim is either supported by a reproducible comparison/citation or rewritten as a bounded statement; no unsupported “no software,” “most popular,” “uniquely,” or “most optimal” wording remains.
- **Confidence:** High.
- **Classification:** Tone, evidence, and overclaiming.

### M4. JSS software and code markup is not implemented

- **Location/evidence:** `main.tex` redefines `\pkg` as bold at line 98, contains 65 `\pkg` uses, 91 `\texttt` uses, and **zero** `\proglang` or `\code` uses. R is mostly plain text; function calls, arguments, file paths, and values are set with `\texttt`. Package names are occasionally plain (`gamlss.longitudinal`, `gamlss`, `gamlss2`).
- **Consequence:** The manuscript fails an explicit JSS checklist item, semantics are inconsistent, and later conversion to the JSS class may create heading/bookmark problems.
- **Remedy:** Use the JSS macros rather than redefining them: `\proglang{R}` for languages, `\pkg{...}` for packages, and `\code{...}` for functions, arguments, literal values, commands, and paths. Supply plain bookmark/title variants where markup occurs in headings.
- **Completion test:** No local redefinition of JSS macros remains; a scripted audit finds consistent semantic markup throughout text, headings, captions, tables, and BibTeX titles; plain package/function names occur only where deliberately exempted.
- **Confidence:** High.
- **Classification:** JSS house style and terminology.

### M5. Code input/output presentation conflicts with JSS guidance

- **Location/evidence:** Twelve `lstlisting` blocks are used (e.g., lines 797–851, 861–920, 934–953, 965–1028, 1599–1608, 2024–2027, and 2073–2143). Printed output is encoded as R comments in the same block (lines 814–850, 888–896, 913–919, 992–999, 1016–1027), and blocks contain prose comments such as “# Standard model outputs” (935). JSS explicitly asks authors not to place comments inside verbatim code and provides separate `CodeInput`/`CodeOutput` environments. The appendix listing spans pages 40–41 without a short executable wrapper or stated package-loading context.
- **Consequence:** Readers cannot distinguish commands from console output; copied examples are not directly executable; the custom styling must be discarded during JSS conversion.
- **Remedy:** Convert to JSS `CodeInput`/`CodeOutput` (or supported Sweave/knitr) environments; move explanations into prose; show only output needed for the argument; use JSS prompts/continuation conventions and width 70; include a compact, copy-pasteable setup and workflow.
- **Completion test:** Inputs and outputs are visually and semantically distinct; no explanatory/output comments remain in verbatim input; every displayed line fits the JSS text width; the minimal example runs exactly as printed in the documented environment.
- **Confidence:** High.
- **Classification:** Code communication and reproducibility.

### M6. Caption, table-note, and artifact-reference practice is not JSS-compliant

- **Location/evidence:** Table captions are placed above tables throughout (for example lines 170–171, 1052–1053, 1121–1122, and 1488–1489), whereas current JSS guidance puts captions below both figures and tables. Table 1 uses superscript footnotes outside the caption (263–269), and other tables use `tablenotes` (1368–1370, 1473–1480, 1557–1561), although JSS asks for such annotations in captions. Eight captions lack final punctuation: lines 1266, 1325, 1585, 2033, 2166, 2173, 2180, and 2187. Ten labelled artifacts are never referred to in prose: `fig:plot_dist`, `fig:plot_copula_diagnostics`, `fig:lipid-copula-by-treat`, `fig:lipid-term-plot`, `fig:cooccurence`, `fig:sim_fixed_recovery_nb`, `fig:sim_smooth_recovery_nb`, `tab:marginal-screen`, `tab:nbi-fit-characteristics`, and `tab:nbi-parameter-recovery`. The label for `tab:marginal-screen` precedes its caption (1324–1325), producing a `\caption@xref` value in `main.aux` rather than table number 12.
- **Consequence:** This violates multiple explicit JSS checklist items; some captions are not self-contained; the broken label will yield an incorrect reference if used later.
- **Remedy:** Put every caption below its artifact; write sentence-style, self-contained captions ending in periods; move table notes into captions; place labels immediately after captions; and explicitly discuss every retained artifact in the text. Remove redundant artifacts rather than adding perfunctory references.
- **Completion test:** A source audit finds caption → label order for every artifact, no table footnote annotations, all captions end with periods, and every retained figure/table label is referenced at least once with `Figure~\ref` or `Table~\ref`.
- **Confidence:** High.
- **Classification:** JSS captions and cross-references.

### M7. Several figures/tables are too small or insufficiently self-explanatory

- **Location/evidence:** Table 1 is explicitly `\tiny` (lines 165–270; PDF page 2) and is extremely dense. Figure 3 (page 16) and especially the nine-panel diagnostic Figure 16 (page 45) use annotations materially smaller than the caption. Figure 17 (page 46) has many miniature term labels and confidence intervals. Tables 8–10 on pages 21–23 and Table 17 on page 26 are difficult to scan at normal page size. Captions such as “Standard marginal distribution fit diagnostics” (2166), “Joint distribution diagnostics” (2173), and “Term plot for final model for cholesterol data” (2180) do not identify the fit, sample, panels, intervals, or intended diagnostic conclusion.
- **Consequence:** Readers cannot inspect evidence without extreme zoom, and graphics fail JSS’s pen-to-paper guidance. Sparse captions make the problem worse.
- **Remedy:** Split or simplify dense artifacts; enlarge annotation text to approximately caption size; remove repeated/low-value panels; move full diagnostic grids and detailed tables to supplementary material; write captions that identify data/model, panel order, metric/interval, and takeaway.
- **Completion test:** At 100% viewing and on an A4/letter print, all labels and values are legible without zoom; captions let a reader interpret each artifact without searching the main text; no table uses `\tiny` or `\scriptsize` as a layout rescue.
- **Confidence:** High.
- **Classification:** Visual communication and accessibility.

### M8. The reference list is interrupted by appendix floats

- **Location/evidence:** “References” begins midway down PDF page 42 under Figure 13. References 9–14 share page 43 with Figure 14; references 15–21 share page 44 with Figure 15; Figures 16–17 occupy pages 45–46; references 22–26 resume below Figure 18 on page 47. This arises because bibliography commands at lines 2192–2193 are reached while appendix floats from lines 2147–2189 remain pending.
- **Consequence:** The bibliography is fragmented across diagnostic figures and is hard to navigate; it visibly signals unfinished float control.
- **Remedy:** Resolve float placement and force all appendix floats before the bibliography (normally a deliberate `\clearpage`/JSS-compatible arrangement), then keep the reference list continuous at the end.
- **Completion test:** The rendered PDF contains one uninterrupted References section after the final figure/table, with no floats inserted between entries.
- **Confidence:** High.
- **Classification:** Document structure and typesetting.

### M9. The clinical narrative needs more cautious, reproducible reporting

- **Location/evidence:** Lines 1294 and 1593 say the analysis establishes “previously unidentified” non-mean relationships and “substantially more nuance.” Model building is described as sequential covariate selection with repeated likelihood comparisons and reduction of non-significant terms (1374–1484), followed by unqualified significance language (1565–1593). The data cannot be redistributed, the `simes` source is unresolved, and no protocol is given for the 10% random sample beyond its size (1301–1307).
- **Consequence:** Exploratory, data-driven findings can be mistaken for confirmatory clinical inference; readers cannot reproduce sampling or selection; novelty relative to the parent clinical literature is unsupported.
- **Remedy:** Label the example explicitly exploratory; document sample seed/eligibility and missing-data handling; distinguish pre-specified from selected terms; report uncertainty and selection caveats; avoid clinical novelty/causal language; repair the trial citations and provide a data-access/provenance statement.
- **Completion test:** The section can be read without inferring confirmatory or causal claims, the sampling/model-selection path is fully documented, and all trial/background claims resolve to complete references.
- **Confidence:** Medium-high.
- **Classification:** Applied-reporting tone and provenance.

## MINOR issues

### N1. Headings and PDF bookmarks need sentence style and a coherent hierarchy

- **Location/evidence:** “Statistical Model” (308), “RS Algorithm” (1639), “CG Algorithm” (1690), “Numerical Hessian Used...” (1791), and several derivative headings (1826–1922) use title-style capitalization. Line 1790 jumps from `\section` directly to `\subsubsection` at 1791. `main.log` reports a bookmark-level jump and two “Token not allowed in a PDF string” warnings (1125–1133).
- **Consequence:** Headings fail JSS sentence-style rules and PDF navigation is degraded.
- **Remedy:** Normalize sentence case, restore a section/subsection/subsubsection hierarchy, and provide plain bookmark text for marked-up headings.
- **Completion test:** Bookmark warnings are zero and the PDF outline mirrors the visible heading hierarchy.
- **Confidence:** High.
- **Classification:** Heading style and accessibility.

### N2. Cross-reference prose is inconsistent

- **Location/evidence:** Lowercase `figure`, `table`, `section`, and `appendix` precede references at lines 551, 577, 768, 789, 1041, 1094, 1114, 1170, 1182, 1188, 1234, 1239, 1241, 1286, 1294, 1612, and 1824. Many omit the nonbreaking `~`. Line 1045 hard-codes “Section 3.”
- **Consequence:** References are typographically inconsistent and vulnerable to renumbering/line breaks.
- **Remedy:** Use `Section~\ref`, `Figure~\ref`, `Table~\ref`, and `Appendix~\ref` consistently; never hard-code a section number.
- **Completion test:** A regex audit finds no hard-coded section numbers or lowercase artifact names before `\ref` and no space in place of `~`.
- **Confidence:** High.
- **Classification:** Cross-reference mechanics.

### N3. Copy-editing is required throughout

- **Location/evidence:** Examples include “pacakge” (148), “sentenes” and “GAMLSSrandom” (277), “paramterized” (285), “reproducibiltiy” (306), “alreday” (330), “simplications” (348), “avialblae” (748), “a nother” (1186), “treament” (1565), “paramater” (1567), “dypnoeoa” (1593), “oragnized” and “pacakage” (1612–1616), and “currnetly” (1694). There are recurrent comma splices, “however” joins, doubled words (“the the,” line 323), and overlong sentences (notably 159, 279, 327, 350, 1103, 1241, and 1591).
- **Consequence:** Errors obscure already technical material and lower confidence in exact mathematical/software claims.
- **Remedy:** Perform a substantive line edit after structural revision, then a separate proofread/spellcheck.
- **Completion test:** Automated spelling/grammar review and two human proofreads find no known errors; sentences longer than roughly 40–50 words are individually justified or split.
- **Confidence:** High.
- **Classification:** Language quality.

### N4. Spelling, hyphenation, and statistical style are mixed

- **Location/evidence:** The paper alternates British and US forms (`modelling`/`modeling`, `optimisation`/`optimization`, `behaviour`/`behavior`, `summarise`/`summarize`), and varies `timepoint`/`time point`, `t copula`/`t-copula`/`Student's t copula`, `mis-specification`/`misspecification`, and `log likelihood`/`log-likelihood`. JSS house style calls for `$p$~value`, but the manuscript uses “p-value” (981) and table header `$p$`/“Sig.” conventions without consistent prose.
- **Consequence:** Terminology looks unstable and makes search/cross-reading harder.
- **Remedy:** Adopt the style sheet below and apply it mechanically after content stabilizes.
- **Completion test:** Search results show one approved form for each listed term and consistent statistical notation.
- **Confidence:** High.
- **Classification:** Terminology consistency.

### N5. Notation and parameter interpretation need a final audit

- **Location/evidence:** `$n$` is defined as subjects and `$T$` as time points at 1090, but line 1177 switches to `$N$` for sample size. Line 364 says time points are `t=1,...,n` before using `T`. The prose alternates dependence, covariance, correlation, rank correlation, and Kendall's `\tau`; Table 15 calls `\zeta` “tail correlation,” while earlier it is the t-copula degrees-of-freedom parameter (793). “Fixed or smooth covariates” confuses covariates with effect forms (143 and elsewhere).
- **Consequence:** Readers can misinterpret what is modeled and what the reported coefficients mean.
- **Remedy:** Reserve `n` for subjects, `T` for scheduled time points, and (if needed) `N` for observations; distinguish covariates from linear/smooth effects; distinguish the copula parameter, Kendall's `\tau`, Pearson correlation, and tail-dependence coefficients; describe `\zeta` on its actual link/parameter scale.
- **Completion test:** A notation table or first-use definitions give each symbol/term one meaning, and all tables/captions use those meanings.
- **Confidence:** High.
- **Classification:** Technical terminology.

### N6. The build is not warning-clean even beyond citations

- **Location/evidence:** `main.log` contains 17 overfull `\hbox` reports (15 associated with the Table 15 alignment at source line 1556), 12 underfull `\hbox`, 12 underfull `\vbox`, float-placement changes, a removed `siunitx` option (`detect-inline-weight`), bookmark warnings, and undefined citations. Packages are duplicated (`amsmath`, `graphicx`, `booktabs`, `float`, `siunitx`) at lines 13–38.
- **Consequence:** Some warnings are benign, but the volume masks genuine layout/reference defects and complicates JSS conversion.
- **Remedy:** Remove duplicate/unneeded packages, update obsolete options, repair Table 15 column sizing, and review every remaining warning after switching to the JSS class.
- **Completion test:** Final pdfLaTeX/BibTeX build has zero undefined references/citations and zero overfull boxes; any residual underfull/float warnings are individually inspected and documented.
- **Confidence:** High.
- **Classification:** TeX hygiene.

## Terminology and style sheet

Use this as a manuscript-wide controlled vocabulary after structural revision.

| Concept | Preferred form | Avoid / notes |
|---|---|---|
| R language | `\proglang{R}` | Plain `R` when referring to the language. |
| Package | `\pkg{gamlss.longitudinal}` | Plain or `\texttt{gamlss.longitudinal}`; preserve exact lowercase name. |
| Other packages | `\pkg{gamlss}`, `\pkg{gamlss2}`, `\pkg{geepack}`, etc. | Cite the official package citation at first substantive use. |
| Functions/arguments/values | `\code{gamlss_longitudinal()}`, `\code{method = "RS"}` | `\texttt` for semantic code. |
| GAMLSS | “generalized additive models for location, scale and shape (GAMLSS)” at first use in abstract and body | Capitalized expansion in running prose; unexplained acronym. |
| Model versus software | “the longitudinal GAMLSS model” versus “the `\pkg{gamlss.longitudinal}` package” | Calling the package a model or the method a package. |
| Effect forms | “linear, parametric, or smooth effects of covariates” | “fixed or smooth covariates.” |
| Serial structure | “first-order dependence” (adjective), “dependence of first order” (noun phrase) | Unhyphenated adjective. |
| Pair copulas | “pair-copula construction (PCC)” and “first-order truncated D-vine” | “pair copula deconstruction” unless formally defined. |
| Time | “time point” (noun), “time-point-specific” (compound adjective) | `timepoint`, mixed forms. |
| Copula | Choose “Student's t copula”; choose “Gaussian copula” consistently | Alternating `t copula`, `t-copula`, normal/Gaussian. |
| Dependence measures | “Kendall's `\tau`,” “Pearson correlation,” “tail-dependence coefficient” | Generic “correlation” when the quantity is different. |
| Misspecification | `misspecification`, `misspecified` | Mixed `mis-specification`, `mis-specified`. |
| Optimization | Use US forms consistently: `optimization`, `optimized`, `modeling`, `behavior` | Current British/US mixture. |
| Likelihood | “log likelihood” (noun), “log-likelihood criterion” (adjective) | Mixed hyphenation. |
| Sample notation | `n` subjects, `T` scheduled time points, `N` observed rows only if explicitly needed | Interchanging `n` and `N`. |
| Statistical style | `$p$~value`, `$t$~statistic`, `95\% confidence interval` | `p-value`; bare `95 percent` in one place and `%` elsewhere. |
| Cross-references | `Section~\ref{...}`, `Figure~\ref{...}`, `Table~\ref{...}`, `Appendix~\ref{...}` | Lowercase labels, normal spaces, hard-coded numbers. |
| Abbreviations | Expand once in the abstract and once in the body where necessary | GAMLSS, GEE, GLMM, PCC, PIT, QQ, IRLS, RMSE, IRMSE, GAIC without consistent first-use expansion. |
| Captions/headings | Sentence style; every caption ends in a period | Title capitalization in headings; fragment captions. |

## Bibliography repair list

1. **Resolve six missing keys.** Add/repair `czado2015`, `Beareseo2015`, `topmodels`, `simes`, and the intended GJRM missingness reference now written as `marra2025?`. Delete the duplicated paragraph at 330–332 and its wrong `lambert_copula` key, or replace that key with the existing `lambert_copula-based_2002` if the duplicate is unexpectedly retained.
2. **Convert incompatible fields.** For every cited `@article`, change or export `journaltitle` → `journal` and `date` → `year` (or use the exact engine/style supported by the JSS template). `main.blg` proves that the present combination loses these fields.
3. **Repair `Rigby2005` from the primary record.** The current raw record says 2015, pages 1–2, and no journal; it is not an acceptable citation for the foundational GAMLSS paper. Verify title, authors, year, journal, volume, issue, full pages, and DOI.
4. **Repair unsupported web/online types.** `mark_padgham_2021_5556756`, `Rpackages2e`, and `noauthor_tidy_nodate` use entry types not defined by `plainnat`. Use JSS-compatible `@Manual`, `@Book`, or `@Misc` records with authors/organizations, year, version/edition, access date where appropriate, and durable URL/DOI.
5. **Deduplicate the R Packages sources.** `Rpackages2e` and `wickham_r_2023` describe the same resource but one has no author/year and the other has malformed publisher/date handling. Retain one authoritative edition with authors, year, edition, publisher, and URL.
6. **Cite R itself and all named packages.** Add the official `citation()`/CRAN citation for R and, as substantively used, `gamlss.longitudinal`, `gamlss`, `gamlss2`, `geepack`, `glmtoolbox`, `VGAM`, `lme4`, `mgcv`, `MASS`, `GJRM`, `VineCopula`, `rvinecopulib`, `gamCopula`, `mvtnorm`, `topmodels`, and any package used to generate reported output. A method article does not automatically replace a package/version citation.
7. **Bring every BibTeX title into JSS title style.** Protect acronyms/proper software markup with braces and use `\pkg`, `\proglang`, and `\code` inside titles where appropriate. Current compiled titles are largely sentence case.
8. **Remove leaked reference-manager metadata.** `_eprint`, “Publisher:” notes, local `file` fields, `Google-Books-ID`, abstracts, keywords, and very long tracking/search URLs should not leak into the compiled bibliography. Prefer DOI URLs or stable publisher/CRAN/project URLs.
9. **Verify article metadata.** The compiled list currently omits journals/years for most items and contains suspicious details such as Smith et al. `105(492)`, Klein and Kneib with only `26:841–860`, and book author “Harry. Joe.” Check every cited record against Crossref/publisher/CRAN and normalize author names, en dashes, issue/article numbers, and DOI capitalization.
10. **Update preprint/software status at submission.** Confirm whether `sareffhibbert2025comparisoncopulabasedmixedmodel` has a newer version or publication; record arXiv version/date if it remains a preprint. Give `gamlss.longitudinal` an archived release/DOI or repository commit and version in its own citation.
11. **Use one bibliography, not the 192-record export.** The manuscript cites only 32 unique keys and currently compiles 26 entries. A curated submission bibliography will make validation and title-style conversion tractable.
12. **QA the compiled reference list, not only `.bib`.** Confirm author–year citation formatting under the JSS class, no raw metadata, correct hyperlinks, and a continuous References section after all floats.

## Top 10 actions, in order

1. Rebase the manuscript on the current JSS article template and confirm the correct front matter, citations, code environments, and PDF appearance.
2. Resolve all substantive draft decisions and delete the revision/todo machinery only after the corresponding prose, notation, and software behavior are settled.
3. Complete and test the reproducibility protocol: archived version/commit, environment, seeds, runtime, manifest, exact commands, hashes, and clinical-data exception.
4. Rebuild the bibliography from primary/official records; resolve all six missing keys and all 72 BibTeX warnings; add R and package citations.
5. Rewrite the title and abstract to JSS style, with an expanded acronym, concrete result, runtime/accuracy trade-off, first-order limitation, and verified availability statement.
6. Restructure and shorten the opening: remove lines 330–332, consolidate competitor discussion, and bring the package contract/minimal example forward.
7. Audit every novelty/comparative/clinical claim for scope, citation, effect size, and exploratory versus confirmatory tone.
8. Apply JSS semantic markup and convert listings to `CodeInput`/`CodeOutput`; ensure the printed workflow is executable and within width.
9. Rebuild artifacts: captions below all figures/tables, sentence style and periods, notes in captions, caption–label order, every retained artifact referenced, and readable fonts/tables.
10. Perform final line editing and a warning-clean pdfLaTeX/BibTeX build, then inspect all pages to verify bookmark hierarchy, float order, continuous references, and no overfull content.

## Final assessment

The paper is promising as a JSS software article because it combines a distinctive package interface, practical diagnostics, sensitivity guidance, and an applied example. The current frozen manuscript is nevertheless a working draft, not a JSS-ready submission. The decisive obstacles are objective and testable: wrong document class, unresolved citations and malformed bibliography, unfinished scientific wording, and incomplete reproduction metadata. Once those are fixed, the next editorial priority should be compression and calibration—state the contribution earlier, remove duplicate theory, use less absolute language, and make the code and visual evidence easier to verify.
