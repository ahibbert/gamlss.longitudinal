# Independent visual-communication review for Journal of Statistical Software

## Scope and recommendation

I independently inspected all 47 supplied page PNGs at their original 918 x 1188 resolution and checked the frozen 47-page PDF directly. The review covers page architecture, margins, hierarchy, author notes, floats, code and output, every caption, Figures 1-18, and Tables 1-21. I did not consult earlier reviews or internal audits.

**Recommendation: major revision; do not advance the current PDF to editorial or production review.** Two production-integrity defects are blocking: visible author/TODO material remains throughout the manuscript, and float placement breaks reading order, including code and the reference list. After those are cleared, several tables and multi-panel figures still require redesign for JSS-scale readability and accessibility.

## Five strengths

1. The underlying page grid is disciplined: body text, equations, section headings, running heads, and page numbers are generally aligned and consistent.
2. Most code blocks stay within the text width, use a distinct monospaced treatment, and visually support the package-workflow narrative.
3. Tables generally use restrained horizontal rules, meaningful spanners, and logical row grouping rather than visually heavy grids.
4. The figure sequence covers the full user journey--exploration, fitting, validation, diagnostics, and application--and most visuals are explicitly called out in the surrounding prose.
5. Rendering is sharp, with no broken glyphs or rasterization artifacts; several graphics already use useful redundant encodings, especially point shape in Figure 13 and line type/markers in Figure 18.

## BLOCKING issues

### B1. Visible drafting notes and margin comments are part of the rendered paper

- **Location:** pages 1, 2, 3, 5, 6, 7, 9, 11, 13, 17, 20, 21, 30, 31, 32, and 33.
- **Evidence:** page 1 contains a large red “Next steps/Key Revisions” box plus an orange margin question; pages 2-7 contain orange editorial comments; pages 9, 11, 17, 20, 21, and 32 contain red “Next steps” boxes; page 13 has an orange optimizer note; page 30 contains a red runtime placeholder; page 31 says to report the operating system, R version, versions, cores, seed policy, and runtime; page 33 contains red “[Draft wording ...]” text. Several orange boxes extend into the outer margin and visually approach the trim edge.
- **Consequence:** this is an unambiguous submission blocker. The material exposes unfinished decisions, interrupts the argument, creates apparent margin overflow, and makes the PDF look like an internal markup copy rather than a manuscript.
- **Remedy:** remove the drafting-note layer/content at source, resolve each note in prose or delete it, and rebuild the frozen PDF from a clean submission mode. Do not merely recolor or hide the boxes behind white objects.
- **Completion test:** a fresh text extraction contains no “Next steps”, “Draft wording”, runtime placeholder, or note-to-author language; a manual pass over all pages shows no red/orange editorial boxes or marginal callouts. Purposeful red/orange data marks in plots may remain.
- **Classification/confidence:** **BLOCKING -- production integrity; confidence: very high.**

### B2. Float order breaks code, appendix hierarchy, and the reference list

- **Location:** Figure 2 on page 14; appendices C-I on pages 39-47; Tables 19-21 on pages 40-41; Figures 13-18 on pages 42-47.
- **Evidence:** Figure 2 appears at the top of page 14 before the remainder of the `margin_screen` output that began on page 13. Appendix G begins a code block with two lines at the foot of page 39, then Tables 19 and 20 occupy page 40 and Table 21 starts page 41 before the code resumes. Appendix headings C, D, and E appear consecutively on page 39 without content beneath them. Headings H and I appear together on page 41 before either appendix’s figures. Figure 13 then appears on page 42, the “References” heading starts immediately below it, and Figures 14-18 are interleaved with references on pages 43-47.
- **Consequence:** the visual reading order is incorrect. Executable code is not reproducible as presented, appendix ownership is ambiguous, and the bibliography appears to begin before five manuscript figures have been presented.
- **Remedy:** impose float barriers at section/appendix boundaries; keep Figure 2 after the complete code/output sequence it illustrates; place Tables 18-21 together before Appendix G or after its complete code block; place Figures 13-14 under Appendix H and Figures 15-18 under Appendix I; start References only after all floats have been flushed. Populate, merge, or remove empty appendix headings C-E.
- **Completion test:** every heading is followed by its own content; no code block is interrupted by a figure or table; Figure 13 precedes heading I; Figures 14-18 precede “References”; and no figure/table appears after the first reference entry.
- **Classification/confidence:** **BLOCKING -- reading order and reproducibility; confidence: very high.**

## MAJOR issues

### M1. Several tables are below practical JSS reading size

- **Location:** Table 1 (page 2), Table 3 (page 8), Table 6 (page 19), Table 10 (page 22), Table 13 (page 26), and Table 14 (page 27). Table 2 (page 7) is borderline because of its mathematical cells.
- **Evidence:** direct PDF span inspection finds text as small as 5-6 pt on pages 2 and 8, 6 pt on page 19, and 6.5-7 pt on pages 26-27. Table 1 compresses 13 columns plus four notes; Table 3 compresses long copula formulas; Table 10 uses a dense multilevel 13-column header; Table 13 adds a very small explanatory footnote; Table 14 renders model formulas in particularly small monospaced type.
- **Consequence:** column headings, symbols, formulas, and notes cannot be read comfortably at 100% screen view or on a normal paper print. Readers will miss definitions or mis-associate values with columns.
- **Remedy:** reduce rather than scale: split comparison tables into coherent panels, abbreviate repeated labels, move detailed formulas/secondary metrics to an appendix or supplement, and use landscape only if it preserves the reading sequence. Table body and notes should be approximately caption-sized, not 5-6 pt.
- **Completion test:** at 100% on a standard screen and at actual-size A4/Letter print, every header, cell, symbol, and footnote is readable without zoom; no table text is materially smaller than its caption.
- **Classification/confidence:** **MAJOR -- legibility; confidence: high.**

### M2. Dense multi-panel figures use undersized labels and inefficient panels

- **Location:** Figure 3 (page 16), Figure 4 (page 18), Figures 10-11 (page 37), Figure 12 (page 38), Figure 14 (page 43), Figure 16 (page 45), and Figure 17 (page 46).
- **Evidence:** these contain between 6 and 18 panels. Panel subtitles, axes, legends, and annotations are markedly smaller than the captions. Figures 10-12 reserve empty or sparse facet positions for parameters unavailable in one model; Figure 11 has a largely empty theta panel; Figures 3 and 16 compress nine distinct diagnostic questions into one page-width graphic; Figure 17 places 16 term panels on one page, with long or rotated labels that are difficult to decode.
- **Consequence:** the central evidence is technically present but not visually available at normal reading size. Dense diagnostic dashboards encourage superficial pattern matching and conceal outliers, uncertainty, and parameter identity.
- **Remedy:** split dashboards into two or more conceptually grouped figures; remove empty facets; directly label curves where possible; enlarge strip, axis, legend, and annotation text to roughly caption size; and move secondary diagnostic panels to the appendix. Use a stable panel order and common axes only where comparisons are valid.
- **Completion test:** all annotations are readable at 100% without zoom; no facet is empty merely to preserve a rectangular grid; each figure supports one clear comparison; and a reader can identify model, parameter, metric, and uncertainty from the figure alone.
- **Classification/confidence:** **MAJOR -- information density; confidence: high.**

### M3. Low contrast and color-only distinctions weaken accessibility

- **Location:** Figure 1 (page 4), Figure 5 (page 24), Figure 7 (page 27), Figure 8 (page 29), and Figure 9 (page 30); smaller concerns also apply to Figures 10-13.
- **Evidence:** Figure 1 uses very pale, tiny blue points. Figures 7-9 rely heavily on cyan versus salmon/red. In Figure 9 the observed histogram bins are so faint that the fitted contours dominate and “fitted versus observed” is not a fair visual comparison. Figure 8 repeats two “Treatment group” legends--one with labels Control/Treatment and one with 0/1--and uses translucent fills that will merge in grayscale. Figure 5 uses separate red and blue scales with very light cells but no cell values or explicit neutral reference.
- **Consequence:** readers with color-vision deficiencies, grayscale print users, and readers on low-contrast displays may not distinguish groups or may fail to see the observed data layer.
- **Remedy:** use a color-vision-safe palette plus redundant line type, shape, outline, or direct labels; raise point/bin contrast; remove the duplicate Figure 8 legend by mapping both layers to one consistently labelled factor; and give heatmaps an explicit neutral reference and/or selected cell annotations.
- **Completion test:** every group and fitted/observed layer remains distinguishable in grayscale and under common red-green and blue-yellow color-vision simulations; Figure 9’s observed bins remain visible without suppressing the contours; Figure 8 has one unambiguous legend.
- **Classification/confidence:** **MAJOR -- accessibility and evidential balance; confidence: high.**

### M4. Three captions contradict or underspecify the plotted data

- **Location:** Figure 4 (page 18), Figure 8 (page 29), Figure 11 (page 37), and Figure 9 (page 30).
- **Evidence:** Figures 4 and 11 call the outcome “absolute bias”, yet their y-axes contain negative values; absolute bias cannot be negative. Figure 8 calls the clinical histograms “true distribution histograms”, whereas the figure title and surrounding clinical-example prose describe observed responses. Figure 9’s one-line caption does not explain that shaded bins are observed adjacent-time dependence and contours are the fitted t-copula, nor does it identify the grouping encoding.
- **Consequence:** readers can misinterpret signed error as magnitude, mistake observed clinical data for known truth, and fail to understand which marks support the fitted-versus-observed comparison.
- **Remedy:** decide whether Figures 4 and 11 show signed bias or absolute bias and make the transformation, y-axis, caption, and prose agree. Replace “true” with “observed” in Figure 8 unless a genuine known generating distribution is plotted. Expand Figure 9’s caption to define panels, bins, contours, and groups.
- **Completion test:** an independent reader can state the plotted estimand and encoding from each caption; no “absolute” quantity has negative plotted values; “true”, “simulated”, “empirical”, and “observed” are used only for their correct data source.
- **Classification/confidence:** **MAJOR -- interpretation integrity; confidence: very high.**

### M5. Page packing is uneven and sometimes suppresses the visual hierarchy

- **Location:** page 26 (Figure 6 plus Tables 12-13), page 37 (Figures 10-11 plus Table 16), pages 42-47 (figures mixed with references), and page 31 (large residual whitespace).
- **Evidence:** pages 26 and 37 stack three numbered objects with limited separation, forcing small table/figure text; page 31 leaves roughly the lower quarter empty; pages 42-47 alternate large visuals and bibliography entries. The imbalance is partly caused by unresolved notes and uncontrolled floats.
- **Consequence:** primary and secondary results receive the same visual weight, captions become cramped, and the document rhythm alternates between overpacked and underfilled pages.
- **Remedy:** reflow only after B1-B2 are fixed; give primary figures sufficient page area, move secondary results to the appropriate appendix, and permit a dedicated figure page where a diagnostic dashboard genuinely warrants it.
- **Completion test:** no page requires reducing labels to fit three major objects; no reference page contains a manuscript float; and isolated whitespace reflects a deliberate section break rather than float failure.
- **Classification/confidence:** **MAJOR -- hierarchy and page composition; confidence: high.**

## MINOR issues

### m1. Caption style and punctuation are inconsistent

- **Location:** examples include Table 8 (page 21), Table 12 (page 26), Table 18 (page 39), Figure 9 (page 30), and Figures 15-18 (pages 44-47).
- **Evidence:** several captions omit terminal punctuation or read as title fragments; capitalization varies (“Total Cholesterol”, “Treatment v Control”); some captions name the function/output but not the visual takeaway.
- **Consequence:** the manuscript departs from sentence-style JSS captions and makes scanning less consistent.
- **Remedy:** use sentence case and terminal punctuation throughout; define abbreviations/encodings on first use; prefer a concise declarative description over a fragment.
- **Completion test:** every numbered object has a grammatical, sentence-style caption with consistent punctuation and capitalization.
- **Classification/confidence:** **MINOR -- editorial consistency; confidence: high.**

### m2. Code output is visually light even where its width is correct

- **Location:** pages 12-14 and 16-17.
- **Evidence:** output is rendered in a small, light gray, slanted monospaced face. The distinction from input is useful, but some numeric columns and comments are low contrast at print size.
- **Consequence:** long summaries are harder to scan, particularly in grayscale or on a projector.
- **Remedy:** retain the input/output distinction but darken output, avoid unnecessary italics for dense numeric results, and excerpt only the rows required by the prose.
- **Completion test:** output is readable in grayscale at actual size, and the narrative does not require the reader to parse an unabridged console dump.
- **Classification/confidence:** **MINOR -- code typography; confidence: medium-high.**

### m3. Reference pages contain visually noisy database fields and raw links

- **Location:** pages 42-47, especially references 6, 7, and 11.
- **Evidence:** raw URLs wrap over multiple lines and labels such as “Publisher:” and `_eprint:` appear in the bibliography. Long Google Books and publisher URLs create uneven word spacing and distracting line breaks.
- **Consequence:** the reference list looks like an unclean BibTeX export and competes visually with the adjacent figures.
- **Remedy:** use the JSS bibliography style consistently, retain DOI links where available, remove database-only fields, and suppress redundant raw URLs.
- **Completion test:** no `_eprint:` or extraneous database field is printed; DOI/URL presentation is consistent; references contain no manuscript figures after B2 is resolved.
- **Classification/confidence:** **MINOR -- bibliography presentation; confidence: high.**

### m4. A few otherwise serviceable visuals need small legend/annotation enlargement

- **Location:** Figures 2 (page 14), 6 (page 26), 13 (page 42), 15 (page 44), and 18 (page 47); Tables 2, 7-9, 11-12, and 15-21.
- **Evidence:** the primary marks and cells are legible, but legends, facet strips, or table notes are smaller than the adjacent caption. Figure 18’s line-key text and Figure 13’s facet strips are particular examples.
- **Consequence:** useful visuals become less accessible than necessary even though their structure can be retained.
- **Remedy:** enlarge labels/notes, shorten text where necessary, and check at actual print size.
- **Completion test:** captions, legends, annotations, and table notes have comparable visual size and remain readable at 100%.
- **Classification/confidence:** **MINOR -- finishing; confidence: high.**

## CSV-ready figure and table audit

Status is one of `keep`, `redesign`, `move`, or `remove`. “Move” means the object can be retained but must be relocated to restore reading order; the action may also request local redesign.

```csv
item,page,status,priority,issue,action
"Figure 1",4,"redesign","major","Pale tiny points and small axes reduce contrast and print readability","Increase point opacity/size and annotation fonts; use a darker accessible treatment"
"Figure 2",14,"move","blocking","Float interrupts the margin_screen output begun on page 13","Place after the complete related code/output; retain graphic with slightly larger labels"
"Figure 3",16,"redesign","major","Nine diagnostics are compressed into caption-subscale panels","Split into grouped diagnostics and enlarge all panel text/legends"
"Figure 4",18,"redesign","major","Six small facets and caption says absolute bias although negative values are plotted","Enlarge/simplify facets and reconcile signed versus absolute bias"
"Figure 5",24,"redesign","major","Very light dual heatmaps have small keys and weak neutral-value guidance","Use accessible scales, explicit zero/neutral reference, and selected numeric annotations"
"Figure 6",26,"keep","minor","Useful seven-panel comparison but legend/strip text is borderline small","Retain structure; enlarge annotations and verify one consistent distribution key"
"Figure 7",27,"redesign","major","Cyan/salmon contours rely on color and do not directly expose model differences","Use redundant line styles/direct labels and an accessible palette"
"Figure 8",29,"redesign","major","Duplicate treatment legends, color-reliant overlays, and incorrect 'true histogram' wording","Unify legend, add redundant encoding, raise fill contrast, and label data as observed"
"Figure 9",30,"redesign","major","Observed bins are nearly invisible and the caption does not define marks","Increase bin contrast, distinguish fit redundantly, and make caption self-contained"
"Figure 10",37,"redesign","major","Small legend/axes and an empty model-parameter facet waste space","Remove empty facet, enlarge labels, and use a compact comparison layout"
"Figure 11",37,"redesign","major","Sparse tiny facets and absolute-bias caption conflicts with negative values","Reformat as a compact coefficient plot and correct the estimand wording"
"Figure 12",38,"redesign","major","Small annotations and empty/sparse facet structure reduce usable area","Remove unavailable facets and enlarge legend, strips, and ribbon explanation"
"Figure 13",42,"move","blocking","Missingness figure appears after both H and I headings, obscuring appendix ownership","Place under H before heading I; retain shape redundancy and enlarge strips"
"Figure 14",43,"move","blocking","Eighteen-panel H figure appears after References begins and text is very small","Move under H before I/References, then split or enlarge the facet grid"
"Figure 15",44,"move","blocking","I figure is embedded in the reference list","Move under I before References; retain plot with consistent caption punctuation"
"Figure 16",45,"move","blocking","Nine-panel I diagnostic appears alone inside References with caption-subscale labels","Move under I before References and split/enlarge diagnostic groups"
"Figure 17",46,"move","blocking","Sixteen-panel term plot is inside References and labels are difficult to decode","Move under I before References; redesign as grouped coefficient/term plots"
"Figure 18",47,"move","blocking","Tail-calibration figure appears before the last references rather than in Appendix I","Move under I before References; enlarge legend and retain line-style redundancy"
"Table 1",2,"redesign","major","Thirteen crowded columns and notes render at roughly 5-6 pt","Split comparisons, abbreviate repeated fields, and restore caption-size body text"
"Table 2",7,"keep","minor","Logical compact summary but mathematical entries and notes are borderline small","Retain; enlarge formula cells/notes and verify symbol definitions"
"Table 3",8,"redesign","major","Long copula formulas are compressed to roughly 5-7 pt","Split formula and tail-dependence content or move full formulas to an appendix"
"Table 4",18,"keep","minor","Clear compact design with a slightly small note","Retain and enlarge the note to match caption scale"
"Table 5",19,"keep","minor","Readable comparison with ample width but dense statistical notation","Retain; standardize decimal alignment and note size"
"Table 6",19,"redesign","major","Large parameter matrix reaches approximately 6 pt and is hard to scan","Split by parameter block or move detailed rows to an appendix; enlarge type"
"Table 7",20,"keep","minor","Useful guidance table but prose-heavy cells are dense","Retain; tighten wording and keep body/notes at readable size"
"Table 8",21,"keep","minor","Strong decision-oriented structure; caption lacks consistent terminal punctuation","Retain and standardize caption/line breaks"
"Table 9",21,"keep","minor","Compact scenario key is visually effective","Retain; enlarge mathematical expression only if final print test fails"
"Table 10",22,"redesign","major","Dense 13-column multilevel header creates high decoding load","Split AR(1) and exchangeable results or separate fit from dependence metrics"
"Table 11",23,"keep","minor","Narrower benchmark table is readable but metric group labels are visually close","Retain; add spacing or subtle rule between metric groups"
"Table 12",26,"keep","minor","Short AIC comparison is clear but caption style is inconsistent","Retain and standardize caption punctuation/decimal alignment"
"Table 13",26,"redesign","major","Body and long explanatory note render around 6.5 pt on an already dense page","Move detail to appendix or simplify columns; enlarge note and body"
"Table 14",27,"redesign","major","Model formulas and sequential-comparison note are too small for reliable reading","Use wrapped formula labels or a staged model table with fewer columns"
"Table 15",28,"keep","minor","Dense but well grouped full-page coefficient table","Retain; enlarge/signpost the final link-function note and verify decimal alignment"
"Table 16",37,"keep","minor","Readable summary but crowded beneath two figures","Retain after page reflow; give it more vertical separation"
"Table 17",38,"keep","minor","Parameter blocks are compact and coherent","Retain; enlarge table note/caption if needed after reflow"
"Table 18",39,"keep","minor","Clear scenario table; caption punctuation is inconsistent","Retain and standardize caption/abbreviation notes"
"Table 19",40,"move","blocking","Float interrupts the Appendix G simulation code","Place with Tables 18-21 before G or after the complete code block"
"Table 20",40,"move","blocking","Float interrupts the Appendix G simulation code","Place with Tables 18-21 before G or after the complete code block"
"Table 21",41,"move","blocking","Float interrupts the Appendix G simulation code immediately before it resumes","Place with Tables 18-21 before G or after the complete code block"
```

## Publication-ready completion gate

The revised PDF should be considered visually ready only when: (1) all drafting notes are absent; (2) code, appendices, floats, and references follow a continuous reading order; (3) no table or figure annotation is materially smaller than its caption; (4) Figures 4, 8, and 11 use correct data/estimand language; (5) every comparison survives grayscale and color-vision checks; and (6) a final 100% screen and actual-size print inspection of all pages finds no unreadable cell, legend, strip, or note.
