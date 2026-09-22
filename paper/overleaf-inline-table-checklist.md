# Overleaf inline-table conversion checklist

The asset publisher intentionally never edits `main.tex`. Before submission, convert each public data-backed inline table to `\\input{tables/...}` and add the path to `paper/manifest.csv`. Keep narrative-only tables inline.

Current exact classifications from canonical `main.tex`:

- Externalize `tab:simulation-design` as `tables/simulation-design.tex` and generate it from the BCPE/t and NBI scenario registry.
- Externalize `tab:sim_jvs_scenarios` as `tables/joint-vs-separate-scenario-design.tex` and generate it from the joint-versus-separate case definitions.
- Keep `tab:r-software-comparison`, `tab:copula_families`, `tab:copula_functions_tail_dependence`, `tab:optimizer-guidance`, and `tab:sensitivity-guidance` inline as non-data explanatory content.
- Do not add private-data tables to the public publisher allowlist.

## ACTG workflow assets

The ACTG workflow now generates and publishes the following manuscript-ready
tables. Replace the corresponding inline `tabular` blocks in `main.tex` with
the listed `\input{...}` calls, while retaining the table environments,
captions, labels, and surrounding prose in `main.tex`:

- Cohort and missingness summary: `\input{tables/actg-cohort.tex}`.
- Joint margin--copula shortlist: `\input{tables/actg-joint-family-screen.tex}`.
- Final parameter-specific model specification: `\input{tables/actg-final-model.tex}`.
- Treatment effects with model-based and sandwich uncertainty: `\input{tables/actg-treatment-effects.tex}`.
- Treatment-by-visit fitted predictions: `\input{tables/actg-fitted-predictions.tex}`.
- Observed-versus-simulated marginal and dependence summaries: `\input{tables/actg-simulation-comparison.tex}`.
- Full covariate-selection path for the appendix: `\input{tables/actg-covariate-selection.tex}`.

The ACTG figure paths already match the filenames used by the current draft:

- `charts/final-model-plot_dist.png`
- `charts/nbi_plot_adjusted.png`
- `charts/zisichel_plot_adjusted.png`
- `charts/gumbel_copula_fit.png`
- `charts/t_copula_fit.png`
- `charts/plot_margin_trt.png`
- `charts/plot_copula_trt.png`
- `charts/final_model_copula_diagnostics.png`
- `charts/final_model_plot.png`
- `charts/plot_terms.png`
- `charts/ACTG no intermittent missingness model.png`

After making the table substitutions, confirm that each table is included only
once and that its existing `\label{...}` remains unique. Do not copy generated
captions or labels into the `.tex` assets: those remain manuscript text.

- Confirm every simulation-design, recovery, convergence, runtime, sensitivity, and benchmark table is externalized.
- Describe only manifest-listed public tables as regenerated.
- Run `Rscript paper/audit-manuscript.R --paper-repo <clone>` after each TeX edit.
- Run the publisher in dry-run mode, inspect its complete allowlisted diff, then use `--apply`.
- Compile the paper clone and visually compare figures, captions, cross-references, and table widths.
