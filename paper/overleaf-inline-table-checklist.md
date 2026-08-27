# Overleaf inline-table conversion checklist

The asset publisher intentionally never edits `main.tex`. Before submission, convert each public data-backed inline table to `\\input{tables/...}` and add the path to `paper/manifest.csv`. Keep narrative-only tables inline.

Current exact classifications from canonical `main.tex`:

- Externalize `tab:simulation-design` as `tables/simulation-design.tex` and generate it from the BCPE/t and NBI scenario registry.
- Externalize `tab:sim_jvs_scenarios` as `tables/joint-vs-separate-scenario-design.tex` and generate it from the joint-versus-separate case definitions.
- Keep `tab:r-software-comparison`, `tab:copula_families`, `tab:copula_functions_tail_dependence`, `tab:optimizer-guidance`, and `tab:sensitivity-guidance` inline as non-data explanatory content.
- Keep `tab:marginal-screen`, `tab:copula-lik-selection`, `tab:likelihood-comparison`, and `tab:covariate-estimates` classified as private LIPID tables. They must not be added to the public publisher allowlist.

- Confirm every simulation-design, recovery, convergence, runtime, sensitivity, and benchmark table is externalized.
- Leave all LIPID tables classified as private-data publication artifacts; do not describe them as regenerated.
- Run `Rscript paper/audit-manuscript.R --paper-repo <clone>` after each TeX edit.
- Run the publisher in dry-run mode, inspect its complete allowlisted diff, then use `--apply`.
- Compile the paper clone and visually compare figures, captions, cross-references, and table widths.
