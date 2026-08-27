# Public analysis-ready inputs

These files are compact, non-clinical inputs for the `paper` profile. They are retained so manuscript assets can be rebuilt without ignored local result directories.

- `bcpe-t/`: fixed-seed per-replicate BCPE/t simulation estimates and audit metrics consumed by `make_bcpe_t_paper_outputs.R`.
- `nbi-clayton/`: fixed-seed per-replicate NBI/Clayton simulation estimates and audit metrics consumed by `make_nbi_paper_outputs.R`.
- `copula-misspecification/`: public Gamma simulation results used to rebuild the manuscript summary heatmap.
- `correlation-misspecification/`: public standard-model benchmark summaries and TeX producers.
- `joint-vs-separate/`: per-replicate delta CSVs and scenario summaries used to regenerate the three approved TeX tables.
- `missingness/`: aggregate public simulation summaries used to redraw the two missingness figures. The `full` profile recreates their upstream checkpointed simulation.

No LIPID or RAND records, synthetic analogues, fitted objects, empirical pseudo-observations, residuals, design matrices, or clinical data-derived intermediates are permitted here. Input hashes are recorded for every replication run.
