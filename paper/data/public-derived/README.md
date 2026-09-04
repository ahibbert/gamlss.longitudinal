# Public analysis-ready inputs

These files are compact, non-clinical inputs for the `paper` profile. They are retained so manuscript assets can be rebuilt without ignored local result directories.

- `bcpe-t/`: fixed-seed per-replicate BCPE/t simulation estimates and audit metrics consumed by `make_bcpe_t_paper_outputs.R`.
- `nbi-clayton/`: fixed-seed per-replicate NBI/Clayton simulation estimates and audit metrics consumed by `make_nbi_paper_outputs.R`.
- `copula-misspecification/`: a signed six-file public Gamma simulation bundle used to rebuild the manuscript summary heatmap and claims. It contains exact fit attempts, same-dataset paired effects with Monte Carlo uncertainty, selection attempts and summaries, a structured warning audit, and the generation-time execution manifest. The registered production design is exactly Gamma margins crossed with Gaussian and Clayton generating/fitted copulas, two Kendall's $\tau$ levels, three subject sample sizes, and 100 replicates: 1,200 generated datasets and 2,400 fit attempts. Only a fresh `ga-nc-v2` checkpoint run is eligible; older wider-copula or `ga-nc-v1` evidence must not be published.
- `correlation-misspecification/`: public standard-model benchmark summaries and TeX producers.
- `joint-vs-separate/`: per-replicate delta CSVs and scenario summaries used to regenerate the three approved TeX tables.
- `missingness/`: aggregate public simulation summaries used to redraw the two missingness figures. The `full` profile recreates their upstream checkpointed simulation.

No private-application or RAND records, synthetic analogues, fitted objects, empirical pseudo-observations, residuals, design matrices, or data-derived intermediates are permitted here. Input hashes are recorded for every replication run.
