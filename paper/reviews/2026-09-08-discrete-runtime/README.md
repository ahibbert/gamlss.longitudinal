# Exact-discrete runtime investigation

This directory contains bounded profiling code for the ACTG application. The
profiler uses fixed, treatment-stratified subject samples, one outer iteration,
one inner iteration per parameter block, no covariance calculation, and a
ten-minute safety budget. It therefore diagnoses optimizer cost without
starting an unrestricted full fit.

```text
Rscript paper/reviews/2026-09-08-discrete-runtime/profile-actg.R --subjects=150 --model=full
```

Use `--run=converge` for a realistic fit capped at ten minutes. The default
`--run=bounded` performs only one inner and one outer iteration.

Generated profiles and summaries are written under the ignored `results/`
directory.

## Findings

The original bounded 50-subject profile took 5.93 seconds. Two general
likelihood bottlenecks were identified: row-by-row Gaussian copula CDF calls
and nested R loops for the 64-node stable rectangle quadrature. Replacing both
with vectorised implementations reduced the same bounded profile to 1.67
seconds.

The full-formula, 300-subject convergence benchmark was then used to assess
the discrete-margin path. All runs used seed 20260908, 1,200 scheduled rows,
857 observed responses, the DEL margin, Gaussian copula, segmented missingness,
and no covariance calculation.

| Implementation | Elapsed seconds |
|---|---:|
| Vectorised copula CDF and rectangle quadrature | 263.61 |
| One marginal CDF plus derived companions | 116.10 |
| Vectorised Delaporte CDF | 76.76 |
| Stable vectorised Delaporte survival | 58.83 |

These timings are local diagnostic measurements, not portable performance
guarantees. The final implementation also provides exact vectorised shortcuts
for every count family in the current capability registry: PO, NBI, DEL, ZIP,
ZAP, and ZINBI. Extreme-tail checks use direct log-survival calculations rather
than subtracting a rounded CDF from one.
