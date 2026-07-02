# Continuous BCPE/t Evidence

## Reviewed Sources

- Scripts moved to `R (testing)/01-continuous-bcpe-t/`.
- Compact results retained under `results/jss-exploratory/01-continuous-bcpe-t/`.
- Raw per-replicate `.rds` files archived under `results/_archive/testing-cleanup-20260612/results/bcpe_t_current_defaults_rep100_*`.

## Retained Evidence

- Current defaults 100-replicate summaries for CG, RS joint, and RS separate.
- Comparison tables and figures from `bcpe_t_current_defaults_rep100_comparison`.
- Root-level BCPE summary CSVs and retained preview plots moved to `loose-root-results`.

## Captured Insights

- The retained summaries preserve fixed-effect recovery, smooth recovery,
  predictive scores, joint distribution diagnostics, fit logs, and publication
  table/figure artifacts.
- Bulky per-replicate `.rds` objects were archived because equivalent by-rep
  CSVs and compact summaries remain active.

## Archive Decision

Keep summaries active for paper review; use the archive only if an unpublished
diagnostic needs raw fit objects.

