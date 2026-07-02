# Testing Cleanup and Joint-vs-Separate Simulation Spec

## Scope

This spec covers cleanup of exploratory testing scripts and generated outputs
that sit outside the formal package test suite. It does not apply to
`tests/testthat/`, which remains the package regression test suite.

The cleanup should preserve paper-relevant evidence before deleting or
archiving old exploratory files.

## Target Buckets

Use the same buckets for exploratory scripts and exploratory results:

1. `01-continuous-bcpe-t`
   - Continuous BCPE and t-copula simulations.
   - Parameter recovery, smooth recovery, fit comparison against standard
     `gamlss`/`gamlss2` style baselines.
2. `02-discrete-delaporte-clayton`
   - Discrete Delaporte margin and Clayton copula simulations.
   - Parameter recovery, fit quality, and inference calibration.
3. `03-joint-vs-separate-optimization`
   - Joint versus separate optimisation investigations.
   - Candidate mining for the final paper simulation.
4. `04-missingness-dropout-sensitivity`
   - Missingness/dropout sensitivity checks.
   - Comparison of complete-case, missingness-aware, and sensitivity designs.
5. `05-standard-model-benchmarking`
   - Placeholder bucket for comparisons with standard GLM, GEE, GLMM, GAM,
     and related comparator models.

Recommended exploratory layout:

```text
R (testing)/
  01-continuous-bcpe-t/
  02-discrete-delaporte-clayton/
  03-joint-vs-separate-optimization/
  04-missingness-dropout-sensitivity/
  05-standard-model-benchmarking/
  _archive-delete-after-insight-capture/

results/jss-exploratory/
  01-continuous-bcpe-t/
  02-discrete-delaporte-clayton/
  03-joint-vs-separate-optimization/
  04-missingness-dropout-sensitivity/
  05-standard-model-benchmarking/
```

The formal paper workflow should continue to write final regenerated artifacts
under `results/jss-replication/<profile>/`.

## Current Inventory

Observed exploratory script categories in `R (testing)/`:

- Continuous BCPE/t: `simulation_bcpe_t_gamlss_comparison.R`,
  `combine_bcpe_t_default_100rep_outputs.R`,
  `make_bcpe_t_comparison_latex_table.R`,
  `make_bcpe_t_paper_outputs.R`,
  `score_bcpe_t_saved_predictive.R`,
  `plot_bcpe_t_defaults_pilot3.R`.
- Discrete Delaporte/Clayton: `simulation_del_clayton_gamlss_comparison.R`,
  `plot_del_clayton_defaults.R`, `make_nbi_paper_outputs.R`,
  `diagnose_rs_clayton_showcase_scenarios.R`,
  `diagnose_rs_clayton_tau_stress_parallel.R`.
- Joint vs separate: `balanced_tau_calibrated_rs_grid.R`,
  `broad_rs_joint_copula_screen.R`,
  `controlled_rs_joint_effect_strength_sweep.R`,
  `stress_rs_joint_effect_strength_search.R`,
  `tau_calibrated_copula_pilot.R`,
  `selected_rs_joint_copula_report.R`,
  `showcase_rs_joint_fit_quality_report.R`,
  `plot_stress_effect_search_trends.R`.
- Missingness: `simulation_bcpe_t_missingness_comparison.R`.
- Benchmarking: current implementation is mostly in package code:
  `R/benchmark-*.R`, plus outputs in `results/adoption_benchmarks_*` and notes
  under `inst/benchmarks/`.

Large exploratory result directories currently include:

- `results/bcpe_t_current_defaults_rep100_*`: about 1.5 GB combined, mostly
  per-replicate `.rds` files plus compact CSV summaries.
- `results/rs_joint_*`: about 180 MB combined, with thousands of cell-level
  `.rds` files plus summary/candidate CSVs and PNG reports.
- `results/adoption_benchmarks_*`: standard-model benchmark pilot and extended
  outputs.
- `results/rs_clayton_showcase_fit_quality`: discrete/joint optimisation
  report tables and repeated timestamped fit summaries.

## Cleanup Rules

Before deleting anything, create a capture note for each bucket:

```text
paper/notes/01-continuous-bcpe-t-evidence.md
paper/notes/02-discrete-delaporte-clayton-evidence.md
paper/notes/03-joint-vs-separate-optimization-evidence.md
paper/notes/04-missingness-dropout-sensitivity-evidence.md
paper/notes/05-standard-model-benchmarking-evidence.md
```

Each capture note should record:

- source scripts reviewed;
- output directories reviewed;
- key tables/figures retained;
- conclusions that are paper-relevant;
- known failed or abandoned designs and why they were abandoned;
- raw files deleted or archived after capture.

Retain:

- final or candidate-summary CSVs;
- by-replicate CSVs needed to reproduce tables;
- final report `.md`, `.tex`, `.pdf`, and PNG figure outputs;
- scripts that are either still executable or contain unreimplemented logic;
- generated manifests and output hashes from `results/jss-replication/`.

Delete only after capture:

- cell-level `.rds` files when equivalent `*_by_rep.csv`,
  `*_summary.csv`, and `*_candidates.csv` files exist;
- duplicate timestamped CSVs after hash or content comparison confirms they are
  redundant;
- smoke output directories once a non-smoke run supersedes them and the smoke
  run is not referenced by tests;
- process test logs such as `_startprocess*.out` and `_startprocess*.err`;
- one-off plot previews that are not referenced by a capture note.

Do not delete:

- `tests/testthat/`;
- private-data placeholders or application workflow stubs;
- any file with uncommitted user edits unless it is explicitly listed in the
  cleanup manifest and reviewed in the same pass.

## Joint-vs-Separate Candidate Mining

Mine candidates from summary files first; avoid opening thousands of cell-level
`.rds` files until the summaries identify specific scenarios or replicates.

Primary candidate sources:

- `results/rs_joint_broad_copula_screen/broad_top_candidates.csv`
- `results/rs_joint_broad_copula_screen/top_cases_loglik_logscore_all_rmse_improved.csv`
- `results/rs_joint_controlled_effect_sweep/controlled_top_candidates.csv`
- `results/rs_joint_stress_effect_search/stress_balanced_candidates.csv`
- `results/rs_joint_stress_effect_search/stress_max_loglik_candidates.csv`
- `results/rs_joint_tau_calibrated_pilot/tau_calibrated_balanced_candidates.csv`
- `results/rs_joint_balanced_tau_grid/balanced_tau_balanced_candidates.csv`
- `results/rs_joint_selected_copula_report/selected_joint_vs_separate_summary.csv`
- `results/rs_clayton_showcase_fit_quality/joint_vs_separate_delta_summary.csv`

Candidate scenarios should be selected by a transparent score:

```text
candidate_score =
  scaled(mean joint - separate held-out log score per observation)
  + scaled(mean joint - separate train joint log-likelihood)
  + scaled(reduction in theta/tau RMSE)
  + scaled(reduction in mu/sigma RMSE)
  - convergence penalty
  - boundary/pathology penalty
  - excessive runtime penalty
```

Selection should include:

- 3 to 5 clean examples where joint optimisation improves held-out fit and
  dependence recovery;
- 1 to 2 tie/control examples where joint and separate optimisation are
  practically equivalent;
- optionally 1 cautionary example where joint optimisation improves likelihood
  but worsens another estimand, if that is important for an honest discussion.

Exclude candidates with:

- frequent non-convergence in either method;
- fitted dependence repeatedly at parameter bounds;
- improvement driven only by train likelihood with no held-out gain;
- severe runtime blow-up unless runtime is the point of the example.

## Final Joint-vs-Separate Simulation

Implement the final simulation in `paper/R/03-joint-vs-separate-optimization.R`
so it replaces the current stub target `module_03_joint_vs_separate`.

Minimum methods:

- separate optimisation;
- joint optimisation with the current production defaults.

Optional methods, only if already supported cleanly:

- joint optimisation with alternative warm-start or optimiser controls;
- `gamlss2` marginal comparator for marginal-only context.

Recommended outputs:

- `results/jss-replication/<profile>/data/03-joint-vs-separate-candidate-selection.csv`
- `results/jss-replication/<profile>/data/03-joint-vs-separate-results.csv`
- `results/jss-replication/<profile>/tables/03-joint-vs-separate-summary.csv`
- `results/jss-replication/<profile>/tables/03-joint-vs-separate-metric-wins.csv`
- `results/jss-replication/<profile>/figures/03-joint-vs-separate-deltas.png`
- `results/jss-replication/<profile>/figures/03-joint-vs-separate-metric-dashboard.png`

Core metrics:

- convergence rate;
- elapsed seconds;
- train joint log-likelihood;
- held-out joint log score per observation;
- marginal log score;
- copula log score;
- AIC where comparable;
- parameter bias and RMSE for fixed effects;
- smooth integrated RMSE and max absolute error;
- theta/tau RMSE and absolute error over time;
- 95 percent interval coverage and interval width;
- PIT mean absolute error and PIT KS p-value;
- lower 5 percent and upper 5 percent tail errors;
- 90th percentile MAE and upper-tail error.

Reuse existing benchmark metric names where possible:

- `benchmark_mean_rmse`
- `benchmark_mean_mae`
- `benchmark_mean_bias`
- `benchmark_neg_log_score`
- `benchmark_q90_mae`
- `benchmark_upper_tail_error_90`
- `benchmark_interval_coverage_95`
- `benchmark_interval_width_95`
- `benchmark_pit_mean_abs_error`
- `benchmark_pit_ks_p_value`
- `benchmark_tail_error_lower_05`
- `benchmark_tail_error_upper_05`
- `benchmark_theta_time_abs_error`
- `smooth_eta_rmse`
- `smooth_eta_max_abs_error`
- `elapsed_sec`

Summarise method comparisons with `summarise_benchmark_results()` where the
result table shape permits it. Add joint-vs-separate specific delta columns for
paper interpretation, for example:

- `delta_test_log_score_per_obs`
- `delta_train_joint_loglik`
- `delta_copula_log_score`
- `delta_theta_rmse`
- `delta_tau_rmse`
- `delta_elapsed_sec`

## Acceptance Criteria

Cleanup is complete when:

- each exploratory script has been moved into one of the five buckets or listed
  as deleted in an evidence note;
- each large result directory has a retain/delete decision in an evidence note;
- no paper-relevant conclusion depends only on deleted raw files;
- the `05-standard-model-benchmarking` placeholder exists in the exploratory
  structure;
- `paper/R/03-joint-vs-separate-optimization.R` no longer writes a stub and
  produces real data, tables, and figures through the targets workflow;
- `paper/manifest.csv` records the final non-stub outputs for module 03.

