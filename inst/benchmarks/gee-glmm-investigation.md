# GEE/GLMM Benchmark Investigation

This note starts the benchmark track for comparing `gamlss.longitudinal` with
the models applied users already reach for: GEEs, GLMMs, GAMMs, and marginal
GAMLSS references. The first package scaffold is
`benchmark_standard_models()`, which fits common comparator models and returns
a tidy result table with availability, success, elapsed time, log-likelihood
where meaningful, AIC where meaningful, and response-scale MAE/RMSE.

## Current Harness

The existing `run_coverage_simulations()` harness already supports:

- internal longitudinal methods: `rs_separate`, `rs_joint`, and `cg`;
- marginal references: `gamlss`, plus optional `gamlss2`;
- common simulation axes: family, copula, sample size, time points,
  dependence strength, mean-covariate and scale-covariate designs,
  missingness, and starts;
- recorded fit outcomes: convergence/failure taxonomy, runtime, marginal
  log-likelihood, joint log-likelihood, parameter error, and runtime summaries.

That makes it the right base for the broader adoption benchmark. The new
`benchmark_standard_models()` helper can now be used inside or alongside that
harness to add GEE/GLMM/GAM comparator rows before expanding into repeated
simulation coverage metrics.

The benchmark track now also has an executable adoption layer:

- `adoption_benchmark_scenarios()` returns named simulation scenarios with the
  family, copula, design, methods, primary metrics, and claim being tested.
- `run_adoption_benchmarks()` repeats those scenarios through
  `run_coverage_simulations()` and attaches a `summarise_benchmark_results()`
  win/tie/loss table.

This gives the investigation a stable unit of work: add or refine scenarios,
run them repeatedly outside CRAN checks, and promote mature scenarios into
vignettes or papers once the results are stable.

## Local Comparator Availability

Checked in the current development environment:

- `geepack`: available; use for GEE baselines.
- `lme4`: available; use for Gaussian, Poisson, and binomial GLMM baselines.
- `mgcv`: available; use for GAM/GAMM-style smooth baselines where relevant.
- `glmmTMB`: not installed locally; keep optional until a benchmark runner
  declares the dependency explicitly.

## Minimal Scaffold

```r
comparators <- benchmark_standard_models(
  data = dat,
  formula = response ~ treatment + time + age_scaled,
  subject_var = "subject",
  family = "gaussian",
  comparators = c("gee", "glmm", "gam", "glmmTMB")
)

comparators$results
```

This is deliberately a first-stage scaffold. It establishes comparable fitting
records before the larger simulation harness adds estimand-specific truth,
coverage, calibration, and tail metrics.

## Repeated Simulation Scaffold

The coverage harness can also include standard comparator methods explicitly:

```r
results <- run_coverage_simulations(
  families = "NO",
  copulas = "N",
  methods = c("rs_separate", "gee", "glmm", "gam"),
  designs = "covariate",
  n = 80,
  times = 1:3,
  write_results = FALSE
)

results[
  results$method %in% c("gee", "glmm", "gam"),
  c(
    "method",
    "success",
    "benchmark_estimator",
    "benchmark_mean_rmse",
    "benchmark_q90_mae",
    "benchmark_upper_tail_error_90",
    "benchmark_interval_coverage_95",
    "benchmark_tail_error_lower_05",
    "benchmark_tail_error_upper_05"
  )
]
```

Comparator methods are opt-in; the default coverage grid remains focused on the
native and marginal GAMLSS methods. Standard comparators currently map only
GAMLSS families with clear mean-model analogues (`NO`/`NO2`, `PO`/`ZIP`/`ZIP2`,
`GA`, and `EXP`). Other families are recorded as unsupported comparator-family
rows rather than treated as numerical failures.

The first benchmark metrics are deliberately separated by estimand:

- `benchmark_mean_bias`, `benchmark_mean_mae`, and `benchmark_mean_rmse` compare
  fitted response means with simulated `true_mu`.
- `benchmark_q90_mae` compares the model-implied 90th percentile with the
  simulated truth, where the comparator distribution is available.
- `benchmark_upper_tail_error_90` compares the model-implied probability above
  the true 90th percentile with the simulated truth. This is available for
  Gaussian, Gamma, and Poisson-style comparator families.
- `benchmark_theta_time_abs_error` compares the estimated time effect in the
  copula parameter with the simulated truth on the link scale.
- `benchmark_mae` and `benchmark_rmse` compare fitted response means with the
  noisy observed response.
- `benchmark_interval_coverage_95`, `benchmark_pit_ks_p_value`,
  `benchmark_pit_mean_abs_error`, and 5% tail errors are residual-calibrated
  Gaussian diagnostics. They are available for Gaussian comparator families and
  should be treated as in-sample diagnostic screens, not final inferential
  coverage claims.

These metrics are now populated for successful native
`gamlss.longitudinal` rows as well as standard comparator rows, so the summary
tables can score the package's own model against GEE, GLMM, and GAM baselines on
the same simulated rows.

Use `summarise_benchmark_results()` after repeated simulations to get a compact
win/tie/loss table:

```r
summary <- summarise_benchmark_results(
  results,
  metrics = c(
    "benchmark_mean_rmse",
    "benchmark_q90_mae",
    "benchmark_upper_tail_error_90",
    "benchmark_interval_coverage_95",
    "benchmark_pit_mean_abs_error",
    "elapsed_sec"
  )
)

summary$summary
```

The summary scores each method inside a simulation case, labels the case result,
then aggregates by method and metric. This gives the benchmark section a stable
language for claims such as "wins on mean RMSE", "ties on interval calibration",
or "loses on runtime".

For the standard adoption benchmark plan, start from the named scenarios:

```r
adoption_benchmark_scenarios()
```

Then run a small smoke pass:

```r
smoke <- run_adoption_benchmarks(
  scenarios = "gaussian_heteroskedastic",
  reps = 2,
  methods = c("rs_separate", "gee", "glmm", "gam"),
  max_elapsed_sec = 20,
  write_results = FALSE
)

smoke$summary$summary
```

And run a full opt-in benchmark outside CRAN checks:

```r
bench <- run_adoption_benchmarks(
  reps = 100,
  methods = c("rs_separate", "gee", "glmm", "gam"),
  max_elapsed_sec = 60,
  write_results = TRUE
)

bench$summary$summary
```

## Opt-In Runner Script

The package also ships a reproducible runner for longer benchmark campaigns:

```r
script <- system.file(
  "benchmarks",
  "run-adoption-benchmarks.R",
  package = "gamlss.longitudinal"
)
source(script)
```

The script writes:

- `adoption_benchmark_results.csv` for all fitted rows;
- `adoption_benchmark_summary.csv` for win/tie/loss summaries;
- `adoption_benchmark_case_results.csv` for per-scenario, per-replicate scores;
- `adoption_benchmark_scenarios.csv` for the exact scenario plan;
- `adoption_benchmark_comparator_status.csv` for optional backend availability;
- `adoption_benchmark_report.md` for a Markdown evidence report;
- `adoption_benchmark_object.rds` for the complete R object.

The Markdown report includes both overall headline tables and scenario-level
tables. Scenario-level tables are restricted to each scenario's declared primary
metrics. Use that section for adoption claims; the aggregate table is mainly a
quick orientation because each scenario answers a different applied question.

Use `write_benchmark_report()` directly when benchmark results are already in
memory:

```r
write_benchmark_report(
  bench,
  path = "results/adoption_benchmarks/adoption_benchmark_report.md"
)
```

Use environment variables to control longer runs without editing package code:

```sh
GAMLSS_LONGITUDINAL_ADOPTION_REPS=100 \
GAMLSS_LONGITUDINAL_ADOPTION_SCENARIOS=gaussian_heteroskedastic,gamma_positive \
GAMLSS_LONGITUDINAL_ADOPTION_METHODS=rs_separate,gee,glmm,gam \
GAMLSS_LONGITUDINAL_ADOPTION_MAX_ELAPSED_SEC=60 \
GAMLSS_LONGITUDINAL_ADOPTION_OUTPUT_DIR=results/adoption_benchmarks \
Rscript -e "source(system.file('benchmarks', 'run-adoption-benchmarks.R', package = 'gamlss.longitudinal'))"
```

## Pilot Sanity Run

A bounded local pilot run on 2026-05-28 used all five named scenarios, three
replicates per scenario, and methods `rs_separate`, `gee`, `glmm`, and `gam`:

```sh
GAMLSS_LONGITUDINAL_ADOPTION_REPS=3 \
GAMLSS_LONGITUDINAL_ADOPTION_METHODS=rs_separate,gee,glmm,gam \
GAMLSS_LONGITUDINAL_ADOPTION_MAX_ELAPSED_SEC=30 \
GAMLSS_LONGITUDINAL_ADOPTION_OUTPUT_DIR=results/adoption_benchmarks_pilot \
Rscript -e "devtools::load_all(quiet = TRUE); source('inst/benchmarks/run-adoption-benchmarks.R')"
```

The pilot is a wiring and plausibility check, not citation-grade benchmark
evidence. It produced 60 successful fits with finite native
`gamlss.longitudinal` metrics for mean, quantile, tail, interval, dependence,
and runtime estimands. The headline pilot summary was:

- dependence recovery: `rs_separate` was the only method with finite
  `benchmark_theta_time_abs_error` because the standard baselines do not model
  time-varying copula dependence;
- interval calibration: `rs_separate` had the highest win-or-tie rate in the
  pilot (`10/15` wins);
- mean RMSE and 90th-percentile MAE: GEE won most often in this small pilot,
  with `rs_separate` close enough to motivate larger scenario-specific runs;
- runtime: GEE was fastest in all pilot cases;
- all pilot claims need larger replication before they should appear in a
  paper, vignette conclusion, or package-site performance claim.

## Extended Local Run

A larger local run on 2026-05-28 used all five named scenarios, 20 replicates
per scenario, and methods `rs_separate`, `gee`, `glmm`, and `gam`:

```sh
GAMLSS_LONGITUDINAL_ADOPTION_REPS=20 \
GAMLSS_LONGITUDINAL_ADOPTION_METHODS=rs_separate,gee,glmm,gam \
GAMLSS_LONGITUDINAL_ADOPTION_MAX_ELAPSED_SEC=30 \
GAMLSS_LONGITUDINAL_ADOPTION_OUTPUT_DIR=results/adoption_benchmarks_extended \
Rscript -e "devtools::load_all(quiet = TRUE); source('inst/benchmarks/run-adoption-benchmarks.R')"
```

This produced 400 successful fits (`100/100` for each method) in the local
development environment. The scenario-level report supports cautious,
scenario-specific lessons:

- time-varying dependence: only `rs_separate` had finite
  `benchmark_theta_time_abs_error`, because the GEE, GLMM, and GAM baselines do
  not estimate the time-varying copula parameter;
- Gaussian heteroskedasticity: `rs_separate` led the declared 90th-percentile
  metric and was competitive on mean RMSE, while GEE remained fastest;
- missing visits: `rs_separate` led the declared interval and mean metrics in
  this run, while GEE remained fastest;
- positive Gamma and Poisson count scenarios: GEE led the simple mean,
  quantile, tail, and runtime headline metrics under these data-generating
  mechanisms.

This is useful evidence for the package positioning rather than a universal
win claim. It says the adoption case is strongest when distributional
calibration or dependence is the estimand, and weakest when the applied question
is only a simple marginal mean that a fast GEE already targets well.

## Citation-Grade Evidence Checklist

The current scaffold is ready for repeated opt-in runs, but the local pilot and
20-replicate extended run should still be treated as development evidence. A
paper, package-site claim, or vignette conclusion should be based on a locked
benchmark plan with:

1. a fixed scenario table exported by `adoption_benchmark_scenarios()`;
2. a fixed method set and package-version record from
   `benchmark_comparator_status()`;
3. enough replicates per scenario to estimate win/tie/loss rates with stable
   Monte Carlo error;
4. scenario-level primary metrics declared before results are inspected;
5. explicit unsupported-method labels for comparators that cannot represent the
   estimand;
6. sensitivity runs for sample size, time-grid length, missing visits, and
   dependence strength;
7. saved CSV, RDS, and Markdown report artifacts from
   `run-adoption-benchmarks.R`;
8. a short "when not to use" interpretation that names where GEE or GLMM remains
   simpler or stronger.

The benchmark report should be read scenario-by-scenario. Aggregate headlines
are useful for orientation, but they should not be used as a universal ranking
because each scenario targets a different estimand.

## First Comparisons To Add

The first comparison families are now represented in the named adoption
scenarios. Continue to refine them narrowly and keep the win/tie/loss cases
interpretable:

- Gaussian longitudinal outcome with heteroskedasticity:
  compare `geepack::geeglm()`, `lme4::lmer()`, `mgcv::gam()`, and
  `gamlss_longitudinal()` with the `scale` design on variance calibration,
  prediction intervals, tail behaviour, residual dependence, and runtime.
- Gamma or log-normal positive outcome:
  compare GEE/GLMM mean fits against GAMLSS-copula margins on mean prediction,
  quantile calibration, tail under/over-coverage, and convergence.
- Count outcome:
  compare Poisson/negative-binomial GEE/GLMM where supported against
  GAMLSS-copula count margins on dispersion, upper-tail probabilities, and
  subject-level trajectory simulation.
- Time-varying dependence:
  compare a GAMLSS-copula model with `theta ~ time` against exchangeable GEE,
  random-intercept GLMM, and GAM baselines. Use
  `benchmark_theta_time_abs_error` for dependence recovery and mean/tail
  metrics for response-scale context.

## Metrics

The benchmark result table should add these columns beyond the current coverage
metrics:

- `comparator_class`: `gamlss_longitudinal`, `gee`, `glmm`, `gamm`, or
  `marginal_gamlss`;
- `estimator`: package/function label such as `geepack::geeglm`;
- `estimand`: marginal mean, conditional mean, quantile, tail probability, or
  dependence;
- `coverage_95`: empirical interval coverage for the estimand;
- `calibration_error`: PIT/quantile calibration error for held-out rows;
- `tail_error`: lower/upper tail probability error;
- `trajectory_score`: proper score for subject trajectories when available.

## Open Implementation Questions

- GEEs target marginal mean parameters; GLMMs target conditional mean
  parameters. Benchmarks must state the estimand before comparing bias.
- Some GLMM competitors cannot represent heteroskedasticity or arbitrary GAMLSS
  shape parameters. Those should be labelled as unsupported rather than failed.
- `glmmTMB` should remain optional because it is valuable but heavy; benchmark
  runners can skip it when unavailable.
- Small-sample comparisons need repeated simulation with scheduled or opt-in
  workflows rather than CRAN-time tests.
