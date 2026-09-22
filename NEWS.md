# gamlss.longitudinal 0.2.0

This release is the stable package baseline before the experimental joint
missingness and full higher-order D-vine extension. Development of that method
is isolated on `codex/experimental-missingness-higher-vine`.

## Count-margin and Student-t stability

- Student-t starting-value screening now evaluates complex GAMLSS margin
  parameters on the distribution scale, broadcasts scalar copula starts across
  every response pair, omits unnecessary margin-score calculations, and picks
  the best finite degrees-of-freedom grid candidate.
- Discrete Student-t rectangle probabilities now use direct vectorized
  Gauss--Legendre density quadrature instead of subtracting four separately
  integrated CDF corners. This removes the ACTG `ZIPIG` worker crash, avoids
  tail cancellation, and substantially reduces fitting time.
- Invalid or excluded response pairs no longer trigger needless copula
  evaluation.
- `select_joint_distribution()` isolates candidates in separate R workers by
  default, supports per-candidate timeouts and resumable atomic checkpoints,
  and retains worker failures as structured result rows rather than allowing a
  native dependency fault to terminate the screening session.

## Conditional forward simulation

- `simulate()` now accepts `condition_on = "fitted"` or a custom history data
  frame to generate future trajectories conditional on each subject's most
  recent usable response. Adjacent future visits retain the fitted copula
  dependence recursively; omitted scheduled visits begin an explicitly
  independent segment. Structured conditioning provenance is attached to the
  result, and gap-induced segment starts produce one aggregate classed warning.
- Existing unconditional simulation remains the default. Analytic conditional
  prediction and integration of dependence across missing visits remain out of
  scope.

## Inference diagnostics

- Hessian validation concerns no longer discard covariance estimates that were
  successfully calculated. The covariance, standard errors, and downstream
  inference are returned with a classed warning and recorded diagnostic codes;
  errors are reserved for cases where usable inference cannot be calculated.

## Copula compatibility

- All six supported copulas (`N`, `C`, `F`, `G`, `J`, and `t`) are now
  available for every supported homogeneous marginal family, including
  dynamically discovered `gamlss.dist` families.

## Interface consistency

- `select_margin()`, `select_copula()`, and `gamlss_longitudinal()` now report
  concise progress by default. Set `verbose = 0` for completely quiet calls;
  higher levels retain detailed fitting diagnostics.
- `plot_margin_fit()` now treats `bins` as a maximum for discrete margins,
  aggregating wide integer supports and fitted masses over shared bins while
  retaining unit-width bins for smaller supports.
- `plot_margin_fit()` gains `separate_zero = TRUE` to reserve a dedicated zero
  bin when wide discrete supports are grouped, making zero inflation visible.
- `plot_margin_fit()` and `plot_copula_fit()` now accept `by` to facet observed
  and fitted overlays by a supplied grouping column, including treatment arms;
  `by_time = TRUE` remains the time-faceting shortcut. Both helpers preserve
  the supplied factor levels so matching groups appear in the same panel order.
- `plot_margin_comparison()` overlays several fitted count margins on shared
  empirical probability bars across time, with an optional boundary-corrected
  KDE reference and consistent aggregation for wide count supports.
- Segmented fits now compute and cache their routed subject-clustered sandwich
  covariance during fitting, so a compatible `summary()` call does not repeat
  the expensive inference calculation.
- Sandwich summaries now show model-based and subject-clustered robust standard
  errors side by side. Tests, p-values, intervals, and the backward-compatible
  `std_error` column continue to use the robust covariance.
- `gamlss_longitudinal()` now uses `data`, `optimization = "joint"` or
  `"separate"`, `optimizer_control`, and `start` as its primary fitting
  interface. The previous `dataset`, `include_dlcopdpar`, `start_from`, and
  flat optimizer arguments remain available with deprecation warnings for one
  transition period.
- Smoothing initialization and selection are now explicit shared controls:
  `lambda_start` and `lambda_penalty` in `gamlss_longitudinal_control()`.
- `verbose = TRUE` is equivalent to the default level `1`; integer levels 2
  and 3 provide detailed tracing.
- Fitted models, optimizer controls, model summaries, distribution selections,
  diagnostics, and model checks have compact, consistent print methods.
- Package ggplot output now shares `theme_gamlss_longitudinal()` and a common,
  colorblind-safe visual vocabulary.
- `plot_copula_diagnostics()` defaults to a focused four-panel overview and
  provides `type = "calibration"`, `"residuals"`, `"tails"`, or `"all"`.
  Its primary controls are now `empirical`, `n_bins`, and `ncol`.
- `plot_dist()` now uses `data`; its previous `dataset` argument is retained as
  a deprecated alias.
- `model_spec()` now uses `$optimization` as its canonical field while retaining
  the equivalent `$optimisation` field for one transition period.
- The compatibility aliases `plotDist()`, `plot.terms()`,
  `plot_copula_overlay()`, and `screen_margin()` now direct users to their
  canonical replacements.
