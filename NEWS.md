# gamlss.longitudinal 0.1.0

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

- `gamlss_longitudinal()` now uses `data`, `optimization = "joint"` or
  `"separate"`, `optimizer_control`, and `start` as its primary fitting
  interface. The previous `dataset`, `include_dlcopdpar`, `start_from`, and
  flat optimizer arguments remain available with deprecation warnings for one
  transition period.
- Smoothing initialization and selection are now explicit shared controls:
  `lambda_start` and `lambda_penalty` in `gamlss_longitudinal_control()`.
- Fitting is quiet by default. `verbose = TRUE` prints a concise progress
  summary; integer levels 2 and 3 provide detailed tracing.
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
