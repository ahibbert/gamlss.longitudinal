# Public medical application candidates

This exploratory analysis compares four public longitudinal clinical-trial
datasets before one is promoted into the package vignette and paper. It is not
part of the manuscript replication graph.

Run the standard margin and copula screens plus marginal-parameter covariate
screens with:

```text
Rscript paper/reviews/2026-09-05-medical-candidates/run-screening.R
```

Individual dataset names (`tlc`, `actg`, `skin`, or `epilepsy`) can be supplied
as trailing arguments. ACTG currently requires the `bcmixed` package. TLC and
skin-cancer data are read from the public Applied Longitudinal Analysis site;
epilepsy data are supplied by `MASS`.

The resulting comparison and recommendation are recorded in
[`findings.md`](findings.md). Generated CSV and RDS files remain under the
ignored `results/` directory.

Missing ACTG follow-up counts remain `NA`. They are omitted only from the
temporary marginal screens. Any later joint fit must use
`missingness = "segment"` and report that explicit composite-likelihood
assumption.

The screening runner also compares prespecified `theta` formulas using the
selected family and a fixed-margin pair likelihood. For discrete responses this
is the same exact rectangle objective used by `select_copula()`; it is a rapid
covariate screen, not a substitute for the final jointly re-estimated model.
It also exports the coefficients from the AIC-preferred marginal model so that
the direction and magnitude of non-mean effects can be reviewed.

After screening, run a deliberately small, checkpointed joint-model grid with:

```text
Rscript paper/reviews/2026-09-05-medical-candidates/run-joint-covariates.R tlc
```

The joint comparison fits three nested specifications: constant non-mean and
copula parameters, preferred marginal-parameter covariates, and the same model
with prespecified non-time covariates in `theta`. Fits are cached individually
under the ignored `results/` directory.
