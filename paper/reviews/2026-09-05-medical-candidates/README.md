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

Missing ACTG follow-up counts remain `NA`. They are omitted only from the
temporary marginal screens. Any later joint fit must use
`missingness = "segment"` and report that explicit composite-likelihood
assumption.

After screening, run a deliberately small, checkpointed joint-model grid with:

```text
Rscript paper/reviews/2026-09-05-medical-candidates/run-joint-covariates.R tlc
```

The joint comparison fits three nested specifications: constant non-mean and
copula parameters, preferred marginal-parameter covariates, and the same model
with prespecified non-time covariates in `theta`. Fits are cached individually
under the ignored `results/` directory.
