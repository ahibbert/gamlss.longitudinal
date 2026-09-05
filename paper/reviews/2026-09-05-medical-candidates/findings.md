# Medical application candidate findings

## Purpose and method

Four public longitudinal clinical-trial datasets were screened as possible
paper and vignette applications. Each analysis used `select_margin()` followed
by `select_copula()`. Prespecified formulas were then compared for non-mean
margin parameters and for the copula parameter.

The dependence-covariate comparison is a fixed-margin screen. It refits the
selected margin using the full mean formula and the AIC-preferred non-mean
formula, then compares `theta` formulas using adjacent observations. Discrete
responses use exact CDF rectangles. These screens rank candidate formulas but
do not replace final joint estimation or valid post-selection inference.

## Comparative result

| Dataset | Subjects | Observed visits | Raw adjacent Spearman | Selected margin | Selected copula | Best non-mean formula; AIC gain | Best theta formula; AIC gain (BIC gain) | Assessment |
|---|---:|---:|---:|---|---|---|---|---|
| TLC lead trial | 100 | 400 | 0.559 | GG | Frank | `sigma + nu ~ treatment`; 106.26 | `theta ~ treatment`; 38.39 (34.69) | Best complete showcase; full joint fit validated |
| ACTG 193A HIV trial | 1,177 | 3,352 | 0.788 | DEL | Gaussian | `sigma ~ treatment + baseline`; 254.91 | `theta ~ treatment`; 31.86 (15.29) | Strongest large medical dataset, but current joint fitting is slow |
| Skin Cancer Prevention Trial | 1,683 | 7,081 | 0.232 | DEL | Gaussian | `sigma ~ age`; 4.15 | `theta ~ sex + skin type`; 3.04 (-10.15) | Weak evidence after complexity penalty |
| Progabide epilepsy trial | 59 | 236 | 0.664 | DEL | Gaussian | `sigma ~ treatment`; 0.19 | constant; 0.00 (0.00) | Does not meet the desired covariate criterion |

An AIC or BIC gain is the constant-parameter criterion minus the criterion for
the listed model; positive values favour the covariate model.

## TLC lead trial

The generalized gamma margin was clearly preferred to Normal, Gamma and
log-normal alternatives (the next-best margin was 24.05 AIC points behind).
The selected Frank copula also beat Gaussian by 11.36 AIC points.

In the marginal screen, treatment affected both generalized-gamma non-mean
parameters. On their link scales:

- `sigma`: intercept -0.823; placebo contrast -0.726 (log link).
- `nu`: intercept 2.115; placebo contrast -3.574 (identity link).

The fixed-margin dependence screen estimated Frank theta as 2.279 for succimer
and 7.722 for placebo, equivalent to Kendall tau of approximately 0.240 and
0.592. The full joint model strengthened this separation: theta 2.134 versus
9.525, or tau approximately 0.226 versus 0.652.

All three TLC joint fits converged quickly. Adding the two marginal-parameter
effects reduced AIC from 2606.52 to 2516.01; adding treatment to theta reduced
it further to 2477.46. This dataset satisfies every desired modelling feature
and is the best immediate paper/vignette choice.

## ACTG 193A HIV trial

This is the largest high-correlation candidate and retains missing follow-up
counts as `NA`: 3,352 observed responses and 1,356 missing responses across
1,177 patients. DEL was preferred to NBI by 165.90 AIC points. Gaussian was
preferred to Clayton in the exact discrete copula screen.

The AIC-preferred marginal model placed treatment and standardized log
baseline CD4 on DEL sigma. On the log link, its estimates were:

- intercept -0.302;
- treatment contrasts 0.254, 0.743 and 0.765 for regimens 2--4 versus 1;
- standardized log-baseline contrast -0.363.

At mean baseline, the treatment-specific sigma values are therefore about
0.74, 0.95, 1.55 and 1.59. A one-standard-deviation increase in log baseline
multiplies sigma by about 0.70.

Treatment also strongly improved the dependence screen. The Gaussian-link
coefficients were 0.592 with contrasts 0.102, 0.263 and 0.305. These correspond
to correlations of approximately 0.532, 0.601, 0.694 and 0.715, or Kendall tau
of 0.357, 0.410, 0.488 and 0.507.

This is a compelling substantive application, but it is not yet the fastest
submission route. A full segmented joint fit progressed normally but exceeded
17 minutes within an optimizer iteration and was stopped. It should be pursued
only if the larger trial is worth a dedicated performance pass.

## Skin Cancer Prevention Trial

DEL was selected, but its margin advantage over NBI was only 6.05 AIC points.
Age on sigma improved marginal AIC by 4.15, with a log-link coefficient of
-0.229 (about a 20% reduction in sigma per standard deviation of age). BIC
preferred constant sigma.

The best dependence screen used sex and skin type. Its estimated Gaussian
correlations ranged from about 0.125 for female/other skin to 0.319 for
male/burns skin (Kendall tau about 0.080 to 0.207). Although AIC improved by
3.04, BIC worsened by 10.15. This is too weak for the principal demonstration.

## Progabide epilepsy trial

DEL and Gaussian were selected and the raw repeated-measure correlation was
high. However, adding treatment to sigma improved the marginal AIC by only
0.19. Its log-link estimate was 0.491, but the full joint model had worse AIC
than the constant-parameter fit (1285.21 versus 1284.49). No theta covariate
improved AIC. The dataset should not be advanced for this paper.

## Recommendation

Use TLC as the immediate replacement application because it combines a
distinctly GAMLSS margin, strong and interpretable treatment effects on both
non-mean marginal parameters and dependence, complete four-visit follow-up,
and fast joint fitting. Keep ACTG as the leading alternative if a performance
pass makes full segmented joint fitting practical. Do not invest further in
the skin or epilepsy candidates for the main manuscript.

