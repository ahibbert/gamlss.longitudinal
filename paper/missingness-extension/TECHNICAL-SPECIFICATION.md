# Longitudinal GAMLSS missingness extension: technical specification

Status: design baseline for review; not yet an implementation claim.

This document specifies the first continuous-outcome implementation. The code
architecture must retain a capability boundary through which discrete margins
can be added later without changing the user-facing meaning of the missingness
model.

## 1. Primary contribution and scope

The primary methodological contribution is an observed-data model that jointly
estimates longitudinal GAMLSS margin parameters, serial dependence parameters,
and a user-specified observation regression while accounting for potentially
informative response missingness. A full, pooled-regression D-vine is supporting
methodology: it prevents a restrictive first-order outcome model from being
confounded with the missingness adjustment.

The first implementation covers:

- scheduled continuous longitudinal responses;
- complete records, leading missing responses, intermittent gaps, monotone
  dropout, and declared mixtures of intermittent missingness and dropout;
- fully observed exogenous covariates and observation-history summaries that
  are defined at every at-risk scheduled visit;
- arbitrary registered continuous GAMLSS margins, with Gaussian and Gamma
  serving as the initial reference families;
- the existing first-tree copula model plus all higher D-vine trees;
- one pooled higher-tree formula per copula parameter;
- a logistic observation regression with fixed and smooth terms;
- outcome-level and conditional-innovation selection connections;
- one penalized joint objective for all estimated components; and
- fixed-association sensitivity analysis as the primary MNAR analysis.

Initially excluded are missing covariates, informative visit scheduling,
outcomes undefined after death or another terminal event, and exact discrete
higher-vine likelihoods. These exclusions are capability boundaries, not claims
that the cases can be treated as ordinary response missingness.

## 2. Scheduled data and risk sets

For subject i, let t = 1,...,T_i index scheduled visits in chronological order,
Y_it be the response that exists at that visit, and R_it indicate whether it is
observed. The numerical visit time is d_it. Vine topology follows scheduled
order; elapsed-time predictors use differences in d_it.

An at-risk indicator A_it determines whether an observation-process factor is
present.

- `process = "intermittent"`: A_it = 1 at every scheduled visit.
- `process = "dropout"`: A_i1 = 1 and later visits remain at risk only until
  the first declared dropout. No Bernoulli factors are added after dropout.
- `process = "mixed"`: an explicit at-risk or missingness-reason variable must
  distinguish temporary non-observation from permanent dropout.

A final run of missing responses cannot, from response values alone, reveal
whether it is dropout or intermittent non-observation through study end. The
software must not silently infer this scientific distinction. Monotone data can
be validated automatically after the user declares dropout.

The scheduled panel must contain observation-model covariates even where Y_it is
missing. Subject-constant covariates may be filled across an inserted panel only
after equality within subject has been verified. Time-varying exogenous
covariates may not be invented by the package.

## 3. Full-data outcome model

Let

    u_it = F_it(y_it; psi_it),

where every component of psi_it has its usual GAMLSS additive predictor. For an
edge joining visits a < b with conditioning set D_ab = {a+1,...,b-1}, define
the conditional arguments u_ia|D and u_ib|D recursively using the lower-tree
h-functions. The full D-vine density is

    f_Yi(y_i) = product_t f_it(y_it)
                * product_{ell=1}^{T_i-1}
                  product_{t=ell+1}^{T_i}
                  c^Y_{i,t-ell,t|D}(u_left|D, u_right|D; theta_e).

Tree 1 retains the package's current copula family and `theta.formula` /
`zeta.formula` semantics. Trees 2,...,T_i-1 use a pooled edge table and one
formula per copula parameter. For example,

    higher.theta.formula = ~ tree + s(time_distance)

shares coefficients across all higher edges without forcing their fitted
dependence to be identical.

Required generated higher-edge variables are:

- `tree`: D-vine tree number;
- `order_distance`: number of scheduled steps between endpoints;
- `condition_size`: number of conditioning visits;
- `time_left`, `time_right`, `time_mid`, and `time_distance`;
- `left_<name>` and `right_<name>` for endpoint covariates; and
- unprefixed subject-constant covariates after within-subject validation.

`time_distance` uses the original numeric visit scale, not the internal ordered
index. A factor visit variable can define topology but cannot define elapsed
distance unless a separate numeric time variable is supplied.

Setting every higher-tree copula to independence must reproduce the existing
first-order outcome likelihood. Optional truncation is a comparator and
diagnostic; the new model defaults to the full available vine.

## 4. User-specified observation model

At each at-risk visit, the user supplies a one-sided formula such as

    observation.formula = ~ treatment + sex + s(time) + last_observed_response

The package creates R_it internally. The formula defines

    logit(pi_it) = z_it^T gamma + sum_j s_j(z_ijt).

All formula variables must be available immediately before the observation
decision. Raw Y_it is therefore not permitted in `observation.formula`: it is
unavailable when R_it = 0. Dependence on the current potentially unobserved
outcome is introduced by the selection connection below.

Useful generated history variables should include:

- `previous_observed`;
- `n_previous_observed`;
- `last_observed_response`;
- `time_since_last_observed`; and
- `previously_missing`.

The exact names will be finalized with the implementation, but their
definitions must depend only on recorded history. Users may also supply their
own fully defined history summaries.

## 5. Selection connections

Let K_alpha(s,v) be a bivariate selection copula joining a response score S_it
to a latent observation uniform V_it. Define R_it = 1 when V_it <= pi_it. Then

    q_it(s) = P(R_it = 1 | S_it = s, history)
            = partial K_alpha(s, pi_it) / partial s.

The independence copula gives q_it(s) = pi_it.

Two response scores remain supported.

### 5.1 Conditional-innovation connection

Define

    W_it = G_it(Y_it | Y_i,1:(t-1), X_i),

where G_it is the sequential conditional CDF implied by the complete full
D-vine outcome model. Under a correctly specified continuous outcome model,
W_it is uniform conditionally on the complete past. Use S_it = W_it.

This connection asks whether observation depends on whether the current outcome
is unexpectedly high or low given the subject's history. It also gives the
calibration

    P(R_it = 1 | recorded history, at risk) = pi_it

under the stated sequential construction.

### 5.2 Outcome-level connection

Use S_it = U_it = F_it(Y_it; psi_it). This asks whether observation depends on
the current outcome's location in its population marginal distribution.

Because U_it need not be uniform conditional on the subject's recorded history,
pi_it is a baseline observation parameter rather than automatically the
conditional observation rate. Calibration and scientific sensitivity ranges
must therefore be reported separately for this connection.

The association alpha is initially supplied on Kendall's tau scale and mapped
to the selected copula's native parameter. The primary analysis refits every
other component jointly over a user-specified tau grid. Estimating alpha from
the same incomplete responses is not the default because the MNAR association
is generally weakly or non-identifiable without additional restrictions.

## 6. Normalized joint law and observed likelihood

For a complete latent outcome vector and an allowed observation pattern,

    p(y_i, r_i) = f_Yi(y_i)
                  * product_{t:A_it=1}
                    q_it(S_it)^r_it {1-q_it(S_it)}^(1-r_it).

Post-dropout statuses are deterministic and contribute no additional Bernoulli
factor. For fixed y_i, the sequential observation-process probabilities sum to
one across all allowed patterns, so summing over patterns preserves f_Yi and its
population-level GAMLSS margins.

If O_i and M_i index observed and missing scheduled responses, the subject's
observed-data contribution is

    L_i = integral p(y_iO, y_iM, r_i) d y_iM.

With u_it = F_it(y_it), changing variables for missing responses gives

    L_i = {product_{t in O_i} f_it(y_it)}
          * integral_[0,1]^|M_i|
            c_Yi(u_iO, u_iM)
            * observation_factors(u_iO, u_iM) d u_iM.

This unit-hypercube form is the numerical reference because missing marginal
densities cancel with the Jacobian. Conditional innovations needed by the
observation factors are recovered from the candidate full u vector using the
same D-vine recursion as the outcome density.

## 7. Analytic reductions and numerical integration

Integration is used only after applying valid reductions.

1. Complete subject: no integration.
2. Higher-tree independence: use the first-order outcome path exactly.
3. Terminal dropout with innovation selection: integrate the first missing
   innovation analytically to 1 - pi at dropout, and integrate the unrestricted
   future outcome tail to one.
4. A terminal intermittent run with innovation selection and predictors based
   only on recorded history: integrate its sequential innovations analytically,
   giving the product of the corresponding 1 - pi terms.
5. Selection independence: observation factors leave the outcome
   marginalization unchanged. Terminal outcome tails integrate to one, but an
   interior gap still requires the correct marginal density for later observed
   visits.
6. Outcome-level selection generally retains at least a one-dimensional
   integral at the first unobserved terminal visit because U_it is not uniform
   conditional on past outcomes.

Remaining low-dimensional integrals use fixed deterministic quadrature nodes
during optimization so the objective is reproducible and sufficiently smooth.
The result is audited at higher resolution. A controlled randomized
quasi-Monte Carlo fallback may be used for larger missing blocks, with repeated
scrambles and an error estimate recorded in the fitted object.

Quadrature nodes, edge topology, model matrices, and missingness-pattern plans
are cached outside the optimizer. Integration diagnostics include dimension,
method, node count, resolution difference or scramble error, and the number of
non-finite integrand evaluations.

## 8. Proposed public API

The API should use specification objects rather than adding many loosely
related scalar arguments to `gamlss_longitudinal()`.

```r
fit <- gamlss_longitudinal(
  data = dat,
  margin_dist = GA(),
  copula_dist = "N",
  time_var = "visit_time",
  subject_var = "id",
  mu.formula = response ~ treatment + s(visit_time),
  sigma.formula = ~ treatment,
  theta.formula = ~ s(visit_time),
  higher_dependence = vine_higher(
    copula_dist = "N",
    theta.formula = ~ tree + s(time_distance),
    max_tree = Inf
  ),
  observation = observation_model(
    formula = ~ treatment + s(visit_time) + previously_missing,
    process = "intermittent",
    connection = "innovation",
    copula_dist = "N",
    association_tau = 0.25
  ),
  optimization = "joint"
)
```

An absent `higher_dependence` specification preserves the current first-order
model. An absent `observation` specification preserves the current missingness
contract until the new observed-data likelihood is explicitly requested.

Required observation-process choices are `"intermittent"`, `"dropout"`, and
`"mixed"`. Required connections are `"innovation"` and `"outcome"`.
`association_tau = 0` is the independence/MAR submodel for copula families in
which zero Kendall association is independence. Family-specific validation must
not assume that a native zero parameter always implies independence.

## 9. Joint estimation contract

The estimator minimizes one penalized negative observed-data log likelihood
over:

- all active GAMLSS margin coefficients;
- first-tree copula coefficients;
- pooled higher-tree copula coefficients;
- observation-regression coefficients; and
- all active smooth coefficients and smoothing components.

Fixed sensitivity associations are inputs, not estimated coefficients. A
sequential or unadjusted fit may supply starting values but is not the reported
estimator. The initial correctness path may use numerical gradients. Analytic
or recursively propagated derivatives can replace them only after comparison
with high-accuracy numerical derivatives.

## 10. Simulation demonstrator

Before the full publication simulation, build a small reviewable demonstrator
from paired generated complete datasets.

Each replicate retains the complete outcomes, then applies a declared
missingness process. Fit:

1. the complete-data oracle;
2. ordinary unadjusted `gamlss_longitudinal()` on observed responses;
3. the adjusted model under independence/MAR;
4. the adjusted model with the correct MNAR connection and association; and
5. the adjusted model with the alternative connection.

The first visual report shows true and estimated trajectories for location,
scale, and a scientifically relevant upper quantile. A second panel shows bias
or error against missingness strength. Every result is paired by generated
complete dataset. The report also displays observation rates, pattern counts,
convergence, integration error, and runtime.

The demonstrator establishes inspectability, not publication-level operating
characteristics. The later study adds Monte Carlo uncertainty, interval
coverage, association sensitivity grids, outcome/dependence misspecification,
multiple sample sizes and visit counts, and bootstrap uncertainty where needed.

## 11. Implementation modules

Proposed new internal modules are:

- `vine-design.R`: scheduled topology and generated edge covariates;
- `vine-recursion.R`: conditional PITs, h-functions, and full density;
- `vine-higher-spec.R`: pooled higher-tree formula contract;
- `observation-model-spec.R`: user-facing observation specification;
- `observation-history.R`: generated recorded-history covariates;
- `observation-risk.R`: intermittent/dropout/mixed risk sets;
- `observation-selection.R`: connection-specific probabilities;
- `observed-likelihood-plan.R`: pattern-specific analytic reductions;
- `observed-likelihood-integrate.R`: deterministic integration and audits; and
- `model-joint-objective.R`: one penalized objective adapter.

Existing margin, copula backend, formula, smoothing, optimizer, prediction,
simulation, and reporting code should be reused behind explicit interfaces.
The existing adjacent-pair cache remains the first-order fast path rather than
being stretched into a full-vine state representation.

## 12. Design questions deliberately left open

These choices do not block the mathematical and indexing tests:

- the first three- or four-parameter margin used in the substantive example;
- the substantive application and scientifically calibrated tau range;
- whether the initial higher-tree family must match tree 1;
- the deterministic quadrature library and high-dimensional fallback; and
- whether estimating association under an exclusion restriction belongs in the
  first paper or a follow-up.

