# Causal Machine Learning for Differential Variance

`cmldiffvar()` infers the differential variance between two treatments
using causal machine learning methods.

## Usage

``` r
cmldiffvar(
  data_tbl,
  estimand_type = "absolute",
  estimator_type = "tmle",
  confidence_level = 0.95,
  propensity_score_adj_var_names,
  cond_exp_outcome_adj_var_names,
  treatment_var_name,
  propensity_score_var_name = NULL,
  outcome_var_name,
  propensity_score_library = c("SL.glm", "SL.earth", "SL.ranger"),
  cond_exp_outcome_library = c("SL.glm", "SL.earth", "SL.ranger"),
  cond_exp_sq_outcome_library = c("SL.glm.interaction", "SL.ranger"),
  num_nuisance_sl_folds = 5,
  cross_fit = FALSE,
  num_cross_fit_folds = 5,
  parallel = FALSE
)
```

## Arguments

- data_tbl:

  A [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`tibble`](https://tibble.tidyverse.org/reference/tibble.html).

- estimand_type:

  A `character` indicating whether to estimate the absolute
  (`"absolute"`) or relative (`"relative"`) differential variance.

- estimator_type:

  A `character` indicating whether to use a one-step (`"one-step"`) or a
  targeted maximum likelihood estimator (`"tmle"`) in the cross-fitting
  procedure.

- confidence_level:

  A `numeric` between \$0.1\$ and \$0.99\$ providing the confidence
  level used to compute confidence intervals. Defaults to `0.95`.

- propensity_score_adj_var_names:

  A `character` vector providing the columns names of the adjustment set
  variables for propensity score estimation stored in `data_tbl`.
  Ignored if `propensity_score_var_name` is `NULL`.

- cond_exp_outcome_adj_var_names:

  A `character` vector providing the columns names of the adjustment set
  variables for conditional expected (squared) outcome estimation stored
  in `data_tbl`.

- treatment_var_name:

  A `character` providing the column name of the treatment assignment
  indicator stored in `data_tbl`.

- propensity_score_var_name:

  An optional `character` providing the column name of the treatment
  assignment indicator stored in `data_tbl`. Defaults to `NULL`. See the
  Details section for more information.

- outcome_var_name:

  A `character` providing the column name of the outcome variable stored
  in `data_tbl`.

- propensity_score_library:

  A `character` vector of candidate learners used by the SuperLearner
  estimator of the propensity score. Defaults to
  `c("SL.glm", "SL.earth", "SL.ranger")`.

- cond_exp_outcome_library:

  A `character` vector of candidate learners used by the SuperLearner
  estimator of the expected outcome conditional on confounders and
  treatment assignment. Defaults to
  `c("SL.glm", "SL.earth", "SL.ranger")`.

- cond_exp_sq_outcome_library:

  A `character` vector of candidate learners used by the SuperLearner
  estimator of the expected squared outcome conditional on confounders
  and treatment assignment. Defaults to
  `c("SL.glm.interaction", "SL.ranger")`.

- num_nuisance_sl_folds:

  A `numeric` indicating the number of folds to use in cross-validated
  SuperLearner estimators. Defaults to `5`.

- cross_fit:

  A `logical flag` determining whether cross-fitted estimators are used.
  Defaults to `TRUE`.

- num_cross_fit_folds:

  A `numeric` setting the number of folds to use by the cross-fitting
  procedures, if used. Defaults to `5`.

- parallel:

  A `logical flag` indicating whether the cross-fitting procedure, if
  used, should be parallelized. Defaults to `FALSE`.

## Value

A one-row [tibble](https://tibble.tidyverse.org/reference/tibble.html)
containing the following columns:

- `estimand`: The scale of the differential variance estimand

- `estimate`: The differential variance estimate

- `se`: The estimator's standard error

- `ci_low`: The lower bound of the Wald-type confidence interval

- `ci_high`: The upper bound of the Wald-type confidence interval

- `p_value`: The p-value of a two-sided test using the z-score of the
  differential variance estimate

## Details

`cmldiffvar()` assumes that the data in `data_tbl` are generated
according to a parallel group study design with a binary treatment
variable, a continuous outcome, and continuous or binary
treatment-outcome confounders.

Under standard causal identifiability conditions — namely, consistency,
positivity, and full exchangeability — `cmldiffvar()` performs inference
on the differential variance of the potential outcomes. Differential
variance is defined on the absolute scale as the difference in potential
outcome standard deviations of two treatments. On the relative scale,
differential variance is defined as the ratio of the potential outcome
variances. The scale of the differential variance estimated by
`cmldiffvar()` is specified by the `estimand_type` parameter.

These differential variance estimands rely on three nuisance parameters:
the propensity score, the expected outcome conditional on confounders
and treatment assignment, and the expected squared outcome conditional
on confounders and treatment assignment.

`cmldiffvar()` implements cross-fitted one-step and targeted maximum
likelihood estimators of the differential variance estimands. These
cross-fitted estimators are used by default; non-cross-fitted estimators
can be used by setting `cross_fit = FALSE`. The number of folds used in
these cross-fitting procedures is set by `num_cross_fit_folds`. The
differential variance estimands' nuisance parameters are flexibly
estimated with
[SuperLearner](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)
estimators.

The estimators implemented in `cmldiffvar()` are consistent if *at
least* one of the following conditions is satisfied: (1) the propensity
score is consistently estimated, and (2) the expected outcome
conditional on confounders and treatment assignment and the expected
squared outcome conditional on confounders and treatment assignment are
consistently estimated.

The estimators implemented in `cmldiffvar()` are asymptotically linear —
meaning their asymptotic sampling distribution is normally distributed
about the true differential variance parameter — if *all* nuisance
parameter estimators are consistently estimated at a rate of
\\o_P(n^{-1/4})\\. The confidence intervals, standard errors, and
p-values reported by `cmldiffvar()` are incorrect if these rate
conditions are not satisfied.

When the data is the product of a randomized study with known propensity
scores, these propensity scores can be provided to the
`propensity_score_var_name` parameter. The conditions required of the
estimators to be consistent and asymptotic linear are automatically
satisfied when know propensity scores are used.

The cross-fitted estimation procedures can be parallelized by setting
`parallel = TRUE`. Parallelization relies on the
[future](https://future.futureverse.org/reference/future.html) package.
Instructions for setting up parallel processing are available in the
[future](https://future.futureverse.org/reference/future.html) package's
vignettes.
