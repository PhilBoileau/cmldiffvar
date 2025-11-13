# One-Step Estimator of Differential Variance

`one_step_diff_var_estimator_fun()` estimates the differential variance
using a one-step estimator.

## Usage

``` r
one_step_diff_var_estimator_fun(
  clean_tbl,
  propensity_score_adj_var_names,
  cond_exp_outcome_adj_var_names,
  treatment_var_name,
  propensity_score_var_name,
  outcome_var_name,
  propensity_score_sl_fit,
  cond_exp_outcome_sl_fit,
  cond_exp_sq_outcome_sl_fit,
  estimand_type
)
```

## Arguments

- clean_tbl:

  A pre-processed
  [tibble](https://tibble.tidyverse.org/reference/tibble.html) that is
  ready for nuisance parameter estimation.

- propensity_score_adj_var_names:

  A `character` vector providing the columns names of the adjustment set
  variables for propensity score estimation stored in `clean_tbl`.

- cond_exp_outcome_adj_var_names:

  A `character` vector providing the columns names of the adjustment set
  variables for outcome regression estimation stored in `clean_tbl`.

- treatment_var_name:

  A `character` providing the column name of the treatment assignment
  indicator stored in `clean_tbl`.

- propensity_score_var_name:

  An optional `character` providing the column name of the treatment
  assignment indicator stored in `data_tbl`.

- outcome_var_name:

  A `character` providing the column name of the outcome variable stored
  in `clean_tbl`.

- propensity_score_sl_fit:

  A
  [SuperLearner::SuperLearner](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)
  object corresponding to the estimated propensity score.

- cond_exp_outcome_sl_fit:

  A
  [SuperLearner::SuperLearner](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)
  object corresponding to the estimated conditional expected outcome.

- cond_exp_sq_outcome_sl_fit:

  A
  [SuperLearner::SuperLearner](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)
  object corresponding to the estimated conditional expected squared
  outcome.

- estimand_type:

  A `character` indicating whether to estimate an absolute or relative
  effect. Set this parameter to `"absolute"` to estimate the difference
  of group-specific standard deviations. Set this parameter to
  `"relative"` to estimate the ratio of the group-specific variances.

## Value

A named list with the following components:

- `estimate`: A `numeric` estimate of the selected estimand.

- `eif`: A `numeric` vector of the efficient influence function of the
  selected estimand.
