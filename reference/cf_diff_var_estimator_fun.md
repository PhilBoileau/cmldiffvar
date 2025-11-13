# Cross-Fitted Differential Variance Estimation

`cf_diff_var_estimator_fun()` produces a cross-fitted estimate of the
select differential variance estimand, using either a one-step or
targeted maximum likelihood estimator.

## Usage

``` r
cf_diff_var_estimator_fun(
  fold,
  clean_tbl,
  propensity_score_adj_var_names,
  cond_exp_outcome_adj_var_names,
  treatment_var_name,
  propensity_score_var_name,
  outcome_var_name,
  propensity_score_library,
  cond_exp_outcome_library,
  cond_exp_sq_outcome_library,
  num_nuisance_sl_folds,
  estimator_type,
  estimand_type
)
```

## Arguments

- fold:

  A [fold](http://tlverse.org/origami/reference/make_folds.md) object
  indicating which observations in `clean_tbl` are members of the
  training and validation sets.

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

- propensity_score_library:

  A `character` vector of candidate learners used by the SuperLearner
  estimator.

- cond_exp_outcome_library:

  A `character` vector of candidate learners used by the SuperLearner
  estimator.

- cond_exp_sq_outcome_library:

  A `character` vector of candidate learners used by the SuperLearner
  estimator.

- num_nuisance_sl_folds:

  A `numeric` indicating the number of folds to use in cross-validated
  SuperLearner estimator.

- estimator_type:

  A `character` indicating whether to use a one-step or a targeted
  maximum likelihood estimator in the cross-fitting procedure.

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
