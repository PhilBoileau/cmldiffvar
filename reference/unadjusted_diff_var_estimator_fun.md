# Estimator of Differential Variance

`unadjusted_diff_var_estimator_fun()` estimates the differential
variance using an unadjusted estimator. This estimator is regular and
asymptotically linear when the treatment assignment mechanism does not
depend on pre-treatment covariates, such as in a randomized controlled
trial.

## Usage

``` r
unadjusted_diff_var_estimator_fun(
  clean_tbl,
  treatment_var_name,
  propensity_score_var_name,
  outcome_var_name,
  estimand_type
)
```

## Arguments

- clean_tbl:

  A pre-processed
  [tibble](https://tibble.tidyverse.org/reference/tibble.html) that is
  ready for nuisance parameter estimation.

- treatment_var_name:

  A `character` providing the column name of the treatment assignment
  indicator stored in `clean_tbl`.

- propensity_score_var_name:

  An optional `character` providing the column name of the treatment
  assignment indicator stored in `data_tbl`.

- outcome_var_name:

  A `character` providing the column name of the outcome variable stored
  in `clean_tbl`.

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
