# Generate Counterfactual Dataset

`generate_counterfactural_tbl_fun()` generates a counterfactual version
of `clean_tbl`.

## Usage

``` r
generate_counterfactural_tbl_fun(
  clean_tbl,
  treatment_group,
  propensity_score_adj_var_names,
  cond_exp_outcome_adj_var_names,
  treatment_var_name,
  propensity_score_var_name,
  propensity_score_sl_fit,
  cond_exp_outcome_sl_fit,
  cond_exp_sq_outcome_sl_fit
)
```

## Arguments

- clean_tbl:

  A pre-processed
  [tibble](https://tibble.tidyverse.org/reference/tibble.html) that is
  ready for nuisance parameter estimation.

- treatment_group:

  A binary `numeric` indicating the treatment group for which the
  variance is being estimated.

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

## Value

A counterfactual
[tibble](https://tibble.tidyverse.org/reference/tibble.html), with the
additional of propensity score, conditional expected outcome, and
conditional expected square outcome estimates.
