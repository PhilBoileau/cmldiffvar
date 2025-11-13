# SuperLearned Conditional Expected Outcome

`estimate_cond_exp_outcome_fun()` estimates the conditional expected
outcome using a SuperLearner estimator implemented in the
[`SuperLearner::SuperLearner()`](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)
function.

## Usage

``` r
estimate_cond_exp_outcome_fun(
  clean_tbl,
  adj_set_var_names,
  treatment_var_name,
  outcome_var_name,
  cond_exp_outcome_library,
  num_nuisance_sl_folds
)
```

## Arguments

- clean_tbl:

  A pre-processed
  [tibble](https://tibble.tidyverse.org/reference/tibble.html) that is
  ready for nuisance parameter estimation.

- adj_set_var_names:

  A `character` vector providing the columns names of the adjustment set
  variables stored in `clean_tbl`.

- treatment_var_name:

  A `character` providing the column name of the treatment assignment
  indicator stored in `clean_tbl`.

- outcome_var_name:

  A `character` providing the column name of the outcome variable stored
  in `clean_tbl`.

- cond_exp_outcome_library:

  A `character` vector of candidate learners used by the SuperLearner
  estimator.

- num_nuisance_sl_folds:

  A `numeric` indicating the number of folds to use in cross-validated
  SuperLearner estimator.

## Value

A
[SuperLearner::SuperLearner](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)
class object containing the estimated conditional expected outcome.
