# Efficient Influence Function of Group-Specific Variance

`efficient_influence_function_fun()` computes the efficient influence
function based on the observed data and nuisance parameter estimates for
the specified treatment group.

## Usage

``` r
efficient_influence_function_fun(
  treatment_group,
  treatment_vec,
  outcome_vec,
  ps_est_vec,
  cond_exp_outcome_est_vec,
  cond_exp_sq_outcome_est_vec,
  mean_est
)
```

## Arguments

- treatment_group:

  A binary `numeric` indicating the treatment group for which the
  variance is being estimated.

- treatment_vec:

  A binary `numeric` vector reporting observations' observed treatment
  groups.

- outcome_vec:

  A `numeric` vector reporting observations' observed outcomes.

- ps_est_vec:

  A `numeric` vector bounded to the unit interval corresponding to the
  predicted probability of being assigned to the treatment group.

- cond_exp_outcome_est_vec:

  A `numeric` corresponding to the predicted expected outcome
  conditional on confounders and treatment equal to `treatment_group`.

- cond_exp_sq_outcome_est_vec:

  A `numeric` corresponding to the predicted expected squared outcome
  conditional on confounders and treatment equal to `treatment_group`.

- mean_est:

  A `numeric` corresponding to the estimated outcome population mean of
  the designated `treatment_group`.

## Value

A `numeric` vector corresponding to the efficient influence function.
