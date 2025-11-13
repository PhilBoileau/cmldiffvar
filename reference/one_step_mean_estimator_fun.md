# One-Step Group-Specific Mean Estimator

`one_step_mean_estimator_fun()` estimates the group-specific mean using
a one-step estimator.

## Usage

``` r
one_step_mean_estimator_fun(
  treatment_group,
  treatment_vec,
  outcome_vec,
  ps_est_vec,
  cond_exp_outcome_est_vec
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

## Value

A `numeric` estimate of the `treatment_group`-specific mean.
