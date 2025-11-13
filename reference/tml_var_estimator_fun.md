# Targeted Minimum Loss-Based Estimator of Group-Specific Variance

`tml_var_estimator_fun()` estimates of the group-specific variance using
a targeted minimum loss-based estimator.

## Usage

``` r
tml_var_estimator_fun(
  treatment_group,
  treatment_vec,
  outcome_vec,
  ps_est_vec,
  cond_exp_outcome_est_vec,
  cond_exp_sq_outcome_est_vec
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

## Value

A list containing the `numeric` estimate of the
`treatment_group`-specific marginal outcome variance and the `numeric`
efficient influence function vector.
