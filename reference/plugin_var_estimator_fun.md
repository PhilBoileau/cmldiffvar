# Plug-In Estimator of the Group-Specific Variance

`plugin_var_estimator_fun()` estimates of the group-specific variance
using a plug-in estimator.

## Usage

``` r
plugin_var_estimator_fun(cond_exp_outcome_est_vec, cond_exp_sq_outcome_est_vec)
```

## Arguments

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
