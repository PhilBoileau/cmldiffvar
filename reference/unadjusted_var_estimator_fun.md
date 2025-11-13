# Unadjsuted Group-Specific Variance Estimator

`unadjusted_var_estimator_fun()` estimates of the group-specific
variance using an unadjusted plug-in estimator. Note that this estimator
is regular and asymptotically linear in a nonparametric model when
treatment assignment is randomized. The efficient influence function
returned by this estimator corresponds to that of the complete data.

## Usage

``` r
unadjusted_var_estimator_fun(
  treatment_group,
  treatment_vec,
  outcome_vec,
  ps_est_vec
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

## Value

A list containing the `numeric` estimate of the
`treatment_group`-specific marginal outcome variance and the `numeric`
efficient influence function vector.
