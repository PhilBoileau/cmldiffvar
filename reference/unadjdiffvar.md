# Unadjusted Differential Variance Inference

`unadjdiffvar()` infers the differential variance of two treatments
using an unadjusted estimator.

## Usage

``` r
unadjdiffvar(
  data_tbl,
  estimand_type = "absolute",
  confidence_level = 0.95,
  treatment_var_name,
  propensity_score_var_name = NULL,
  outcome_var_name
)
```

## Arguments

- data_tbl:

  A [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`tibble`](https://tibble.tidyverse.org/reference/tibble.html).

- estimand_type:

  A `character` indicating whether to estimate the absolute
  (`"absolute"`) or relative (`"relative"`) differential variance.

- confidence_level:

  A `numeric` between \$0.1\$ and \$0.99\$ providing the confidence
  level used to compute confidence intervals. Defaults to `0.95`.

- treatment_var_name:

  A `character` providing the column name of the treatment assignment
  indicator stored in `data_tbl`.

- propensity_score_var_name:

  An optional `character` providing the column name of the treatment
  assignment indicator stored in `data_tbl`. Defaults to `NULL`. See the
  Details section for more information.

- outcome_var_name:

  A `character` providing the column name of the outcome variable stored
  in `data_tbl`.

## Value

A one-row [tibble](https://tibble.tidyverse.org/reference/tibble.html)
containing the following columns:

- `estimand`: The scale of the differential variance estimand

- `estimate`: The differential variance estimate

- `se`: The estimator's standard error

- `ci_low`: The lower bound of the Wald-type confidence interval

- `ci_high`: The upper bound of the Wald-type confidence interval

- `p_value`: The p-value of a two-sided test using the z-score of the
  differential variance estimate

## Details

Under standard causal identifiability conditions — namely, consistency,
positivity, and full exchangeability — `unadjdiffvar()` performs
inference on the differential variance of the potential outcomes.
Differential variance is defined on the absolute scale as the difference
in potential outcome variances of two treatments. On the relative scale,
differential variance is defined as the ratio of the potential outcome
variances. The scale of the differential variance estimated by
`unadjdiffvar()` is specified by the `estimand_type` parameter.

`unadjdiffvar()` assumes that the data in `data_tbl` are generated
according to a parallel group study design with a binary treatment
variable and a continuous outcome. Treatment is assumed to be assigned
at random. If the treatment assignment probabilities are known, such as
in a randomized controlled trial, then these may be passed to the
estimator using the `propensity_score_var_name` argument. If there are
(un)measured confounders, this estimator will be biased.
