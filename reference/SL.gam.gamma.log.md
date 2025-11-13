# SuperLearner wrapper for GAM with Gamma family and log link

A SuperLearner wrapper that implements generalized additive models using
the Gamma family with a log link.

## Usage

``` r
SL.gam.gamma.log(Y, X, newX, obsWeights, deg.gam = 2, cts.num = 4, ...)
```

## Arguments

- Y:

  A numeric `vector` of outcome values.

- X:

  A numeric `matrix` or `data.frame` of covariates and treatment.

- newX:

  A numeric `matrix` or `data.frame` of predictors.

- obsWeights:

  Not used.

- deg.gam:

  A numeric representing the degrees of the GAM. Defaults to 2.

- cts.num:

  A numeric indicating the minimum number of unique values a numeric
  covariate requires to be considered as a continuous variable. Defaults
  to 4.

- ...:

  Any additional arguments.

## Value

A list with components:

- `pred`: A numeric vector of predictions on `newX`.

- `fit`: A list containing the fitted model object.

## Details

A copy of `[SuperLearner::SL.gam]` but using the Gammma family.
