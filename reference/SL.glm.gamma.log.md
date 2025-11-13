# SuperLearner wrapper for GLM with Gamma family and log link

A SuperLearner wrapper that implements generalized linear models using
the Gamma family with a log link.

## Usage

``` r
SL.glm.gamma.log(Y, X, newX, ...)
```

## Arguments

- Y:

  A numeric `vector` of outcome values.

- X:

  A numeric `matrix` or `data.frame` of covariates and treatment.

- newX:

  A numeric `matrix` or `data.frame` of predictors.

- ...:

  Any additional arguments.

## Value

A list with components:

- `pred`: A numeric vector of predictions on `newX`.

- `fit`: A list containing the fitted model object.

## Details

The predictor matrix `X` is assumed to have a binary treatment column as
its last column. Thus, we don't consider the last column when squared
terms are defined. In addition to the main terms, coefficients for the
squared covariates are included in the model.
