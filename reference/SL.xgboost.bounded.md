# SuperLearner wrapper of SL.xgboost that ensures positive predictions

A SuperLearner wrapper for SL.xgboost that ensures positive predictions.
Predictions below 1e-4 are clipped to a minimum value of 1e-4.

## Usage

``` r
SL.xgboost.bounded(..., ntrees = 100, lower = 1e-04)
```

## Arguments

- ...:

  SL.xgboost arguments.

- ntrees:

  Number of trees, default is 100.

- lower:

  Lower bound to clip predictions to, default is 1e-4.

## Value

A list with components:

- `pred`: A numeric vector of predictions.

- `fit`: A list containing the fitted model object.
