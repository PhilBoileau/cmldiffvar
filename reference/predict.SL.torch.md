# Prediction method wrapper for base learners of class SL.torch

This function provides the prediction method for objects of class
SL.torch. It converts `newdata` into a torch tensor, applies the fitted
torch model, and returns predictions as a numeric vector

## Usage

``` r
# S3 method for class 'SL.torch'
predict(object, newdata, ...)
```

## Arguments

- object:

  An object of class SL.torch.

- newdata:

  A numeric `matrix` or `data.frame` of predictors.

- ...:

  Additional arguments (not currently used).

## Value

A numeric vector of predictions.
