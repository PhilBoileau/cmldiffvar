# SuperLearner wrapper for Multivariate Adaptive Regression Splines

A SuperLearner wrapper that implements Multivariate Adaptive Regression
Splines using the Gamma family with a log link.

## Usage

``` r
SL.earth.gamma.log(
  Y,
  X,
  newX,
  degree = 2,
  penalty = 3,
  nk = max(21, 2 * ncol(X) + 1),
  pmethod = "backward",
  nfold = 0,
  ncross = 1,
  minspan = 0,
  endspan = 0,
  ...
)
```

## Arguments

- Y:

  A numeric `vector` of outcome values.

- X:

  A numeric `matrix` or `data.frame` of covariates and treatment.

- newX:

  A numeric `matrix` or `data.frame` of predictors.

- degree:

  Maximum degree of interaction. Default is 3, meaning build a model
  with interaction terms.

- penalty:

  Generalized Cross Validation (GCV) penalty per knot. Defaults is 3.

- nk:

  Maximum number of model terms before pruning, i.e., the maximum number
  of terms created by the forward pass. Includes the intercept. The
  default is semi-automatically calculated from the number of predictors
  but may need adjusting.

- pmethod:

  Pruning method. One of: "backward", "none", "exhaustive", "forward",
  "seqrep", "cv". Default is "backward". Specify pmethod="cv" to use
  cross-validation to select the number of terms. This selects the
  number of terms that gives the maximum mean out-of-fold RSq on the
  fold models. Requires the nfold argument. Use "none" to retain all the
  terms created by the forward pass.

- nfold:

  Number of cross-validation folds. Default is 0, no cross validation.

- ncross:

  Only applies if nfold\>1. Number of cross-validations. Each
  cross-validation has nfold folds. Default 1.

- minspan:

  Minimum number of observations between knots. (This increases
  resistance to runs of correlated noise in the input data.) The default
  minspan=0 is treated specially and means calculate the minspan
  internally, as per Friedman's MARS paper section 3.8 with alpha =
  0.05.

- endspan:

  Minimum number of observations before the first and after the final
  knot. The default endspan=0 is treated specially and means calculate
  the endspan internally, as per the MARS paper equation 45 with alpha =
  0.05.

- ...:

  Any additional arguments.

## Value

A list with components:

- `pred`: A numeric vector of predictions on `newX`.

- `fit`: A list containing the fitted model object.

## Details

A copy of `[SuperLearner::SL.earth]()` but using the Gammma family,
which is itself a wrapper of `[earth::earth]()`. See `[earth::earth]()`
for more details.
