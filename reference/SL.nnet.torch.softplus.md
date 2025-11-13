# SuperLearner wrapper for neural network using torch

A SuperLearner wrapper that implements a neural network using torch,
with a softplus activation to encourage strictly positive predictions.
By default, the neural network uses a single hidden layer with 10 hidden
units, an Adam optimizer with a learning rate of 0.01, an MSE loss
function, and 100 epochs.

## Usage

``` r
SL.nnet.torch.softplus(
  Y,
  X,
  newX,
  n_hidden = 10,
  learning_rate = 0.01,
  loss_fn = torch::nn_mse_loss(),
  epochs = 100,
  lower = 1e-04,
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

- n_hidden:

  Number of hidden units in the hidden layer, default is 10.

- learning_rate:

  Learning rate for the Adam optimizer, default is 0.01.

- loss_fn:

  The loss function used during training, default is MSE.

- epochs:

  The number of epochs for training, default is 100.

- lower:

  Lower bound to clip predictions to, default is 1e-4

- ...:

  Any additional arguments.

## Value

A list with components:

- `pred`: A numeric vector of predictions on `newX`.

- `fit`: A list containing the fitted model object.
