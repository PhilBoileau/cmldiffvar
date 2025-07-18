# Set global variable `self
utils::globalVariables(c("self"))

#' SuperLearner wrapper for GLM with Gamma family and identity link
#'
#' @description A SuperLearner wrapper that implements generalized linear models
#'   using the Gamma family with an identity link.
#'
#' @details The predictor matrix `X` is assumed to have a binary treatment
#'   column as its last column. Thus, we don't consider the last column when
#'   squared terms are defined. In addition to the main terms, coefficients for
#'   the squared covariates are included in the model.
#'
#' @param Y A numeric `vector` of outcome values.
#' @param X A numeric `matrix` or `data.frame` of covariates and treatment.
#' @param newX A numeric `matrix` or `data.frame` of predictors.
#' @param ... Any additional arguments.
#'
#' @return A list with components:
#' * `pred`: A numeric vector of predictions on `newX`.
#' * `fit`: A list containing the fitted model object.
#'
#' @export
SL.glm.gamma.identity <- function(Y, X, newX, ...) {
  # Get column names of predictor matrix X
  col_names <- colnames(X)

  # Get treatment column name
  treatment_col_name <- col_names[length(col_names)]

  # Define terms for the GLM formula
  main_terms <- paste(col_names, collapse="+")
  sq_terms <- paste0("I(", col_names[-length(col_names)], "^2)", collapse="+")
  prod_terms <- paste(
    utils::combn(col_names, 2, function(x) {
      sorted <- sort(x)
      paste0("I(", sorted[1], "*", sorted[2], ")")
    }),
    collapse = "+"
  )

  # Define the full formula as a string
  formula_str <- paste0("Y ~ ", main_terms, "+", prod_terms, "+", sq_terms)

  # Define the formula
  formula <- stats::as.formula(formula_str)

  # Define the data.frame for training
  data_train <- data.frame(Y = Y, X)

  # Make a matrix out of the formula
  formula_mat <- stats::model.matrix(formula, data = data_train)

  # Create starting values vector
  start_vals <- stats::setNames(
    rep(0, ncol(formula_mat)), colnames(formula_mat)
  )

  # Set coefficient for treatment and intercept to 1
  start_vals[treatment_col_name] <- 1
  start_vals[1] <- 1

  # Fit the glm
  suppressWarnings(
    fit_glm <-stats::glm(
      formula,
      data = data_train,
      family = stats::Gamma(link = "identity"),
      start = start_vals,
      maxit = 1000
    )
  )

  # Compute predictions using newX
  pred <- stats::predict(fit_glm, newdata = newX, type = "response")

  # Wrap and return
  fit = list(model = fit_glm)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- "SL.glm.gamma"

  return(out)
}


#' SuperLearner wrapper for GLM with Gamma family and log link
#'
#' @description
#' A SuperLearner wrapper that implements GLM using Gamma family
#' with a log link.
#'
#' @details
#' The predictor matrix `X` is assumed to have a binary
#' treatment column as its last column. Thus, we don't consider the last column
#' when squared terms are defined.
#'
#' @param Y A numeric `vector` of outcome values.
#' @param X A numeric `matrix` or `data.frame` of covariates and treatment.
#' @param newX A numeric `matrix` or `data.frame` of predictors.
#' @param ... Any additional arguments.
#' @return A list with components:
#' * `pred`: A numeric vector of predictions on `newX`.
#' * `fit`: A list containing the fitted model object.
SL.glm.gamma.log <- function(Y, X, newX, ...) {
  # Get column names of predictor matrix X
  col_names <- colnames(X)

  # Get treatment column name
  treatment_col_name <- col_names[length(col_names)]

  # Define terms for the GLM formula
  main_terms <- paste(col_names, collapse="+")
  sq_terms <- paste0("I(", col_names[-length(col_names)], "^2)", collapse="+")
  prod_terms <- paste(
    utils::combn(col_names, 2, function(x) {
      sorted <- sort(x)
      paste0("I(", sorted[1], "*", sorted[2], ")")
    }),
    collapse = "+"
  )
  # Define the full formula as a string
  formula_str <- paste0("Y ~ ", main_terms, "+", prod_terms, "+", sq_terms)

  # Define the formula
  formula <- stats::as.formula(formula_str)

  # Define the data.frame for training
  data_train <- data.frame(Y = Y, X)

  # Make a matrix out of the formula
  formula_mat <- stats::model.matrix(formula, data = data_train)

  # Create starting values vector
  start_vals <- stats::setNames(rep(0, ncol(formula_mat)), colnames(formula_mat))

  # Set coefficient for treatment and intercept to 1
  start_vals[treatment_col_name] <- 1
  start_vals[1] <- 1

  # Fit the glm
  fit_glm <-
    stats::glm(
      formula,
      data = data_train,
      family = stats::Gamma(link = "log"),
      start = start_vals,
      maxit = 100
    )

  # Compute predictions using newX
  pred <- stats::predict(fit_glm, newdata = newX, type = "response")

  # Wrap and return
  fit = list(model = fit_glm)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- "SL.glm.gamma"
  return(out)
}


#' SuperLearner wrapper for neural network using torch
#'
#' @description A SuperLearner wrapper that implements a neural network using
#'   torch, with a softplus activation to encourage strictly positive
#'   predictions. By default, the neural network uses a single hidden layer with
#'   10 hidden units, an Adam optimizer with a learning rate of 0.01, an MSE
#'   loss function, and 100 epochs.
#'
#' @param Y A numeric `vector` of outcome values.
#' @param X A numeric `matrix` or `data.frame` of covariates and treatment.
#' @param newX A numeric `matrix` or `data.frame` of predictors.
#' @param n_hidden Number of hidden units in the hidden layer, default is 10.
#' @param learning_rate Learning rate for the Adam optimizer, default is 0.01.
#' @param loss_fn The loss function used during training, default is MSE.
#' @param epochs The number of epochs for training, default is 100.
#' @param ... Any additional arguments.
#'
#' @return A list with components:
#' * `pred`: A numeric vector of predictions on `newX`.
#' * `fit`: A list containing the fitted model object.
#'
#' @export
SL.torch.softplus <- function(
    Y,
    X,
    newX,
    n_hidden=10,
    learning_rate=0.01,
    loss_fn=torch::nn_mse_loss(),
    epochs=100,
    ...
  ){

  # Use correct data structures
  X_mat <- as.matrix(X)
  newX_mat <- as.matrix(newX)
  Y_vec <- as.numeric(Y)

  # Define predictors and outcome for training in correct format
  x_train <- torch::torch_tensor(X_mat, dtype = torch::torch_float())
  y_train <- torch::torch_tensor(
    Y_vec, dtype = torch::torch_float()
  )$unsqueeze(2)

  # Define held out predictors for predictions in correct format
  x_new <- torch::torch_tensor(newX_mat, dtype = torch::torch_float())

  # Get input features
  n_input_features <- ncol(X_mat)

  # Define the neural network with one hidden layer and a softplus activation
  nnet <- torch::nn_module(
    initialize = function() {
      self$input_layer <- torch::nn_linear(n_input_features, n_hidden)
      self$outcome_layer <- torch::nn_linear(n_hidden, 1)
    },
    forward = function(x) {
      x |>
        self$input_layer() |>
        torch::nnf_softplus() |>
        self$outcome_layer()
    }
  )

  # Get neural network
  model <- nnet()

  # Define optimizer
  optimizer <- torch::optim_adam(model$parameters, lr = learning_rate)

  # Train the neural network
  for (epoch in 1:epochs) {
    model$train()
    optimizer$zero_grad()
    output <- model(x_train)
    loss <- loss_fn(output, y_train)
    loss$backward()
    optimizer$step()
  }

  # Switch to evaluation mode
  model$eval()

  # Compute predictions on held out predictors
  torch::with_no_grad({
    preds <- as.numeric(model(x_new)$squeeze())
  })

  # Wrap and return
  fit <- list(model = model)
  class(fit) <- "SL.torch.softplus"

  return(list(pred = preds, fit = fit))
}

#' SuperLearner wrapper of SL.xgboost that ensures positive predictions
#'
#' @description
#' A SuperLearner wrapper for SL.xgboost that ensures positive predictions. Any
#' predictions below 1e-4 are clipped to 1e-4.
#'
#' @param ... SL.xgboost arguments.
#' @param ntrees Number of trees, default is 100.
#' @param lower Lower bound to clip predictions to, default is 1e-4.
#' @return A list with components:
#' * `pred`: A numeric vector of predictions.
#' * `fit`: A list containing the fitted model object.
SL.xgboost.wrapper <- function(..., ntrees = 100, lower = 1e-4) {
  out <- SL.xgboost(..., ntrees=ntrees)
  out$pred <- pmax(out$pred, lower)
  return(out)
}
