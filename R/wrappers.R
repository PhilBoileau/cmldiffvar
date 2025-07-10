#' SuperLearner wrapper for GLM with Gamma family
#'
#' A SuperLearner wrapper that implements GLM using Gamma family.
#'
#' @param Y A numeric `vector` of outcome values.
#' @param X A numeric `matrix` or `data.frame` of covariates and treatment.
#' @param newX A numeric `matrix` or `data.frame` of predictors.
#' @return A list with components:
#' \describe{
#'   \item{pred}{A numeric vector of predictions on \code{newX}.}
#'   \item{fit}{A list containing the fitted model object.}}
SL.glm.gamma <- function(Y, X, newX, ...) {
  # Getting column names of predictor matrix X
  col_names <- colnames(X)

  # Defining terms for the GLM formula
  main_terms <- paste(col_names, collapse="+")
  sq_terms <- paste0("I(", col_names[-length(col_names)], "^2)", collapse="+")
  prod_terms <- paste(
    combn(col_names, 2, function(x) {
      sorted <- sort(x)
      paste0("I(", sorted[1], "*", sorted[2], ")")
    }),
    collapse = "+"
  )
  formula_str <- paste0("Y ~ ", main_terms, "+", prod_terms, "+", sq_terms)

  # Defining the formula
  formula <- as.formula(formula_str)

  # Defining the data.frame for training
  df <- data.frame(Y = Y, X)

  # Computing starting values using a linear model
  lm_fit <- lm(formula, data = df)
  start_vals <- coef(lm_fit)
  pred_vals <- predict(lm_fit)

  # If starting values are not valid, generate valid ones
  if(any(pred_vals <= 0)){
    start_vals <- get_valid_glm_gamma_start(formula, data=df)
  }

  # Fit the glm
  fit_glm <-
    glm2::glm2(
      formula,
      data = df,
      family = Gamma(link = "identity"),
      start = start_vals,
      maxit = 100
    )

  # Compute predictions using newX
  pred <- predict(fit_glm, newdata = newX, type = "response")

  # Wrap and return
  fit = list(model = fit_glm)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- "SL.glm.gamma"
  return(out)
}

#' SuperLearner wrapper for GAM with Gamma family
#'
#' A SuperLearner wrapper that implements GAM using Gamma family and log link.
#'
#' @param Y A numeric `vector` of outcome values.
#' @param X A numeric `matrix` or `data.frame` of covariates and treatment.
#' @param newX A numeric `matrix` or `data.frame` of predictors.
#' @return A list with components:
#' \describe{
#'   \item{pred}{A numeric vector of predictions on \code{newX}.}
#'   \item{fit}{A list containing the fitted model object.}}
SL.gam.gamma <- function(Y, X, newX, ...) {
  # Define data.frame for training
  data_train <- data.frame(Y = Y, X)

  # Fit the GAM
  fit_gam <- gam(
    Y ~ .,
    data = data_train,
    family = Gamma(link = "log"),
    maxit=100
  )
  # Compute predictions
  pred <- predict(fit_gam, newdata = newX, type = "response")

  # Wrap and return
  fit = list(model = fit_gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- "SL.gam.gamma"
  return(out)
}

#' SuperLearner wrapper for neural network
#'
#' A SuperLearner wrapper that implements a neural network using torch and using
#' a softplus activation to encourage strictly positive predictons.
#'
#' @param Y A numeric `vector` of outcome values.
#' @param X A numeric `matrix` or `data.frame` of covariates and treatment.
#' @param newX A numeric `matrix` or `data.frame` of predictors.
#' @return A list with components:
#' \describe{
#'   \item{pred}{A numeric vector of predictions on \code{newX}.}
#'   \item{fit}{A list containing the fitted model object.}}
SL.torch.softplus <- function(Y, X, newX, ...) {
  # Load required library
  library(torch)

  # Use correct data structures
  X_mat <- as.matrix(X)
  newX_mat <- as.matrix(newX)
  Y_vec <- as.numeric(Y)

  # Define predictors and outcome for training in correct format
  x_train <- torch_tensor(X_mat, dtype = torch_float())
  y_train <- torch_tensor(Y_vec, dtype = torch_float())$unsqueeze(2)

  # Define held out predictors for predictions in correct format
  x_new <- torch_tensor(newX_mat, dtype = torch_float())

  # Get input features
  n_input_features <- ncol(X_mat)

  # Define the number of hidden neurons
  n_hidden <- 10

  # Define the neural network with one hidden layer and a softplus activation
  nnet <- nn_module(
    initialize = function() {
      self$input_layer <- nn_linear(n_input_features, n_hidden)
      self$outcome_layer <- nn_linear(n_hidden, 1)
    },
    forward = function(x) {
      x |>
        self$input_layer() |>
        nnf_softplus() |>
        self$outcome_layer()
    }
  )

  # Get neural network
  model <- nnet()

  # Define optimizer
  optimizer <- optim_adam(model$parameters, lr = 0.01)

  # Define loss function
  loss_fn <- nn_mse_loss()

  # Train the neural network
  epochs <- 100
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
  with_no_grad({
    preds <- as.numeric(model(x_new)$squeeze())
  })

  # Wrap and return
  fit <- list(model = model)
  class(fit) <- "SL.torch.softplus"

  return(list(pred = preds, fit = fit))
}
