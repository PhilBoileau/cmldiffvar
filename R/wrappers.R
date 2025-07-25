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
  fit <-list(object = fit_glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}


#' SuperLearner wrapper for GLM with Gamma family and log link
#'
#' @description A SuperLearner wrapper that implements generalized linear models
#'   using the Gamma family with a log link.
#'
#'
#' @details The predictor matrix `X` is assumed to have a binary treatment
#'   column as its last column. Thus, we don't consider the last column when
#'   squared terms are defined. In addition to the main terms, coefficients for
#'   the squared covariates are included in the model.
#'
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
      maxit = 1000
    )

  # Compute predictions using newX
  pred <- stats::predict(fit_glm, newdata = newX, type = "response")

  # Wrap and return
  fit <-list(object = fit_glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
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
#' @param lower Lower bound to clip predictions to, default is 1e-4
#' @param ... Any additional arguments.
#'
#' @return A list with components:
#' * `pred`: A numeric vector of predictions on `newX`.
#' * `fit`: A list containing the fitted model object.
#'
#' @export
SL.nnet.torch.softplus <- function(
    Y,
    X,
    newX,
    n_hidden=10,
    learning_rate=0.01,
    loss_fn=torch::nn_mse_loss(),
    epochs=100,
    lower = 1e-4,
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
    preds <- pmax(as.numeric(model(x_new)$squeeze()), lower)
  })

  # Wrap and return
  fit <- list(object = model)
  class(fit) <- "SL.nnet"
  out <- list(pred = preds, fit = fit)
  return(out)
}

#' SuperLearner wrapper of SL.xgboost that ensures positive predictions
#'
#' @description A SuperLearner wrapper for SL.xgboost that ensures positive
#'   predictions. Predictions below 1e-4 are clipped to a minimum value of 1e-4.
#'
#'
#' @param ... SL.xgboost arguments.
#' @param ntrees Number of trees, default is 100.
#' @param lower Lower bound to clip predictions to, default is 1e-4.
#'
#' @return A list with components:
#' * `pred`: A numeric vector of predictions.
#' * `fit`: A list containing the fitted model object.
#'
#' @export
SL.xgboost.bounded <- function(..., ntrees = 100, lower = 1e-4) {
  out <- SuperLearner::SL.xgboost(..., ntrees=ntrees)
  out$pred <- pmax(out$pred, lower)
  return(out)
}

#' SuperLearner wrapper for GAM with Gamma family and log link
#'
#' @description A SuperLearner wrapper that implements generalized additive
#'   models using the Gamma family with a log link.
#'
#' @details A copy of `[SuperLearner::SL.gam]` but using the Gammma family.
#'
#' @param Y A numeric `vector` of outcome values.
#' @param X A numeric `matrix` or `data.frame` of covariates and treatment.
#' @param obsWeights An optional vector of weights to be used in the fitting
#'   process.
#' @param deg.gam A numeric representing the degrees of the GAM. Defaults to 2.
#' @param cts.num A numeric indicating the minimum number of unique values a
#'   numeric covariate requires to be considered as a continuous variable.
#'   Defaults to 4.
#'
#' @return A list with components:
#' * `pred`: A numeric vector of predictions on `newX`.
#' * `fit`: A list containing the fitted model object.
#'
#' @export
SL.gam.gamma.log <- function(
    Y, X, newX, obsWeights, deg.gam = 2, cts.num = 4, ...
) {

  # requireNamespace() alone does not work. requireNamespace, unlike require(),
  # does not attached the package and allow the formula to parse correctly with
  # s(), gam::s() doesn't work, is not recognized as a special function
  if (!requireNamespace('gam'))
    stop("SL.gam requires the gam package, but it isn't available")

  # check if gam attached, if not, then attached
  if (!"package:gam" %in% search()) attachNamespace('gam')

  if ("mgcv" %in% loadedNamespaces())
    warning(paste0(
      "mgcv and gam packages are both in use. You might see an error because ",
      "both packages use the same function names."
    ))

  # create the formula for gam with a spline for each continuous variable
  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
  if (sum(!cts.x) > 0) {
    gam.model <- as.formula(
      paste(
        "Y~",
        paste(
          paste(
            "s(", colnames(X[, cts.x, drop = FALSE]), ",", deg.gam,")", sep=""
          ),
          collapse = "+"
        ),
        "+",
        paste(colnames(X[, !cts.x, drop=FALSE]), collapse = "+")
      )
    )
  } else {
    gam.model <- as.formula(
      paste(
        "Y~",
        paste(
          paste(
            "s(", colnames(X[, cts.x, drop = FALSE]), ",", deg.gam, ")", sep=""
          ),
          collapse = "+")
        )
      )
  }

  # fix for when all variables are binomial
  if (sum(!cts.x) == length(cts.x)) {
    gam.model <- as.formula(
      paste("Y~", paste(colnames(X), collapse = "+"), sep = "")
    )
  }

  fit.gam <- gam::gam(
    gam.model,
    data = X,
    family = Gamma(link = "log"),
    control = gam::gam.control(maxit = 50, bf.maxit = 50),
    weights = obsWeights
  )

  if(packageVersion('gam') >= "1.15") {
    # updated gam class in version 1.15
    pred <- gam::predict.Gam(fit.gam, newdata = newX, type = "response")
  } else {
    stop(
      paste0(
        "This SL.gam wrapper requires gam version >= 1.15, please update the",
        "gam package with 'update.packages('gam')'"
      )
    )
  }

  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gam")
  return(out)
}

#' SuperLearner wrapper for Multivariate Adaptive Regression Splines
#'
#' @description A SuperLearner wrapper that implements Multivariate Adaptive
#'   Regression Splines using the Gamma family with a log link.
#'
#' @details A copy of `[SuperLearner::SL.earth]()` but using the Gammma family,
#' which is itself a wrapper of `[earth::earth]()`. See `[earth::earth]()` for
#' more details.
#'
#' @param Y A numeric `vector` of outcome values.
#' @param X A numeric `matrix` or `data.frame` of covariates and treatment.
#' @param obsWeights An optional vector of weights to be used in the
#'   SuperLearner fitting process. Not used.
#' @param degree Maximum degree of interaction. Default is 3, meaning build a
#'   model with interaction terms.
#' @param penalty Generalized Cross Validation (GCV) penalty per knot. Defaults
#'   is 3.
#' @param nk Maximum number of model terms before pruning, i.e., the maximum
#'   number of terms created by the forward pass. Includes the intercept. The
#'   default is semi-automatically calculated from the number of predictors but
#'   may need adjusting.
#' @param pmethod Pruning method. One of: "backward", "none", "exhaustive",
#'   "forward", "seqrep", "cv". Default is "backward". Specify pmethod="cv" to
#'   use cross-validation to select the number of terms. This selects the number
#'   of terms that gives the maximum mean out-of-fold RSq on the fold models.
#'   Requires the nfold argument. Use "none" to retain all the terms created by
#'   the forward pass.
#' @param nfold Number of cross-validation folds. Default is 0, no cross
#'   validation.
#' @param ncross Only applies if nfold>1. Number of cross-validations. Each
#'   cross-validation has nfold folds. Default 1.
#' @param minspan Minimum number of observations between knots. (This increases
#'   resistance to runs of correlated noise in the input data.) The default
#'   minspan=0 is treated specially and means calculate the minspan internally,
#'   as per Friedman's MARS paper section 3.8 with alpha = 0.05.
#' @param endspan Minimum number of observations before the first and after the
#'   final knot. The default endspan=0 is treated specially and means calculate
#'   the endspan internally, as per the MARS paper equation 45 with  alpha =
#'   0.05.
#'
#' @return A list with components:
#' * `pred`: A numeric vector of predictions on `newX`.
#' * `fit`: A list containing the fitted model object.
#'
#' @export
SL.earth.gamma <- function(
    Y, X, newX, obsWeights, id, degree = 2, penalty = 3,
    nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", nfold = 0,
    ncross = 1, minspan = 0, endspan = 0, ...
) {

  fit.earth <- earth::earth(
    x = X, y = Y, degree = degree, nk = nk, penalty = penalty,
    pmethod = pmethod, nfold = nfold, ncross = ncross, minspan = minspan,
    endspan = endspan, glm = list(family = Gamma(link = "log"))
  )
  pred <- predict(fit.earth, newdata = newX, type = "response")
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)

}
