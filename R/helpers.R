#' Valid starting value searcher
#'
#' Automatically look for valid starting values for Gamma GLMs, since Gamma
#' assumes that \eqn{X\beta > 0}.
#'
#' @param formula A `formula` object describing the fit.
#' @param data The data used for the formula.
#' @param max_tries Maximum number of iterations to find valid starting values,
#' default is 100.
#' @param verbose Flag to determine if the function should print the attempt at
#' which it found valid starting values, default is FALSE.
#' @return `start_vals` containing valid starting values.
get_valid_glm_gamma_start <- function(
    formula, data,
    max_tries = 100,
    verbose = FALSE
){
  # Define design matrix
  X <- stats::model.matrix(formula, data)

  # Get numver of coefficients
  num_coef <- ncol(X)

  # Generate good starting values
  attempt <- 1
  while (attempt <= max_tries) {
    start_vals <- stats::runif(num_coef, min = 0.1, max = 1)

    pred_vals <- X %*% start_vals

    if (all(pred_vals > 0)) {
      if (verbose)
        print(paste("Found valid starting values on attempt: ", attempt))
      return(start_vals)
    }

    attempt <- attempt + 1
  }
}
