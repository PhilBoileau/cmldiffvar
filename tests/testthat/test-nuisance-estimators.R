test_that("propensity score estimator wrapper learns from data", {

  # load required libraries
  library(SuperLearner)
  library(dplyr)

  # set seed for reproducibility
  set.seed(234642)

  # compute negative log loss risks of a correctly and misspecified SL
  num_iters <- 100
  neg_log_loss_risk_tbl <- lapply(
    seq_len(num_iters),
    function(iter_idx) {

      # grab a sample of the population
      sample_tbl <- slice_sample(toy_population_tbl, n = 100)

      # estimate propensity score witch correctly specified model
      ps_sl_correct_fit <- estimate_propensity_score_fun(
        sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        propensity_score_library = c("SL.glm", "SL.mean"),
        num_nuisance_sl_folds = 5
      )

      # estimate propensity score with misspecified model
      ps_sl_misspecified_fit <- estimate_propensity_score_fun(
        sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        propensity_score_library = c("SL.mean"),
        num_nuisance_sl_folds = 5
      )

      # compute negative log likelihood risks
      ps_sl_correct_fit_risk <- -mean(
        sample_tbl$treatment * log(ps_sl_correct_fit$SL.predict) +
          (1 - sample_tbl$treatment) * log(1 - ps_sl_correct_fit$SL.predict)
      )
      ps_sl_misspecified_fit_risk <- -mean(
        sample_tbl$treatment * log(ps_sl_misspecified_fit$SL.predict) +
          (1 - sample_tbl$treatment) *
          log(1 - ps_sl_misspecified_fit$SL.predict)
      )

      # return a tibble row with results
      return(tibble(
        "ps_sl_correct_fit_risk" = ps_sl_correct_fit_risk,
        "ps_sl_misspecified_fit_risk" = ps_sl_misspecified_fit_risk)
      )
    }
  ) |>
    bind_rows()

  # ensure that correctly specified model is outperforming misspecified model
  expect_true(
    mean(neg_log_loss_risk_tbl$ps_sl_correct_fit_risk) <
      mean(neg_log_loss_risk_tbl$ps_sl_misspecified_fit_risk)
  )

})


test_that("conditional expected outcome estimator wrapper learns from data", {

  # load required libraries
  library(SuperLearner)
  library(dplyr)

  # set seed for reproducibility
  set.seed(234642)

  # compute negative log loss risks of a correctly and misspecified SL
  num_iters <- 100
  mse_tbl <- lapply(
    seq_len(num_iters),
    function(iter_idx) {

      # grab a sample of the population
      sample_tbl <- slice_sample(toy_population_tbl, n = 100)

      # estimate propensity score witch correctly specified model
      cond_exp_outcome_correct_fit <- estimate_cond_exp_outcome_fun(
        sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        cond_exp_outcome_library = c("SL.glm", "SL.mean"),
        num_nuisance_sl_folds = 5
      )

      # estimate propensity score with misspecified model
      cond_exp_outcome_misspecified_fit <- estimate_cond_exp_outcome_fun(
        sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        cond_exp_outcome_library = c("SL.mean"),
        num_nuisance_sl_folds = 5
      )

      # compute mean squared errors
      cond_exp_outcome_sl_correct_fit_risk <- mean(
        (sample_tbl$outcome - cond_exp_outcome_correct_fit$SL.predict)^2
      )
      cond_exp_outcome_sl_misspecified_fit_risk <- mean(
        (sample_tbl$outcome - cond_exp_outcome_misspecified_fit$SL.predict)^2
      )

      # return a tibble row with results
      return(tibble(
        "cond_exp_outcome_sl_correct_fit_risk" =
          cond_exp_outcome_sl_correct_fit_risk,
        "cond_exp_outcome_sl_misspecified_fit_risk" =
          cond_exp_outcome_sl_misspecified_fit_risk)
      )
    }
  ) |>
    bind_rows()

  # ensure that correctly specified model is outperforming misspecified model
  expect_true(
    mean(mse_tbl$cond_exp_outcome_sl_correct_fit_risk) <
      mean(mse_tbl$cond_exp_outcome_sl_misspecified_fit_risk)
  )

})


test_that("conditional expected outcome^2 estimator wrapper learns from data", {

  # load required libraries
  library(SuperLearner)
  library(earth)
  library(dplyr)

  # set seed for reproducibility
  set.seed(234642)

  # compute negative log loss risks of a correctly and misspecified SL
  num_iters <- 100
  mse_tbl <- lapply(
    seq_len(num_iters),
    function(iter_idx) {

      # grab a sample of the population
      sample_tbl <- slice_sample(toy_population_tbl, n = 100)

      # estimate propensity score witch correctly specified model
      cond_exp_sq_outcome_correct_fit <- estimate_cond_exp_sq_outcome_fun(
        sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        cond_exp_sq_outcome_library = c("SL.glm", "SL.earth"),
        num_nuisance_sl_folds = 5
      )

      # estimate propensity score with misspecified model
      cond_exp_sq_outcome_misspecified_fit <- estimate_cond_exp_sq_outcome_fun(
        sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        cond_exp_sq_outcome_library = c("SL.mean"),
        num_nuisance_sl_folds = 5
      )

      # compute mean squared errors
      cond_exp_sq_outcome_sl_correct_fit_risk <- mean(
        (sample_tbl$outcome^2 - cond_exp_sq_outcome_correct_fit$SL.predict)^2
      )
      cond_exp_sq_outcome_sl_misspecified_fit_risk <- mean(
        (sample_tbl$outcome^2 -
           cond_exp_sq_outcome_misspecified_fit$SL.predict)^2
      )

      # return a tibble row with results
      return(tibble(
        "cond_exp_sq_outcome_sl_correct_fit_risk" =
          cond_exp_sq_outcome_sl_correct_fit_risk,
        "cond_exp_sq_outcome_sl_misspecified_fit_risk" =
          cond_exp_sq_outcome_sl_misspecified_fit_risk)
      )
    }
  ) |>
    bind_rows()

  # ensure that correctly specified model is outperforming misspecified model
  expect_true(
    mean(mse_tbl$cond_exp_sq_outcome_sl_correct_fit_risk) <
      mean(mse_tbl$cond_exp_sq_outcome_sl_misspecified_fit_risk)
  )

})


test_that("conditional expected outcome^2 estimator SL.glm.gamma.identity
          wrapper predicts > 0", {

  # Load required libraries
  library(SuperLearner)
  library(dplyr)

  # Set seed for reproducibility
  set.seed(234642)

  # Compute SL predictions
  num_iters <- 100
  pred_tbl <- lapply(
    seq_len(num_iters),
    function(iter_idx) {

      # Grab a sample of the population
      sample_tbl <- slice_sample(toy_population_tbl, n = 100)

      # Estimate propensity score witch correctly specified model
      cond_exp_sq_outcome_correct_fit <- estimate_cond_exp_sq_outcome_fun(
        sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        cond_exp_sq_outcome_library = c("SL.glm.gamma.identity"),
        num_nuisance_sl_folds = 5
      )

      # Estimate propensity score with misspecified model
      cond_exp_sq_outcome_misspecified_fit <- estimate_cond_exp_sq_outcome_fun(
        sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        cond_exp_sq_outcome_library = c("SL.glm"),
        num_nuisance_sl_folds = 5
      )

      # Return a tibble row with results
      return(tibble(
        "EY2_hat_correct_fit" =
          cond_exp_sq_outcome_correct_fit$SL.predict,
        "EY2_hat_misspecified_fit" =
          cond_exp_sq_outcome_misspecified_fit$SL.predict,
        "treatment" = sample_tbl$treatment)
      )
    }
  ) |>
    bind_rows()

  # Ensure that correctly specified model predictions are strictly positive
  expect_true(all(pred_tbl$EY2_hat_correct_fit > 0))

  # Ensure that misspecified model predictions are not all strictly positive
  expect_true(min(pred_tbl$EY2_hat_misspecified_fit) <= 0)
})


test_that("conditional expected outcome^2 estimator SL.glm.gamma.log
          wrapper predicts > 0", {

  # Load required libraries
  library(SuperLearner)
  library(dplyr)

  # Set seed for reproducibility
  set.seed(234642)

  # Compute SL predictions
  num_iters <- 100
  pred_tbl <- lapply(
    seq_len(num_iters),
    function(iter_idx) {

      # Grab a sample of the population
      sample_tbl <- slice_sample(toy_population_tbl, n = 100)

      # Estimate propensity score
      cond_exp_sq_outcome_fit <- estimate_cond_exp_sq_outcome_fun(
        sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        cond_exp_sq_outcome_library = c("SL.glm.gamma.log"),
        num_nuisance_sl_folds = 5
      )

      # Return a tibble row with results
      return(tibble(
        "EY2_hat_fit" =
          cond_exp_sq_outcome_fit$SL.predict,
        "treatment" = sample_tbl$treatment)
      )
    }
  ) |>
    bind_rows()

  # Ensure that model predictions are strictly positive
  expect_true(all(pred_tbl$EY2_hat_fit > 0))
})


test_that("conditional expected outcome^2 estimator SL.torch.softplus wrapper
          predicts > 0", {

  # Load required libraries
  library(SuperLearner)
  library(dplyr)
  library(torch)

  # Set seed for reproducibility
  set.seed(234642)

  # Compute SL predictions
  num_iters <- 100
  pred_tbl <- lapply(
    seq_len(num_iters),
    function(iter_idx) {

      # Grab a sample of the population
      sample_tbl <- slice_sample(toy_population_tbl, n = 100)

      # Estimate propensity score
      cond_exp_sq_outcome_fit <- estimate_cond_exp_sq_outcome_fun(
        sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        cond_exp_sq_outcome_library = c("SL.torch.softplus"),
        num_nuisance_sl_folds = 5
      )

      # Return a tibble row with results
      return(tibble(
        "EY2_hat_fit" =
          cond_exp_sq_outcome_fit$SL.predict,
        "treatment" = sample_tbl$treatment)
      )
    }
  ) |>
    bind_rows()

  # Ensure that model predictions are strictly positive
  expect_true(all(pred_tbl$EY2_hat_fit > 0))
})


test_that("conditional expected outcome^2 estimator SL.xgboost.wrapper wrapper
          predicts > 0", {

  # Load required libraries
  library(SuperLearner)
  library(dplyr)
  library(xgboost)

  # Set seed for reproducibility
  set.seed(234642)

  # Compute SL predictions
  num_iters <- 100
  pred_tbl <- lapply(
    seq_len(num_iters),
    function(iter_idx) {

      # Grab a sample of the population
      sample_tbl <- slice_sample(toy_population_tbl, n = 100)

      # Estimate propensity score
      cond_exp_sq_outcome_fit <- estimate_cond_exp_sq_outcome_fun(
        sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        cond_exp_sq_outcome_library = c("SL.xgboost.wrapper"),
        num_nuisance_sl_folds = 5
      )

      # Return a tibble row with results
      return(tibble(
        "EY2_hat_fit" =
          cond_exp_sq_outcome_fit$SL.predict,
        "treatment" = sample_tbl$treatment)
      )
    }
  ) |>
    bind_rows()

  # Ensure that model predictions are strictly positive
  expect_true(all(pred_tbl$EY2_hat_fit > 0))
})


test_that("conditional expected outcome^2 estimator SL.gam.gamma.log wrapper
          predicts > 0", {

  # Load required libraries
  library(SuperLearner)
  library(dplyr)
  library(mgcv)

  # Set seed for reproducibility
  set.seed(234642)

  # Compute SL predictions
  num_iters <- 100
  pred_tbl <- lapply(
    seq_len(num_iters),
    function(iter_idx) {

      # Grab a sample of the population
      sample_tbl <- slice_sample(toy_population_tbl, n = 100)

      # Estimate propensity score
      cond_exp_sq_outcome_fit <- estimate_cond_exp_sq_outcome_fun(
        sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        cond_exp_sq_outcome_library = c("SL.gam.gamma.log"),
        num_nuisance_sl_folds = 5
      )

      # Return a tibble row with results
      return(tibble(
        "EY2_hat_fit" =
          cond_exp_sq_outcome_fit$SL.predict,
        "treatment" = sample_tbl$treatment)
      )
    }
  ) |>
    bind_rows()

  # Ensure that model predictions are strictly positive
  expect_true(all(pred_tbl$EY2_hat_fit > 0))
})


test_that("conditional expected outcome^2 estimator SL.glm.gamma wrapper
          has consistent coefficients under the assumption of normally
          distributed outcomes", {

  # Load required libraries
  library(SuperLearner)
  library(dplyr)

  # Set seed for reproducibility
  set.seed(234642)

  # Compute SL predictions
  num_iters <- 1e3
  coef_list <- list()
  for(iter_idx in 1:num_iters){
    # Grab a sample of the population
    sample_tbl <- slice_sample(toy_population_tbl, n = 100)

    # Estimate propensity score with correctly specified model
    cond_exp_sq_outcome_correct_fit <- estimate_cond_exp_sq_outcome_fun(
      sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      outcome_var_name = "outcome",
      cond_exp_sq_outcome_library = c("SL.glm.gamma.identity"),
      num_nuisance_sl_folds = 5
    )

    # Grab the fit
    fit <- cond_exp_sq_outcome_correct_fit$fitLibrary$SL.glm.gamma_All$model

    # Save the coefficients of the fit
    coef_list[[iter_idx]] <- coef(fit)
  }

  # Compute the mean and the standard error of the coefficients
  coef_mat <- matrix(unlist(coef_list), ncol = length(coef_list[[1]]), byrow=T)
  coef_mean <- colMeans(coef_mat)
  coef_sd <- apply(coef_mat, 2, sd)

  # Compute lower and upper bound of each coefficient
  lower <- coef_mean - 1.96 * coef_sd
  upper <- coef_mean + 1.96 * coef_sd

  # Define true values.
  # Functional form is: E[Y^2|A,W] = 2 + 2W + 16A + 4AW + W^2
  true_values <- c(2, 2, 16, 4, 1)

  # ensure true values are within their respective bounds
  expect_true(all(true_values >= lower & true_values <= upper))
})


test_that("SL.glm.gamma wrapper adapts to multiple confounders", {

  # Load required libraries
  library(SuperLearner)
  library(dplyr)

  # Set seed for reproducibility
  set.seed(234642)

  # Grab a sample of the population
  sample_tbl <- slice_sample(multiple_confounders_toy_population_tbl, n = 100)

  # Estimate propensity score with correctly specified model
  cond_exp_sq_outcome_correct_fit <- estimate_cond_exp_sq_outcome_fun(
    sample_tbl,
    confounder_var_names = c("confounder_1", "confounder_2", "confounder_3"),
    treatment_var_name = "treatment",
    outcome_var_name = "outcome",
    cond_exp_sq_outcome_library = c("SL.glm.gamma.identity"),
    num_nuisance_sl_folds = 5
  )

  # Grab the fit
  fit <- cond_exp_sq_outcome_correct_fit$fitLibrary$SL.glm.gamma_All$model

  # Get coefficient names
  colnames <- names(fit$coefficients)
  real_colnames <- c(
    "(Intercept)",
    "treatment", "confounder_1", "confounder_2", "confounder_3",
    "I(confounder_1 * treatment)", "I(confounder_2 * treatment)",
    "I(confounder_3 * treatment)",
    "I(confounder_1 * confounder_2)", "I(confounder_1 * confounder_3)",
    "I(confounder_2 * confounder_3)",
    "I(confounder_1^2)", "I(confounder_2^2)", "I(confounder_3^2)"
  )

  # Expect all factors
  expect_true(setequal(colnames, real_colnames))
})
