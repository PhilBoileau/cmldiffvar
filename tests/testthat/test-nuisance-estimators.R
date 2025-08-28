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
        adj_set_var_names = "confounder",
        treatment_var_name = "treatment",
        propensity_score_library = c("SL.glm", "SL.mean"),
        num_nuisance_sl_folds = 5
      )

      # estimate propensity score with misspecified model
      ps_sl_misspecified_fit <- estimate_propensity_score_fun(
        sample_tbl,
        adj_set_var_names = "confounder",
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
        adj_set_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        cond_exp_outcome_library = c("SL.glm", "SL.mean"),
        num_nuisance_sl_folds = 5
      )

      # estimate propensity score with misspecified model
      cond_exp_outcome_misspecified_fit <- estimate_cond_exp_outcome_fun(
        sample_tbl,
        adj_set_var_names = "confounder",
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
        adj_set_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        cond_exp_sq_outcome_library = c("SL.glm", "SL.earth"),
        num_nuisance_sl_folds = 5
      )

      # estimate propensity score with misspecified model
      cond_exp_sq_outcome_misspecified_fit <- estimate_cond_exp_sq_outcome_fun(
        sample_tbl,
        adj_set_var_names = "confounder",
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
