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


test_that("conditional event hazard estimator wrapper learns from data", {

  # load required libraries
  library(SuperLearner)
  library(dplyr)

  # set seed for reproducibility
  set.seed(86181)

  # compute negative log loss risks of an adjusted and unadjusted estimator
  num_iters <- 50
  neg_log_loss_risk_tbl <- lapply(
    seq_len(num_iters),
    function(iter_idx) {

      # grab a sample of the population
      sample_tbl <- generate_test_data(n_obs = 50) |>
        melt_tte_data_fun(
          baseline_var_names = "w",
          treatment_var_name = "a",
          outcome_var_name = "time",
          censoring_var_name = "censoring",
          time_cutoff = 50
        )

      # estimate propensity score witch adjusted estimator
      cond_event_haz_sl_adj_fit <- estimate_cond_event_haz_fun(
        clean_long_tbl = sample_tbl,
        adj_set_var_names = "w",
        treatment_var_name = "a",
        outcome_var_name = "L",
        cond_event_haz_library = c("SL.ranger"),
        num_nuisance_sl_folds = 5
      )

      # estimate propensity score with unadjusted estimator
      cond_event_haz_sl_unadj_fit <- estimate_cond_event_haz_fun(
        clean_long_tbl = sample_tbl,
        adj_set_var_names = "w",
        treatment_var_name = "a",
        outcome_var_name = "L",
        cond_event_haz_library = c("SL.mean"),
        num_nuisance_sl_folds = 5
      )

      # compute negative log likelihood risks
      outcome_vec <- sample_tbl |>
        filter(I == 1) |>
        pull(L)
      adj_fit_pred_vec <- cond_event_haz_sl_adj_fit$SL.predict |>
        as.vector()
      cond_event_haz_sl_adj_fit_risk <- -mean(
        outcome_vec * log(adj_fit_pred_vec) +
          (1 - outcome_vec) * log(1 - adj_fit_pred_vec)
      )
      unadj_fit_pred_vec <-
        cond_event_haz_sl_unadj_fit$SL.predict |>
        as.vector()
      cond_event_haz_sl_unadj_fit_risk <- -mean(
        outcome_vec * log(unadj_fit_pred_vec) +
          (1 - outcome_vec) * log(1 - unadj_fit_pred_vec)
      )

      # return a tibble row with results
      return(
        tibble(
          "cond_event_haz_sl_adj_fit_risk" =
            cond_event_haz_sl_adj_fit_risk,
          "cond_event_haz_sl_unadj_fit_risk" =
            cond_event_haz_sl_unadj_fit_risk
        )
      )
    }
  ) |>
    bind_rows()

  # ensure that adjusted estimator is outperforming unadjusted
  expect_true(
    mean(neg_log_loss_risk_tbl$cond_event_haz_sl_adj_fit_risk) <
      mean(neg_log_loss_risk_tbl$cond_event_haz_sl_unadj_fit_risk)
  )

})


test_that("conditional censoring hazard estimator wrapper learns from data", {

  # load required libraries
  library(SuperLearner)
  library(dplyr)

  # set seed for reproducibility
  set.seed(81345)

  # compute negative log loss risks of an adjusted and unadjusted estimator
  num_iters <- 50
  neg_log_loss_risk_tbl <- lapply(
    seq_len(num_iters),
    function(iter_idx) {

      # grab a sample of the population
      sample_tbl <- generate_test_data(n_obs = 50) |>
        melt_tte_data_fun(
          baseline_var_names = "w",
          treatment_var_name = "a",
          outcome_var_name = "time",
          censoring_var_name = "censoring",
          time_cutoff = 50
        )

      # estimate propensity score witch adjusted estimator
      cond_censoring_haz_sl_adj_fit <- estimate_cond_censoring_haz_fun(
        clean_long_tbl = sample_tbl,
        adj_set_var_names = "w",
        treatment_var_name = "a",
        outcome_var_name = "R",
        cond_censoring_haz_library = c("SL.ranger"),
        num_nuisance_sl_folds = 5
      )

      # estimate propensity score with unadjusted estimator
      cond_censoring_haz_sl_unadj_fit <- estimate_cond_censoring_haz_fun(
        clean_long_tbl = sample_tbl,
        adj_set_var_names = "w",
        treatment_var_name = "a",
        outcome_var_name = "R",
        cond_censoring_haz_library = c("SL.mean"),
        num_nuisance_sl_folds = 5
      )

      # compute negative log likelihood risks
      outcome_vec <- sample_tbl |>
        filter(J == 1) |>
        pull(R)
      adj_fit_pred_vec <- cond_censoring_haz_sl_adj_fit$SL.predict |>
        as.vector()
      cond_censoring_haz_sl_adj_fit_risk <- -mean(
        outcome_vec * log(adj_fit_pred_vec) +
          (1 - outcome_vec) * log(1 - adj_fit_pred_vec)
      )
      unadj_fit_pred_vec <-
        cond_censoring_haz_sl_unadj_fit$SL.predict |>
        as.vector()
      cond_censoring_haz_sl_unadj_fit_risk <- -mean(
        outcome_vec * log(unadj_fit_pred_vec) +
          (1 - outcome_vec) * log(1 - unadj_fit_pred_vec)
      )

      # return a tibble row with results
      return(
        tibble(
          "cond_censoring_haz_sl_adj_fit_risk" =
            cond_censoring_haz_sl_adj_fit_risk,
          "cond_censoring_haz_sl_unadj_fit_risk" =
            cond_censoring_haz_sl_unadj_fit_risk
        )
      )
    }
  ) |>
    bind_rows()

  # ensure that adjusted estimator is outperforming unadjusted
  expect_true(
    mean(neg_log_loss_risk_tbl$cond_censoring_haz_sl_adj_fit_risk) <
      mean(neg_log_loss_risk_tbl$cond_censoring_haz_sl_unadj_fit_risk)
  )

})
