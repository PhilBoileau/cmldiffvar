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


test_that("one-step marginal survival estimator is consistent", {

  # load required libraries
  library(dplyr)
  library(SuperLearner)

  # set seed for reproducibility
  set.seed(1681615)

  # grab a sample of the population and melt their data
  survival_estimate_vec <- sapply(
    seq_len(100),
    function(iter_idx) {
      sample_tbl <- generate_test_data(n_obs = 500, null_marginal = FALSE)
      long_sample_tbl <- sample_tbl |>
        melt_tte_data_fun(
          baseline_var_names = "w",
          treatment_var_name = "a",
          outcome_var_name = "time",
          censoring_var_name = "censoring",
          time_cutoff = 50
        )

      # fit nuisance parameter estimators
      ps_sl_fit <- estimate_propensity_score_fun(
        sample_tbl,
        adj_set_var_names = "w",
        treatment_var_name = "a",
        propensity_score_library = c("SL.glm", "SL.mean"),
        num_nuisance_sl_folds = 5
      )
      cond_event_haz_sl_fit <- estimate_cond_event_haz_fun(
        clean_long_tbl = long_sample_tbl,
        adj_set_var_names = "w",
        treatment_var_name = "a",
        outcome_var_name = "L",
        cond_event_haz_library = c("SL.glm"),
        num_nuisance_sl_folds = 5
      )
      cond_censoring_haz_sl_fit <- estimate_cond_censoring_haz_fun(
        clean_long_tbl = long_sample_tbl,
        adj_set_var_names = "w",
        treatment_var_name = "a",
        outcome_var_name = "R",
        cond_censoring_haz_library = c("SL.glm"),
        num_nuisance_sl_folds = 5
      )

      # construct the counterfactual dataset under the treatment group
      sample_long_treatment_tbl <- generate_long_counterfactural_tbl_fun(
        clean_long_tbl = long_sample_tbl,
        treatment_group = 1,
        propensity_score_adj_var_names = "w",
        cond_event_haz_adj_var_names = "w",
        cond_censoring_haz_adj_var_names = "w",
        treatment_var_name = "a",
        propensity_score_var_name = NULL,
        propensity_score_sl_fit = ps_sl_fit,
        cond_event_haz_sl_fit = cond_event_haz_sl_fit,
        cond_censoring_haz_sl_fit = cond_censoring_haz_sl_fit
      )

      # extract the conditional survival time estimate at time = 50 in long
      # format
      num_unique_times <- sample_long_treatment_tbl |>
        pull(cmldiffvar_long_time) |>
        unique() |>
        length()
      cond_surv_at_trunc_est_vec <- sample_long_treatment_tbl |>
        filter(cmldiffvar_long_time <= 50) |>
        group_by(cmldiffvar_id) |>
        slice_tail(n = 1) |>
        select(cmldiffvar_id, pred_cond_event_survival) |>
        ungroup() |>
        slice(rep(1:n(), each = num_unique_times)) |>
        pull(pred_cond_event_survival)

      # estimate the survival probability among the treated at time = 50
      one_step_marginal_survival_estimator_fun(
        treatment_group = 1,
        cmldiffvar_id = sample_long_treatment_tbl$cmldiffvar_id,
        cmldiffvar_time_long = sample_long_treatment_tbl$cmldiffvar_long_time,
        time_cutoff = 50,
        treatment_vec = sample_long_treatment_tbl$a,
        ps_est_vec = sample_long_treatment_tbl$pred_propensity_score,
        cond_survival_est_vec =
          sample_long_treatment_tbl$pred_cond_event_survival,
        cond_surv_at_trunc_est_vec = cond_surv_at_trunc_est_vec,
        cond_censoring_survival_est_vec =
          sample_long_treatment_tbl$pred_cond_censoring_survival,
        I_vec = sample_long_treatment_tbl$I,
        L_vec = sample_long_treatment_tbl$L,
        cond_event_haz_est_vec = sample_long_treatment_tbl$pred_cond_event_haz
      )
    }
  )

  # approximate true survival probability under treatment at time=50
  population_tbl <- generate_test_data(n_obs = 10000, null_marginal = FALSE)
  pop_surv_at_50 <- mean(population_tbl$potential_time_1 >= 50)

  # make sure the empirical bias is small
  expect_lte(abs(mean(survival_estimate_vec) - pop_surv_at_50), 0.025)

})
