test_that("counterfactual table generator produces counterfactual tables", {

  # load required libraries
  library(dplyr)
  library(SuperLearner)

  # set seed for reproducibility
  set.seed(234566)

  # grab a sample of the population
  sample_tbl <- slice_sample(toy_population_tbl, n = 100)

  # fit nuisance parameter estimators
  ps_sl_fit <- estimate_propensity_score_fun(
    sample_tbl,
    adj_set_var_names = "confounder",
    treatment_var_name = "treatment",
    propensity_score_library = c("SL.glm", "SL.mean"),
    num_nuisance_sl_folds = 5
  )
  cond_exp_outcome_fit <- estimate_cond_exp_outcome_fun(
    sample_tbl,
    adj_set_var_names = "confounder",
    treatment_var_name = "treatment",
    outcome_var_name = "outcome",
    cond_exp_outcome_library = c("SL.glm", "SL.mean"),
    num_nuisance_sl_folds = 5
  )
  cond_exp_sq_outcome_fit <- estimate_cond_exp_sq_outcome_fun(
    sample_tbl,
    adj_set_var_names = "confounder",
    treatment_var_name = "treatment",
    outcome_var_name = "outcome",
    cond_exp_sq_outcome_library = c("SL.glm", "SL.glm.gamma.log"),
    num_nuisance_sl_folds = 5
  )

  # construct the counterfactual dataset under the treatment group
  sample_treatment_tbl <- generate_counterfactural_tbl_fun(
    sample_tbl,
    treatment_group = 1,
    propensity_score_adj_var_names = "confounder",
    cond_exp_outcome_adj_var_names = "confounder",
    treatment_var_name = "treatment",
    propensity_score_var_name = NULL,
    ps_sl_fit,
    cond_exp_outcome_fit,
    cond_exp_sq_outcome_fit
  )

  # all treatment indicators are set to treatment (1)
  expect(sum(sample_treatment_tbl$treatment), nrow(sample_treatment_tbl))

  # predictions from nuisance parameters match manual calculations
  expect_equal(
    as.numeric(sample_treatment_tbl$pred_propensity_score),
    as.numeric(
      predict(ps_sl_fit, newdata = sample_tbl |> select(confounder))$pred
    )
  )
  expect_equal(
    as.numeric(sample_treatment_tbl$pred_cond_exp_outcome),
    as.numeric(
      predict(
        cond_exp_outcome_fit,
        newdata = sample_tbl |>
          mutate(treatment = 1) |>
          select(treatment, confounder)
      )$pred
    )
  )
  expect_equal(
    as.numeric(sample_treatment_tbl$pred_cond_exp_sq_outcome),
    as.numeric(
      predict(
        cond_exp_sq_outcome_fit,
        newdata = sample_tbl |>
          mutate(treatment = 1) |>
          select(treatment, confounder)
      )$pred
    )
  )

})

test_that(
  "longitudinal counterfactual table generator produces counterfactual tables",
  {

    # load required libraries
    library(dplyr)
    library(SuperLearner)

    # set seed for reproducibility
    set.seed(841651)

    # grab a sample of the population and melt their data
    sample_tbl <- generate_test_data(n_obs = 50)
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

    # all treatment indicators are set to treatment (1)
    expect(
      sum(sample_long_treatment_tbl$counterfactual_treatment),
      nrow(sample_long_treatment_tbl)
    )

    # predictions from nuisance parameters match manual calculations
    expect_equal(
      as.numeric(sample_long_treatment_tbl$pred_propensity_score),
      as.numeric(
        predict(ps_sl_fit, newdata = long_sample_tbl |> select(w))$pred
      )
    )
    expect_equal(
      as.numeric(sample_long_treatment_tbl$pred_cond_event_haz),
      as.numeric(
        predict(
          cond_event_haz_sl_fit,
          newdata = long_sample_tbl |>
            mutate(a = 1) |>
            select(cmldiffvar_long_time, a, w)
        )$pred
      )
    )
    expect_equal(
      as.numeric(sample_long_treatment_tbl$pred_cond_censoring_haz),
      as.numeric(
        predict(
          cond_censoring_haz_sl_fit,
          newdata = long_sample_tbl |>
            mutate(a = 1) |>
            select(cmldiffvar_long_time, a, w)
        )$pred
      )
    )
  }
)
