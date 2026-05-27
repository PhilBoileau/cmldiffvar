test_that(
  "tte_cmldiffvar runs without issue in null homogeneous effect setting", {

  # load required libraries
  library(dplyr)
  library(SuperLearner)

  # set seed for reproducibility
  set.seed(61168)

  # sample a dataset
  sample_tbl <- generate_test_data(n_obs = 1000, null_marginal = TRUE)

  # produce differential variance inference
  tte_diff_var_results <- tte_cmldiffvar(
    data_tbl = sample_tbl,
    homogeneity_type = "additive",
    restriction_time = 100,
    homogeneous_effect_range_min = -5,
    homogeneous_effect_range_max = 5,
    propensity_score_adj_var_names = "w",
    cond_event_haz_adj_var_names = "w",
    cond_censoring_haz_adj_var_names = "w",
    treatment_var_name = "a",
    outcome_var_name = "time",
    censoring_var_name = "censoring",
    propensity_score_sl_library = c("SL.glm", "SL.mean", "SL.ranger"),
    cond_event_haz_sl_library = c("SL.glm", "SL.mean", "SL.ranger"),
    cond_censoring_haz_sl_library = c("SL.glm", "SL.mean", "SL.ranger")
  )

  # check that null hypothesis is not rejected
  expect_true(tte_diff_var_results$hypothesis_test_tbl$p_value > 0.05)
})

test_that(
  paste("tte_cmldiffvar runs without issue in null homogeneous multiplicative",
        "effect setting"),
{

  # load required libraries
  library(dplyr)
  library(SuperLearner)

  # set seed for reproducibility
  set.seed(3516813)

  # sample a dataset
  sample_tbl <- generate_test_data(
    n_obs = 600, null_marginal = TRUE, hazard_model = FALSE
  )

  # produce differential variance inference
  tte_diff_var_results <- tte_cmldiffvar(
    data_tbl = sample_tbl,
    homogeneity_type = "multiplicative",
    restriction_time = 90,
    homogeneous_effect_range_min = 0.95,
    homogeneous_effect_range_max = 1.2,
    propensity_score_adj_var_names = "w",
    cond_event_haz_adj_var_names = "w",
    cond_censoring_haz_adj_var_names = "w",
    treatment_var_name = "a",
    outcome_var_name = "time",
    censoring_var_name = "censoring",
    propensity_score_sl_library = c("SL.glm", "SL.mean", "SL.ranger"),
    cond_event_haz_sl_library = c("SL.glm", "SL.mean", "SL.ranger"),
    cond_censoring_haz_sl_library = c("SL.glm", "SL.mean", "SL.ranger")
  )

  # check that null hypothesis is not rejected
  expect_true(tte_diff_var_results$hypothesis_test_tbl$p_value > 0.05)
})

test_that(
  paste("tte_cmldiffvar runs without issue in non-null heterogeneous",
        "multiplicative effect setting"),
{

  # load required libraries
  library(dplyr)
  library(SuperLearner)

  # set seed for reproducibility
  set.seed(71345)

  # sample a dataset
  sample_tbl <- generate_test_data(
    n_obs = 500, null_marginal = FALSE, hazard_model = FALSE
  )

  # produce differential variance inference
  tte_diff_var_results <- tte_cmldiffvar(
    data_tbl = sample_tbl,
    homogeneity_type = "multiplicative",
    restriction_time = 120,
    homogeneous_effect_range_min = 0.9,
    homogeneous_effect_range_max = 1.1,
    propensity_score_adj_var_names = "w",
    cond_event_haz_adj_var_names = "w",
    cond_censoring_haz_adj_var_names = "w",
    treatment_var_name = "a",
    outcome_var_name = "time",
    censoring_var_name = "censoring",
    propensity_score_sl_library = c("SL.glm"),
    cond_event_haz_sl_library = c("SL.glm"),
    cond_censoring_haz_sl_library = c("SL.glm")
  )

  # check that null hypothesis is not rejected
  expect_true(tte_diff_var_results$hypothesis_test_tbl$p_value < 0.05)
})
