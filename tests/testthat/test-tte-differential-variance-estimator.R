test_that(
  paste("one_step_tte_diff_var_estimator_fun performs as expected in the null",
        "additive setting"),
{

  # load required libraries
  library(dplyr)
  library(SuperLearner)

  # set seed for reproducibility
  set.seed(154884)

  # sample a dataset
  sample_tbl <- generate_test_data(n_obs = 500, null_marginal = TRUE)

  # produce a clean, longitudinal version
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
    propensity_score_library = c("SL.glm"),
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

  # produce the differential variance inference table under a homogeneous
  # additive model
  inference_tbl <- one_step_tte_diff_var_estimator_fun(
    clean_long_tbl = long_sample_tbl,
    propensity_score_adj_var_names = "w",
    cond_event_haz_adj_var_names = "w",
    cond_censoring_haz_adj_var_names = "w",
    treatment_var_name = "a",
    propensity_score_var_name = NULL,
    propensity_score_sl_fit = ps_sl_fit,
    cond_event_haz_sl_fit = cond_event_haz_sl_fit,
    cond_censoring_haz_sl_fit = cond_censoring_haz_sl_fit,
    max_time_cutoff = 50,
    effect_range_vec = c(-3, 5),
    homogeneity_type = "additive"
  )

  # make sure that the homogeneity test is not rejected
  max_p_value <- inference_tbl |>
    pull(p_value) |>
    max()
  expect_true(max_p_value > 0.05)

})

test_that(
  paste("one_step_tte_diff_var_estimator_fun performs as expected in the",
        "non-null additive setting"),
  {

    # load required libraries
    library(dplyr)
    library(SuperLearner)

    # set seed for reproducibility
    set.seed(2356239)

    # sample a dataset
    sample_tbl <- generate_test_data(n_obs = 500, null_marginal = FALSE)

    # produce a clean, longitudinal version
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
      propensity_score_library = c("SL.glm"),
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

    # produce the differential variance inference table under a homogeneous
    # additive model
    inference_tbl <- one_step_tte_diff_var_estimator_fun(
      clean_long_tbl = long_sample_tbl,
      propensity_score_adj_var_names = "w",
      cond_event_haz_adj_var_names = "w",
      cond_censoring_haz_adj_var_names = "w",
      treatment_var_name = "a",
      propensity_score_var_name = NULL,
      propensity_score_sl_fit = ps_sl_fit,
      cond_event_haz_sl_fit = cond_event_haz_sl_fit,
      cond_censoring_haz_sl_fit = cond_censoring_haz_sl_fit,
      max_time_cutoff = 40,
      effect_range_vec = c(-3, 5),
      homogeneity_type = "additive"
    )

    # make sure that the homogeneity test is rejected
    max_p_value <- inference_tbl |>
      pull(p_value) |>
      max()
    expect_true(max_p_value < 0.05)

  })


test_that(
  paste("one_step_tte_diff_var_estimator_fun performs as expected in the null",
        "multiplicative setting"),
  {

    # load required libraries
    library(dplyr)
    library(SuperLearner)

    # set seed for reproducibility
    set.seed(45642)

    # sample a dataset
    sample_tbl <- generate_test_data(
      n_obs = 500, null_marginal = TRUE, hazard_model = FALSE
    )

    # produce a clean, longitudinal version
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
      propensity_score_library = c("SL.glm"),
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

    # produce the differential variance inference table under a homogeneous
    # additive model
    inference_tbl <- one_step_tte_diff_var_estimator_fun(
      clean_long_tbl = long_sample_tbl,
      propensity_score_adj_var_names = "w",
      cond_event_haz_adj_var_names = "w",
      cond_censoring_haz_adj_var_names = "w",
      treatment_var_name = "a",
      propensity_score_var_name = NULL,
      propensity_score_sl_fit = ps_sl_fit,
      cond_event_haz_sl_fit = cond_event_haz_sl_fit,
      cond_censoring_haz_sl_fit = cond_censoring_haz_sl_fit,
      max_time_cutoff = 50,
      effect_range_vec = c(0.95, 1.2),
      homogeneity_type = "multiplicative"
    )

    # make sure that the homogeneity test is not rejected
    max_p_value <- inference_tbl |>
      pull(p_value) |>
      max()
    expect_true(max_p_value > 0.05)

})

test_that(
  paste("one_step_tte_diff_var_estimator_fun performs as expected in the",
        "non-null multiplicative setting"),
  {

    # load required libraries
    library(dplyr)
    library(SuperLearner)

    # set seed for reproducibility
    set.seed(182593)

    # sample a dataset
    sample_tbl <- generate_test_data(
      n_obs = 500, null_marginal = FALSE, hazard_model = FALSE
    )

    # produce a clean, longitudinal version
    long_sample_tbl <- sample_tbl |>
      melt_tte_data_fun(
        baseline_var_names = "w",
        treatment_var_name = "a",
        outcome_var_name = "time",
        censoring_var_name = "censoring",
        time_cutoff = 100
      )

    # fit nuisance parameter estimators
    ps_sl_fit <- estimate_propensity_score_fun(
      sample_tbl,
      adj_set_var_names = "w",
      treatment_var_name = "a",
      propensity_score_library = c("SL.glm"),
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

    # produce the differential variance inference table under a homogeneous
    # additive model
    inference_tbl <- one_step_tte_diff_var_estimator_fun(
      clean_long_tbl = long_sample_tbl,
      propensity_score_adj_var_names = "w",
      cond_event_haz_adj_var_names = "w",
      cond_censoring_haz_adj_var_names = "w",
      treatment_var_name = "a",
      propensity_score_var_name = NULL,
      propensity_score_sl_fit = ps_sl_fit,
      cond_event_haz_sl_fit = cond_event_haz_sl_fit,
      cond_censoring_haz_sl_fit = cond_censoring_haz_sl_fit,
      max_time_cutoff = 100,
      effect_range_vec = c(0.9, 1.1),
      homogeneity_type = "multiplicative"
    )

    # make sure that the homogeneity test is rejected
    max_p_value <- inference_tbl |>
      pull(p_value) |>
      max()
    expect_true(max_p_value < 0.05)

})
