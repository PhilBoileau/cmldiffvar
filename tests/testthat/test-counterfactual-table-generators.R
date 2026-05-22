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

  # predictions from nuisance parameters match manula calculations
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
