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
    confounder_var_names = "confounder",
    treatment_var_name = "treatment",
    propensity_score_library = c("SL.glm", "SL.mean"),
    num_folds = 5
  )
  cond_exp_outcome_fit <- estimate_cond_exp_outcome_fun(
    sample_tbl,
    confounder_var_names = "confounder",
    treatment_var_name = "treatment",
    outcome_var_name = "outcome",
    cond_exp_outcome_library = c("SL.glm", "SL.mean"),
    num_folds = 5
  )
  cond_exp_sq_outcome_fit <- estimate_cond_exp_sq_outcome_fun(
    sample_tbl,
    confounder_var_names = "confounder",
    treatment_var_name = "treatment",
    outcome_var_name = "outcome",
    cond_exp_sq_outcome_library = c("SL.glm", "SL.earth"),
    num_folds = 5
  )

  # construct the counterfactual dataset under the treatment group
  sample_treatment_tbl <- generate_counterfactural_tbl_fun(
    sample_tbl,
    treatment_group = 1,
    confounder_var_names = "confounder",
    treatment_var_name = "treatment",
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

test_that("differential variance estimators are consistent", {

  # load required libraries
  library(dplyr)
  library(SuperLearner)
  library(earth)

  set.seed(713452)

  # calculate estimand
  var_treatment <- var(toy_population_tbl$potential_outcome_treatment)
  var_control <- var(toy_population_tbl$potential_outcome_control)
  estimand <- var_treatment - var_control

  # compute bias
  num_iters <- 100
  one_step_diff_var_ests <- rep(NA, num_iters)
  tml_diff_var_ests <- rep(NA, num_iters)
  for (iter in seq_len(num_iters)) {

    # grab a sample of the population
    sample_tbl <- slice_sample(toy_population_tbl, n = 1000) |>
      mutate(sq_outcome = outcome^2)

    # fit nuisance parameter estimators
    ps_sl_fit <- estimate_propensity_score_fun(
      sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      propensity_score_library = c("SL.glm", "SL.mean"),
      num_folds = 5
    )
    cond_exp_outcome_fit <- estimate_cond_exp_outcome_fun(
      sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      outcome_var_name = "outcome",
      cond_exp_outcome_library = c("SL.glm", "SL.mean"),
      num_folds = 5
    )
    cond_exp_sq_outcome_fit <- estimate_cond_exp_sq_outcome_fun(
      sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      outcome_var_name = "outcome",
      cond_exp_sq_outcome_library = c("SL.glm", "SL.earth"),
      num_folds = 5
    )

    # one-step estimate
    one_step_diff_var_ests[iter] <- one_step_diff_var_estimator_fun(
      sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      outcome_var_name = "outcome",
      ps_sl_fit,
      cond_exp_outcome_fit,
      cond_exp_sq_outcome_fit
    )

    # TML estimate
    tml_diff_var_ests[iter] <- tml_diff_var_estimator_fun(
      sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      outcome_var_name = "outcome",
      ps_sl_fit,
      cond_exp_outcome_fit,
      cond_exp_sq_outcome_fit
    )

  }

  # expect negligible bias in large sample sizes
  # that is < 5% relative bias (0.05 * estimand = 0.4)
  expect_lt(abs(mean(one_step_diff_var_ests - estimand)), 0.4)
  expect_lt(abs(mean(tml_diff_var_ests - estimand)), 0.4)

})
