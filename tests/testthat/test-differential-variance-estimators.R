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
    num_nuisance_sl_folds = 5
  )
  cond_exp_outcome_fit <- estimate_cond_exp_outcome_fun(
    sample_tbl,
    confounder_var_names = "confounder",
    treatment_var_name = "treatment",
    outcome_var_name = "outcome",
    cond_exp_outcome_library = c("SL.glm", "SL.mean"),
    num_nuisance_sl_folds = 5
  )
  cond_exp_sq_outcome_fit <- estimate_cond_exp_sq_outcome_fun(
    sample_tbl,
    confounder_var_names = "confounder",
    treatment_var_name = "treatment",
    outcome_var_name = "outcome",
    cond_exp_sq_outcome_library = c("SL.glm", "SL.earth"),
    num_nuisance_sl_folds = 5
  )

  # construct the counterfactual dataset under the treatment group
  sample_treatment_tbl <- generate_counterfactural_tbl_fun(
    sample_tbl,
    treatment_group = 1,
    confounder_var_names = "confounder",
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

test_that("group-specific one-step mean estimator is consistent", {

  # load required libraries
  library(dplyr)
  library(SuperLearner)

  # set seed for reproducibility
  set.seed(723424)

  # calculate estimands
  mean_treatment <- mean(toy_population_tbl$potential_outcome_treatment)

  # compute bias
  num_iters <- 100
  os_est_vec <- sapply(
    seq_len(num_iters),
    function(iter) {

      # grab a sample of the population
      sample_tbl <- slice_sample(toy_population_tbl, n = 1000)

      # fit nuisance parameter estimators
      ps_sl_fit <- estimate_propensity_score_fun(
        sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        propensity_score_library = c("SL.glm", "SL.mean"),
        num_nuisance_sl_folds = 5
      )
      cond_exp_outcome_fit <- estimate_cond_exp_outcome_fun(
        sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        cond_exp_outcome_library = c("SL.glm", "SL.mean"),
        num_nuisance_sl_folds = 5
      )
      cond_exp_sq_outcome_fit <- estimate_cond_exp_sq_outcome_fun(
        sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        cond_exp_sq_outcome_library = c("SL.glm", "SL.earth"),
        num_nuisance_sl_folds = 5
      )

      # create the counterfactual dataset for the treatment group
      sample_treatment_tbl <- generate_counterfactural_tbl_fun(
        sample_tbl,
        treatment_group = 1,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        propensity_score_var_name = NULL,
        ps_sl_fit,
        cond_exp_outcome_fit,
        cond_exp_sq_outcome_fit
      )

      # one-step estimate
      one_step_mean_estimator_fun(
        treatment_group = 1,
        sample_treatment_tbl$treatment,
        sample_treatment_tbl$outcome,
        sample_treatment_tbl$pred_propensity_score,
        sample_treatment_tbl$pred_cond_exp_outcome
      )
    }
  )

  # expect negligible bias in large sample sizes
  # that is < 5% relative bias (0.05 * estimand = 0.15)
  expect_lt(abs(mean(os_est_vec - mean_treatment)), 0.15)

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
  abs_estimand <- var_treatment - var_control
  rel_estimand <- var_treatment / var_control

  # compute bias
  num_iters <- 100
  abs_one_step_diff_var_ests <- rep(NA, num_iters)
  abs_tml_diff_var_ests <- rep(NA, num_iters)
  rel_one_step_diff_var_ests <- rep(NA, num_iters)
  rel_tml_diff_var_ests <- rep(NA, num_iters)
  for (iter in seq_len(num_iters)) {

    # grab a sample of the population
    sample_tbl <- slice_sample(toy_population_tbl, n = 1000)

    # fit nuisance parameter estimators
    ps_sl_fit <- estimate_propensity_score_fun(
      sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      propensity_score_library = c("SL.glm", "SL.mean"),
      num_nuisance_sl_folds = 5
    )
    cond_exp_outcome_fit <- estimate_cond_exp_outcome_fun(
      sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      outcome_var_name = "outcome",
      cond_exp_outcome_library = c("SL.glm", "SL.mean"),
      num_nuisance_sl_folds = 5
    )
    cond_exp_sq_outcome_fit <- estimate_cond_exp_sq_outcome_fun(
      sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      outcome_var_name = "outcome",
      cond_exp_sq_outcome_library = c("SL.mean", "SL.earth", "SL.glm"),
      num_nuisance_sl_folds = 5
    )

    # absolute one-step estimate
    abs_one_step_diff_var_ests[iter] <- one_step_diff_var_estimator_fun(
      sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      propensity_score_var_name = NULL,
      outcome_var_name = "outcome",
      ps_sl_fit,
      cond_exp_outcome_fit,
      cond_exp_sq_outcome_fit,
      estimand_type = "absolute"
    )$estimate

    # absolute TML estimate
    abs_tml_diff_var_ests[iter] <- tml_diff_var_estimator_fun(
      sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      propensity_score_var_name = NULL,
      outcome_var_name = "outcome",
      ps_sl_fit,
      cond_exp_outcome_fit,
      cond_exp_sq_outcome_fit,
      estimand_type = "absolute"
    )$estimate

    # relative one-step estimate
    rel_one_step_diff_var_ests[iter] <- one_step_diff_var_estimator_fun(
      sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      propensity_score_var_name = NULL,
      outcome_var_name = "outcome",
      ps_sl_fit,
      cond_exp_outcome_fit,
      cond_exp_sq_outcome_fit,
      estimand_type = "relative"
    )$estimate

    # relative TML estimate
    rel_tml_diff_var_ests[iter] <- tml_diff_var_estimator_fun(
      sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      propensity_score_var_name = NULL,
      outcome_var_name = "outcome",
      ps_sl_fit,
      cond_exp_outcome_fit,
      cond_exp_sq_outcome_fit,
      estimand_type = "relative"
    )$estimate

  }

  # expect negligible bias in large sample sizes for absolute estimate
  # that is < 5% relative bias (0.05 * estimand = 0.4)
  expect_lt(abs(mean(abs_one_step_diff_var_ests - abs_estimand)), 0.4)
  expect_lt(abs(mean(abs_tml_diff_var_ests - abs_estimand)), 0.4)

  # expect negligible bias in large sample sizes for relative estimate
  # that is < 5% relative bias (0.05 * estimand = 0.4)
  expect_lt(abs(mean(rel_one_step_diff_var_ests - rel_estimand)), 0.25)
  expect_lt(abs(mean(rel_tml_diff_var_ests - rel_estimand)), 0.25)

})

test_that("TMLE approximately solves the EIF", {

  # load required libraries
  library(dplyr)
  library(SuperLearner)
  library(earth)

  set.seed(713452)

  # grab a sample of the population
  sample_tbl <- slice_sample(toy_population_tbl, n = 1000)

  # fit nuisance parameter estimators
  ps_sl_fit <- estimate_propensity_score_fun(
    sample_tbl,
    confounder_var_names = "confounder",
    treatment_var_name = "treatment",
    propensity_score_library = c("SL.glm", "SL.mean"),
    num_nuisance_sl_folds = 5
  )
  cond_exp_outcome_fit <- estimate_cond_exp_outcome_fun(
    sample_tbl,
    confounder_var_names = "confounder",
    treatment_var_name = "treatment",
    outcome_var_name = "outcome",
    cond_exp_outcome_library = c("SL.glm", "SL.mean"),
    num_nuisance_sl_folds = 5
  )
  cond_exp_sq_outcome_fit <- estimate_cond_exp_sq_outcome_fun(
    sample_tbl,
    confounder_var_names = "confounder",
    treatment_var_name = "treatment",
    outcome_var_name = "outcome",
    cond_exp_sq_outcome_library = c("SL.mean", "SL.earth", "SL.glm"),
    num_nuisance_sl_folds = 5
  )

  # absolute TML estimator's eif
  abs_tml_diff_var_eif <- tml_diff_var_estimator_fun(
    sample_tbl,
    confounder_var_names = "confounder",
    treatment_var_name = "treatment",
    propensity_score_var_name = NULL,
    outcome_var_name = "outcome",
    ps_sl_fit,
    cond_exp_outcome_fit,
    cond_exp_sq_outcome_fit,
    estimand_type = "absolute"
  )$eif

  # relative TML estimator's eif
  rel_tml_diff_var_eif <- tml_diff_var_estimator_fun(
    sample_tbl,
    confounder_var_names = "confounder",
    treatment_var_name = "treatment",
    propensity_score_var_name = NULL,
    outcome_var_name = "outcome",
    ps_sl_fit,
    cond_exp_outcome_fit,
    cond_exp_sq_outcome_fit,
    estimand_type = "relative"
  )$eif

  expect_lt(mean(abs_tml_diff_var_eif), .Machine$double.eps * 10e3)
  expect_lt(mean(rel_tml_diff_var_eif), .Machine$double.eps * 10e3)

})

test_that("cross-fitted differential variance estimators are consistent", {

  # load required libraries
  library(dplyr)
  library(SuperLearner)
  library(earth)
  library(origami)

  set.seed(2484322)

  # calculate estimand
  var_treatment <- var(toy_population_tbl$potential_outcome_treatment)
  var_control <- var(toy_population_tbl$potential_outcome_control)
  abs_estimand <- var_treatment - var_control

  # compute bias
  num_iters <- 100
  abs_one_step_diff_var_ests <- rep(NA, num_iters)
  abs_tmle_diff_var_ests <- rep(NA, num_iters)
  for (iter in seq_len(num_iters)) {

    # grab a sample of the population
    sample_tbl <- slice_sample(toy_population_tbl, n = 1000)

    # split the sample data into folds
    folds <- make_folds(sample_tbl, fold_fun = folds_vfold, V = 5L)

    # cross-validate the diff var estimators
    cf_os_diff_var_ests <- cross_validate(
      cv_fun = cf_diff_var_estimator_fun,
      folds = folds,
      clean_tbl = sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      propensity_score_var_name = NULL,
      outcome_var_name = "outcome",
      propensity_score_library = c("SL.mean", "SL.glm"),
      cond_exp_outcome_library = c("SL.mean", "SL.glm"),
      cond_exp_sq_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
      num_nuisance_sl_folds = 5,
      estimator_type = "one-step",
      estimand_type = "absolute"
    )
    cf_tmle_diff_var_ests <- cross_validate(
      cv_fun = cf_diff_var_estimator_fun,
      folds = folds,
      clean_tbl = sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      propensity_score_var_name = NULL,
      outcome_var_name = "outcome",
      propensity_score_library = c("SL.mean", "SL.glm"),
      cond_exp_outcome_library = c("SL.mean", "SL.glm"),
      cond_exp_sq_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
      num_nuisance_sl_folds = 5,
      estimator_type = "tmle",
      estimand_type = "absolute"
    )

    # absolute one-step estimate
    abs_one_step_diff_var_ests[iter] <- mean(cf_os_diff_var_ests$estimates)
    abs_tmle_diff_var_ests[iter] <- mean(cf_tmle_diff_var_ests$estimates)

  }

  # expect negligible bias in large sample sizes for absolute estimate
  # that is < 5% relative bias (0.05 * estimand = 0.4)
  expect_lt(abs(mean(abs_one_step_diff_var_ests - abs_estimand)), 0.4)
  expect_lt(abs(mean(abs_tmle_diff_var_ests - abs_estimand)), 0.4)

})

test_that("cross-fitted TMLE approximately solves the EIF", {

  # load required libraries
  library(dplyr)
  library(SuperLearner)
  library(earth)
  library(origami)

  # grab a sample of the population
  sample_tbl <- slice_sample(toy_population_tbl, n = 1000) |>
    mutate(sq_outcome = outcome^2)

  # split the sample data into folds
  folds <- make_folds(sample_tbl, fold_fun = folds_vfold, V = 5L)

  # compute the eif of the cross-fitted procedure
  cf_tmle_diff_var_ests <- cross_validate(
    cv_fun = cf_diff_var_estimator_fun,
    folds = folds,
    clean_tbl = sample_tbl,
    confounder_var_names = "confounder",
    treatment_var_name = "treatment",
    propensity_score_var_name = NULL,
    outcome_var_name = "outcome",
    propensity_score_library = c("SL.mean", "SL.glm"),
    cond_exp_outcome_library = c("SL.mean", "SL.glm"),
    cond_exp_sq_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
    num_nuisance_sl_folds = 5,
    estimator_type = "tmle",
    estimand_type = "absolute"
  )

  expect_lt(mean(cf_tmle_diff_var_ests$eif), .Machine$double.eps * 10e3)

})

test_that(
  "cross-fitted differential variance estimators are asymptotically linear",
{

  # load required libraries
  library(dplyr)
  library(SuperLearner)
  library(earth)
  library(origami)

  set.seed(18349321)

  # calculate estimand
  var_treatment <- var(toy_population_tbl$potential_outcome_treatment)
  var_control <- var(toy_population_tbl$potential_outcome_control)
  abs_estimand <- var_treatment - var_control

  # approximate coverage
  num_iters <- 100
  abs_one_step_covered <- rep(NA, num_iters)
  abs_tmle_covered <- rep(NA, num_iters)
  for (iter in seq_len(num_iters)) {

    # grab a sample of the population
    sample_tbl <- slice_sample(toy_population_tbl, n = 1000)

    # split the sample data into folds
    folds <- make_folds(sample_tbl, fold_fun = folds_vfold, V = 5L)

    # cross-validate the diff var estimators
    cf_os_diff_var_ests <- cross_validate(
      cv_fun = cf_diff_var_estimator_fun,
      folds = folds,
      clean_tbl = sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      propensity_score_var_name = NULL,
      outcome_var_name = "outcome",
      propensity_score_library = c("SL.mean", "SL.glm"),
      cond_exp_outcome_library = c("SL.mean", "SL.glm"),
      cond_exp_sq_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
      num_nuisance_sl_folds = 5,
      estimator_type = "one-step",
      estimand_type = "absolute"
    )
    cf_tmle_diff_var_ests <- cross_validate(
      cv_fun = cf_diff_var_estimator_fun,
      folds = folds,
      clean_tbl = sample_tbl,
      confounder_var_names = "confounder",
      treatment_var_name = "treatment",
      propensity_score_var_name = NULL,
      outcome_var_name = "outcome",
      propensity_score_library = c("SL.mean", "SL.glm"),
      cond_exp_outcome_library = c("SL.mean", "SL.glm"),
      cond_exp_sq_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
      num_nuisance_sl_folds = 5,
      estimator_type = "tmle",
      estimand_type = "absolute"
    )

    # compute CIs for one-step
    cf_os_estimate <- mean(cf_os_diff_var_ests$estimates)
    cf_os_se <- sqrt(var(cf_os_diff_var_ests$eif) / nrow(sample_tbl))
    cf_os_ci_lower <- cf_os_estimate - 1.96 * cf_os_se
    cf_os_ci_upper <- cf_os_estimate + 1.96 * cf_os_se

    # compute CIs for TMLE
    cf_tmle_estimate <- mean(cf_tmle_diff_var_ests$estimates)
    cf_tmle_se <- sqrt(var(cf_tmle_diff_var_ests$eif) / nrow(sample_tbl))
    cf_tmle_ci_lower <- cf_tmle_estimate - 1.96 * cf_tmle_se
    cf_tmle_ci_upper <- cf_tmle_estimate + 1.96 * cf_tmle_se


    # absolute one-step estimate
    abs_one_step_covered[iter] <- ifelse(
      cf_os_ci_lower < abs_estimand && abs_estimand < cf_os_ci_upper, 1, 0
    )
    abs_tmle_covered[iter] <- ifelse(
      cf_tmle_ci_lower < abs_estimand && abs_estimand < cf_tmle_ci_upper, 1, 0
    )

  }

  # expect approximate 95% coverage, with some wiggle room
  expect_lt(abs(mean(abs_one_step_covered) - 0.95), 0.025)
  expect_lt(abs(mean(abs_tmle_covered) - 0.95), 0.025)

})


test_that(
  paste0(
    "cross-fitted differential variance estimators are asymptotically ",
    "linear when propensity score is known"
  ),
{

    # load required libraries
    library(dplyr)
    library(SuperLearner)
    library(earth)
    library(origami)

    set.seed(6162342)

    # calculate estimand
    var_treatment <- var(toy_population_tbl$potential_outcome_treatment)
    var_control <- var(toy_population_tbl$potential_outcome_control)
    abs_estimand <- var_treatment - var_control

    # approximate coverage
    num_iters <- 100
    abs_one_step_covered <- rep(NA, num_iters)
    abs_tmle_covered <- rep(NA, num_iters)
    for (iter in seq_len(num_iters)) {

      # grab a sample of the population
      sample_tbl <- slice_sample(toy_population_tbl, n = 1000)

      # split the sample data into folds
      folds <- make_folds(sample_tbl, fold_fun = folds_vfold, V = 5L)

      # cross-validate the diff var estimators
      cf_os_diff_var_ests <- cross_validate(
        cv_fun = cf_diff_var_estimator_fun,
        folds = folds,
        clean_tbl = sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        propensity_score_var_name = "propensity_score",
        outcome_var_name = "outcome",
        cond_exp_outcome_library = c("SL.mean", "SL.glm"),
        cond_exp_sq_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        num_nuisance_sl_folds = 5,
        estimator_type = "one-step",
        estimand_type = "absolute"
      )
      cf_tmle_diff_var_ests <- cross_validate(
        cv_fun = cf_diff_var_estimator_fun,
        folds = folds,
        clean_tbl = sample_tbl,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        propensity_score_var_name = "propensity_score",
        outcome_var_name = "outcome",
        cond_exp_outcome_library = c("SL.mean", "SL.glm"),
        cond_exp_sq_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        num_nuisance_sl_folds = 5,
        estimator_type = "tmle",
        estimand_type = "absolute"
      )

      # compute CIs for one-step
      cf_os_estimate <- mean(cf_os_diff_var_ests$estimates)
      cf_os_se <- sqrt(var(cf_os_diff_var_ests$eif) / nrow(sample_tbl))
      cf_os_ci_lower <- cf_os_estimate - 1.96 * cf_os_se
      cf_os_ci_upper <- cf_os_estimate + 1.96 * cf_os_se

      # compute CIs for TMLE
      cf_tmle_estimate <- mean(cf_tmle_diff_var_ests$estimates)
      cf_tmle_se <- sqrt(var(cf_tmle_diff_var_ests$eif) / nrow(sample_tbl))
      cf_tmle_ci_lower <- cf_tmle_estimate - 1.96 * cf_tmle_se
      cf_tmle_ci_upper <- cf_tmle_estimate + 1.96 * cf_tmle_se


      # absolute one-step estimate
      abs_one_step_covered[iter] <- ifelse(
        cf_os_ci_lower < abs_estimand && abs_estimand < cf_os_ci_upper, 1, 0
      )
      abs_tmle_covered[iter] <- ifelse(
        cf_tmle_ci_lower < abs_estimand && abs_estimand < cf_tmle_ci_upper, 1, 0
      )

    }

    # expect approximate 95% coverage, with some wiggle room
    expect_lt(abs(mean(abs_one_step_covered) - 0.95), 0.025)
    expect_lt(abs(mean(abs_tmle_covered) - 0.95), 0.025)

  })
