test_that("estimators of treatment group-specific variance are consistent", {

  # calculate estimand
  var_treatment <- var(toy_population_tbl$potential_outcome_treatment)
  var_control <- var(toy_population_tbl$potential_outcome_control)

  # compute bias
  num_iters <- 100
  one_step_var_treatment_est_vec <- rep(NA, num_iters)
  one_step_var_control_est_vec <- rep(NA, num_iters)
  tml_var_treatment_est_vec <- rep(NA, num_iters)
  tml_var_control_est_vec <- rep(NA, num_iters)
  for (iter in seq_len(num_iters)) {

    # grab a sample of the population
    sample_tbl <- dplyr::slice_sample(toy_population_tbl, n = 1000) |>
      dplyr::mutate(sq_outcome = outcome^2)

    # estimate nuisance parameters
    ps_fit <- glm(
      treatment ~ confounder, family = "binomial", data = sample_tbl
    )
    cond_outcome_fit <- earth::earth(
      outcome ~ treatment + confounder,
      data = sample_tbl
    )
    cond_sq_outcome_fit <- earth::earth(
      sq_outcome ~ treatment + confounder + treatment^2 + confounder^2,
      data = sample_tbl,
      degree = 2
    )

    # predict conditional outcomes under treatment
    sample_treatment_tbl <- sample_tbl |>
      dplyr::mutate(treatment = 1)
    cond_outcome_treatment_est <- predict(
      cond_outcome_fit,
      sample_treatment_tbl
    )
    cond_sq_outcome_treatment_est <- predict(
      cond_sq_outcome_fit, sample_treatment_tbl
    )

    # predict conditional outcomes under control
    sample_control_tbl <- sample_tbl |>
      dplyr::mutate(treatment = 0)
    cond_outcome_control_est <- predict(cond_outcome_fit, sample_control_tbl)
    cond_sq_outcome_control_est <- predict(
      cond_sq_outcome_fit, sample_control_tbl
    )

    # predict propensity score
    ps_est <- predict(ps_fit, type = "response")

    # predicted mean under treatment and control (one-step estimator)
    one_step_mean_treatment_est <- mean(
      (sample_tbl$treatment == 1) *
        (sample_tbl$outcome - cond_outcome_treatment_est) /
        ps_est + cond_outcome_treatment_est
    )
    one_step_mean_control_est <- mean(
      (sample_tbl$treatment == 0) *
        (sample_tbl$outcome - cond_outcome_control_est) /
        ps_est + cond_outcome_control_est
    )

    # calculate variance of treatment one-step estimate
    one_step_var_treatment_est_vec[iter] <- one_step_var_estimator_fun(
      treatment_group = 1,
      sample_tbl$treatment,
      sample_treatment_tbl$outcome,
      ps_est,
      cond_outcome_treatment_est,
      cond_sq_outcome_treatment_est,
      one_step_mean_treatment_est
    )$estimate

    # calculate variance of treatment TML estimate
    tml_var_treatment_est_vec[iter] <- tml_var_estimator_fun(
      treatment_group = 1,
      treatment_vec = sample_tbl$treatment,
      outcome_vec = sample_treatment_tbl$outcome,
      ps_est_vec = ps_est,
      cond_exp_outcome_est_vec = cond_outcome_treatment_est,
      cond_exp_sq_outcome_est_vec = cond_sq_outcome_treatment_est
    )$estimate

    # calculate variance of control one-step estimate
    one_step_var_control_est_vec[iter] <- one_step_var_estimator_fun(
      treatment_group = 0,
      sample_tbl$treatment,
      sample_control_tbl$outcome,
      ps_est,
      cond_outcome_control_est,
      cond_sq_outcome_control_est,
      one_step_mean_control_est
    )$estimate

    # calculate variance of control TML estimate
    tml_var_control_est_vec[iter] <- tml_var_estimator_fun(
      treatment_group = 0,
      sample_tbl$treatment,
      sample_control_tbl$outcome,
      ps_est,
      cond_outcome_control_est,
      cond_sq_outcome_control_est
    )$estimate

  }

  # check one-step estimator empirical bias
  os_empirical_bias_var_treatment <- mean(
    one_step_var_treatment_est_vec - var_treatment
  )
  expect_lt(abs(os_empirical_bias_var_treatment), 1)
  os_empirical_bias_var_control <- mean(
    one_step_var_control_est_vec - var_control
  )
  expect_lt(abs(os_empirical_bias_var_control), 1)

  # check TML estimator empirical bias
  tml_empirical_bias_var_treatment <- mean(
    tml_var_treatment_est_vec - var_treatment
  )
  expect_lt(abs(tml_empirical_bias_var_treatment), 1)
  tml_empirical_bias_var_control <- mean(
    tml_var_control_est_vec - var_control
  )
  expect_lt(abs(tml_empirical_bias_var_control), 1)

})


test_that("TMLE of treatment group-specific variance solves the EIF", {

  # grab a sample of the population
  sample_tbl <- dplyr::slice_sample(toy_population_tbl, n = 1000) |>
    dplyr::mutate(sq_outcome = outcome^2)

  # estimate nuisance parameters
  ps_fit <- glm(
    treatment ~ confounder, family = "binomial", data = sample_tbl
  )
  cond_outcome_fit <- earth::earth(
    outcome ~ treatment + confounder,
    data = sample_tbl
  )
  cond_sq_outcome_fit <- earth::earth(
    sq_outcome ~ treatment + confounder + treatment^2 + confounder^2,
    data = sample_tbl,
    degree = 2
  )

  # predict conditional outcomes under treatment
  sample_treatment_tbl <- sample_tbl |>
    dplyr::mutate(treatment = 1)
  cond_outcome_treatment_est <- predict(
    cond_outcome_fit,
    sample_treatment_tbl
  )
  cond_sq_outcome_treatment_est <- predict(
    cond_sq_outcome_fit, sample_treatment_tbl
  )

  # predict conditional outcomes under control
  sample_control_tbl <- sample_tbl |>
    dplyr::mutate(treatment = 0)
  cond_outcome_control_est <- predict(cond_outcome_fit, sample_control_tbl)
  cond_sq_outcome_control_est <- predict(
    cond_sq_outcome_fit, sample_control_tbl
  )

  # predict propensity score
  ps_est <- predict(ps_fit, type = "response")

  # calculate EIF of variance of treatment TML estimate
  tml_var_treatment_eif <- tml_var_estimator_fun(
    treatment_group = 1,
    treatment_vec = sample_tbl$treatment,
    outcome_vec = sample_treatment_tbl$outcome,
    ps_est_vec = ps_est,
    cond_exp_outcome_est_vec = cond_outcome_treatment_est,
    cond_exp_sq_outcome_est_vec = cond_sq_outcome_treatment_est
  )$eif

  # calculate EIF of variance of control TML estimate
  tml_var_control_eif <- tml_var_estimator_fun(
    treatment_group = 0,
    sample_tbl$treatment,
    sample_control_tbl$outcome,
    ps_est,
    cond_outcome_control_est,
    cond_sq_outcome_control_est
  )$eif

  # make sure mean EIF is approximately equal to zero
  # that is, that the TMLE tilting procedure solves the estimating equation
  # defined by the empirical bias term in the von Mises expansion
  expect_lt(mean(tml_var_treatment_eif), 1e-6)
  expect_lt(mean(tml_var_control_eif), 1e-6)

})

test_that("estimators are consistent when regressions are estimated by means", {

  library(dplyr)
  library(SuperLearner)

  set.seed(621341)

  # a DGP function for assessing multiple robusteness
  # the true variance of Y1 and Y0 are 7.72 and 3.11, respectively
  var_y1 <- 7.72
  var_y0 <- 3.11
  multiple_robustness_dgp_fun <- function(n) {

    # Define confounder variables
    W1 <- stats::rbinom(n = n, size = 1, p = 0.3)
    W2 <- stats::rnorm(n = n, mean = 0, sd = 1)

    # Compute propensity score
    pA <- stats::plogis(q = (1 + W1 + W2)/4)

    # Define treatment variable
    A <- stats::rbinom(n = n, size = 1, p = pA)

    # Compute mean of potential outcomes
    # We have set our DGP such that: E[Y(a)] = 1 + a + W1 + W2 + a*W2 + W1*W2
    meanYA1 <- 1 + 1 + W1 + W2 + 1*W2 + W1*W2
    meanYA0 <- 1 + 0 + W1 + W2 + 0*W2 + W1*W2

    # Define potential outcomes
    YA1 <- stats::rnorm(n = n, mean = meanYA1, sd = sqrt(2))
    YA0 <- stats::rnorm(n = n, mean = meanYA0, sd = 1)

    # Compute outcome variable
    Y <- A * YA1 + (1 - A) * YA0

    # return a tibble of the observed data
    tibble(
      "Y" = Y,
      "A" = A,
      "W1" = W1,
      "W2" = W2
    )
  }

  # generate a dataset
  sample_tbl <- multiple_robustness_dgp_fun(n = 1000)

  # estimate the propensity score with a well-specified estimator
  ps_fit <- glm(
    A ~ W1 + W2, family = "binomial", data = sample_tbl
  )

  # estimate the outcome regressions by their group means
  cond_outcome_fit <- glm(Y ~ 1, data = sample_tbl)
  cond_sq_outcome_fit <- glm(Y^2 ~ 1, data = sample_tbl)

  # predict conditional outcomes under treatment
  sample_treatment_tbl <- sample_tbl |> mutate(A = 1)
  cond_outcome_treatment_est <- predict(
    cond_outcome_fit,
    sample_treatment_tbl
  )
  cond_sq_outcome_treatment_est <- predict(
    cond_sq_outcome_fit, sample_treatment_tbl
  )

  # predict conditional outcomes under control
  sample_control_tbl <- sample_tbl |> mutate(A = 0)
  cond_outcome_control_est <- predict(cond_outcome_fit, sample_control_tbl)
  cond_sq_outcome_control_est <- predict(
    cond_sq_outcome_fit, sample_control_tbl
  )

  # predict propensity score
  ps_est <- predict(ps_fit, type = "response")

  # predicted mean under treatment and control (one-step estimator)
  one_step_mean_treatment_est <- mean(
    (sample_tbl$A == 1) *
      (sample_tbl$Y - cond_outcome_treatment_est) /
      ps_est + cond_outcome_treatment_est
  )
  one_step_mean_control_est <- mean(
    (sample_tbl$A == 0) *
      (sample_tbl$Y - cond_outcome_control_est) /
      ps_est + cond_outcome_control_est
  )

  # calculate variance of treatment one-step estimate
  one_step_var_treatment_est <- one_step_var_estimator_fun(
    treatment_group = 1,
    treatment_vec = sample_tbl$A,
    outcome_vec = sample_treatment_tbl$Y,
    ps_est_vec = ps_est,
    cond_exp_outcome_est_vec = cond_outcome_treatment_est,
    cond_exp_sq_outcome_est_vec = cond_sq_outcome_treatment_est,
    mean_est = one_step_mean_treatment_est
  )$estimate

  # calculate variance of treatment TML estimate
  tml_var_treatment_est <- tml_var_estimator_fun(
    treatment_group = 1,
    treatment_vec = sample_tbl$A,
    outcome_vec = sample_treatment_tbl$Y,
    ps_est_vec = ps_est,
    cond_exp_outcome_est_vec = cond_outcome_treatment_est,
    cond_exp_sq_outcome_est_vec = cond_sq_outcome_treatment_est
  )$estimate

  # calculate variance of control one-step estimate
  one_step_var_control_est <- one_step_var_estimator_fun(
    treatment_group = 0,
    treatment_vec = sample_tbl$A,
    outcome_vec = sample_treatment_tbl$Y,
    ps_est_vec = ps_est,
    cond_exp_outcome_est_vec = cond_outcome_control_est,
    cond_exp_sq_outcome_est_vec = cond_sq_outcome_control_est,
    mean_est = one_step_mean_control_est
  )$estimate

  # calculate variance of control TML estimate
  tml_var_control_est <- tml_var_estimator_fun(
    treatment_group = 0,
    treatment_vec = sample_tbl$A,
    outcome_vec = sample_treatment_tbl$Y,
    ps_est_vec = ps_est,
    cond_exp_outcome_est_vec = cond_outcome_control_est,
    cond_exp_sq_outcome_est_vec = cond_sq_outcome_control_est
  )$estimate

  # make sure that estimates are close to the true parameter values
  expect_equal(
    one_step_var_treatment_est, var_y1, tolerance = 0.2
  )
  expect_equal(
    tml_var_treatment_est, var_y1, tolerance = 0.2
  )
  expect_equal(
    one_step_var_control_est, var_y0, tolerance = 0.2
  )
  expect_equal(
    tml_var_control_est, var_y0, tolerance = 0.2
  )

})

test_that("unadjusted estimators of group-specific are consistent", {

  library(dplyr)

  # generate randomized study data
  set.seed(7234135)
  n_pop <- 100000
  propensity_score <- 0.5
  propensity_score_vec <- rep(propensity_score, n_pop)
  treatment_vec <- rbinom(n_pop, 1, propensity_score)
  outcome_treatment_vec <- rnorm(n = n_pop, mean = 3, sd = 3)
  outcome_control_vec <- rnorm(n = n_pop, mean = 1, sd = 1)
  outcome_vec <- sapply(
    seq_len(n_pop),
    function(pop_idx) {
      if (treatment_vec[pop_idx] == 1) {
        outcome_treatment_vec[pop_idx]
      } else {
        outcome_control_vec[pop_idx]
      }
    }
  )
  population_tbl <- tibble(
    propensity_score = propensity_score_vec,
    treatment = treatment_vec,
    outcome = outcome_vec
  )

  # estimate bias
  num_iters <- 100
  var_estimates_tbl <- lapply(
    seq_len(num_iters),
    function(iter) {

      # sample from population
      sample_tbl <- slice_sample(population_tbl, n = 500)

      # estimate treatment group variance
      treatment_var_est <- unadjusted_var_estimator_fun(
        treatment_group = 1,
        treatment_vec = sample_tbl$treatment,
        outcome_vec = sample_tbl$outcome,
        ps_est_vec = sample_tbl$propensity_score
      )

      # estimate control group variance
      control_var_est <- unadjusted_var_estimator_fun(
        treatment_group = 0,
        treatment_vec = sample_tbl$treatment,
        outcome_vec = sample_tbl$outcome,
        ps_est_vec = sample_tbl$propensity_score
      )

      # return estimates
      data.frame(
        treatment_var_est = treatment_var_est$estimate,
        control_var_est = control_var_est$estimate,
        treatment_eif_mean = mean(treatment_var_est$eif),
        control_eif_mean = mean(control_var_est$eif)
      )
    }
  ) |>
    bind_rows()

  # assess empirical absolute bias
  abs_bias_treatment <- abs(mean(
    var_estimates_tbl$treatment_var_est - var(outcome_treatment_vec)
  ))
  abs_bias_control <- abs(mean(
    var_estimates_tbl$control_var_est - var(outcome_control_vec)
  ))
  expect_lt(abs_bias_treatment, 1)
  expect_lt(abs_bias_control, 1)

  # ensure that the EIF is solved for each estimator
  expect_true(all(abs(var_estimates_tbl$treatment_eif_mean) < 1e-10))
  expect_true(all(abs(var_estimates_tbl$control_eif_mean) < 1e-10))

})


test_that("One-step TTE variance estimator is consistent", {

  # load required libraries
  library(dplyr)
  library(SuperLearner)

  # set seed for reproducibility
  set.seed(8134613)

  # grab a sample of the population and melt their data
  var_estimate_vec <- sapply(
    seq_len(50),
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

      # estimate the survival probability among the treated at time = 50
      one_step_tte_var_estimator_fun(
        counterfactual_long_tbl = sample_long_treatment_tbl,
        treatment_group = 1,
        treatment_var_name = "a",
        time_cutoff = 40
      )
    }
  )

  # approximate true RMST under control at time=40
  population_tbl <- generate_test_data(n_obs = 10000, null_marginal = FALSE)
  var_at_40 <- population_tbl |>
    mutate(
      trunc_potential_time_1 = if_else(
        potential_time_1 > 40, 40, potential_time_1
      )
    ) |>
    summarise(var_at_40 = var(trunc_potential_time_1)) |>
    pull(var_at_40)

  # make sure the empirical bias is small for the SD
  expect_lte(abs(mean(sqrt(var_estimate_vec)) - sqrt(var_at_40)), 0.5)
})
