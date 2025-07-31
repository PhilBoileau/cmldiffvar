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
