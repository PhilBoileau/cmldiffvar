test_that("estimators of treatment group-specific variance is consistent", {

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
    )

    # calculate variance of treatment TML estimate
    tml_var_treatment_est_vec[iter] <- tml_var_estimator_fun(
      treatment_group = 1,
      treatment_vec = sample_tbl$treatment,
      outcome_vec = sample_treatment_tbl$outcome,
      ps_est_vec = ps_est,
      cond_exp_outcome_est_vec = cond_outcome_treatment_est,
      cond_exp_sq_outcome_est_vec = cond_sq_outcome_treatment_est
    )

    # calculate variance of control one-step estimate
    one_step_var_control_est_vec[iter] <- one_step_var_estimator_fun(
      treatment_group = 0,
      sample_tbl$treatment,
      sample_control_tbl$outcome,
      ps_est,
      cond_outcome_control_est,
      cond_sq_outcome_control_est,
      one_step_mean_control_est
    )

    # calculate variance of control TML estimate
    tml_var_control_est_vec[iter] <- tml_var_estimator_fun(
      treatment_group = 0,
      sample_tbl$treatment,
      sample_control_tbl$outcome,
      ps_est,
      cond_outcome_control_est,
      cond_sq_outcome_control_est
    )

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
