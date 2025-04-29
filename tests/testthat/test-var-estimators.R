test_that("estimator of treatment group-specific variance is consistent", {

  # generate a large sample from a simple DGP
  set.seed(83452235)
  pop_size <- 1000000
  confounder_vec <- rnorm(n = pop_size)
  ps_vec <- plogis(0.1 * confounder_vec)
  treatment_vec <- sapply(
    ps_vec,
    function(ps_ind) {rbinom(n = 1, size = 1, prob = ps_ind)}
  )
  outcome_treatment_vec <- sapply(
    seq_len(pop_size),
    function(obs_idx) {
      rnorm(
        n = 1,
        mean = 1 + confounder_vec[obs_idx],
        sd = 3
      )
    }
  )
  outcome_control_vec <- sapply(
    seq_len(pop_size),
    function(obs_idx) {
      rnorm(
        n = 1,
        mean = confounder_vec[obs_idx],
        sd = 1
      )
    }
  )
  outcome_obs_vec <- treatment_vec * outcome_treatment_vec +
    (1 - treatment_vec) * outcome_control_vec

  population_tbl <- dplyr::tibble(
    confounder = confounder_vec,
    treatment = treatment_vec,
    outcome = outcome_obs_vec
  )

  # calculate estimand
  var_treatment <- var(outcome_treatment_vec)
  var_control <- var(outcome_control_vec)

  # grab a sample of the population
  sample_tbl <- dplyr::slice_sample(population_tbl, n = 1000) |>
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
    sq_outcome ~ treatment + confounder,
    data = sample_tbl,
    degree = 2
  )

  # predict conditional outcomes under treatment
  sample_treatment_tbl <- sample_tbl |>
    dplyr::mutate(treatment = 1)
  cond_outcome_treatment_est <- predict(cond_outcome_fit, sample_treatment_tbl)
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

  # calculate variance of treatment estimate
  one_step_var_treatment_est <- one_step_var_estimator_fun(
    treatment_group = 1,
    sample_tbl$treatment,
    sample_treatment_tbl$outcome,
    ps_est,
    cond_outcome_treatment_est,
    cond_sq_outcome_treatment_est,
    one_step_mean_treatment_est
  )

  # calculate variance of control estimate
  one_step_var_control_est <- one_step_var_estimator_fun(
    treatment_group = 0,
    sample_tbl$treatment,
    sample_control_tbl$outcome,
    ps_est,
    cond_outcome_control_est,
    cond_sq_outcome_control_est,
    one_step_mean_control_est
  )

  # check error
  expect_lt(abs(one_step_var_treatment_est - var_treatment), 1)
  expect_lt(abs(one_step_var_control_est - var_control), 1)

})
