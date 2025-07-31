test_that("serial cmldiffvarr estimates diff vars without errors", {

  library(dplyr)

  # generate randomized study data
  set.seed(3976468)
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

  # sample tibble
  sample_tbl <- slice_sample(population_tbl, n = 500)

  # unadjusted estimator for absolute diff var with propensity score
  expect_no_error(
    sample_tbl |>
      unadjdiffvar(
        estimand_type = "absolute",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_var_name = "propensity_score"
      )
  )

  # unadjusted estimator for absolute diff var without propensity score
  expect_no_error(
    sample_tbl |>
      unadjdiffvar(
        estimand_type = "absolute",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_var_name = NULL
      )
  )

  # unadjusted estimator for relative diff var with propensity score
  expect_no_error(
    sample_tbl |>
      unadjdiffvar(
        estimand_type = "relative",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_var_name = "propensity_score"
      )
  )

  # unadjusted estimator for relative diff var without propensity score
  expect_no_error(
    sample_tbl |>
      unadjdiffvar(
        estimand_type = "relative",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_var_name = NULL
      )
  )

})
