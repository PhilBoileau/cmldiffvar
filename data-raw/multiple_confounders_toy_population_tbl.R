# generate population-level data from a simple DGP with multiple confounders

# set seed for reproducibility
set.seed(42)

# Set population size
pop_size <- 100000

# Generate multiple confounders (W1, W2, W3)
confounders <- matrix(rnorm(pop_size * 3), ncol = 3)
colnames(confounders) <- c("W1", "W2", "W3")

# Generate propensity scores
confounder_score <-
  0.5 * confounders[,1] - 0.3 * confounders[,2] + 0.2 * confounders[,3]
ps_vec <- plogis(confounder_score)

# Generate treatment assignment based on propensity score
treatment_vec <- rbinom(pop_size, 1, ps_vec)

# Generate potential outcomes under treatment
outcome_treatment_vec <- sapply(
  seq_len(pop_size),
  function(obs_idx) {
    rnorm(
      n = 1,
      mean =
        3 + 0.5 * confounders[obs_idx,1] + 0.3 * confounders[obs_idx,2],
      sd = 3
    )
  }
)

# Generate potential outcomes under control
outcome_control_vec <- sapply(
  seq_len(pop_size),
  function(obs_idx) {
    rnorm(
      n = 1,
      mean = 1 + 0.5 * confounders[obs_idx,1] - 0.2 * confounders[obs_idx,2],
      sd = 1
    )
  }
)

# Compute observed outcome
outcome_obs_vec <- treatment_vec * outcome_treatment_vec +
  (1 - treatment_vec) * outcome_control_vec

# Assemble into a tibble
multiple_confounders_toy_population_tbl <- tibble(
  W1 = confounders[,1],
  W2 = confounders[,2],
  W3 = confounders[,3],
  treatment = treatment_vec,
  outcome = outcome_obs_vec,
  potential_outcome_treatment = outcome_treatment_vec,
  potential_outcome_control = outcome_control_vec
)

# add population data to package
usethis::use_data(multiple_confounders_toy_population_tbl,
                  overwrite = TRUE,
                  compress = "xz")
