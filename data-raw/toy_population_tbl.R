# generate population-level data from a simple DGP

# set seed for reproducibility
set.seed(83452235)

# set population size
pop_size <- 100000

# generate the confounder variable
confounder_vec <- rnorm(n = pop_size)

# generate the propensity score
ps_vec <- plogis(confounder_vec)

# generate treatment assignment based on propensity score
treatment_vec <- sapply(
  ps_vec,
  function(ps_ind) {rbinom(n = 1, size = 1, prob = ps_ind)}
)

# generate the potential outcome under treatment
outcome_treatment_vec <- sapply(
  seq_len(pop_size),
  function(obs_idx) {
    rnorm(
      n = 1,
      mean = 3 + confounder_vec[obs_idx],
      sd = 3
    )
  }
)

# generate the potential outcome under control
outcome_control_vec <- sapply(
  seq_len(pop_size),
  function(obs_idx) {
    rnorm(
      n = 1,
      mean = 1 + confounder_vec[obs_idx],
      sd = 1
    )
  }
)

# specify the observed outcome based on treatment assignment
outcome_obs_vec <- treatment_vec * outcome_treatment_vec +
  (1 - treatment_vec) * outcome_control_vec

# assemble the population tibble
toy_population_tbl <- dplyr::tibble(
  confounder = confounder_vec,
  treatment = treatment_vec,
  outcome = outcome_obs_vec,
  potential_outcome_treatment = outcome_treatment_vec,
  potential_outcome_control = outcome_control_vec
)

# add population data to package
usethis::use_data(toy_population_tbl, overwrite = TRUE, compress = "xz")
