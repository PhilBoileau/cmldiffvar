estimate_propensity_score_fun <- function(
  clean_tbl,
  confounder_var_names,
  treatment_var_name,
  propensity_score_library,
  v_folds
) {

  # extract dependent and independent variables
  treatment_vec <- clean_tbl |>
    dplyr::select(dplyr::all_of(treatment_var_name)) |>
    as.numeric()
  confounders_tbl <- clean_tbl |>
    dplyr::select(dplyr::all_of(confounder_var_names))

  # estimate propensity score
  propensity_score_sl_fit <- SuperLearner::CV.SuperLearner(
    Y = treatment_vec,
    X = confounders_tbl,
    family = binomial(),
    SL.library = propensity_score_library,
    V = v_folds
  )

  return(propensity_score_sl_fit)

}

estimate_cond_exp_outcome_fun <- function(
    clean_tbl,
    confounder_var_names,
    treatment_var_name,
    outcome_var_name,
    cond_exp_outcome_library,
    v_folds
) {

  # extract dependent and independent variables
  outcome_vec <- clean_tbl |>
    dplyr::select(dplyr::all_of(outcome_var_name)) |>
    as.numeric()
  confounders_and_treatment_tbl <- clean_tbl |>
    dplyr::select(dplyr::all_of(confounder_var_names, outcome_var_name))

  # estimate conditional expected outcome
  cond_exp_outcome_sl_fit <- SuperLearner::CV.SuperLearner(
    Y = outcome_vec,
    X = confounders_and_treatment_tbl,
    family = gaussian(),
    SL.library = cond_exp_outcome_library,
    V = v_folds
  )

  return(cond_exp_outcome_sl_fit)

}

estimate_cond_exp_sq_outcome_fun <- function(
    clean_tbl,
    confounder_var_names,
    treatment_var_name,
    outcome_var_name,
    cond_exp_sq_outcome_library,
    v_folds
) {

  # extract dependent and independent variables
  outcome_vec <- clean_tbl |>
    dplyr::select(dplyr::all_of(outcome_var_name)) |>
    as.numeric()
  sq_outcome_vec <- outcome_vec^2
  confounders_and_treatment_tbl <- clean_tbl |>
    dplyr::select(dplyr::all_of(confounder_var_names, outcome_var_name))

  cond_exp_sq_outcome_sl_fit <- SuperLearner::CV.SuperLearner(
    Y = sq_outcome_vec,
    X = confounders_and_treatment_tbl,
    family = gaussian(),
    SL.library = cond_exp_sq_outcome_library,
    V = v_folds
  )

  return(cond_exp_sq_outcome_sl_fit)

}
