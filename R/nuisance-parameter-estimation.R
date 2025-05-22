#' SuperLearned Propensity Score
#'
#' `estimate_propensity_score_fun()` estimates the propensity score using a
#' cross-validated SuperLearner estimator implemented in the
#' [SuperLearner::CV.SuperLearner()] function.
#'
#' @param clean_tbl A pre-processed [tibble] that is ready for nuisance
#'   parameter estimation.
#' @param confounder_var_names A `character` vector providing the columns names
#'   of the treatment--outcome confounders stored in `clean_tbl`.
#' @param treatment_var_name A `character` providing the column name of the
#'   treatment assignment indicator stored in `clean_tbl`.
#' @param propensity_score_library A `character` vector of candidate learners
#'   used by the SuperLearner estimator.
#' @param num_folds A `numeric` indicating the number of folds to use in
#'   cross-validated SuperLearner estimator.
#'
#' @keywords internal
#'
#' @returns A [SuperLearner::CV.SuperLearner] class object containing the
#'   estimated propensity score.
#'
estimate_propensity_score_fun <- function(
  clean_tbl,
  confounder_var_names,
  treatment_var_name,
  propensity_score_library,
  num_folds
) {

  # extract dependent and independent variables
  treatment_vec <- clean_tbl |>
    dplyr::select(dplyr::all_of(treatment_var_name)) |>
    dplyr::pull()
  confounders_tbl <- clean_tbl |>
    dplyr::select(dplyr::all_of(confounder_var_names))

  # estimate propensity score
  propensity_score_sl_fit <- SuperLearner::SuperLearner(
    Y = treatment_vec,
    X = confounders_tbl,
    family = binomial(),
    SL.library = propensity_score_library,
    cvControl = list("V" = num_folds)
  )

  return(propensity_score_sl_fit)

}

#' SuperLearned Conditional Expected Outcome
#'
#' `estimate_cond_exp_outcome_fun()` estimates the conditional expected outcome
#' using a cross-validated SuperLearner estimator implemented in the
#' [SuperLearner::CV.SuperLearner()] function.
#'
#' @inheritParams estimate_propensity_score_fun
#' @param outcome_var_name A `character` providing the column name of the
#'   outcome variable stored in `clean_tbl`.
#' @param cond_exp_outcome_library A `character` vector of candidate learners
#'   used by the SuperLearner estimator.
#'
#' @keywords internal
#'
#' @returns A [SuperLearner::CV.SuperLearner] class object containing the
#'   estimated conditional expected outcome.
estimate_cond_exp_outcome_fun <- function(
    clean_tbl,
    confounder_var_names,
    treatment_var_name,
    outcome_var_name,
    cond_exp_outcome_library,
    num_folds
) {

  # extract dependent and independent variables
  outcome_vec <- clean_tbl |>
    dplyr::select(dplyr::all_of(outcome_var_name)) |>
    dplyr::pull()
  confounders_and_treatment_tbl <- clean_tbl |>
    dplyr::select(dplyr::all_of(c(confounder_var_names, treatment_var_name)))

  # estimate conditional expected outcome
  cond_exp_outcome_sl_fit <- SuperLearner::SuperLearner(
    Y = outcome_vec,
    X = confounders_and_treatment_tbl,
    family = gaussian(),
    SL.library = cond_exp_outcome_library,
    cvControl = list("V" = num_folds)
  )

  return(cond_exp_outcome_sl_fit)

}

#' SuperLearned Conditional Expected Outcome
#'
#' `estimate_cond_exp_sq_outcome_fun()` estimates the conditional expected
#' squared outcome using a cross-validated SuperLearner estimator implemented in
#' the [SuperLearner::CV.SuperLearner()] function.
#'
#' @inheritParams estimate_cond_exp_outcome_fun
#' @param cond_exp_sq_outcome_library A `character` vector of candidate learners
#'   used by the SuperLearner estimator.
#'
#' @keywords internal
#'
#' @returns A [SuperLearner::CV.SuperLearner] class object containing the
#'   estimated conditional expected squared outcome.
#'
estimate_cond_exp_sq_outcome_fun <- function(
    clean_tbl,
    confounder_var_names,
    treatment_var_name,
    outcome_var_name,
    cond_exp_sq_outcome_library,
    num_folds
) {

  # extract dependent and independent variables
  outcome_vec <- clean_tbl |>
    dplyr::select(dplyr::all_of(outcome_var_name)) |>
    dplyr::pull()
  sq_outcome_vec <- outcome_vec^2
  confounders_and_treatment_tbl <- clean_tbl |>
    dplyr::select(dplyr::all_of(c(confounder_var_names, treatment_var_name)))

  cond_exp_sq_outcome_sl_fit <- SuperLearner::SuperLearner(
    Y = sq_outcome_vec,
    X = confounders_and_treatment_tbl,
    family = gaussian(),
    SL.library = cond_exp_sq_outcome_library,
    cvControl = list("V" = num_folds)
  )

  return(cond_exp_sq_outcome_sl_fit)

}
