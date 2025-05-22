#' Generate Counterfactual Dataset
#'
#' `generate_counterfactural_tbl_fun()` generates a counterfactual version of
#' `clean_tbl`.
#'
#' @inheritParams estimate_propensity_score_fun
#' @inheritParams one_step_var_estimator_fun
#' @param propensity_score_sl_fit A [SuperLearner::SuperLearner] object
#'   corresponding to the estimated propensity score.
#' @param cond_exp_outcome_sl_fit A [SuperLearner::SuperLearner] object
#'   corresponding to the estimated conditional expected outcome.
#' @param cond_exp_sq_outcome_sl_fit A [SuperLearner::SuperLearner] object
#'   corresponding to the estimated conditional expected squared outcome.
#'
#' @keywords internal
#'
#' @returns A counterfactual  [tibble], with the additional of propensity score,
#'   conditional expected outcome, and conditional expected square outcome
#'   estimates.
#'
generate_counterfactural_tbl_fun <- function(
  clean_tbl,
  treatment_group,
  confounder_var_names,
  treatment_var_name,
  propensity_score_sl_fit,
  cond_exp_outcome_sl_fit,
  cond_exp_sq_outcome_sl_fit
) {

  # create counterfactual dataset: all units assigned to treatment_group
  clean_counterfactual_tbl <- clean_tbl |>
    dplyr::mutate(
      !!treatment_var_name := treatment_group,
      counterfactual_treatment = treatment_group
    )

  # predict nuisance parameters under counterfactual scenarios
  pred_propensity_score <- predict(
    propensity_score_sl_fit,
    newdata = clean_counterfactual_tbl |>
      dplyr::select(dplyr::all_of(confounder_var_names)),
    onlySL = TRUE
  )$pred
  pred_cond_exp_outcome <- predict(
    cond_exp_outcome_sl_fit,
    newdata = clean_counterfactual_tbl |>
      dplyr::select(
        dplyr::all_of(
          c(confounder_var_names, treatment_var_name)
        )
      ),
    onlySL = TRUE
  )$pred
  pred_cond_exp_sq_outcome <- predict(
    cond_exp_sq_outcome_sl_fit,
    newdata = clean_counterfactual_tbl |>
      dplyr::select(
        dplyr::all_of(
          c(confounder_var_names, treatment_var_name)
        )
      ),
    onlySL = TRUE
  )$pred

  # add predictions to counterfactual dataset
  clean_counterfactual_tbl <- clean_counterfactual_tbl |>
    dplyr::mutate(
      pred_propensity_score = pred_propensity_score,
      pred_cond_exp_outcome = pred_cond_exp_outcome,
      pred_cond_exp_sq_outcome = pred_cond_exp_sq_outcome
    )

  # add back true treatment assignment vector for record keeping
  clean_counterfactual_tbl[[treatment_var_name]] <-
    clean_tbl[[treatment_var_name]]

  return(clean_counterfactual_tbl)

}


#' One-Step Estimator of Differential Variance
#'
#' `one_step_diff_var_estimator_fun()` estimates the differential variance using
#' a one-step estimator.
#'
#' @inheritParams generate_counterfactural_tbl_fun
#' @inheritParams estimate_cond_exp_outcome_fun
#'
#' @returns A `numeric` estimate of the difference in treatment group variances.
#'
#' @keywords internal
#'
one_step_diff_var_estimator_fun <- function(
  clean_tbl,
  confounder_var_names,
  treatment_var_name,
  outcome_var_name,
  propensity_score_sl_fit,
  cond_exp_outcome_sl_fit,
  cond_exp_sq_outcome_sl_fit
) {

  # generate the counterfactual datasets
  clean_treatment_tbl <- generate_counterfactural_tbl_fun(
    clean_tbl,
    treatment_group = 1,
    confounder_var_names,
    treatment_var_name,
    propensity_score_sl_fit,
    cond_exp_outcome_sl_fit,
    cond_exp_sq_outcome_sl_fit
  )
  clean_control_tbl <- generate_counterfactural_tbl_fun(
    clean_tbl,
    treatment_group = 0,
    confounder_var_names,
    treatment_var_name,
    propensity_score_sl_fit,
    cond_exp_outcome_sl_fit,
    cond_exp_sq_outcome_sl_fit
  )

  # estimate group-specific means
  treatment_group_mean_est <- one_step_mean_estimator_fun(
    treatment_group = 1,
    treatment_vec = clean_treatment_tbl[[treatment_var_name]],
    outcome_vec = clean_treatment_tbl[[outcome_var_name]],
    ps_est_vec = clean_treatment_tbl$pred_propensity_score,
    cond_exp_outcome_est_vec = clean_treatment_tbl$pred_cond_exp_outcome
  )
  control_group_mean_est <- one_step_mean_estimator_fun(
    treatment_group = 0,
    treatment_vec = clean_control_tbl[[treatment_var_name]],
    outcome_vec = clean_control_tbl[[outcome_var_name]],
    ps_est_vec = clean_control_tbl$pred_propensity_score,
    cond_exp_outcome_est_vec = clean_control_tbl$pred_cond_exp_outcome
  )

  # estimate group-specific variances
  treatment_group_var_est <- one_step_var_estimator_fun(
    treatment_group = 1,
    treatment_vec = clean_treatment_tbl[[treatment_var_name]],
    outcome_vec = clean_treatment_tbl[[outcome_var_name]],
    ps_est_vec = clean_treatment_tbl$pred_propensity_score,
    cond_exp_outcome_est_vec = clean_treatment_tbl$pred_cond_exp_outcome,
    cond_exp_sq_outcome_est_vec = clean_treatment_tbl$pred_cond_exp_sq_outcome,
    mean_est = treatment_group_mean_est
  )
  control_group_var_est <- one_step_var_estimator_fun(
    treatment_group = 0,
    treatment_vec = clean_control_tbl[[treatment_var_name]],
    outcome_vec = clean_control_tbl[[outcome_var_name]],
    ps_est_vec = clean_control_tbl$pred_propensity_score,
    cond_exp_outcome_est_vec = clean_control_tbl$pred_cond_exp_outcome,
    cond_exp_sq_outcome_est_vec = clean_control_tbl$pred_cond_exp_sq_outcome,
    mean_est = control_group_mean_est
  )

  # estimate differential variance
  os_diff_var_est <- treatment_group_var_est - control_group_var_est

  return(os_diff_var_est)
}

#' Targeted Minimum Loss-Based Estimator Estimator of Differential Variance
#'
#' `tml_diff_var_estimator_fun()` estimates the differential variance using
#' a targeted minimum loss-based estimator.
#'
#' @inheritParams one_step_diff_var_estimator_fun
#'
#' @inherit one_step_diff_var_estimator_fun return
#'
#' @keywords internal
#'
tml_diff_var_estimator_fun <- function(
  clean_tbl,
  confounder_var_names,
  treatment_var_name,
  outcome_var_name,
  propensity_score_sl_fit,
  cond_exp_outcome_sl_fit,
  cond_exp_sq_outcome_sl_fit
) {

  # generate the counterfactual datasets
  clean_treatment_tbl <- generate_counterfactural_tbl_fun(
    clean_tbl,
    treatment_group = 1,
    confounder_var_names,
    treatment_var_name,
    propensity_score_sl_fit,
    cond_exp_outcome_sl_fit,
    cond_exp_sq_outcome_sl_fit
  )
  clean_control_tbl <- generate_counterfactural_tbl_fun(
    clean_tbl,
    treatment_group = 0,
    confounder_var_names,
    treatment_var_name,
    propensity_score_sl_fit,
    cond_exp_outcome_sl_fit,
    cond_exp_sq_outcome_sl_fit
  )

  # estimate group-specific variances
  treatment_group_var_est <- tml_var_estimator_fun(
    treatment_group = 1,
    treatment_vec = clean_treatment_tbl[[treatment_var_name]],
    outcome_vec = clean_treatment_tbl[[outcome_var_name]],
    ps_est_vec = clean_treatment_tbl$pred_propensity_score,
    cond_exp_outcome_est_vec = clean_treatment_tbl$pred_cond_exp_outcome,
    cond_exp_sq_outcome_est_vec = clean_treatment_tbl$pred_cond_exp_sq_outcome
  )
  control_group_var_est <- tml_var_estimator_fun(
    treatment_group = 0,
    treatment_vec = clean_control_tbl[[treatment_var_name]],
    outcome_vec = clean_control_tbl[[outcome_var_name]],
    ps_est_vec = clean_control_tbl$pred_propensity_score,
    cond_exp_outcome_est_vec = clean_control_tbl$pred_cond_exp_outcome,
    cond_exp_sq_outcome_est_vec = clean_control_tbl$pred_cond_exp_sq_outcome
  )

  # estimate differential variance
  tml_diff_var_est <- treatment_group_var_est - control_group_var_est

  return(tml_diff_var_est)
}
