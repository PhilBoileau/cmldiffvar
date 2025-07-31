#' Generate Counterfactual Dataset
#'
#' `generate_counterfactural_tbl_fun()` generates a counterfactual version of
#' `clean_tbl`.
#'
#' @inheritParams estimate_propensity_score_fun
#' @inheritParams one_step_var_estimator_fun
#' @param propensity_score_var_name An optional `character` providing the column
#'   name of the treatment assignment indicator stored in `data_tbl`.
#' @param propensity_score_sl_fit A [SuperLearner::SuperLearner] object
#'   corresponding to the estimated propensity score.
#' @param cond_exp_outcome_sl_fit A [SuperLearner::SuperLearner] object
#'   corresponding to the estimated conditional expected outcome.
#' @param cond_exp_sq_outcome_sl_fit A [SuperLearner::SuperLearner] object
#'   corresponding to the estimated conditional expected squared outcome.
#'
#' @importFrom rlang :=
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
  propensity_score_var_name,
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
  if (is.null(propensity_score_var_name)) {
    pred_propensity_score <- SuperLearner::predict.SuperLearner(
      propensity_score_sl_fit,
      newdata = clean_counterfactual_tbl |>
        dplyr::select(dplyr::all_of(confounder_var_names)),
      onlySL = TRUE
    )$pred
  } else {
    pred_propensity_score <- clean_tbl[[propensity_score_var_name]]
  }

  pred_cond_exp_outcome <- SuperLearner::predict.SuperLearner(
    cond_exp_outcome_sl_fit,
    newdata = clean_counterfactual_tbl |>
      dplyr::select(
        dplyr::all_of(
          c(confounder_var_names, treatment_var_name)
        )
      ),
    onlySL = TRUE
  )$pred
  pred_cond_exp_sq_outcome <- SuperLearner::predict.SuperLearner(
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
#' @param estimand_type A `character` indicating whether to estimate an absolute
#'   or relative effect. Set this parameter to `"absolute"` to estimate the
#'   difference of group-specific variances. Set this parameter to `"relative"`
#'   to estimate the ratio of the group-specific variances.
#'
#' @returns A `numeric` estimate of the difference in treatment group variances.
#'
#' @keywords internal
#'
one_step_diff_var_estimator_fun <- function(
  clean_tbl,
  confounder_var_names,
  treatment_var_name,
  propensity_score_var_name,
  outcome_var_name,
  propensity_score_sl_fit,
  cond_exp_outcome_sl_fit,
  cond_exp_sq_outcome_sl_fit,
  estimand_type
) {

  # generate the counterfactual datasets
  clean_treatment_tbl <- generate_counterfactural_tbl_fun(
    clean_tbl,
    treatment_group = 1,
    confounder_var_names,
    treatment_var_name,
    propensity_score_var_name,
    propensity_score_sl_fit,
    cond_exp_outcome_sl_fit,
    cond_exp_sq_outcome_sl_fit
  )
  clean_control_tbl <- generate_counterfactural_tbl_fun(
    clean_tbl,
    treatment_group = 0,
    confounder_var_names,
    treatment_var_name,
    propensity_score_var_name,
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

  # estimate group-specific variances and compute EIFs
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

  # estimate differential variance and compute EIFs
  if (estimand_type == "absolute") {

    estimate <- treatment_group_var_est$estimate -
      control_group_var_est$estimate
    eif <- treatment_group_var_est$eif - control_group_var_est$eif

  } else if (estimand_type == "relative") {

    estimate <- treatment_group_var_est$estimate /
      control_group_var_est$estimate
    eif <- (treatment_group_var_est$eif / control_group_var_est$estimate) -
      (treatment_group_var_est$estimate / (control_group_var_est$estimate^2)) *
        control_group_var_est$eif

  }

  # prepare list of results
  results_ls <- list(
    estimate = estimate,
    eif = eif
  )

  return(results_ls)
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
  propensity_score_var_name,
  outcome_var_name,
  propensity_score_sl_fit,
  cond_exp_outcome_sl_fit,
  cond_exp_sq_outcome_sl_fit,
  estimand_type
) {

  # generate the counterfactual datasets
  clean_treatment_tbl <- generate_counterfactural_tbl_fun(
    clean_tbl,
    treatment_group = 1,
    confounder_var_names,
    treatment_var_name,
    propensity_score_var_name,
    propensity_score_sl_fit,
    cond_exp_outcome_sl_fit,
    cond_exp_sq_outcome_sl_fit
  )
  clean_control_tbl <- generate_counterfactural_tbl_fun(
    clean_tbl,
    treatment_group = 0,
    confounder_var_names,
    treatment_var_name,
    propensity_score_var_name,
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

  # estimate differential variance and compute EIFs
  if (estimand_type == "absolute") {

    estimate <- treatment_group_var_est$estimate -
      control_group_var_est$estimate
    eif <- treatment_group_var_est$eif - control_group_var_est$eif

  } else if (estimand_type == "relative") {

    estimate <- treatment_group_var_est$estimate /
      control_group_var_est$estimate
    eif <- (treatment_group_var_est$eif / control_group_var_est$estimate) -
      (treatment_group_var_est$estimate / (control_group_var_est$estimate^2)) *
      control_group_var_est$eif

  }

  # prepare list of results
  results_ls <- list(
    estimate = estimate,
    eif = eif
  )

  return(results_ls)
}


#' Cross-Fitted Differential Variance Estimation
#'
#' `cf_diff_var_estimator_fun()` produces a cross-fitted estimate of the select
#' differential variance estimand, using either a one-step or targeted maximum
#' likelihood estimator.
#'
#' @param fold A [fold][origami::make_folds] object indicating which
#'   observations in `clean_tbl` are members of the training and validation
#'   sets.
#' @inheritParams one_step_diff_var_estimator_fun
#' @inheritParams estimate_propensity_score_fun
#' @inheritParams estimate_cond_exp_outcome_fun
#' @inheritParams estimate_cond_exp_sq_outcome_fun
#' @param estimator_type A `character` indicating whether to use a one-step or a
#'   targeted maximum likelihood estimator in the cross-fitting procedure.
#'
#' @inherit one_step_diff_var_estimator_fun return
#'
#' @keywords internal
cf_diff_var_estimator_fun <- function(
  fold,
  clean_tbl,
  confounder_var_names,
  treatment_var_name,
  propensity_score_var_name,
  outcome_var_name,
  propensity_score_library,
  cond_exp_outcome_library,
  cond_exp_sq_outcome_library,
  num_nuisance_sl_folds,
  estimator_type,
  estimand_type
) {

  # split the data into training and validation
  train_tbl <- origami::training(clean_tbl)
  valid_tbl <- origami::validation(clean_tbl)

  # estimate the nuisance parameters on the training data
  if (is.null(propensity_score_var_name)) {
    propensity_score_sl_fit <- estimate_propensity_score_fun(
      train_tbl,
      confounder_var_names = confounder_var_names,
      treatment_var_name = treatment_var_name,
      propensity_score_library = propensity_score_library,
      num_nuisance_sl_folds = num_nuisance_sl_folds
    )
  } else {
    propensity_score_sl_fit <- NULL
  }
  cond_exp_outcome_sl_fit <- estimate_cond_exp_outcome_fun(
    train_tbl,
    confounder_var_names = confounder_var_names,
    treatment_var_name = treatment_var_name,
    outcome_var_name = outcome_var_name,
    cond_exp_outcome_library = cond_exp_outcome_library,
    num_nuisance_sl_folds = num_nuisance_sl_folds
  )
  cond_exp_sq_outcome_sl_fit <- estimate_cond_exp_sq_outcome_fun(
    train_tbl,
    confounder_var_names = confounder_var_names,
    treatment_var_name = treatment_var_name,
    outcome_var_name = outcome_var_name,
    cond_exp_sq_outcome_library = cond_exp_sq_outcome_library,
    num_nuisance_sl_folds = num_nuisance_sl_folds
  )

  # estimate the differential variances on the validation data
  if (estimator_type == "one-step") {

    diff_var_est <- one_step_diff_var_estimator_fun(
      valid_tbl,
      confounder_var_names = confounder_var_names,
      treatment_var_name = treatment_var_name,
      propensity_score_var_name = propensity_score_var_name,
      outcome_var_name = outcome_var_name,
      propensity_score_sl_fit = propensity_score_sl_fit,
      cond_exp_outcome_sl_fit = cond_exp_outcome_sl_fit,
      cond_exp_sq_outcome_sl_fit = cond_exp_sq_outcome_sl_fit,
      estimand_type = estimand_type
    )

  } else if (estimator_type == "tmle") {

    diff_var_est <- tml_diff_var_estimator_fun(
      valid_tbl,
      confounder_var_names = confounder_var_names,
      treatment_var_name = treatment_var_name,
      propensity_score_var_name = propensity_score_var_name,
      outcome_var_name = outcome_var_name,
      propensity_score_sl_fit = propensity_score_sl_fit,
      cond_exp_outcome_sl_fit = cond_exp_outcome_sl_fit,
      cond_exp_sq_outcome_sl_fit = cond_exp_sq_outcome_sl_fit,
      estimand_type = estimand_type
    )

  }

  # return the differential variance estimate and EIF
  valid_out_ls <- list(
    estimates = diff_var_est$estimate,
    eif = diff_var_est$eif
  )

  return(valid_out_ls)
}


unadjusted_diff_var_estimator_fun <- function(
  clean_tbl,
  treatment_var_name,
  propensity_score_var_name,
  outcome_var_name,
  estimand_type
) {

  # check if propensity score values are provided
  # use unadjusted propensity score estimator otherwise
  if (is.null(propensity_score_var_name)) {
    propensity_score_vec <- mean(clean_tbl[[treatment_var_name]])
  } else {
    propensity_score_vec <- clean_tbl[[propensity_score_var_name]]
  }

  # estimate the treatment-group specific variance
  treatment_group_var <- unadjusted_var_estimator_fun(
    treatment_group = 1,
    treatment_vec = clean_tbl[[treatment_var_name]],
    outcome_vec = clean_tbl[[outcome_var_name]],
    ps_est_vec = propensity_score_vec
  )

  # estimate the control-group specific variance
  control_group_var <- unadjusted_var_estimator_fun(
    treatment_group = 0,
    treatment_vec = clean_tbl[[treatment_var_name]],
    outcome_vec = clean_tbl[[outcome_var_name]],
    ps_est_vec = propensity_score_vec
  )

  # assemble the point estimate and eif of the estimate based on the estimand
  if (estimand_type == "absolute") {
    estimate <- treatment_group_var$estimate - control_group_var$estimate
    eif <- treatment_group_var$eif - control_group_var$eif
  } else if (estimand_type == "relative") {
    estimate <- treatment_group_var$estimate / control_group_var$estimate
    eif <- treatment_group_var$eif / control_group_var$estimate -
      treatment_group_var$estimate * control_group_var$eif /
      control_group_var$estimate^2
  }

  # prepare list of results
  results_ls <- list(
    estimate = estimate,
    eif = eif
  )

  return(results_ls)

}
