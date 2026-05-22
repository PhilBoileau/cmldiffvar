#' One-Step Estimator of Differential Variance
#'
#' `one_step_diff_var_estimator_fun()` estimates the differential variance using
#' a one-step estimator.
#'
#' @inheritParams generate_counterfactural_tbl_fun
#' @inheritParams estimate_cond_exp_outcome_fun
#' @param estimand_type A `character` indicating whether to estimate an absolute
#'   or relative effect. Set this parameter to `"absolute"` to estimate the
#'   difference of group-specific standard deviations. Set this parameter to
#'   `"relative"` to estimate the ratio of the group-specific variances.
#'
#' @returns A named list with the following components:
#'  * `estimate`: A `numeric` estimate of the selected estimand.
#'  * `eif`: A `numeric` vector of the efficient influence function of the
#'    selected estimand.
#'
#' @keywords internal
#'
one_step_diff_var_estimator_fun <- function(
  clean_tbl,
  propensity_score_adj_var_names,
  cond_exp_outcome_adj_var_names,
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
    propensity_score_adj_var_names,
    cond_exp_outcome_adj_var_names,
    treatment_var_name,
    propensity_score_var_name,
    propensity_score_sl_fit,
    cond_exp_outcome_sl_fit,
    cond_exp_sq_outcome_sl_fit
  )
  clean_control_tbl <- generate_counterfactural_tbl_fun(
    clean_tbl,
    treatment_group = 0,
    propensity_score_adj_var_names,
    cond_exp_outcome_adj_var_names,
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

    # check if variance estimates are approximately non-positive
    # set to sqrt(.Machine$double.eps) and output a warning if so
    estimate_bound <- sqrt(.Machine$double.eps)
    if (treatment_group_var_est$estimate < estimate_bound) {
      treatment_group_var_est$estimate <- estimate_bound
      warning(
        paste(
          "The treatment group's variance estimate is negative. Its value has",
          "been set to sqrt(.Machine$double.eps)."
        ),
        call. = FALSE
      )
    }
    if (control_group_var_est$estimate <= estimate_bound) {
      control_group_var_est$estimate <- estimate_bound
      warning(
        paste(
          "The control group's variance estimate is negative. Its value has",
          "been set to sqrt(.Machine$double.eps)."
        ),
        call. = FALSE
      )
    }

    # estimate the absolute difference variance estimand
    estimate <- sqrt(treatment_group_var_est$estimate) -
      sqrt(control_group_var_est$estimate)
    eif <- treatment_group_var_est$eif /
      (2 * sqrt(treatment_group_var_est$estimate)) -
      control_group_var_est$eif /
        (2 * sqrt(control_group_var_est$estimate))

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
  propensity_score_adj_var_names,
  cond_exp_outcome_adj_var_names,
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
    propensity_score_adj_var_names,
    cond_exp_outcome_adj_var_names,
    treatment_var_name,
    propensity_score_var_name,
    propensity_score_sl_fit,
    cond_exp_outcome_sl_fit,
    cond_exp_sq_outcome_sl_fit
  )
  clean_control_tbl <- generate_counterfactural_tbl_fun(
    clean_tbl,
    treatment_group = 0,
    propensity_score_adj_var_names,
    cond_exp_outcome_adj_var_names,
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

    # check if variance estimates are approximately non-positive
    # set to sqrt(.Machine$double.eps) and output a warning if so
    estimate_bound <- sqrt(.Machine$double.eps)
    if (treatment_group_var_est$estimate < estimate_bound) {
      treatment_group_var_est$estimate <- estimate_bound
      warning(
        paste(
          "The treatment group's variance estimate is negative. Its value has",
          "been set to sqrt(.Machine$double.eps)."
        ),
        call. = FALSE
      )
    }
    if (control_group_var_est$estimate <= estimate_bound) {
      control_group_var_est$estimate <- estimate_bound
      warning(
        paste(
          "The control group's variance estimate is negative. Its value has",
          "been set to sqrt(.Machine$double.eps)."
        ),
        call. = FALSE
      )
    }

    # estimate the absolute difference variance estimand
    estimate <- sqrt(treatment_group_var_est$estimate) -
      sqrt(control_group_var_est$estimate)
    eif <- treatment_group_var_est$eif /
      (2 * sqrt(treatment_group_var_est$estimate)) -
      control_group_var_est$eif / (2 * sqrt(control_group_var_est$estimate))

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
  propensity_score_adj_var_names,
  cond_exp_outcome_adj_var_names,
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
      adj_set_var_names = propensity_score_adj_var_names,
      treatment_var_name = treatment_var_name,
      propensity_score_library = propensity_score_library,
      num_nuisance_sl_folds = num_nuisance_sl_folds
    )
  } else {
    propensity_score_sl_fit <- NULL
  }
  cond_exp_outcome_sl_fit <- estimate_cond_exp_outcome_fun(
    train_tbl,
    adj_set_var_names = cond_exp_outcome_adj_var_names,
    treatment_var_name = treatment_var_name,
    outcome_var_name = outcome_var_name,
    cond_exp_outcome_library = cond_exp_outcome_library,
    num_nuisance_sl_folds = num_nuisance_sl_folds
  )
  cond_exp_sq_outcome_sl_fit <- estimate_cond_exp_sq_outcome_fun(
    train_tbl,
    adj_set_var_names = cond_exp_outcome_adj_var_names,
    treatment_var_name = treatment_var_name,
    outcome_var_name = outcome_var_name,
    cond_exp_sq_outcome_library = cond_exp_sq_outcome_library,
    num_nuisance_sl_folds = num_nuisance_sl_folds
  )

  # estimate the differential variances on the validation data
  if (estimator_type == "one-step") {

    diff_var_est <- one_step_diff_var_estimator_fun(
      valid_tbl,
      propensity_score_adj_var_names = propensity_score_adj_var_names,
      cond_exp_outcome_adj_var_names = cond_exp_outcome_adj_var_names,
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
      propensity_score_adj_var_names = propensity_score_adj_var_names,
      cond_exp_outcome_adj_var_names = cond_exp_outcome_adj_var_names,
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


#' Estimator of Differential Variance
#'
#' `unadjusted_diff_var_estimator_fun()` estimates the differential variance
#' using an unadjusted estimator. This estimator is regular and asymptotically
#' linear when the treatment assignment mechanism does not depend on
#' pre-treatment covariates, such as in a randomized controlled trial.
#'
#' @inheritParams one_step_diff_var_estimator_fun
#'
#' @inherit one_step_diff_var_estimator_fun return
#'
#' @keywords internal
#'
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
    estimate <- sqrt(treatment_group_var$estimate) -
      sqrt(control_group_var$estimate)
    eif <- treatment_group_var$eif / (2 * sqrt(treatment_group_var$estimate)) -
      control_group_var$eif / (2 * sqrt(control_group_var$estimate))
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
