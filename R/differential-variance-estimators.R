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
#'   corresponding to the estimated conditional expexted squared outcome.
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
    dplyr::mutate(!!treatment_var_name := treatment_group)

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

  return(clean_counterfactual_tbl)

}
