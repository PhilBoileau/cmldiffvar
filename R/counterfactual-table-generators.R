#' Generate Counterfactual Dataset For a Continuous Outcome
#'
#' `generate_counterfactural_tbl_fun()` generates a counterfactual version of
#' `clean_tbl`.
#'
#' @inheritParams estimate_propensity_score_fun
#' @inheritParams one_step_var_estimator_fun
#' @param propensity_score_adj_var_names A `character` vector providing the
#'   columns names of the adjustment set variables for propensity score
#'   estimation stored in `clean_tbl`.
#' @param cond_exp_outcome_adj_var_names A `character` vector providing the
#'   columns names of the adjustment set variables for outcome regression
#'   estimation stored in `clean_tbl`.
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
    propensity_score_adj_var_names,
    cond_exp_outcome_adj_var_names,
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
        dplyr::select(dplyr::all_of(propensity_score_adj_var_names)),
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
          c(cond_exp_outcome_adj_var_names, treatment_var_name)
        )
      ),
    onlySL = TRUE
  )$pred
  pred_cond_exp_sq_outcome <- SuperLearner::predict.SuperLearner(
    cond_exp_sq_outcome_sl_fit,
    newdata = clean_counterfactual_tbl |>
      dplyr::select(
        dplyr::all_of(
          c(cond_exp_outcome_adj_var_names, treatment_var_name)
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
