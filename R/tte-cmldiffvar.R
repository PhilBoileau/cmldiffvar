#' Causal Machine Learning for Differential Variance with Time-to-Event Outcomes
#'
#' @description
#' `tte_cmldiffvar()` performs a homogeneous treatment effect hypothesis testing
#' procedure based on the relative differential variance in datasets with
#' right-censored, restricted time-to-event outcomes. Two types of homogeneous
#' treatment hypotheses can be considered: a homogeneous additive effect, and a
#' homogeneous multiplicative effect.
#'
#' @param data_tbl A [`data.frame`] or [`tibble`][tibble::tibble].
#' @param restriction_time A `numeric` integer corresponding to the
#'   time-to-event outcomes restriction time.
#' @param homogeneity_type A `character` indicating the type of homogeneous
#'   treatment effect under the null hypothesis: `"additive"` or
#'   `"multiplicative"`.
#' @param homogeneous_effect_range_min A `numeric` providing the minimum
#'   suspected possible effect size. The scaling of the effect size is dictated
#'   by the `homogeneity_type` argument.
#' @param homogeneous_effect_range_max A `numeric` providing the maximum
#'   suspected possible effect size. The scaling of the effect size is dictated
#'   by the `homogeneity_type` argument.
#' @param propensity_score_adj_var_names A `character` vector providing the
#'   columns names of the adjustment set variables for propensity score
#'   estimation stored in `data_tbl`. Ignored if `propensity_score_var_name` is
#'   `NULL`.
#' @param cond_event_haz_adj_var_names A `character` vector providing the
#'   columns names of the adjustment set variables for conditional event hazard
#'   rate estimation stored in `data_tbl`.
#' @param cond_censoring_haz_adj_var_names A `character` vector providing the
#'   columns names of the adjustment set variables for conditional censoring
#'   hazard rate estimation stored in `data_tbl`.
#' @param treatment_var_name A `character` providing the column name of the
#'   treatment assignment indicator stored in `data_tbl`.
#' @param propensity_score_var_name An optional `character` providing the column
#'   name of the treatment assignment indicator stored in `data_tbl`. Defaults
#'   to `NULL`. See the Details section for more information.
#' @param outcome_var_name A `character` providing the column name of the
#'   time-to-event outcome variable stored in `data_tbl`.
#' @param censoring_var_name A `character` providing the column name of the
#'   censoring indicator stored in `data_tbl`.
#' @param propensity_score_sl_library A `character` vector of candidate learners
#'   used by the SuperLearner estimator of the propensity score. Defaults to
#'   `c("SL.glm", "SL.earth", "SL.ranger")`.
#' @param cond_event_haz_sl_library A `character` vector of candidate learners
#'   used by the SuperLearner estimator of the conditional event hazard rate.
#'   Defaults to `c("SL.glm", "SL.earth", "SL.ranger")`.
#' @param cond_censoring_haz_sl_library A `character` vector of candidate learners
#'   used by the SuperLearner estimator of the conditional event hazard rate.
#'   Defaults to `c("SL.glm", "SL.earth", "SL.ranger")`.
#' @param num_nuisance_sl_folds A `numeric` indicating the number of folds to
#'   use in cross-validated SuperLearner estimators. Defaults to `5`.
#'
#' @returns A `list` containing two items:
#'   - `hypothesis_test_tbl`: A one-row [tibble][tibble::tibble] reporting the
#'     result of the homogeneous treatment effect hypothesis testing procedure.
#'   - `inference_tbl`: A [tibble][tibble::tibble] reporting the differential
#'     variance estimates and associated p-values considered for the homogeneous
#'     treatment effect hypothesis testing procedure.
#'
#' @export
#'
tte_cmldiffvar <- function(
  data_tbl,
  restriction_time,
  homogeneity_type,
  homogeneous_effect_range_min,
  homogeneous_effect_range_max,
  propensity_score_adj_var_names,
  cond_event_haz_adj_var_names,
  cond_censoring_haz_adj_var_names,
  treatment_var_name,
  propensity_score_var_name = NULL,
  outcome_var_name,
  censoring_var_name,
  propensity_score_sl_library = c("SL.glm", "SL.ranger"),
  cond_event_haz_sl_library = c("SL.glm", "SL.ranger"),
  cond_censoring_haz_sl_library = c("SL.glm", "SL.ranger"),
  num_nuisance_sl_folds = 5
) {

  # checks ----
  checkmate::assert_choice(homogeneity_type, c("additive", "multiplicative"))
  if (is.null(propensity_score_var_name)) {
    checkmate::assert_character(propensity_score_adj_var_names)
  }
  checkmate::assert_character(cond_event_haz_adj_var_names)
  checkmate::assert_character(cond_censoring_haz_adj_var_names)
  checkmate::assert_character(treatment_var_name)
  checkmate::assert_character(outcome_var_name)
  checkmate::assert_int(num_nuisance_sl_folds, lower = 2, upper = 20)

  # check dataset contains appropriate variables
  if (is.null(propensity_score_var_name)) {
    clean_tbl_var_names <- c(
      propensity_score_adj_var_names, cond_event_haz_adj_var_names,
      cond_censoring_haz_adj_var_names, treatment_var_name, censoring_var_name,
      outcome_var_name
    )
  } else {
    checkmate::assert_character(propensity_score_var_name)
    checkmate::assert_numeric(
      data_tbl[[propensity_score_var_name]], lower = 0.001, upper = 0.999
    )
    clean_tbl_var_names <- c(
      propensity_score_adj_var_names, cond_event_haz_adj_var_names,
      cond_censoring_haz_adj_var_names, treatment_var_name, censoring_var_name,
      propensity_score_var_name, outcome_var_name
    )
  }
  checkmate::assert_names(clean_tbl_var_names, subset.of = colnames(data_tbl))

  # only retain relevant columns in data_tbl
  clean_tbl <- data_tbl |> dplyr::select(dplyr::all_of(clean_tbl_var_names))

  # check dataset contains only numerics and factors
  checkmate::assert_data_frame(
    clean_tbl,
    types = c("numeric", "factor"),
    any.missing = FALSE,
    min.rows = 50
  )

  # transform dataset into long format ----
  clean_long_tbl <- clean_tbl |>
    melt_tte_data_fun(
      baseline_var_names = unique(
          c(propensity_score_adj_var_names,
            cond_event_haz_adj_var_names,
            cond_censoring_haz_adj_var_names)
        ),
      treatment_var_name = treatment_var_name,
      outcome_var_name = outcome_var_name,
      censoring_var_name = censoring_var_name,
      propensity_score_var_name = propensity_score_var_name,
      time_cutoff = restriction_time
    )

  # compute nuisance parameter estimates ----

  if (is.null(propensity_score_var_name)) {
    propensity_score_sl_fit <- estimate_propensity_score_fun(
      data_tbl,
      adj_set_var_names = propensity_score_adj_var_names,
      treatment_var_name = treatment_var_name,
      propensity_score_library = propensity_score_sl_library,
      num_nuisance_sl_folds = num_nuisance_sl_folds
    )
  } else {
    propensity_score_sl_fit <- NULL
  }

  cond_event_haz_sl_fit <- estimate_cond_event_haz_fun(
    clean_long_tbl = clean_long_tbl,
    adj_set_var_names = cond_event_haz_adj_var_names,
    treatment_var_name = treatment_var_name,
    outcome_var_name = "L",
    cond_event_haz_library = cond_event_haz_sl_library,
    num_nuisance_sl_folds = num_nuisance_sl_folds
  )

  cond_censoring_haz_sl_fit <- estimate_cond_censoring_haz_fun(
    clean_long_tbl = clean_long_tbl,
    adj_set_var_names = cond_censoring_haz_adj_var_names,
    treatment_var_name = treatment_var_name,
    outcome_var_name = "R",
    cond_censoring_haz_library = cond_censoring_haz_sl_library,
    num_nuisance_sl_folds = num_nuisance_sl_folds
  )

  # produce the differential variance inference data ----
  inference_tbl <- one_step_tte_diff_var_estimator_fun(
    clean_long_tbl = clean_long_tbl,
    propensity_score_adj_var_names = propensity_score_adj_var_names,
    cond_event_haz_adj_var_names = cond_event_haz_adj_var_names,
    cond_censoring_haz_adj_var_names = cond_censoring_haz_adj_var_names,
    treatment_var_name = treatment_var_name,
    propensity_score_var_name = propensity_score_var_name,
    propensity_score_sl_fit = propensity_score_sl_fit,
    cond_event_haz_sl_fit = cond_event_haz_sl_fit,
    cond_censoring_haz_sl_fit = cond_censoring_haz_sl_fit,
    restriction_time = restriction_time,
    effect_range_vec =
      c(homogeneous_effect_range_min, homogeneous_effect_range_max),
    homogeneity_type = homogeneity_type
  )

  # assemble and output results ----
  list(
    "hypothesis_test_tbl" = dplyr::tibble(
      hypothesis_type = paste(
        "homogeneous", homogeneity_type, "treatment effect"
        ),
      p_value = inference_tbl |> dplyr::pull(.data$p_value) |> max()
    ),
    "inference_tbl" = inference_tbl
  )
}
