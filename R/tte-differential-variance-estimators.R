one_step_tte_diff_var_estimator_fun <- function(
  clean_long_tbl,
  propensity_score_adj_var_names,
  cond_event_haz_adj_var_names,
  cond_censoring_haz_adj_var_names,
  treatment_var_name,
  propensity_score_var_name,
  propensity_score_sl_fit,
  cond_event_haz_sl_fit,
  cond_censoring_haz_sl_fit,
  max_time_cutoff,
  effect_range_vec,
  homogeneity_type
) {

  # generate the counterfactual datasets
  clean_long_treatment_tbl <- generate_long_counterfactural_tbl_fun(
    clean_long_tbl = clean_long_tbl,
    treatment_group = 1,
    propensity_score_adj_var_names = propensity_score_adj_var_names,
    cond_event_haz_adj_var_names = cond_event_haz_adj_var_names,
    cond_censoring_haz_adj_var_names = cond_censoring_haz_adj_var_names,
    treatment_var_name = treatment_var_name,
    propensity_score_var_name = propensity_score_var_name,
    propensity_score_sl_fit = propensity_score_sl_fit,
    cond_event_haz_sl_fit = cond_event_haz_sl_fit,
    cond_censoring_haz_sl_fit = cond_censoring_haz_sl_fit
  )
  clean_long_control_tbl <- generate_long_counterfactural_tbl_fun(
    clean_long_tbl = clean_long_tbl,
    treatment_group = 0,
    propensity_score_adj_var_names = propensity_score_adj_var_names,
    cond_event_haz_adj_var_names = cond_event_haz_adj_var_names,
    cond_censoring_haz_adj_var_names = cond_censoring_haz_adj_var_names,
    treatment_var_name = treatment_var_name,
    propensity_score_var_name = propensity_score_var_name,
    propensity_score_sl_fit = propensity_score_sl_fit,
    cond_event_haz_sl_fit = cond_event_haz_sl_fit,
    cond_censoring_haz_sl_fit = cond_censoring_haz_sl_fit
  )

  # compute the uncentered marginal survival EIF values until the max time
  full_uncentered_marginal_survival_eif_treatment_tbl <-
    full_uncentered_marginal_survival_eif_tbl_fun(
      counterfactual_long_tbl = clean_long_treatment_tbl,
      treatment_group = 1,
      treatment_var_name = treatment_var_name,
      time_cutoff = max_time_cutoff
    )
  full_uncentered_marginal_survival_eif_control_tbl <-
    full_uncentered_marginal_survival_eif_tbl_fun(
      counterfactual_long_tbl = clean_long_treatment_tbl,
      treatment_group = 0,
      treatment_var_name = treatment_var_name,
      time_cutoff = max_time_cutoff
    )

  # compute the range of possible truncation times based on the provided range
  # of effects and the selected homogeneity model
  if (homogeneity_type == "additive") {
    time_range_vec <- max_time_cutoff -
      seq(effect_range_vec[1], effect_range_vec[2])
  } else if (homogeneity_type == "multiplicative") {
    time_range_vec <- seq(
      floor(1 / effect_range_vec[2] * max_time_cutoff),
      ceiling(1 / effect_range_vec[1] * max_time_cutoff)
    )
  }

  # create a tibble to keep track of cutoffs for the secondary truncation
  trunc_tbl <- dplyr::tibble(time_range = time_range_vec) |>
    dplyr::mutate(
      trunc_time_treatment = dplyr::if_else(
        max_time_cutoff - time_range > 0, max_time_cutoff,
        2 * max_time_cutoff - time_range
      ),
      trunc_time_control = dplyr::if_else(
        max_time_cutoff - time_range > 0, time_range,
        max_time_cutoff
      )
    )

  # estimate the marginal tte variance under treatment at all the required times
  unique_trunc_time_treament_vec <- trunc_tbl |>
    dplyr::pull(trunc_time_treatment) |>
    unique()
  full_uncentered_marginal_survival_eif_treatment_tbl <-
    full_uncentered_marginal_survival_eif_tbl_fun(
      counterfactual_long_tbl = clean_long_treatment_tbl,
      treatment_group = 1,
      treatment_var_name = treatment_var_name,
      time_cutoff = max_time_cutoff
    )
  tte_var_treatment_tbl <- lapply(
    unique_trunc_time_treament_vec,
    function(trunc_time) {

      # estimate the RMST under treatment at trunc_time
      rmst_est_at_trunc_time <- one_step_rmst_estimator_fun(
        counterfactual_long_tbl = clean_long_treatment_tbl,
        treatment_group = 1,
        treatment_var_name = treatment_var_name,
        time_cutoff = trunc_time
      )

      # compute the uncentered marginal tte var EIFs
      uncentered_marginal_tte_var_eif_treatment_tbl <-
        uncentered_marginal_tte_var_eif_fun(
          full_uncentered_marginal_survival_eif_treatment_tbl,
          time_cutoff = trunc_time,
          rmst_est = rmst_est_at_trunc_time
        )

      # compute the plug-in estimate
      plugin_est <- tte_var_plugin_estimator_fun(
        counterfactual_long_tbl = clean_long_treatment_tbl,
        time_cutoff = trunc_time
      )

      # compute the one-step estimate
      one_step_est <- one_step_marginal_tte_var_estimator_fun(
        uncentered_marginal_tte_var_eif_treatment_tbl
      )

      # return the estimate and uncentered eif
      dplyr::tibble(
        trunc_time_treatment = trunc_time,
        one_step_est_treatment = one_step_est,
        eif_treatment =
          list(
            uncentered_marginal_tte_var_eif_treatment_tbl$uncentered_eif_val -
              plugin_est
          )
      )
    }
  ) |>
    dplyr::bind_rows()

  # estimate the marginal tte variance under control at all the required times
  unique_trunc_time_control_vec <- trunc_tbl |>
    dplyr::pull(trunc_time_control) |>
    unique()
  full_uncentered_marginal_survival_eif_control_tbl <-
    full_uncentered_marginal_survival_eif_tbl_fun(
      counterfactual_long_tbl = clean_long_control_tbl,
      treatment_group = 0,
      treatment_var_name = treatment_var_name,
      time_cutoff = max_time_cutoff
    )
  tte_var_control_tbl <- lapply(
    unique_trunc_time_control_vec,
    function(trunc_time) {

      # estimate the RMST under control at trunc_time
      rmst_est_at_trunc_time <- one_step_rmst_estimator_fun(
        counterfactual_long_tbl = clean_long_control_tbl,
        treatment_group = 0,
        treatment_var_name = treatment_var_name,
        time_cutoff = trunc_time
      )

      # compute the uncentered marginal tte var EIFs
      uncentered_marginal_tte_var_eif_control_tbl <-
        uncentered_marginal_tte_var_eif_fun(
          full_uncentered_marginal_survival_eif_control_tbl,
          time_cutoff = trunc_time,
          rmst_est = rmst_est_at_trunc_time
        )

      # compute the plug-in estimate
      plugin_est <- tte_var_plugin_estimator_fun(
        counterfactual_long_tbl = clean_long_control_tbl,
        time_cutoff = trunc_time
      )

      # compute the one-step estimate
      one_step_est <- one_step_marginal_tte_var_estimator_fun(
        uncentered_marginal_tte_var_eif_control_tbl
      )

      # return the estimate and uncentered eif
      dplyr::tibble(
        trunc_time_control = trunc_time,
        one_step_est_control = one_step_est,
        eif_control =
          list(
            uncentered_marginal_tte_var_eif_control_tbl$uncentered_eif_val -
              plugin_est
          )
      )
    }
  ) |>
    dplyr::bind_rows()

  # add the marginal variance estimates and uncentered EIFs to the trunc_tbl
  estimates_tbl <- trunc_tbl |>
    dplyr::left_join(tte_var_treatment_tbl, by = "trunc_time_treatment") |>
    dplyr::left_join(tte_var_control_tbl, by = "trunc_time_control")

  # compute the differential variance estimates, standard errors,
  # test statistics, and p-values
  num_obs <- clean_long_tbl |>
    dplyr::pull(cmldiffvar_id) |>
    unique() |>
    length()
  inference_tbl <- estimates_tbl |>
    dplyr::mutate(
      diff_var_est = one_step_est_treatment / one_step_est_control,
      eif_diff_var = purrr::pmap(
        list(
          eif_treatment, one_step_est_treatment, eif_control,
          one_step_est_control
        ),
        ~ ..1 / ..4 - ..3 * ..2 / ..4^2
      ),
      var_eif_diff_var = purrr::map_dbl(eif_diff_var, var),
      se_diff_var = sqrt(var_eif_diff_var / num_obs),
      test_statistic = dplyr::case_when(
        homogeneity_type == "additive" ~ (diff_var_est - 1) / se_diff_var,
        homogeneity_type == "multiplicative" ~
          (diff_var_est - (max_time_cutoff / time_range)^2) / se_diff_var
      ),
      p_value = 2 * pmin(
        pnorm(test_statistic, lower.tail = TRUE),
        pnorm(test_statistic, lower.tail = FALSE)
      ),
      hypothesized_homogeneous_effect = dplyr::case_when(
        homogeneity_type == "additive" ~ max_time_cutoff - time_range,
        homogeneity_type == "additive" ~ max_time_cutoff / time_range
      )
    ) |>
    dplyr::select(
      hypothesized_homogeneous_effect,
      one_step_est_treatment,
      one_step_est_control,
      diff_var_est,
      se_diff_var,
      test_statistic,
      p_value
    ) |>
    dplyr::arrange(hypothesized_homogeneous_effect)


  return(inference_tbl)
}
