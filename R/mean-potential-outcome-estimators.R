#' One-Step Group-Specific Mean Estimator
#'
#' `one_step_mean_estimator_fun()` estimates the group-specific mean using a
#' one-step estimator.
#'
#' @inheritParams one_step_var_estimator_fun
#'
#' @returns A `numeric` estimate of the `treatment_group`-specific mean.
#'
#' @keywords internal
#'
one_step_mean_estimator_fun <- function(
    treatment_group,
    treatment_vec,
    outcome_vec,
    ps_est_vec,
    cond_exp_outcome_est_vec
) {

  # uncentered efficient influence function of the potential outcome mean
  uncentered_eif_vec <-
    (treatment_vec == treatment_group) *
      (outcome_vec - cond_exp_outcome_est_vec) /
        ((treatment_group == 1) * ps_est_vec +
           (treatment_group == 0) * (1 - ps_est_vec)) +
    cond_exp_outcome_est_vec

  # one-step estimate
  mean_est <- mean(uncentered_eif_vec)

  return(mean_est)
}


#' Uncentered Efficient Influence Function of the Marginal Survival Probability
#'
#' @param treatment_group A `numeric` indicating the counterfactual treatment
#'   group assignment.
#' @param cmldiffvar_id A `numeric` vector corresponding to the observations'
#'   assigned identifiers.
#' @param cmldiffvar_time_long A `numeric` vector corresponding to the
#'   observations' longitudinal times.
#' @param time_cutoff A `numeric` indicating the time at which to evaluate the
#'   marginal survival probability.
#' @param treatment_vec A `numeric` vector corresponding to the observations'
#'   assigned treatment groups.
#' @param ps_est_vec A `numeric` vector of the observations' predicted
#'   propensity scores.
#' @param cond_survival_est_vec A `numeric` vector of the observations'
#'   predicted conditional event survival probabilities.
#' @param cond_surv_at_trunc_est_vec A `numeric` vector of the observations'
#'   predicted conditional event survival probabilities at `time_cutoff`.
#' @param cond_censoring_survival_est_vec A `numeric` vector of the
#'   observations' predicted conditional censoring survival probabilities.
#' @param I_vec A `numeric` vector of the I indicators.
#' @param L_vec A `numeric` vector of the J indicators.
#' @param cond_event_haz_est_vec A `numeric` vector of the observations'
#'   predicted conditional event hazards.
#'
#' @returns A [tibble] containing the uncentered efficient influence function of
#'   the marginal survival probabilities at time `time_cutoff`.
#'
#' @keywords internal
#'
uncentered_marginal_survival_eif_fun <- function(
  treatment_group,
  cmldiffvar_id,
  cmldiffvar_time_long,
  time_cutoff,
  treatment_vec,
  ps_est_vec,
  cond_survival_est_vec,
  cond_surv_at_trunc_est_vec,
  cond_censoring_survival_est_vec,
  I_vec,
  L_vec,
  cond_event_haz_est_vec
) {
  # uncentered efficient influence function integrand of the marginal survival
  # in long format, without multiplicative factor
  eif_integrand_vec <- (I_vec * (L_vec - cond_event_haz_est_vec)) /
    (cond_survival_est_vec * cond_censoring_survival_est_vec)

  # IPWs
  ipws <- (treatment_vec == treatment_group) /
    ((treatment_group == 1) * ps_est_vec +
       (treatment_group == 0) * (1 - ps_est_vec))

  # collapse the uncentered efficient influence function by id
  uncentered_eif_tbl <- dplyr::tibble(
    cmldiffvar_id = cmldiffvar_id,
    cmldiffvar_time_long = cmldiffvar_time_long,
    cond_surv_at_trunc_est_vec = cond_surv_at_trunc_est_vec,
    ipws = ipws,
    eif_integrand_vec = eif_integrand_vec
  ) |>
    dplyr::group_by(.data$cmldiffvar_id) |>
    dplyr::mutate(
      int_weight = dplyr::lead(
        .data$cmldiffvar_time_long, default = time_cutoff
      ) - .data$cmldiffvar_time_long
    ) |>
    dplyr::summarize(
      uncentered_eif_integral_val = sum(
        .data$int_weight * .data$eif_integrand_vec
      ),
      # just to retain one value for final addition at end of calculation
      cond_surv_at_trunc_est_vec = min(.data$cond_surv_at_trunc_est_vec),
      ipws = min(.data$ipws),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      uncentered_eif_val = .data$ipws * .data$cond_surv_at_trunc_est_vec *
        .data$uncentered_eif_integral_val + .data$cond_surv_at_trunc_est_vec
    ) |>
    dplyr::select(.data$cmldiffvar_id, .data$uncentered_eif_val)

  return(uncentered_eif_tbl)
}

one_step_marginal_survival_estimator_fun <- function(uncentered_eif_tbl) {

  # compute the one-step marginal survival estimate
  mean(uncentered_eif_tbl$uncentered_eif_val)

}


#' One-Step Estimator of Restricted Mean Survival Time
#'
#' @inheritParams generate_long_counterfactural_tbl_fun
#'
#' @returns A `numeric` corresponding to the one-step estimate of the restricted
#'   mean survival time at the specified `time_cutoff`.
#'
#' @keywords internal
#'
one_step_rmst_estimator_fun <- function(
    counterfactual_long_tbl,
    treatment_group,
    treatment_var_name,
    time_cutoff
) {

  # only consider observations until time cutoff
  counterfactual_long_tbl <- counterfactual_long_tbl |>
    dplyr::filter(.data$cmldiffvar_long_time <= time_cutoff)

  # extract all of the unique observed times
  unique_times <- counterfactual_long_tbl |>
    dplyr::pull(.data$cmldiffvar_long_time) |>
    unique() |>
    sort()

  # compute the marginal survival estimates at all observed times until the
  # cutoff
  marginal_survival_estimates_vec <- sapply(
    unique_times,
    function(iter_time) {

      # extract the conditional survival times truncated at the iter_time
      num_times <- sum(unique_times <= iter_time)
      cond_surv_at_trunc_est_vec <- counterfactual_long_tbl |>
        dplyr::filter(.data$cmldiffvar_long_time <= iter_time) |>
        dplyr::group_by(.data$cmldiffvar_id) |>
        dplyr::slice_tail(n = 1) |>
        dplyr::select(.data$cmldiffvar_id, .data$pred_cond_event_survival) |>
        dplyr::ungroup() |>
        dplyr::slice(rep(1:dplyr::n(), each = num_times)) |>
        dplyr::pull(.data$pred_cond_event_survival)

      # restrict the longitudinal counterfactual dataset to the timeframe
      trunc_counterfactual_long_tbl <- counterfactual_long_tbl |>
        dplyr::filter(.data$cmldiffvar_long_time <= iter_time)

      # one-step estimator of the marginal survival probabilities
      uncentered_eif_tbl <- uncentered_marginal_survival_eif_fun(
        treatment_group = treatment_group,
        cmldiffvar_id = trunc_counterfactual_long_tbl$cmldiffvar_id,
        cmldiffvar_time_long =
          trunc_counterfactual_long_tbl$cmldiffvar_long_time,
        time_cutoff = iter_time,
        treatment_vec = trunc_counterfactual_long_tbl[[treatment_var_name]],
        ps_est_vec = trunc_counterfactual_long_tbl$pred_propensity_score,
        cond_survival_est_vec =
          trunc_counterfactual_long_tbl$pred_cond_event_survival,
        cond_surv_at_trunc_est_vec = cond_surv_at_trunc_est_vec,
        cond_censoring_survival_est_vec =
          trunc_counterfactual_long_tbl$pred_cond_censoring_survival,
        I_vec = trunc_counterfactual_long_tbl$I,
        L_vec = trunc_counterfactual_long_tbl$L,
        cond_event_haz_est_vec =
          trunc_counterfactual_long_tbl$pred_cond_event_haz
      )
      one_step_marginal_survival_estimator_fun(uncentered_eif_tbl)

    }
  )

  # compute a weighted sum of the one-step marginal survival estimates
  rmst_est <- dplyr::tibble(
    cmldiffvar_time_long = unique_times,
    marginal_survival_estimates_vec = marginal_survival_estimates_vec
  ) |>
    dplyr::mutate(
      int_weight = dplyr::lead(
        .data$cmldiffvar_time_long, default = time_cutoff
      ) - .data$cmldiffvar_time_long
    ) |>
    # NOTE: add min(cmldiffvar_time_long) to include survival times between
    # time zero and first recorded event
    dplyr::summarize(
      rmst_est = sum(.data$int_weight * .data$marginal_survival_estimates_vec) +
        min(.data$cmldiffvar_time_long)
    ) |>
    dplyr::pull(.data$rmst_est)

  return(rmst_est)
}
