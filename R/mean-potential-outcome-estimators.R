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


one_step_marginal_survival_estimator_fun <- function(
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
  eif_tbl <- dplyr::tibble(
    cmldiffvar_id = cmldiffvar_id,
    cmldiffvar_time_long = cmldiffvar_time_long,
    cond_surv_at_trunc_est_vec = cond_surv_at_trunc_est_vec,
    ipws = ipws,
    eif_integrand_vec = eif_integrand_vec
  ) |>
    dplyr::group_by(cmldiffvar_id) |>
    dplyr::mutate(
      int_weight = dplyr::lead(
        cmldiffvar_time_long, default = time_cutoff
      ) - cmldiffvar_time_long
    ) |>
    dplyr::summarize(
      uncentered_eif_integral_val = sum(int_weight * eif_integrand_vec),
      # just to retain one value for final addition at end of calculation
      cond_surv_at_trunc_est_vec = min(cond_surv_at_trunc_est_vec),
      ipws = min(ipws),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      uncentered_eif_val = ipws * cond_surv_at_trunc_est_vec *
        uncentered_eif_integral_val + cond_surv_at_trunc_est_vec
    )

  # compute the one-step marginal survival estimate
  mean(eif_tbl$uncentered_eif_val)

}

