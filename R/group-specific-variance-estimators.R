#' One-Step Estimator of Group-Specific Variance
#'
#' @description `one_step_var_estimator_fun()` estimates of the group-specific
#'   variance using a one-step estimator.
#'
#' @param treatment_group A binary `numeric` indicating the treatment group for
#'   which the variance is being estimated.
#' @param treatment_vec A binary `numeric` vector reporting observations'
#'   observed treatment groups.
#' @param outcome_vec A `numeric` vector reporting observations' observed
#'   outcomes.
#' @param ps_est_vec A `numeric` vector bounded to the unit interval
#'   corresponding to the predicted probability of being assigned to the
#'   treatment group.
#' @param cond_exp_outcome_est_vec A `numeric` corresponding to the predicted
#'   expected outcome conditional on confounders and treatment equal to
#'   `treatment_group`.
#' @param cond_exp_sq_outcome_est_vec A `numeric` corresponding to the predicted
#'   expected squared outcome conditional on confounders and treatment equal to
#'   `treatment_group`.
#' @param mean_est A `numeric` corresponding to the estimated outcome population
#'   mean of the designated `treatment_group`.
#'
#' @returns A list containing the `numeric` estimate of the
#'   `treatment_group`-specific marginal outcome variance and the `numeric`
#'   efficient influence function vector.
#'
#' @keywords internal
#'
one_step_var_estimator_fun <- function(
  treatment_group,
  treatment_vec,
  outcome_vec,
  ps_est_vec,
  cond_exp_outcome_est_vec,
  cond_exp_sq_outcome_est_vec,
  mean_est
) {

  # compute the plug-in estimate
  plugin_var_est <- plugin_var_estimator_fun(
    cond_exp_outcome_est_vec,
    cond_exp_sq_outcome_est_vec
  )

  # compute the eif
  eif_vec <- efficient_influence_function_fun(
    treatment_group,
    treatment_vec,
    outcome_vec,
    ps_est_vec,
    cond_exp_outcome_est_vec,
    cond_exp_sq_outcome_est_vec,
    mean_est
  )

  # compute the one-step estimate
  os_est <- mean(eif_vec) + plugin_var_est

  # return the estimate and the EIF
  result_ls <- list(
    estimate = os_est,
    eif = eif_vec
  )

  return(result_ls)
}


#' Targeted Minimum Loss-Based Estimator of Group-Specific Variance
#'
#' @description `tml_var_estimator_fun()` estimates of the group-specific
#'   variance using a targeted minimum loss-based estimator.
#'
#' @inheritParams one_step_var_estimator_fun
#'
#' @inherit one_step_var_estimator_fun return
#'
#' @keywords internal
#'
tml_var_estimator_fun <- function(
  treatment_group,
  treatment_vec,
  outcome_vec,
  ps_est_vec,
  cond_exp_outcome_est_vec,
  cond_exp_sq_outcome_est_vec
) {

  # bound the outcome vector
  min_outcome_vec <- min(outcome_vec)
  max_outcome_vec <- max(outcome_vec)
  bounded_outcome_vec <-
    (outcome_vec - min_outcome_vec) / (max_outcome_vec - min_outcome_vec)

  # bound the conditional expected outcome estimates, and don't let estimator
  # extrapolate outside bounds of observed outcome values
  cond_exp_outcome_est_vec[cond_exp_outcome_est_vec < min_outcome_vec] <-
    min_outcome_vec
  cond_exp_outcome_est_vec[cond_exp_outcome_est_vec > max_outcome_vec] <-
    max_outcome_vec
  bounded_cond_exp_outcome_est_vec <-
    (cond_exp_outcome_est_vec - min_outcome_vec) /
    (max_outcome_vec - min_outcome_vec)

  # bound the conditional expected squared outcome estimates, and don't let
  # estimator extrapolate outside bounds of observed outcome values
  sq_outcome_vec <- outcome_vec^2
  min_sq_outcome_vec <- min(sq_outcome_vec)
  max_sq_outcome_vec <- max(sq_outcome_vec)
  bounded_sq_outcome_vec <- (sq_outcome_vec - min_sq_outcome_vec) /
    (max_sq_outcome_vec - min_sq_outcome_vec)
  cond_exp_sq_outcome_est_vec[
    cond_exp_sq_outcome_est_vec < min_sq_outcome_vec
  ] <- min_sq_outcome_vec
  cond_exp_sq_outcome_est_vec[
    cond_exp_sq_outcome_est_vec > max_sq_outcome_vec
  ] <- max_sq_outcome_vec
  bounded_cond_exp_sq_outcome_est_vec <-
    (cond_exp_sq_outcome_est_vec - min_sq_outcome_vec) /
    (max_sq_outcome_vec - min_sq_outcome_vec)

  # define the clever covariate
  clever_covariate_num <- as.numeric(treatment_vec == treatment_group)
  clever_covariate_denom <- 1 /
    ((treatment_group == 1) * ps_est_vec +
       (treatment_group == 0) * (1 - ps_est_vec))

  # tilt the conditional expected outcome estimator
  # NOTE: Make sure offset terms aren't NaN
  logit_bounded_cond_exp_outcome_est_vec <- stats::qlogis(
    bound_away_from_0_and_1_fun(bounded_cond_exp_outcome_est_vec)
  )
  tilted_bounded_cond_exp_outcome_fit <- suppressWarnings(
    stats::glm(
      bounded_outcome_vec ~ -1 + clever_covariate_num,
      family = "binomial",
      offset = logit_bounded_cond_exp_outcome_est_vec,
      weights = clever_covariate_denom
    )
  )
  if (is.na(stats::coef(tilted_bounded_cond_exp_outcome_fit))) {
    tilted_bounded_cond_exp_outcome_fit$coefficients <- 0
  }
  epsilon_cond_exp_outcome <- unname(
    stats::coef(tilted_bounded_cond_exp_outcome_fit)
  )
  tilted_bounded_cond_exp_outcome_est_vec <- stats::plogis(
    logit_bounded_cond_exp_outcome_est_vec + epsilon_cond_exp_outcome
  )
  tilted_cond_exp_outcome_est_vec <-
    tilted_bounded_cond_exp_outcome_est_vec *
    (max_outcome_vec - min_outcome_vec) + min_outcome_vec

  # tilt the conditional expected squared outcome estimator
  # NOTE: Make sure offset terms aren't NaN
  logit_bounded_cond_exp_sq_outcome_est_vec <- stats::qlogis(
    bound_away_from_0_and_1_fun(bounded_cond_exp_sq_outcome_est_vec)
  )
  tilted_bounded_cond_exp_sq_outcome_fit <- suppressWarnings(
    stats::glm(
      bounded_sq_outcome_vec ~ -1 + clever_covariate_num,
      family = "binomial",
      offset = logit_bounded_cond_exp_sq_outcome_est_vec,
      weights = clever_covariate_denom
    )
  )
  if (is.na(stats::coef(tilted_bounded_cond_exp_sq_outcome_fit))) {
    tilted_bounded_cond_exp_sq_outcome_fit$coefficients <- 0
  }
  epsilon_cond_exp_outcome_sq <- unname(
    stats::coef(tilted_bounded_cond_exp_sq_outcome_fit)
  )
  tilted_bounded_cond_exp_sq_outcome_est_vec <- stats::plogis(
    logit_bounded_cond_exp_sq_outcome_est_vec + epsilon_cond_exp_outcome_sq
  )
  tilted_cond_exp_sq_outcome_est_vec <-
    tilted_bounded_cond_exp_sq_outcome_est_vec *
    (max_sq_outcome_vec - min_sq_outcome_vec) + min_sq_outcome_vec

  # compute the the TMLE
  tml_est <- plugin_var_estimator_fun(
    tilted_cond_exp_outcome_est_vec, tilted_cond_exp_sq_outcome_est_vec
  )

  # compute the EIF
  eif_vec <- efficient_influence_function_fun(
    treatment_group,
    treatment_vec,
    outcome_vec,
    ps_est_vec,
    tilted_cond_exp_outcome_est_vec,
    tilted_cond_exp_sq_outcome_est_vec,
    mean(tilted_cond_exp_outcome_est_vec)
  )

  # return the estimate and the EIF
  result_ls <- list(
    estimate = tml_est,
    eif = eif_vec
  )

  return(result_ls)

}


#' Efficient Influence Function of Group-Specific Variance
#'
#' @description `efficient_influence_function_fun()` computes the efficient
#'   influence function based on the observed data and nuisance parameter
#'   estimates for the specified treatment group.
#'
#' @inheritParams one_step_var_estimator_fun
#'
#' @returns A `numeric` vector corresponding to the efficient influence
#'   function.
#'
#' @keywords internal
efficient_influence_function_fun <- function(
  treatment_group,
  treatment_vec,
  outcome_vec,
  ps_est_vec,
  cond_exp_outcome_est_vec,
  cond_exp_sq_outcome_est_vec,
  mean_est
) {

  # compute the plug-in estimate
  plugin_var_est <- plugin_var_estimator_fun(
    cond_exp_outcome_est_vec,
    cond_exp_sq_outcome_est_vec
  )

  # compute the efficient influence function
  eif_vec <- (treatment_vec == treatment_group) *
    (outcome_vec^2 - 2 * outcome_vec * mean_est +
       2 * cond_exp_outcome_est_vec * mean_est - cond_exp_sq_outcome_est_vec) /
    ((treatment_group == 1) * ps_est_vec +
       (treatment_group == 0) * (1 - ps_est_vec)) +
    cond_exp_sq_outcome_est_vec - 2 * cond_exp_outcome_est_vec * mean_est +
    mean_est^2 - plugin_var_est

  return(eif_vec)

}


#' Plug-In Estimator of the Group-Specific Variance
#'
#' @description `plugin_var_estimator_fun()` estimates of the group-specific
#'   variance using a plug-in estimator.
#'
#' @inheritParams one_step_var_estimator_fun
#'
#' @inherit one_step_var_estimator_fun return
#'
#' @keywords internal
#'
plugin_var_estimator_fun <- function(
  cond_exp_outcome_est_vec,
  cond_exp_sq_outcome_est_vec
) {

  # definition of plug-in estimator
  mean(cond_exp_sq_outcome_est_vec) - mean(cond_exp_outcome_est_vec)^2
}


#' Bound Vectors on Unit Interval Away from Zero and One
#'
#' @param unit_interval_vec A numeric interval with values between zero and one.
#'
#' @returns A `numeric` vector with values bounded away from zero and one by
#'  some small value.
#'
#' @keywords internal
bound_away_from_0_and_1_fun <- function(unit_interval_vec) {

  # define the epsilon
  epsilon <- 2 * .Machine$double.eps

  # shift values in unit_interval_vec away from 0 and 1 by epsilon
  unit_interval_vec[unit_interval_vec < epsilon] <- epsilon
  unit_interval_vec[unit_interval_vec > 1 - epsilon] <- 1 - epsilon

  return(unit_interval_vec)
}


#' Unadjsuted Group-Specific Variance Estimator
#'
#' @description `unadjusted_var_estimator_fun()` estimates of the group-specific
#'   variance using an unadjusted plug-in estimator. Note that this estimator is
#'   regular and asymptotically linear in a nonparametric model when treatment
#'   assignment is randomized. The efficient influence function returned by this
#'   estimator corresponds to that of the complete data.
#'
#' @inheritParams one_step_var_estimator_fun
#'
#' @inherit one_step_var_estimator_fun return
#'
#' @keywords internal
unadjusted_var_estimator_fun <- function(
  treatment_group,
  treatment_vec,
  outcome_vec,
  ps_est_vec
) {

  # compute the estimate
  outcome_vec_subset <- outcome_vec[which(treatment_vec == treatment_group)]
  group_mean_est <- mean(outcome_vec_subset)
  group_var_est <- mean((outcome_vec_subset - group_mean_est)^2)

  # compute the efficient influence function for the full data
  eif_vec <- (treatment_vec == treatment_group) /
    ((treatment_group == 1) * ps_est_vec +
       (treatment_group == 0) * (1 - ps_est_vec)) *
    ((outcome_vec - group_mean_est)^2 - group_var_est)

  # return the estimate and the EIF
  result_ls <- list(
    estimate = group_var_est,
    eif = eif_vec
  )

  return(result_ls)

}


one_step_tte_var_estimator_fun <- function(
    counterfactual_long_tbl,
    treatment_group,
    treatment_var_name,
    time_cutoff
) {

  # only consider observations until time cutoff
  counterfactual_long_tbl <- counterfactual_long_tbl |>
    dplyr::filter(cmldiffvar_long_time <= time_cutoff)

  # extract all of the unique observed times
  unique_times <- counterfactual_long_tbl |>
    dplyr::pull(cmldiffvar_long_time) |>
    unique() |>
    sort()

  # compute the uncentered marginal survival EIF values at all observed times
  # until the cutoff
  uncentered_eif_marginal_survival_values_tbl <- lapply(
    unique_times,
    function(iter_time) {

      # extract the conditional survival times truncated at the iter_time
      num_times <- sum(unique_times <= iter_time)
      cond_surv_at_trunc_est_vec <- counterfactual_long_tbl |>
        dplyr::filter(cmldiffvar_long_time <= iter_time) |>
        dplyr::group_by(cmldiffvar_id) |>
        dplyr::slice_tail(n = 1) |>
        dplyr::select(cmldiffvar_id, pred_cond_event_survival) |>
        dplyr::ungroup() |>
        dplyr::slice(rep(1:dplyr::n(), each = num_times)) |>
        dplyr::pull(pred_cond_event_survival)

      # restrict the longitudinal counterfactual dataset to the time frame
      trunc_counterfactual_long_tbl <- counterfactual_long_tbl |>
        dplyr::filter(cmldiffvar_long_time <= iter_time)

      # compute the uncentered marginal survival EIF values
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

      # add the time point
      uncentered_eif_tbl |>
        dplyr::mutate(cmldiffvar_long_time = iter_time)

    }
  ) |>
    dplyr::bind_rows()

  # compute the weighted sum of the uncentered marginal survival EIF values
  # which equals the integral term in the EIF of the marginal TTE variance
  eif_integral_term_tbl <- uncentered_eif_marginal_survival_values_tbl |>
    dplyr::arrange(cmldiffvar_id, cmldiffvar_long_time) |>
    dplyr::group_by(cmldiffvar_id) |>
    dplyr::mutate(
      partial_weight = (cmldiffvar_long_time - 1) *
        (cmldiffvar_long_time - 2) / 2,
      lead_partial_weight = dplyr::lead(
        partial_weight, default = time_cutoff * (time_cutoff - 1) / 2
      ),
      int_weight = lead_partial_weight - partial_weight
    ) |>
    dplyr::summarize(
      integral_term = 2 * sum(int_weight * uncentered_eif_val),
      .groups = "drop"
    )

  # estimate the RMST using one-step estimator
  rmst_est <- one_step_rmst_estimator_fun(
    counterfactual_long_tbl,
    treatment_group,
    treatment_var_name,
    time_cutoff
  )

  # construct the uncentered EIF tibble
  uncentered_eif_tbl <- eif_integral_term_tbl |>
    dplyr::mutate(
      uncentered_eif_val = integral_term + rmst_est * (2 - rmst_est)
    )

  # compute the one-step marginal variance estimate
  mean(uncentered_eif_tbl$uncentered_eif_val)


}
