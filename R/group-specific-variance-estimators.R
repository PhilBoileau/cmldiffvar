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
#' @returns A `numeric` estimate of the `treatment_group`-specific variance.
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

  return(os_est)
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
  clever_covariate <- (treatment_vec == treatment_group) /
    ((treatment_group == 1) * ps_est_vec +
       (treatment_group == 0) * (1 - ps_est_vec))

  # tilt the conditional expected outcome estimator
  # NOTE: Make sure offset terms aren't NaN
  logit_bounded_cond_exp_outcome_est_vec <-stats::qlogis(
    bound_away_from_0_and_1_fun(bounded_cond_exp_outcome_est_vec)
  )
  tilted_bounded_cond_exp_outcome_fit <- suppressWarnings(
    glm(
      bounded_outcome_vec ~ -1 + clever_covariate,
      family = "binomial",
      offset = logit_bounded_cond_exp_outcome_est_vec
    )
  )
  tilted_bounded_cond_exp_outcome_est_vec <- predict(
    tilted_bounded_cond_exp_outcome_fit,
    type = "response"
  )
  tilted_cond_exp_outcome_est_vec <-
    tilted_bounded_cond_exp_outcome_est_vec *
      (max_outcome_vec - min_outcome_vec) + min_outcome_vec

  # tilt the conditional expected squared outcome estimator
  # NOTE: Make sure offset terms aren't NaN
  logit_bounded_cond_exp_sq_outcome_est_vec <-stats::qlogis(
    bound_away_from_0_and_1_fun(bounded_cond_exp_sq_outcome_est_vec)
  )
  tilted_bounded_cond_exp_sq_outcome_fit <- suppressWarnings(
    glm(
      bounded_sq_outcome_vec ~ -1 + clever_covariate,
      family = "binomial",
      offset = logit_bounded_cond_exp_sq_outcome_est_vec
    )
  )
  tilted_bounded_cond_exp_sq_outcome_est_vec <- predict(
    tilted_bounded_cond_exp_sq_outcome_fit,
    type = "response"
  )
  tilted_cond_exp_sq_outcome_est_vec <-
    tilted_bounded_cond_exp_sq_outcome_est_vec *
      (max_sq_outcome_vec - min_sq_outcome_vec) + min_sq_outcome_vec

  # compute the targeted maximum likelihood estimate
  tml_est <- plugin_var_estimator_fun(
    cond_exp_outcome_est_vec = tilted_cond_exp_outcome_est_vec,
    cond_exp_sq_outcome_est_vec = tilted_cond_exp_sq_outcome_est_vec
  )

  return(tml_est)

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
  epsilon <- .Machine$double.eps

  # shift values in unit_interval_vec away from 0 and 1 by epsilon
  unit_interval_vec[unit_interval_vec < epsilon] <- epsilon
  unit_interval_vec[unit_interval_vec > 1 - epsilon] <- 1 - epsilon

  return(unit_interval_vec)
}
