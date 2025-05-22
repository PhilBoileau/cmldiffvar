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
