#' One-Step Estimator of Treatment-Group-Specific Variance
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
#'   expected outcome conditional on confounders and treatment equal to
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

  os_var_est <- mean(
    (treatment_vec == treatment_group) *
      (outcome_vec^2 - 2 * outcome_vec * mean_est +
        2 * cond_exp_outcome_est_vec * mean_est - cond_exp_sq_outcome_est_vec) /
      ((treatment_group == 1) * ps_est_vec +
        (treatment_group == 0) * (1 - ps_est_vec)) +
    cond_exp_sq_outcome_est_vec - 2 * cond_exp_outcome_est_vec * mean_est +
      mean_est^2
  )

}
