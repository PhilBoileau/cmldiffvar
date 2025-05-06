#' One-Step Estimator of Group-Specific Variance
#'
#' @description `one_step_var_estimator_fun()` estimates of the
#'   treatment-group-specific variance using a one-step estimator.
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

  mean(
    (treatment_vec == treatment_group) *
      (outcome_vec^2 - 2 * outcome_vec * mean_est +
        2 * cond_exp_outcome_est_vec * mean_est - cond_exp_sq_outcome_est_vec) /
      ((treatment_group == 1) * ps_est_vec +
        (treatment_group == 0) * (1 - ps_est_vec)) +
    cond_exp_sq_outcome_est_vec - 2 * cond_exp_outcome_est_vec * mean_est +
      mean_est^2
  )

}


#' Targeted Maximum Likelihood Estimator of Group-Specific Variance
#'
#' @description `tml_var_estimator_fun()` estimates of the
#'   treatment-group-specific variance using a targeted maximum likelihood
#'   estimator.
#'
#' @inheritParams one_step_var_estimator_fun
#'
#' @returns A `numeric` estimate of the `treatment_group`-specific variance.
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

  # tilt the conditional expected outcome estimator
  clever_covariate <- (treatment_vec == treatment_group) /
    ((treatment_group == 1) * ps_est_vec +
       (treatment_group == 0) * (1 - ps_est_vec))
  tilted_bounded_cond_exp_outcome_fit <- suppressWarnings(
    glm(
      bounded_outcome_vec ~ -1 + clever_covariate,
      family = "binomial",
      offset = stats::qlogis(bounded_cond_exp_outcome_est_vec)
    )
  )
  tilted_bounded_cond_exp_outcome_vec <- predict(
    tilted_bounded_cond_exp_outcome_fit,
    type = "response"
  )
  tilted_cond_exp_outcome_vec <-
    tilted_bounded_cond_exp_outcome_vec * (max_outcome_vec - min_outcome_vec) +
      min_outcome_vec

  # tilt the conditional expected squared outcome estimator
  # NOTE: add epsilon to avoid qlogis of zero error
  tilted_bounded_sq_cond_exp_outcome_fit <- suppressWarnings(
    glm(
      bounded_sq_outcome_vec ~ -1 + clever_covariate,
      family = "binomial",
      offset = stats::qlogis(
        bounded_cond_exp_sq_outcome_est_vec + .Machine$double.eps
      )
    )
  )
  tilted_bounded_cond_exp_sq_outcome_vec <- predict(
    tilted_bounded_sq_cond_exp_outcome_fit,
    type = "response"
  )
  tilted_cond_exp_sq_outcome_vec <-
    tilted_bounded_cond_exp_sq_outcome_vec *
      (max_sq_outcome_vec - min_sq_outcome_vec) + min_sq_outcome_vec

  # compute the targeted maximum likelihood estimate
  mean(tilted_cond_exp_sq_outcome_vec) - mean(tilted_cond_exp_outcome_vec)^2

}
