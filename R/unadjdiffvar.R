#' Unadjusted Differential Variance Inference
#'
#' @description `unadjdiffvar()` infers the differential variance of two
#'   treatments using an unadjusted estimator.
#'
#' @details Under standard causal identifiability conditions --- namely,
#'   consistency, positivity, and full exchangeability --- `unadjdiffvar()`
#'   performs inference on the differential variance of the potential outcomes.
#'   Differential variance is defined on the absolute scale as the difference in
#'   potential outcome variances of two treatments. On the relative scale,
#'   differential variance is defined as the ratio of the potential outcome
#'   variances. The scale of the differential variance estimated by
#'   `unadjdiffvar()` is specified by the `estimand_type` parameter.
#'
#'   `unadjdiffvar()` assumes that the data in `data_tbl` are generated
#'   according to a parallel group study design with a binary treatment variable
#'   and a continuous outcome. Treatment is assumed to be assigned at random. If
#'   the treatment assignment probabilities are known, such as in a randomized
#'   controlled trial, then these may be passed to the estimator using the
#'   `propensity_score_var_name` argument. If there are (un)measured
#'   confounders, this estimator will be biased.
#'
#' @inheritParams cmldiffvar
#'
#' @inherit cmldiffvar return
#'
#' @export

unadjdiffvar <- function(
  data_tbl,
  estimand_type = "absolute",
  confidence_level = 0.95,
  treatment_var_name,
  propensity_score_var_name = NULL,
  outcome_var_name
) {

  # checks ----

  checkmate::assert_choice(estimand_type, c("absolute", "relative"))
  checkmate::assert_number(confidence_level, lower = 0.01, upper = 0.99)
  checkmate::assert_character(treatment_var_name)
  checkmate::assert_character(outcome_var_name)

  # check dataset contains appropriate variables
  if (is.null(propensity_score_var_name)) {
    clean_tbl_var_names <- c(treatment_var_name, outcome_var_name)
  } else {
    checkmate::assert_character(propensity_score_var_name)
    checkmate::assert_numeric(
      data_tbl[[propensity_score_var_name]], lower = 0.001, upper = 0.999
    )
    clean_tbl_var_names <- c(
      treatment_var_name, propensity_score_var_name, outcome_var_name
    )
  }
  checkmate::assert_names(clean_tbl_var_names, subset.of = colnames(data_tbl))

  # only retain relevant columns in data_tbl
  clean_tbl <- data_tbl |> dplyr::select(dplyr::all_of(clean_tbl_var_names))

  # check dataset contains only numeric variables
  checkmate::assert_data_frame(
    clean_tbl,
    types = c("numeric"),
    any.missing = FALSE,
    min.rows = 50
  )

  # compute estimates and EIF ----

  # compute the unadjusted estimates and EIF
  diff_var_inference <- unadjusted_diff_var_estimator_fun(
    clean_tbl = clean_tbl,
    treatment_var_name = treatment_var_name,
    propensity_score_var_name = propensity_score_var_name,
    outcome_var_name = outcome_var_name,
    estimand_type = estimand_type
  )
  estimate <- diff_var_inference$estimate
  eif <- diff_var_inference$eif

  # compute confidence interval ----

  # Wald-type confidence intervals relying on asymptotic linearity
  se <- as.numeric(sqrt(stats::var(eif) / nrow(clean_tbl)))
  critical_value <- stats::qnorm(1 - (1 - confidence_level) / 2)
  ci_low <- estimate - critical_value * se
  ci_high <- estimate + critical_value * se


  # compute p-value ----

  # p-value relying on asymptotic linearity
  if (estimand_type == "absolute") {
    z_score <- estimate / se
  } else if (estimand_type == "relative") {
    z_score <- (estimate - 1) / se
  }
  p_value <- 2 * min(
    stats::pnorm(z_score), stats::pnorm(z_score, lower.tail = FALSE)
  )

  # assemble and output results ----

  dplyr::tibble(
    estimand = paste(estimand_type, "differential variance"),
    estimate = estimate,
    se = se,
    ci_low = ci_low,
    ci_high = ci_high,
    p_value = p_value,
    eif = eif
  )

}
