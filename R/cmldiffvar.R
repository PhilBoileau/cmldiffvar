#' Causal Machine Learning for Differential Variance
#'
#' @description `cmldiffvar()` infers the differential variance between two
#'   treatments using causal machine learning methods.
#'
#' @details `cmldiffvar()` assumes that the data in `data_tbl` are generated
#'   according to a parallel group study design with a binary treatment
#'   variable, a continuous outcome, and continuous or binary treatment-outcome
#'   confounders.
#'
#'   Under standard causal identifiability conditions --- namely, consistency,
#'   positivity, and full exchangeability --- `cmldiffvar()` performs inference
#'   on the differential variance of the potential outcomes. Differential
#'   variance is defined on the absolute scale as the difference in potential
#'   outcome standard deviations of two treatments. On the relative scale,
#'   differential variance is defined as the ratio of the potential outcome
#'   variances. The scale of the differential variance estimated by
#'   `cmldiffvar()` is specified by the `estimand_type` parameter.
#'
#'   These differential variance estimands rely on three nuisance parameters:
#'   the propensity score, the expected outcome conditional on confounders and
#'   treatment assignment, and the expected squared outcome conditional on
#'   confounders and treatment assignment.
#'
#'   `cmldiffvar()` implements cross-fitted one-step and targeted maximum
#'   likelihood estimators of the differential variance estimands. These
#'   cross-fitted estimators are used by default; non-cross-fitted estimators
#'   can be used by setting `cross_fit = FALSE`. The number of folds used in
#'   these cross-fitting procedures is set by `num_cross_fit_folds`. The
#'   differential variance estimands' nuisance parameters are flexibly estimated
#'   with [SuperLearner][SuperLearner::SuperLearner] estimators.
#'
#'   The estimators implemented in `cmldiffvar()` are consistent if *at least*
#'   one of the following conditions is satisfied: (1) the propensity score is
#'   consistently estimated, and  (2) the expected outcome conditional on
#'   confounders and treatment assignment and the expected squared outcome
#'   conditional on confounders and treatment assignment are consistently
#'   estimated.
#'
#'   The estimators implemented in `cmldiffvar()` are asymptotically linear ---
#'   meaning their asymptotic sampling distribution is normally distributed
#'   about the true differential variance parameter --- if *all* nuisance
#'   parameter estimators are consistently estimated at a rate of
#'   \eqn{o_P(n^{-1/4})}. The confidence intervals, standard errors, and
#'   p-values reported by `cmldiffvar()` are incorrect if these rate conditions
#'   are not satisfied.
#'
#'   When the data is the product of a randomized study with known propensity
#'   scores, these propensity scores can be provided to the
#'   `propensity_score_var_name` parameter. The conditions required of the
#'   estimators to be consistent and asymptotic linear are automatically
#'   satisfied when know propensity scores are used.
#'
#'   The cross-fitted estimation procedures can be parallelized by setting
#'   `parallel = TRUE`. Parallelization relies on the [future][future::future]
#'   package. Instructions for setting up parallel processing are available in
#'   the [future][future::future] package's vignettes.
#'
#' @param data_tbl A [`data.frame`] or [`tibble`][tibble::tibble].
#' @param estimand_type A `character` indicating whether to estimate the
#'   absolute (`"absolute"`) or relative (`"relative"`) differential variance.
#' @param estimator_type A `character` indicating whether to use a one-step
#'   (`"one-step"`) or a targeted maximum likelihood estimator (`"tmle"`) in the
#'   cross-fitting procedure.
#' @param confidence_level A `numeric` between $0.1$ and $0.99$ providing the
#'   confidence level used to compute confidence intervals. Defaults to `0.95`.
#' @param propensity_score_adj_var_names A `character` vector providing the
#'   columns names of the adjustment set variables for propensity score
#'   estimation stored in `data_tbl`. Ignored if `propensity_score_var_name` is
#'   `NULL`.
#' @param cond_exp_outcome_adj_var_names A `character` vector providing the
#'   columns names of the adjustment set variables for conditional expected
#'   (squared) outcome estimation stored in `data_tbl`.
#' @param treatment_var_name A `character` providing the column name of the
#'   treatment assignment indicator stored in `data_tbl`.
#' @param propensity_score_var_name An optional `character` providing the column
#'   name of the treatment assignment indicator stored in `data_tbl`. Defaults
#'   to `NULL`. See the Details section for more information.
#' @param outcome_var_name A `character` providing the column name of the
#'   outcome variable stored in `data_tbl`.
#' @param propensity_score_library A `character` vector of candidate learners
#'   used by the SuperLearner estimator of the propensity score. Defaults to
#'   `c("SL.glm", "SL.earth", "SL.ranger")`.
#' @param cond_exp_outcome_library A `character` vector of candidate learners
#'   used by the SuperLearner estimator of the expected outcome conditional on
#'   confounders and treatment assignment. Defaults to `c("SL.glm",
#'   "SL.earth", "SL.ranger")`.
#' @param cond_exp_sq_outcome_library A `character` vector of candidate learners
#'   used by the SuperLearner estimator of the expected squared outcome
#'   conditional on confounders and treatment assignment. Defaults to
#'   `c("SL.glm.interaction", "SL.ranger")`.
#' @param num_nuisance_sl_folds A `numeric` indicating the number of folds to
#'   use in cross-validated SuperLearner estimators. Defaults to `5`.
#' @param cross_fit A `logical flag` determining whether cross-fitted estimators
#'   are used. Defaults to `TRUE`.
#' @param num_cross_fit_folds A `numeric` setting the number of folds to use by
#'   the cross-fitting procedures, if used. Defaults to `5`.
#' @param parallel A `logical flag` indicating whether the cross-fitting
#'   procedure, if used, should be parallelized. Defaults to `FALSE`.
#'
#' @returns A one-row [tibble][tibble::tibble] containing the following columns:
#'  - `estimand`: The scale of the differential variance estimand
#'  - `estimate`: The differential variance estimate
#'  - `se`: The estimator's standard error
#'  - `ci_low`: The lower bound of the Wald-type confidence interval
#'  - `ci_high`: The upper bound of the Wald-type confidence interval
#'  - `p_value`: The p-value of a two-sided test using the z-score of the
#'   differential variance estimate
#'
#' @export
#'
cmldiffvar <- function(
  data_tbl,
  estimand_type = "absolute",
  estimator_type = "tmle",
  confidence_level = 0.95,
  propensity_score_adj_var_names,
  cond_exp_outcome_adj_var_names,
  treatment_var_name,
  propensity_score_var_name = NULL,
  outcome_var_name,
  propensity_score_library = c("SL.glm", "SL.earth", "SL.ranger"),
  cond_exp_outcome_library = c("SL.glm", "SL.earth", "SL.ranger"),
  cond_exp_sq_outcome_library = c("SL.glm.interaction", "SL.ranger"),
  num_nuisance_sl_folds = 5,
  cross_fit = FALSE,
  num_cross_fit_folds = 5,
  parallel = FALSE
) {

  # checks ----

  checkmate::assert_choice(estimand_type, c("absolute", "relative"))
  checkmate::assert_choice(estimator_type, c("tmle", "one-step"))
  checkmate::assert_number(confidence_level, lower = 0.01, upper = 0.99)
  if (is.null(propensity_score_var_name)) {
    checkmate::assert_character(propensity_score_adj_var_names)
  }
  checkmate::assert_character(cond_exp_outcome_adj_var_names)
  checkmate::assert_character(treatment_var_name)
  checkmate::assert_character(outcome_var_name)
  checkmate::assert_int(num_nuisance_sl_folds, lower = 2, upper = 20)
  checkmate::assert_flag(cross_fit)
  checkmate::assert_int(num_cross_fit_folds, lower = 2, upper = 20)
  checkmate::assert_flag(parallel)

  # check dataset contains appropriate variables
  if (is.null(propensity_score_var_name)) {
    clean_tbl_var_names <- c(
      propensity_score_adj_var_names, cond_exp_outcome_adj_var_names,
      treatment_var_name, outcome_var_name
    )
  } else {
    checkmate::assert_character(propensity_score_var_name)
    checkmate::assert_numeric(
      data_tbl[[propensity_score_var_name]], lower = 0.001, upper = 0.999
    )
    clean_tbl_var_names <- c(
      propensity_score_adj_var_names, cond_exp_outcome_adj_var_names,
      treatment_var_name, propensity_score_var_name, outcome_var_name
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

  # compute estimates and EIF ----

  # estimate the differential variance estimand
  if (!cross_fit) {

    # estimate the nuisance parameters
    if (is.null(propensity_score_var_name)) {
      propensity_score_sl_fit <- estimate_propensity_score_fun(
        clean_tbl,
        adj_set_var_names = propensity_score_adj_var_names,
        treatment_var_name = treatment_var_name,
        propensity_score_library = propensity_score_library,
        num_nuisance_sl_folds = num_nuisance_sl_folds
      )
    } else {
      propensity_score_sl_fit <- NULL
    }
    cond_exp_outcome_sl_fit <- estimate_cond_exp_outcome_fun(
      clean_tbl,
      adj_set_var_names =  cond_exp_outcome_adj_var_names,
      treatment_var_name = treatment_var_name,
      outcome_var_name = outcome_var_name,
      cond_exp_outcome_library = cond_exp_outcome_library,
      num_nuisance_sl_folds = num_nuisance_sl_folds
    )
    cond_exp_sq_outcome_sl_fit <- estimate_cond_exp_sq_outcome_fun(
      clean_tbl,
      adj_set_var_names =  cond_exp_outcome_adj_var_names,
      treatment_var_name = treatment_var_name,
      outcome_var_name = outcome_var_name,
      cond_exp_sq_outcome_library = cond_exp_sq_outcome_library,
      num_nuisance_sl_folds = num_nuisance_sl_folds
    )

    # estimate the differential variance
    if (estimator_type == "one-step") {

      diff_var_est <- one_step_diff_var_estimator_fun(
        clean_tbl,
        propensity_score_adj_var_names = propensity_score_adj_var_names,
        cond_exp_outcome_adj_var_names = cond_exp_outcome_adj_var_names,
        treatment_var_name = treatment_var_name,
        propensity_score_var_name = propensity_score_var_name,
        outcome_var_name = outcome_var_name,
        propensity_score_sl_fit = propensity_score_sl_fit,
        cond_exp_outcome_sl_fit = cond_exp_outcome_sl_fit,
        cond_exp_sq_outcome_sl_fit = cond_exp_sq_outcome_sl_fit,
        estimand_type = estimand_type
      )

    } else if (estimator_type == "tmle") {

      diff_var_est <- tml_diff_var_estimator_fun(
        clean_tbl,
        propensity_score_adj_var_names = propensity_score_adj_var_names,
        cond_exp_outcome_adj_var_names = cond_exp_outcome_adj_var_names,
        treatment_var_name = treatment_var_name,
        propensity_score_var_name = propensity_score_var_name,
        outcome_var_name = outcome_var_name,
        propensity_score_sl_fit = propensity_score_sl_fit,
        cond_exp_outcome_sl_fit = cond_exp_outcome_sl_fit,
        cond_exp_sq_outcome_sl_fit = cond_exp_sq_outcome_sl_fit,
        estimand_type = estimand_type
      )

    }

    # extract estimate and EIF
    estimate <- diff_var_est$estimate
    eif <- diff_var_est$eif

  } else {

    # create folds for cross-fitting
    folds <- origami::make_folds(
      clean_tbl,
      fold_fun = origami::folds_vfold,
      V = num_cross_fit_folds
    )

    # cross-validate the differential variance estimators
    cf_diff_var_est <- origami::cross_validate(
      cv_fun = cf_diff_var_estimator_fun,
      folds = folds,
      clean_tbl = clean_tbl,
      propensity_score_adj_var_names = propensity_score_adj_var_names,
      cond_exp_outcome_adj_var_names = cond_exp_outcome_adj_var_names,
      treatment_var_name = treatment_var_name,
      propensity_score_var_name = propensity_score_var_name,
      outcome_var_name = outcome_var_name,
      propensity_score_library = propensity_score_library,
      cond_exp_outcome_library = cond_exp_outcome_library,
      cond_exp_sq_outcome_library= cond_exp_sq_outcome_library,
      num_nuisance_sl_folds = num_nuisance_sl_folds,
      estimator_type = estimator_type,
      estimand_type = estimand_type,
      use_future = parallel
    )

    # estimate the differential variance and assemble the EIF
    estimate <- mean(cf_diff_var_est$estimates)
    eif <- cf_diff_var_est$eif

  }


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

  estimator <- paste0(
    ifelse(cross_fit, "cross-fit ", ""),
    estimator_type
  )
  dplyr::tibble(
    estimand = paste(estimand_type, "differential variance"),
    estimate = estimate,
    se = se,
    ci_low = ci_low,
    ci_high = ci_high,
    p_value = p_value
  )

}
